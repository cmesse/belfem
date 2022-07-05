//
// Created by Christian Messe on 01.09.19.
//

#include "constants.hpp"
#include "assert.hpp"
#include "fn_sum.hpp"
#include "fn_dot.hpp"
#include "fn_gesv.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"

#include "cl_Gas.hpp"
#include "GT_globals.hpp"
#include "cl_GT_RefGas.hpp"
#include "cl_GT_RefGasFactory.hpp"

#include "cl_GM_Statevals.hpp"
#include "cl_GM_EoS.hpp"
#include "cl_GM_EoS_Idgas.hpp"
#include "cl_GM_EoS_Cubic.hpp"
#include "cl_GM_EoS_Hydrogen.hpp"
#include "cl_GM_EoS_Oxygen.hpp"
#include "cl_GM_EoS_Methane.hpp"
#include "cl_GM_HelmholtzTransport.hpp"
#include "cl_GM_HelmholtzTransport_Methane.hpp"
#include "cl_Timer.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Gas::Gas() :
            mHeatSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mViscositySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mConductivitySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax )
    {
        // assume air as default gas
        this->initialize( {
                                  "Ar", "CH4", "CO", "CO2", "H2", "He", "Kr",
                                  "N2", "N2O", "Ne", "NO2", "O2", "O3", "Xe" },
                          {
                            9.34e-3,   // Ar
                            2.0e-6,    // CH4
                            2.5e-7,    // CO
                            3.14e-4,   // CO2
                            5.0e-7,    // H2
                            5.24e-6,   // He
                            1.14e-6,   // Kr
                            7.8084e-1, // N2
                            3.1e-7,    // N2O
                            1.818e-5,  // Ne
                            2.0e-8,    // NO2
                            2.0948e-1, // O2
                            8.0e-6,    // O3
                            8.7e-8     // Xe
                          },
                          GasModel::IDGAS );

        /*this->initialize( {"Ar", "N2", "O2"},
        { 0.00934, 0.78084, 0.20942 },
                          GasModel::IDGAS ); */

        //this->initialize( {"Ar", "N2", "O2", "Ar", "CO2"}
        // { 0.78084, 0.209476, 0.009365, 0.000319 } );

    }

//----------------------------------------------------------------------------

    Gas::Gas( const string & aLabel, const GasModel aGasModel ):
            mHeatSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mViscositySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mConductivitySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax )
    {
        Cell< string > tSpecies( 1, aLabel );
        Vector< real > tMolarFractions( 1, 1.0 );

        if( aGasModel == GasModel::HELMHOLTZ )
        {
            if( aLabel == "H2" )
            {
                mHelmholzModel = HelmholtzModel::NormalHydrogen ;
            }
            else if ( aLabel == "O2" )
            {
                mHelmholzModel = HelmholtzModel::Oxygen ;
            }
            else if ( aLabel == "CH4" )
            {
                mHelmholzModel = HelmholtzModel::Methane ;
            }
        }

        this->initialize( tSpecies, tMolarFractions, aGasModel );
    }

//------------------------------------------------------------------------------

    Gas::Gas( const HelmholtzModel aHelmholtzModel ) :
            mHeatSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mViscositySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mConductivitySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mHelmholzModel( aHelmholtzModel )
    {
        string tLabel = "" ;
        switch( aHelmholtzModel )
        {
            case( HelmholtzModel::OrthoHydrogen ):
            case( HelmholtzModel::NormalHydrogen ):
            case( HelmholtzModel::ParaHydrogen ) :
            {
                tLabel = "H2" ;
                break ;
            }
            case( HelmholtzModel::Oxygen ) :
            {
                tLabel = "O2" ;
                break ;
            }
            case( HelmholtzModel::Methane ) :
            {
                tLabel = "CH4" ;
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid Helmholz model");
                break ;
            }
        }

        Cell< string > tSpecies( 1, tLabel );
        Vector< real > tMolarFractions( 1, 1.0 );
        this->initialize( tSpecies, tMolarFractions, GasModel::HELMHOLTZ );
    }

//------------------------------------------------------------------------------

    Gas::Gas(
            const Cell<string> & aSpecies,
            const Vector<real> & aMolarFractions,
            const GasModel aGasModel ) :
            mHeatSpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mViscositySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax ),
            mConductivitySpline( gastables::gNumberOfSplinePoints, 0.0, gastables::gTmax )
    {
        // this catches a bug or wrong usage of one-entry vector
        if( aSpecies.size() == 1 )
        {
            Vector< real > tMolarFractions( 1, 1.0 );
            this->initialize( aSpecies, tMolarFractions, aGasModel );
        }
        else
        {
            this->initialize( aSpecies, aMolarFractions, aGasModel );
        }
    }

//------------------------------------------------------------------------------

    Gas::~Gas()
    {
        // delete components
        for ( gastables::RefGas * tRefGas : mComponents )
        {
            delete tRefGas;
        }

        // delete reference gases that are not in component list
        for ( gastables::RefGas * tRefGas : mExtra )
        {
            delete tRefGas;
        }

        // delete viscosity table
        for ( gastables::RefGas * tRefGas : mViscosityInteractionRefgas )
        {
            delete tRefGas;
        }

        if( mTransport != nullptr )
        {
            delete mTransport ;
        }

        // delete the equation of state
        delete mEoS;
    }

//------------------------------------------------------------------------------

    void
    Gas::initialize(
            const Cell<string> & aSpecies,
            const Vector<real> & aMolarFractions,
            const GasModel aGasModel )
    {
            // remember gas model
            mGasModel = aGasModel ;

            // number of components
            mNumberOfComponents = aSpecies.size();

            // create the reference gas factory
            gastables::RefGasFactory tFactory;

            // make sure that molar fractions fit
            BELFEM_ASSERT( mNumberOfComponents == aMolarFractions.length(),
                          "Length of Gas names and molar fractions does not match ( %u and %u )",
                          ( unsigned int ) mNumberOfComponents,
                          ( unsigned int ) aMolarFractions.length() );

            // allocate work matrix
            mWorkMatrix.set_size( mNumberOfComponents, mNumberOfComponents, 0.0 );
            mWorkVector.set_size( mNumberOfComponents, 0.0 );
            mWorkVector2.set_size( mNumberOfComponents, 0.0 );

            // create the reference gases
            this->create_reference_gases( tFactory, aSpecies );

            // create the viscosity table
            this->create_viscosity_table( tFactory );

            // set number of mass fractions
            mMassFractions.set_size( mNumberOfComponents );

            // allocate mass fractions and molar masses and populate latter
            this->create_mass_properties( aMolarFractions );

            // create the temperature steps for spline creation
            tFactory.create_temperature_steps( mTemperatureSteps );

            // create the help matrix for remixing
            tFactory.create_helpmatrix( mHelpMatrix );

            // reserve memory for heat spline
            mHeatSpline.matrix_data().set_size( 5, gastables::gNumberOfSplinePoints );

            // allocate memory for working vectors
            mWorkMu.set_size( gastables::gNumberOfSplinePoints );
            mWorkLambda.set_size( gastables::gNumberOfSplinePoints );

            // create equation of state
            this->create_eos( aGasModel );

            // make sure that thermo exisis for all components
            this->check_thermo_exists();

            // create the mixture
            this->remix( aMolarFractions, true, true );

            // remember initialization values
            mMolarFractions0 = aMolarFractions ;
    }

//------------------------------------------------------------------------------

    void
    Gas::create_reference_gases(
            gastables::RefGasFactory & aFactory,
            const Cell<string> & aLables )
    {

        // temporary map to remember which gases have been created
        Map<string, uint> tRefgasMap;

        // create the reference gases
        mComponents.set_size( mNumberOfComponents, nullptr );

        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            const string & tLabel = aLables( k );

            // test if entry exists already ( some gases may exist twice, eg reacting and intert H2 )
            if ( !tRefgasMap.key_exists( tLabel ))
            {
                // add entry to temporary map
                tRefgasMap[ tLabel ] = k;
            }

            // create the new component
            mComponents( k ) = aFactory.create_refgas( tLabel );
        }

        // create list of components for formation enthalpy
        Cell<string> tFormation;

        for ( gastables::RefGas * tRefGas : mComponents )
        {
            const Cell<string> & tElements = tRefGas->data()->elements();
            for ( uint k = 0; k < tElements.size(); ++k )
            {
                tFormation.push( tElements( k ));
            }
        }

        unique( tFormation );

        mElements.set_size( tFormation.size(), nullptr );
        mElementNames.set_size( tFormation.size(), "" );

        // count extra elements
        uint tExtraCount = 0;
        for ( string tLabel : tFormation )
        {
            if ( !tRefgasMap.key_exists( reference_element( element_to_molecule( tLabel ))))
            {
                tExtraCount++;
            }
        }

        mExtra.set_size( tExtraCount, nullptr );
        tExtraCount = 0;
        uint tCount = 0;

        for ( string tLabel : tFormation )
        {
            // get name of reference element
            string tReferenceName = reference_element( element_to_molecule( tLabel ));

            // remember original name of element
            mElementNames( tCount ) = tLabel;

            // test if gas has already been created
            if ( tRefgasMap.key_exists( tReferenceName ))
            {
                mElements( tCount++ ) = mComponents( tRefgasMap( tReferenceName ));
            }
            else
            {
                // create a new gas
                gastables::RefGas * tRefGas = aFactory.create_refgas( tReferenceName );

                // add new gas to elements list
                mElements( tCount++ ) = tRefGas;

                // add gas to extra list
                mExtra( tExtraCount++ ) = tRefGas;
            }

        }

        // create the lookup table for formation enthalpy and entropy
        this->create_formation_table();

    }

//------------------------------------------------------------------------------

    string
    Gas::element_to_molecule( const string & aElement )
    {
        // special case for electron
        if ( aElement == "E" )
        {
            return string( "e-" );
        }
        // there are only 7 two atomic gases. Otherwise multiplicity is 1
        else if ( aElement == "Br" )
        {
            return string( "Br2" );
        }
        else if ( aElement == "Cl" )
        {
            return string( "Cl2" );
        }
        else if ( aElement == "F" )
        {
            return string( "F2" );
        }
        else if ( aElement == "I" )
        {
            return string( "I2" );
        }
        else if ( aElement == "O" )
        {
            return string( "O2" );
        }
        else if ( aElement == "N" )
        {
            return string( "N2" );
        }
        else if ( aElement == "H" )
        {
            return string( "H2" );
        }
        else
        {
            return aElement;
        }
    }

//------------------------------------------------------------------------------

    string
    Gas::reference_element( const string & aElement )
    {
        if ( aElement == "Ag" )
        {
            return string( "Ag(cr)" );
        }
        else if ( aElement == "Al" )
        {
            return string( "Al(cr)" );
        }
        else if ( aElement == "B" )
        {
            return string( "B(b)" );
        }
        else if ( aElement == "Ba" )
        {
            return string( "Ba(cr)" );
        }
        else if ( aElement == "Be" )
        {
            return string( "Be(cr)" );
        }
        else if ( aElement == "Br2" )
        {
            return string( "Br2(cr)" );
        }
        else if ( aElement == "C" )
        {
            return string( "C(gr)" );
        }
        else if ( aElement == "Cd" )
        {
            return string( "Cd(cr)" );
        }
        else if ( aElement == "Co" )
        {
            return string( "Co(cr)" );
        }
        else if ( aElement == "Cr" )
        {
            return string( "Cr(cr)" );
        }
        else if ( aElement == "Cs" )
        {
            return string( "Cs(cr)" );
        }
        else if ( aElement == "Cu" )
        {
            return string( "Cu(cr)" );
        }
        else if ( aElement == "Ga" )
        {
            return string( "Ga(cr)" );
        }
        else if ( aElement == "Ge" )
        {
            return string( "Ge(cr)" );
        }
        else if ( aElement == "Hg" )
        {
            return string( "Hg(cr)" );
        }
        else if ( aElement == "I2" )
        {
            return string( "I2(cr)" );
        }
        else if ( aElement == "In" )
        {
            return string( "In(cr)" );
        }
        else if ( aElement == "K" )
        {
            return string( "K(cr)" );
        }
        else if ( aElement == "Li" )
        {
            return string( "Li(cr)" );
        }
        else if ( aElement == "Mg" )
        {
            return string( "Mg(cr)" );
        }
        else if ( aElement == "Mo" )
        {
            return string( "Mo(cr)" );
        }
        else if ( aElement == "Na" )
        {
            return string( "Na(cr)" );
        }
        else if ( aElement == "Nb" )
        {
            return string( "Nb(cr)" );
        }
        else if ( aElement == "Ni" )
        {
            return string( "Ni(cr)" );
        }
        else if ( aElement == "P" )
        {
            return string( "P(cr)" );
        }
        else if ( aElement == "Pb" )
        {
            return string( "Pb(cr)" );
        }
        else if ( aElement == "Rb" )
        {
            return string( "Rb(cr)" );
        }
        else if ( aElement == "Si" )
        {
            return string( "Si(cr)" );
        }
        else if ( aElement == "Sn" )
        {
            return string( "Sn(cr)" );
        }
        else if ( aElement == "Ta" )
        {
            return string( "Ta(cr)" );
        }
        else if ( aElement == "V" )
        {
            return string( "V(cr)" );
        }
        else if ( aElement == "W" )
        {
            return string( "W(cr)" );
        }
        else if ( aElement == "Zn" )
        {
            return string( "Zn(cr)" );
        }
        else
        {
            return aElement;
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::create_viscosity_table( gastables::RefGasFactory & aFactory )
    {
        // allocate matrix and fill it with mNumberOfComponents
        mViscosityInteractionTable.set_size(
                mNumberOfComponents,
                mNumberOfComponents,
                mNumberOfComponents );


        // reset counter
        uint tCount = 0;

        // loop over all gases
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            const string & tA = mComponents( k )->label();

            for ( uint i = k + 1; i < mNumberOfComponents; ++i )
            {
                const string & tB = mComponents( i )->label();

                // test if interaction table exists
                if ( aFactory.interaction_viscosity_exists( tA, tB ))
                {
                    mViscosityInteractionRefgas.push(
                            aFactory.create_interaction_viscosity( tA, tB ));

                    // set counter into table
                    mViscosityInteractionTable( i, k ) = tCount;
                    mViscosityInteractionTable( k, i ) = tCount;

                    // increment counter
                    ++tCount;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::create_mass_properties( const Vector<real> & aMolarFractions )
    {
        mMassFractions.set_size( mNumberOfComponents );
        mMolarMasses.set_size( mNumberOfComponents );

        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            mMolarMasses( k ) = mComponents( k )->data()->M();
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::remix( const Vector<real> & aMolarFractions,
                bool aRemixHeat,
                bool aRemixTransport )
    {

        // reset work temperature value ( needed for viscosity calculation )
        mWorkTemperature = BELFEM_REAL_MAX;

        // recalculate R and M
        this->remix_R( aMolarFractions );

        if ( aRemixHeat )
        {
            mEoS->remix();

            // remix heat spline
            this->remix_heat();

            // update critical point
            this->remix_critical_point();
        }

        if ( aRemixTransport )
        {
            // remix transport spline
            this->remix_transport();
        }
    }
//------------------------------------------------------------------------------

    void
    Gas::remix_mass( const Vector<real> & aMassFractions,
                     bool aRemixHeat,
                     bool aRemixTransport )
    {
        // make sure that input is OK
        BELFEM_ASSERT( aMassFractions.length() == mNumberOfComponents,
                      "size of mass vector does not match ( %u vs %u )",
                      ( unsigned int ) aMassFractions.length(),
                      ( unsigned int ) mNumberOfComponents );


        // reset work temperature value ( needed for viscosity calculation )
        mWorkTemperature = BELFEM_REAL_MAX;

        mMassFractions = aMassFractions / sum( aMassFractions );

        // compute new molar fractions
        for( uint k=0; k<mNumberOfComponents; ++k )
        {
            mMolarFractions( k ) = mMassFractions( k ) / mMolarMasses( k );
        }

        // make molar fractions partition of unity
        mMolarFractions /= sum( mMolarFractions );

        // fixme: use inline multiplication instead
        mMassFractions = mMolarFractions % mMolarMasses;

        // molar mass
        mStatevals.set( BELFEM_STATEVAL_M, sum( mMassFractions ) );

        // mass fractions, part 2
        mMassFractions /= mM;

        // gas constant
        mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / mM );

        if( aRemixHeat )
        {
            // remix equation of state
            mEoS->remix();

            // remix heat spline
            this->remix_heat();

            // update critical point
            this->remix_critical_point();
        }

        if( aRemixTransport )
        {
            this->remix_transport();
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::reset_mixture()
    {
        this->remix( mMolarFractions0 );
    }

//------------------------------------------------------------------------------

    const real &
    Gas::R( const real & aT, const real & aP )
    {
        return mR;
    }


//------------------------------------------------------------------------------

    const real &
    Gas::M( const real & aT, const real & aP )
    {
        return mM;
    }

//------------------------------------------------------------------------------

    void
    Gas::remix_R( const Vector<real> & aMolarFractions )
    {
        // make sure that input is OK
        BELFEM_ASSERT( aMolarFractions.length() == mNumberOfComponents,
                      "size of molar vector does not match ( %u vs %u )",
                      ( unsigned int ) aMolarFractions.length(),
                      ( unsigned int ) mNumberOfComponents );


        // copy molar fractions
        mMolarFractions = aMolarFractions;

        // make molar fractions partition of unity
        mMolarFractions /= sum( mMolarFractions );

        // fixme: use inline multiplication instead
        mMassFractions = mMolarFractions % mMolarMasses;

        // molar mass
        mStatevals.set( BELFEM_STATEVAL_M, sum( mMassFractions ) );

        // mass fractions, part 2
        mMassFractions /= mM;

        // gas constant
        mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / mM );
    }

//------------------------------------------------------------------------------

    // create eos
    void
    Gas::create_eos( const GasModel & aGasModel )
    {
        switch ( aGasModel )
        {
            case ( GasModel::IDGAS ) :
            {
                mEoS = new gasmodels::EoS_Idgas( *this );
                this->link_to_idgas_property_functions();
                break;
            }
            case ( GasModel::SRK ) :
            case ( GasModel::PR ) :
            {
                mEoS = new gasmodels::EoS_Cubic( *this, aGasModel );
                this->link_to_realgas_property_functions();
                break;
            }
            case( GasModel::HELMHOLTZ ) :
            {
                BELFEM_ERROR( mNumberOfComponents == 1,
                             "a Gas can only have one component if the Helmholtz EoS is used" );

                // we need to write these values manually, because remix has not been called yet
                mStatevals.set( BELFEM_STATEVAL_M, mComponents( 0 )->M() ) ;
                mStatevals.set( BELFEM_STATEVAL_R, constant::Rm / mComponents( 0 )->M() ) ;

                // allocate mass and molar fractions
                mMassFractions.set_size( 1, 1.0 );
                mMolarFractions.set_size( 1, 1.0 );
                mMolarFractions0.set_size( 1, 1.0 );

                switch( mHelmholzModel )
                {
                    case( HelmholtzModel::NormalHydrogen ) :
                    case( HelmholtzModel::ParaHydrogen   ) :
                    case( HelmholtzModel::OrthoHydrogen  ) :
                    {
                        mEoS = new gasmodels::EoS_Hydrogen( *this, mHelmholzModel );

                        // dummy object, real thing not implemented yet
                        mTransport = new gasmodels::HelmholtzTransport( *this );


                        break ;
                    }
                    case( HelmholtzModel::Oxygen ) :
                    {
                        mEoS = new gasmodels::EoS_Oxygen( *this );

                        // dummy object, real thing not implemented yet
                        mTransport = new gasmodels::HelmholtzTransport( *this );

                        break ;
                    }
                    case( HelmholtzModel::Methane ) :
                    {
                        mEoS = new gasmodels::EoS_Methane( *this );
                        mTransport = new gasmodels::HelmholtzTransport_Methane( *this );

                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "unknown helmholtz model" );
                        break;
                    }
                }
                this->link_to_helmholtz_property_functions() ;
                break ;
            }
            default :
            {
                BELFEM_ERROR( false, "unknown gas model" );
                break;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::remix_heat()
    {
        // get data object
        Matrix<real> & tData = mHeatSpline.matrix_data();

        // reset heat spline
        tData.fill( 0.0 );

        // loop over all gases
        for ( size_t k = 0; k < mNumberOfComponents; ++k )
        {
            // get component data
            Matrix<real> & tRefGasData = mComponents( k )->heat_spline()->matrix_data();

            // loop over all temperature steps
            for ( size_t j = 0; j < gastables::gNumberOfSplinePoints; ++j )
            {
                for ( size_t i = 0; i < 5; ++i )
                {
                    tData( i, j ) += mMolarFractions( k ) * tRefGasData( i, j );
                }
            }
        }

        // scale unit to kJ/kg
        tData /= mM;

        // additional term for idgas_s
        this->update_mixture_entropy() ;

        mSref = this->idgas_s( gastables::gTref, gastables::gPref );
    }
//------------------------------------------------------------------------------

    void
    Gas::remix_critical_point()
    {
        mEoS->eval_critical_point( mTcrit, mPcrit, mVcrit );

        // stiel thodos parameter
        mGamma = std::pow( mTcrit, 1.0 / 6.0 )
                 * std::pow( mM * 1000, 0.5 )               // scale: kg->g
                 * std::pow( mPcrit / 1.01325e5, -2.0 / 3.0 ) // scale: bar->atm
                 * std::pow( mPcrit * mVcrit / ( mR * mTcrit ), 5 ) // Z_crit
                 / 4184; // scale:J -> cal

        // lucas parameter ( VDI Da 22 - 91 )
        mXi = 0.176 * std::pow( mTcrit, 1.0 / 6.0 )
              * std::pow( mM * 1000, -0.5 )
              * std::pow( mPcrit / 1e5, -2.0 / 3.0 );

    }

//------------------------------------------------------------------------------

    void
    Gas::remix_transport()
    {
        real tT = 0;

        // loop over all points
        for ( uint k = 0; k < gastables::gNumberOfSplinePoints; ++k )
        {
            // evaluate viscosity
            mWorkMu( k ) = this->cea_mu( tT );

            // evaluate thermal conductivity
            mWorkLambda( k ) = this->cea_lambda( tT );

            // increment t
            tT += gastables::gDeltaT;
        }

        // update splines
        mViscositySpline.update_data( mHelpMatrix, mWorkMu );
        mConductivitySpline.update_data( mHelpMatrix, mWorkLambda );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_cp( const real & aT, const real & aP )
    {
        return mHeatSpline.deval( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dcpdT( const real & aT, const real & aP )
    {
        return mHeatSpline.ddeval( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_cv( const real & aT, const real & aP )
    {
        return this->cp( aT, aP ) - this->R( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_gamma( const real & aT, const real & aP )
    {
        real tCp = this->cp( aT, aP );
        return tCp / ( tCp - this->R( aT, aP ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_c( const real & aT, const real & aP )
    {
        return std::sqrt( this->idgas_gamma( aT, aP ) * this->R( aT, aP ) * aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_h( const real & aT, const real & aP )
    {
        return mHeatSpline.eval( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_s( const real & aT, const real & aP )
    {
        /*return mHeatSpline.entropy( aT )
            + this->R( aT, aP ) * std::log( gastables::gPref / aP ); */

        return ( std::log( gastables::gPref / aP ) + mMixtureEntropy )
            * this->R( aT, aP )  + mHeatSpline.entropy( aT ) ;
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dsdT( const real & aT, const real & aP )
    {
        return mHeatSpline.dentropy( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dsdp( const real & aT, const real & aP )
    {
        return -this->R( aT, aP ) / aP;
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_mu( const real & aT, const real & aP )
    {
        return mViscositySpline.eval( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_lambda( const real & aT, const real & aP )
    {
        return mConductivitySpline.eval( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_cp( const real & aT, const real & aP )
    {
        return mHeatSpline.deval( aT )
               + mEoS->cpdep( aT, aP )
               - mEoS->cpdep0( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dcpdT( const real & aT, const real & aP )
    {
        real tCPDEP2 = mEoS->cpdep( aT+1.0, aP )  - mEoS->cpdep0( aT+1.0 );
        real tCPDEP1 = mEoS->cpdep( aT-1.0, aP )  - mEoS->cpdep0( aT-1.0 );

        return mHeatSpline.ddeval( aT ) + 0.5 * ( tCPDEP2 - tCPDEP1 );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_cv( const real & aT, const real & aP )
    {
        // ( 2.11 )
        return this->cp( aT, aP ) - aP * this->v( aT, aP ) * aT
                                    * this->alpha( aT, aP ) * this->beta( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_gamma( const real & aT, const real & aP )
    {
        // ( 2.20 )
        return this->cp( aT, aP ) * this->beta( aT, aP ) /
               ( this->cv( aT, aP ) * this->alpha( aT, aP ));
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_c( const real & aT, const real & aP )
    {
        // ( 2. 151 )
        return std::sqrt( aP * this->v( aT, aP ) * this->beta( aT, aP ) /
                          ( this->alpha( aT, aP ) * this->cv( aT, aP )));
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_h( const real & aT, const real & aP )
    {
        return mHeatSpline.eval( aT )
               + mEoS->hdep( aT, aP )
               - mEoS->hdep0( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_s( const real & aT, const real & aP )
    {
        return mHeatSpline.entropy( aT )
               + mR * ( std::log( gastables::gPref / aP ) )
               + mEoS->sdep( aT, aP )
               - mEoS->sdep0( aT )
               + mSref;
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dsdT( const real & aT, const real & aP )
    {
        return mHeatSpline.dentropy( aT )
               + mEoS->dsdepdT( aT, aP )
               - mEoS->dsdepdT0( aT );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dsdp( const real & aT, const real & aP )
    {
        return -mR / aP + mEoS->dsdepdp( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_mu( const real & aT, const real & aP )
    {
        real aMu = mViscositySpline.eval( aT );

        aMu += this->mu_dep( aMu, aT, aP );

        BELFEM_ASSERT( aMu > 0.0 && aMu < 1.0,
                      "Error in mu" );

        return aMu;
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_lambda( const real & aT, const real & aP )
    {

        real aLambda = mConductivitySpline.eval( aT )
                       + this->lambda_dep( aT, aP );

        BELFEM_ASSERT( aLambda > 0.0 && aLambda < 1.0,
                      "Error in lambda" );

        return aLambda;
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_cp( const real & aT, const real & aP )
    {
        return mEoS->cp( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dcpdT( const real & aT, const real & aP )
    {
        real tT1 = aT * 0.9999 ;
        real tT2 = aT * 1.0001 ;

        return ( mEoS->cp( tT2, aP ) - mEoS->cp( tT1, aP ) ) / ( tT2 - tT1 );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_cv( const real & aT, const real & aP )
    {
        return mEoS->cv( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_gamma( const real & aT, const real & aP )
    {
        return this->cp( aT, aP ) * this->beta( aT, aP ) /
           ( this->cv( aT, aP ) * this->alpha( aT, aP ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_c( const real & aT, const real & aP )
    {
        return mEoS->w( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_h( const real & aT, const real & aP )
    {
        return mEoS->h( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_s( const real & aT, const real & aP )
    {
        return mEoS->s( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dsdT( const real & aT, const real & aP )
    {
        return mEoS->dsdT( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dsdp( const real & aT, const real & aP )
    {
        return mEoS->dsdp( aT, aP ) ;
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_mu( const real & aT, const real & aP )
    {
        return mTransport->mu( aT, aP );
    }
//------------------------------------------------------------------------------

    real
    Gas::helmholtz_lambda( const real & aT, const real & aP )
    {
        return mTransport->lambda( aT, aP );
    }

//------------------------------------------------------------------------------
    void
    Gas::evaluate_viscosity_interaction( const real & aT )
    {
        // - - - - - - - - - - - - - - - - - - - -
        // Step 1: Update Viscosity of components
        // - - - - - - - - - - - - - - - - - - - -
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            mWorkVector( k ) = mComponents( k )->mu( aT );
        }

        uint k;

        // - - - - - - - - - - - - - - - - - - - - - - -
        // Step 2: Calculate interaction parameter Phi
        // - - - - - - - - - - - - - - - - - - - - - - -
        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            // evaluate Interaction
            if ( mComponents( i )->has_viscosity() )
            {
                const double & tM_i = mComponents( i )->data()->M();
                const double & tMu_i = mWorkVector( i );

                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    if ( i != j )
                    {
                        if ( mComponents( j )->has_viscosity() )
                        {
                            const double & tM_j = mComponents( j )->data()->M();
                            const double & tMu_j = mWorkVector( j );

                            k = mViscosityInteractionTable( i, j );
                            // test if interaction parameter exists
                            if ( k < mNumberOfComponents )
                            {
                                // NASA RP-1311 ( 5.7 )
                                mWorkMatrix( i, j ) = tMu_i /
                                                      mViscosityInteractionRefgas( k )->mu( aT )
                                                      * 2.0 * tM_j / ( tM_i + tM_j );
                            }
                            else
                            {
                                // NASA RP-1311 ( 5.5 )
                                mWorkMatrix( i, j ) = 0.25 *
                                                      std::pow(
                                                              1.0 +
                                                              std::sqrt(( tMu_i / tMu_j )
                                                                        * std::sqrt( tM_j / tM_i )), 2 )
                                                      * std::sqrt( 2.0 * tM_j / ( tM_i + tM_j ));
                            }
                        }
                        else
                        {
                            mWorkMatrix( i, j ) = 0.0;
                        }
                    }
                    else
                    {
                        mWorkMatrix( i, j ) = 0.0;
                    }
                } // end j-loop
            }
            else
            {
                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    mWorkMatrix( i, j ) = 0.0;
                }
            }
        } // end i-loop

        // remember temperature
        mWorkTemperature = aT;
    }

//------------------------------------------------------------------------------

    void
    Gas::evaluate_conductivity_interaction( const real & aT )
    {
        if ( aT != mWorkTemperature )
        {
            this->evaluate_viscosity_interaction( aT );
        }

        // - - - - - - - - - - - - - - - - - - - -
        // Step 1: Updatem Conductiviies
        // - - - - - - - - - - - - - - - - - - - -
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            mWorkVector2( k ) = mComponents( k )->lambda( aT );
        }

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            if ( mComponents( i )->has_conductivity() && mComponents( i )->has_viscosity())
            {
                //const double & tMu_i = mWorkVector( i );
                const double & tM_i = mComponents( i )->data()->M();

                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    if ( i != j )
                    {
                        if ( mComponents( j )->has_conductivity() && mComponents( i )->has_viscosity())
                        {
                            //const double & tMu_j = mWorkVector( j );
                            const double & tM_j = mComponents( j )->data()->M();
                            // NASA RP-1311 ( 5.6 )
                            mWorkMatrix( i, j ) *= 1.0 + 2.41 * ( tM_i - tM_j )
                                                         * ( tM_i - 0.142 * tM_j ) /
                                                         std::pow( tM_i + tM_j, 2 );

                            // alternative approach, but worse that Gordon/McBride
                            // Wassiljeva, Mason, Saxena, see VDI D1 ( 108a )
                            //mWorkMatrix( i, j ) = std::pow( 1.0 + std::sqrt( tMu_i / tMu_j
                            // * std::sqrt( tM_j / tM_i) ), 2 ) /
                            //        std::sqrt( 8.0 *( 1.0 + tM_i / tM_j ) );
                        }
                        else
                        {
                            mWorkMatrix( i, j ) = 0.0;
                        }

                    }
                    else
                    {
                        mWorkMatrix( i, j ) = 0.0;
                    }
                } // end j-loop
            }
            else
            {
                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    mWorkMatrix( i, j ) = 0.0;
                }
            }
        }

        // overwrite work temperature
        mWorkTemperature = BELFEM_REAL_MAX;
    }

//------------------------------------------------------------------------------

    real
    Gas::cea_mu( const real & aT )
    {

        this->evaluate_viscosity_interaction( aT );

        const Vector<real> & tX = mMolarFractions;

        // calculate value
        real aMu = 0.0;

        // temporary vector for mixing
        mWorkVector2 = mWorkMatrix * tX;

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            if ( tX( i ) > BELFEM_EPSILON_X && mComponents( i )->has_viscosity() )
            {

                aMu += tX( i ) * mWorkVector( i ) /
                       ( tX( i ) + mWorkVector2( i ));

            }

        }

        return aMu;
    }

//------------------------------------------------------------------------------

    real
    Gas::cea_lambda( const real & aT )
    {

        this->evaluate_conductivity_interaction( aT );

        const Vector<real> & tX = mMolarFractions;

        // calculate value
        real aLambda = 0.0;

        // temporary vector for mixing
        mWorkVector = mWorkMatrix * tX;

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            if ( tX( i ) > BELFEM_EPSILON_X && mComponents( i )->has_conductivity())
            {
                aLambda += tX( i ) * mWorkVector2( i ) /
                           ( tX( i ) + mWorkVector( i ));
            }

        }

        return aLambda;
    }

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Thermodynamic States
//------------------------------------------------------------------------------

    real
    Gas::p( const real & aT, const real & aV )
    {
        mStatevals.update_Tv( aT, aV );

        // check if values are up to date
        if ( !mStatevals.test( BELFEM_STATEVAL_P ))
        {
            mStatevals.set(
                    BELFEM_STATEVAL_P,
                    mEoS->p( aT, aV ));
        }

        return mStatevals.get( BELFEM_STATEVAL_P );
    }

//------------------------------------------------------------------------------

    real
    Gas::v( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_V ))
        {
            mStatevals.set( BELFEM_STATEVAL_V,
                            mEoS->v( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_V );
    }

//------------------------------------------------------------------------------

    real
    Gas::rho( const real & aT, const real & aP )
    {
        return 1.0 / this->v( aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::T( const real & aP, const real & aV )
    {
        mStatevals.update_pv( aP, aV );

        if ( !mStatevals.test( BELFEM_STATEVAL_T ))
        {
            mStatevals.set( BELFEM_STATEVAL_T,
                            mEoS->T( aP, aV ));
        }

        return mStatevals.get( BELFEM_STATEVAL_T );
    }

//------------------------------------------------------------------------------

    real
    Gas::cp( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_CP ))
        {
            mStatevals.set( BELFEM_STATEVAL_CP,
                            ( this->*mFunctionCp )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_CP );
    }

//------------------------------------------------------------------------------

    // dissociation enthalpy ( only for tablegas at this time )
    real
    Gas::hd( const real & aT, const real & aP )
    {
        return 0.0 ;
    }

//------------------------------------------------------------------------------

    real
    Gas::dcpdT( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_DCPDT ))
        {
            mStatevals.set( BELFEM_STATEVAL_DCPDT,
                            ( this->*mFunctiondCpdT )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DCPDT );
    }

//------------------------------------------------------------------------------

    real
    Gas::cv( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_CV ))
        {
            mStatevals.set( BELFEM_STATEVAL_CV,
                            ( this->*mFunctionCv )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_CV );
    }

//------------------------------------------------------------------------------

    real
    Gas::gamma( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_GAMMA ))
        {
            mStatevals.set( BELFEM_STATEVAL_GAMMA,
                            ( this->*mFunctionGamma )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_GAMMA );
    }

//------------------------------------------------------------------------------

    real
    Gas::c( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_C ))
        {
            mStatevals.set( BELFEM_STATEVAL_C,
                            ( this->*mFunctionC )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_C );
    }

//------------------------------------------------------------------------------

    real
    Gas::u( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_U ))
        {
            mStatevals.set( BELFEM_STATEVAL_U,
                            ( this->*mFunctionH )( aT, aP )
                            - aP * this->v( aT, aP ) );
        }

        return mStatevals.get( BELFEM_STATEVAL_U );
    }

//------------------------------------------------------------------------------

    real
    Gas::h( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_H ))
        {
            mStatevals.set( BELFEM_STATEVAL_H,
                            ( this->*mFunctionH )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_H );
    }

//------------------------------------------------------------------------------

    real
    Gas::s( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_S ))
        {
            mStatevals.set( BELFEM_STATEVAL_S,
                            ( this->*mFunctionS )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_S );
    }

//------------------------------------------------------------------------------

    real
    Gas::dsdT( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_DSDT ))
        {
            mStatevals.set( BELFEM_STATEVAL_DSDT,
                            ( this->*mFunctionDSDT )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DSDT );
    }

//------------------------------------------------------------------------------

    real
    Gas::dsdp( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_DSDP ))
        {
            mStatevals.set( BELFEM_STATEVAL_DSDP,
                            ( this->*mFunctionDSDP )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DSDP );
    }

//------------------------------------------------------------------------------

    real
    Gas::mu( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_MU ))
        {
            mStatevals.set( BELFEM_STATEVAL_MU,
                            ( this->*mFunctionMU )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_MU );
    }

//------------------------------------------------------------------------------

    real
    Gas::lambda( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_LAMBDA ))
        {
            mStatevals.set( BELFEM_STATEVAL_LAMBDA,
                            ( this->*mFunctionLAMBDA )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_LAMBDA );
    }

//------------------------------------------------------------------------------

    real
    Gas::Pr( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_PR ))
        {
            mStatevals.set( BELFEM_STATEVAL_PR,
                            this->cp( aT, aP ) * this->mu( aT, aP ) /
                            this->lambda( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_PR );
    }

//------------------------------------------------------------------------------

    // create the table needed for formation enthalpy
    void
    Gas::create_formation_table()
    {
        const uint tNumElements = mElements.size();
        mFormationTable.set_size( mNumberOfComponents, tNumElements, 0.0 );
        mFormationWork.set_size( tNumElements );

        for ( uint j = 0; j < tNumElements; ++j )
        {
            mFormationWork( j ) = mElements( j )->data()->component_multiplicity(
                    mElementNames( j ) );
        }

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            gastables::RefGas * tSpecie = mComponents( i );

            for ( uint j = 0; j < tNumElements; ++j )
            {
                mFormationTable( i, j ) = tSpecie->data()->component_multiplicity(
                        mElementNames( j )) / mFormationWork( j );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::Hf( const real & aT, Vector<real> & aHf )
    {
        // reset vector
        aHf.fill( 0.0 );

        // add species enthalpies to vector
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aHf( tCount++ ) +=  tSpecie->H( aT ) -  tSpecie->H_ref() + tSpecie->data()->Hf();
        }

        // calculate element enthalpies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->H( aT ) - tElement->H_ref();
        }

        aHf.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::Gibbs( const real & aT, Vector<real> & aGibbs )
    {
        // write species entropies into vector
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) = -tSpecie->S( aT );
        }

        // calculate element entropies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->S( aT );
        }

        // add entropies to vector
        aGibbs.vector_data() += mFormationTable * mFormationWork;

        // multiply vector with temperatures
        aGibbs *= aT;

        // add species enthalpies to vector
        tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) +=  tSpecie->H( aT ) -  tSpecie->H_ref() + tSpecie->data()->Hf();
        }

        // calculate element enthalpies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->H( aT ) - tElement->H_ref();
        }

        aGibbs.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::dGibbsdT( const real & aT, Vector< real > & aGibbs )
    {

        // calculate entropy term for species
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) = -( tSpecie->S( aT ) + aT * tSpecie->dSdT( aT ) );
        }

        // calculate element entropies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->S( aT ) + aT * tElement->dSdT( aT );
        }

        aGibbs.vector_data() += mFormationTable * mFormationWork;

        // add specific heats to vector
        tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) += tSpecie->Cp( aT );
        }

        // calculate specific heats
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->Cp( aT );
        }

        aGibbs.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::compute_equilibrium( const real & aT, const real & aP, Vector< real > & aX )
    {
        if( mNumberOfComponents > 1 )
        {
            // minuimum composition a gas may have
            real tEpsilon = 1E-9;

            // maximum number of iterations
            uint tMaxNumIterations = 1000;

            // relaxation factor
            real tOmega0 = 0.9;

            // we reset the work temperature, so that we can use the work vectors safely
            mWorkTemperature = BELFEM_REAL_MAX;

            // step 0: allocate memory
            const real tNumElements = mElements.size();

            if ( mPivotRAND.length() != mNumberOfComponents + 1 )
            {
                mWorkVectorRAND0.set_size( tNumElements + 1 );
                mWorkVectorRAND1.set_size( tNumElements );
                mWorkVectorRAND2.set_size( tNumElements );
                mWorkMatrixRAND.set_size( tNumElements + 1,
                                          tNumElements + 1 );
                mPivotRAND.set_size( tNumElements + 1 );
            }

            // Step 1 : Link variables

            // Gibbs potential at reference pressure
            Vector< real > & tMu0 = mWorkVector;

            // Gibbs potential at given pressure
            Vector< real > & tMu = mWorkVector2;

            // chemical potential
            const Vector< real > & tPsi = mWorkVectorRAND0;

            // value for last equation
            const real & tU = mWorkVectorRAND0( tNumElements );

            // formation table
            Matrix< real > & tA = mFormationTable;

            // System of Equations
            Matrix< real > & tM = mWorkMatrixRAND;

            // Right hand side of system
            Vector< real > & tRHS = mWorkVectorRAND0;

            // Element abundance vectors
            Vector< real > & tB0 = mWorkVectorRAND1;
            Vector< real > & tB = mWorkVectorRAND2;

            // change vector for X
            Vector< real > tDeltaX = mWorkVector2;

            // Step 2 : initial computations
            // Compute Gibbs potential at reference pressure
            this->Gibbs( aT, tMu0 );

            // compute mass balance constraint
            tB0 = trans( tA ) * aX;

            // start loop
            uint tCount = 0;


            //real tSumX ;
            real tNorm = 1.0;

            // avoid having zero components
            for ( uint k = 0; k < mNumberOfComponents; ++k )
            {
                if ( aX( k ) < tEpsilon )
                {
                    aX( k ) = tEpsilon;
                }
            }

            while ( tNorm > tEpsilon )
            {
                // compute sum of X ( may be > 1 )
                //tSumX = sum( tX );

                // compute Gibbs potential at reference pressure
                for ( uint k = 0; k < mNumberOfComponents; ++k )
                {
                    tMu( k ) = tMu0( k ) + constant::Rm * aT * std::log( aP / gastables::gPref * aX( k ));
                }


                // compute mass balance constraint
                tB = trans( tA ) * aX;

                // compute matrix to be solved
                tM.fill( 0.0 );

                for ( uint i = 0; i < tNumElements; ++i )
                {
                    for ( uint j = 0; j < tNumElements; ++j )
                    {
                        for ( uint k = 0; k < mNumberOfComponents; ++k )
                        {
                            tM( i, j ) += tA( k, i ) * tA( k, j ) * aX( k );
                        }
                        tM( tNumElements, j ) = tB( j );
                    }
                    tM( i, tNumElements ) = tB( i );
                }

                tM( tNumElements, tNumElements ) = 0.0;

                // compute RHS
                for ( uint i = 0; i < tNumElements; ++i )
                {
                    tRHS( i ) = tB0( i ) - tB( i );
                    for ( uint k = 0; k < mNumberOfComponents; ++k )
                    {
                        tRHS( i ) += tA( k, i ) * aX( k ) * tMu( k ) / ( constant::Rm * aT );
                    }
                }
                tRHS( tNumElements ) = dot( aX, tMu ) / ( constant::Rm * aT );

                // solve system
                gesv( mWorkMatrixRAND, mWorkVectorRAND0, mPivotRAND );

                // compute change of Mols
                for ( uint k = 0; k < mNumberOfComponents; ++k )
                {
                    tDeltaX( k ) = tU - tMu( k ) / ( constant::Rm * aT );
                    for ( uint i = 0; i < tNumElements; ++i )
                    {
                        tDeltaX( k ) += tA( k, i ) * tPsi( i );
                    }
                }

                real tMinDeltaX = std::abs( min( tDeltaX ));
                real tMaxDeltaX = std::abs( max( tDeltaX ));

                // relaxation factor
                real tOmega = tMaxDeltaX > tMinDeltaX ? tMaxDeltaX : tMinDeltaX;
                if ( tOmega > tOmega0 )
                {
                    tOmega = 1.0 / tOmega * tOmega0;
                }
                else
                {
                    tOmega = tOmega0;
                }

                // adapt X
                for ( uint k = 0; k < mNumberOfComponents; ++k )
                {
                    tDeltaX( k ) *= aX( k );
                }

                tNorm = norm( tDeltaX );

                aX += tOmega * tDeltaX;


                // check for infitite loop
                BELFEM_ERROR( tCount++ < tMaxNumIterations,
                             "To many iterations while trying to find chemical equilibrium." );
            }

            // cleanup
            for ( uint k = 0; k < mNumberOfComponents; ++k )
            {
                if ( aX( k ) < tEpsilon )
                {
                    aX( k ) = 0.0;
                }
            }

            // remix values
            aX /= sum( aX );
        }
        else
        {
            aX.set_size( 1, 1.0 );
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::remix_to_equilibrium( const real & aT, const real & aP,
                               bool aRemixHeat,
                               bool aRemixTransport )
    {
        if( mNumberOfComponents > 1 )
        {
            this->compute_equilibrium( aT, aP, mMolarFractions );
            this->remix( mMolarFractions, aRemixHeat, aRemixTransport );
        }
    }

//------------------------------------------------------------------------------

    real
    Gas::mu_dep( const real & aMu, const real & aT, const real & aP )
    {
        // cutoff value
        real tMuMax = 10.0 * aMu;

        // reduced temperature
        real tTr = aT / mTcrit;
        real tPr = aP / mPcrit;

        // help magnitude ( see Eq. 88 )
        real tX = ( 0.807 * std::pow( tTr, 0.618 )
                    - 0.357 * std::exp( -0.449 * tTr )
                    + 0.34 * std::exp( -4.058 * tTr )
                    + 0.018 ) * 1e-7;

        // correction factor
        real tFid = aMu * mXi / tX;

        if ( tTr < 1.0 )
        {
            real tA = 3.262 + 14.98 * std::pow( tPr, 5.508 );
            real tB = 1.39 + 5.746 * tPr;
            real tZ2 = 0.6 + 0.76 * std::pow( tPr, tA )
                       + ( 6.99 * std::pow( tPr, tB ) - 0.6 )
                         * ( 1.0 - tTr );

            if ( tZ2 > BELFEM_REAL_MAX)
            {
                return tMuMax;
            }
            else
            {
                real tFp = 1.0 + ( tFid - 1.0 ) * std::pow(( tZ2 * tFid ) /
                                                           ( mXi * aMu ), -3 );
                return std::min( 1.0e-7 * tZ2 * tFp / mXi - aMu, tMuMax );
            }
        }
        else
        {
            real tA = 0.001245 / tTr * std::exp( 5.1726 * std::pow( tTr, -0.3286 ));
            real tB = tA * ( 1.6553 * tTr - 1.2723 );
            real tC = 0.4489 / tTr * std::exp( 3.0578 * std::pow( tTr, -37.7332 ));
            real tD = 1.7368 / tTr * std::exp( 2.231 * std::pow( tTr, -7.6351 ));
            real tE = 1.3088;
            real tF = 0.9425 * std::exp( -0.1853 * std::pow( tTr, 0.4489 ));

            real tZ2 = 1.0 + ( tA * std::pow( tPr, tE )) /
                             ( tB * std::pow( tPr, tF ) + 1.0 / ( 1.0 + tC * std::pow( tPr, tD )));

            real tFp = ( 1.0 + ( tFid - 1 ) * std::pow( tZ2, 3 )) / tFid;

            return std::min( aMu * tZ2 * tFp - aMu, tMuMax );

        }
    }

//----------------------------------------------------------------------------

    real
    Gas::lambda_dep( const real & aT, const real & aP )
    {
        real tX = mVcrit / mEoS->v( aT, aP );

        // extrapolated from
        // 10.1002/aic.690100114 Fig 2
        real tC1 = 6.6224e-9;
        real tC2 = 4.079e-9;
        real tC3 = 1.992e-9;
        real tC4 = -1.5328e-9;
        real tC5 = 0.47344e-9;

        return tX * ( tC1 + tX * ( tC2 + tX * ( tC3 + tX * ( tC4 + tX * tC5 )))) / mGamma;
    }


//------------------------------------------------------------------------------

    real
    Gas::T_from_h( const real & aH, const real & aP )
    {
        // reference temperature
        real tTref ;

        if( this->is_idgas() )
        {
            tTref = gastables::gTref ;
        }
        else
        {
            tTref = mTcrit ;
        }

        // reference enthalpy
        real tHref = this->h( tTref, aP );

        real tCp = this->cp( tTref, aP );

        // initial guess
        real aT = ( aH - tHref ) / tCp + tTref;

        // for biseciton fallback
        real tT0 = 0.5 * aT;
        real tT1 = std::min( 2.0 * aT, gastables::gTmax );

        uint tCount = 0;

        real tT = 1;

        real tOmega = 0.95;

        aT = std::max( std::min( aT, 6000.0 ), 200.0 );

        while( std::abs( tT - aT ) > BELFEM_EPSILON_T  )
        {
            // shift t
            tT = aT;

            // newton step
            aT -= tOmega * ( this->h( aT, aP ) - aH ) / this->cp( aT, aP );

            // increment counter
            if( tCount++ == 100 )
            {
                // fallback to bisection
                real tF0 = this->h( tT0, aP ) - aH ;
                real tF;

                aT = ( aH - tHref ) / tCp + gastables::gTref;

                while( std::abs( tT0 - tT1 ) > BELFEM_EPSILON_T )
                {
                    aT = 0.5 * ( tT0 + tT1 );
                    tF =  this->h( tT1, aP ) - aH;

                    if( tF0 * tF >= 0 )
                    {
                        tT0 = aT;
                        tF0 = tF;
                    }
                    else
                    {
                        tT1 = aT;
                    }

                    BELFEM_ERROR( tCount++ < 1000,
                                 "T_from_h did not converge for h=%12.3f, p=%12.3f",
                                 aH, aP );

                }

                break;
            }


           BELFEM_ERROR( tCount < 1000,
                    "T_from_h did not converge for h=%12.3f, p=%12.3f",
                    aH, aP );
        }

        return aT;
    }

//------------------------------------------------------------------------------

    real
    Gas::isen_T( const real & aT0, const real & aP0, const real & aP1 )
    {
        // guess value for new temperature
        real aT1 = aT0 * std::pow( aP1 / aP0, this->R( aT0, aP0 ) / this->cp( aT0, aP0 ));

        // entropy at this state
        real tS = this->s( aT0, aP0 );

        real tT1 = 0.0;

        uint tCount = 0;

        real tDeltaT ;

        while ( std::abs( tT1 - aT1 ) > BELFEM_EPSILON_T )
        {
            tT1 = aT1;
            tDeltaT = ( this->s( aT1, aP1 ) - tS ) / this->dsdT( aT1, aP1 ) ;

            if( aT1 - tDeltaT < 15.0 )
            {
                tDeltaT = aT1 - 15.0 ;
            }
            aT1 -= 0.9 * tDeltaT ;

            ++tCount;


            BELFEM_ERROR( tCount < 1000,
                         "Too many iterations for isen_T ( T0=%f, p0=%f, p1=%f)",
                         ( float ) aT0, ( float ) aP0, ( float ) aP1 );
        }

        return aT1;
    }


// -----------------------------------------------------------------------------

    real
    Gas::isen_p( const real & aT0, const real & aP0, const real & aT1 )
    {
        // guess value for new pressure
        real aP1 = aP0 * std::pow( aT1 / aT0, this->idgas_cp( aT0, aP0 ) / this->R( aT0, aP0 ) );

        // entropy at this state
        real tS = this->s( aT0, aP0 );

        real tP1 = 0.0;

        uint tCount = 0;

        while ( std::abs( tP1 - aP1 ) > BELFEM_EPSILON_P )
        {
            tP1 = aP1;
            aP1 -= ( this->s( aT1, aP1 ) - tS ) / this->dsdp( aT1, aP1 );

            ++tCount;

            BELFEM_ERROR( tCount < 1000,
                         "Too many iterations for isen_p ( T0=%f, p0=%f, T1=%f)",
                         ( float ) aT0, ( float ) aP0, ( float ) aT1 );
        }

        return aP1;
    }

// -----------------------------------------------------------------------------

    void
    Gas::total( const real & aT, const real & aP, const real & aU,
                real & aTt, real & aPt )
    {
        // maximum temperature
        real tTmax = gastables::gTmax - BELFEM_EPSILON_T ;

        // Initial guesses
        real tCp = this->idgas_cp( aT, aP );

        // gas constant
        real tR = this->R( aT, aP );

        // guess ratio of specific heats
        real tGamma = tCp / ( tCp - tR );

        // guess speed of sound
        real tC = std::sqrt( tGamma * tR * aT );

        // guess mach number
        real tMa = aU / tC;

        // solution vector
        Vector<real> tX( 2 );

        real & tTt = tX( 0 );
        real & tPt = tX( 1 );

        // RHS
        Vector<real> tF( 2, 1.0 );

        // Jacobian
        Matrix<real> tJ( 2, 2 );

        // Pivot
        Vector< int > tPivot( 2 );

        // guess total temperature
        tTt = std::min( aT * ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa * tMa ), tTmax -100.0 );

        // guess total preassure
        tPt = aP * std::pow( tTt / aT, tCp / tR );

        // calculate entropy
        real tS = this->s( aT, aP );

        // calculate enthalpy
        real tH = this->h( aT, aP ) + 0.5 * aU * aU;

        // initialize loop counter
        uint tCount = 0;

        aTt = 0.0;
        aPt = 0.0;

        real tOmega0 = 0.99 ;
        real tOmega1;
        real tOmega ;

        while ( true )
        {
            // shift result
            aTt = tX( 0 );
            aPt = tX( 1 );

            // calculate right hand side
            tF( 0 ) = ( this->h( tTt, tPt ) - tH );
            tF( 1 ) = this->s( tTt, tPt ) - tS;

            // calculate Jacobian
            tJ( 0, 0 ) = this->cp( tTt, tPt );
            tJ( 1, 0 ) = this->dsdT( tTt, tPt );
            tJ( 0, 1 ) = this->dhdp( tTt, tPt );
            tJ( 1, 1 ) = this->dsdp( tTt, tPt );

            // solve system
            gesv( tJ, tF, tPivot );

            tOmega1 = 0.9 * std::abs( ( tX( 0 ) - tTmax ) / tF( 0 ) );

            tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;

            // correct result
            tX -= tOmega * tF;

            ++tCount;

            BELFEM_ERROR( tCount < 1000, "Infinite loop at total state calculation" );

            tOmega0 *= 0.99 ;

            // check abort condition
            if ( std::abs( tTt - aTt ) < BELFEM_EPSILON_T &&
                 std::abs( tPt - aPt ) < BELFEM_EPSILON_P )
            {
                // copy result into output
                aTt = tX( 0 );
                aPt = tX( 1 );
                break;
            }
        }


    }

//------------------------------------------------------------------------------
    void
    Gas::expand(
            const real & aA1,
            const real & aT1,
            const real & aP1,
            const real & aU1,
            const real & aA2,
                  real & aT2,
                  real & aP2,
                  real & aU2 )
    {

        aT2 = aT1;
        aP2 = aP1;
        aU2 = aU1;

        if( std::abs( aA1 - aA2 ) < BELFEM_EPSILON )
        {
            return;
        }
        else
        {
            real tOmega = 0.9;

            // mass
            real tMass = this->rho( aT1, aP1 ) * aU1 * aA1;

            // momentum
            real tMomentum = tMass * aU1 + aP1 * aA2; // <-- aA2, not aA1 !

            real tEntropy = this->s( aT1, aP1 );

            real tEnergy = this->h( aT1, aP1 ) + 0.5 * aU1 * aU1;

            // simplified constans for ideal gas


            real tMa1 = aU1 / this->c( aT1, aP1 );

            // RHS
            Vector<real> tF( 3, 0.0 );

            // Jacobian
            Matrix<real> tJ( 3, 3 );

            // Pivot
            Vector<int> tPivot( 3 );

            real tV2;

            uint tCount = 0;

            real tError = 1.0;

            // initial guess
            if( tMa1 < 1.0 )
            {
                real tK = this->gamma( aT1, aP1 );
                real tCp = this->cp( aT1, aP1 );

                // calculate a Borda-Carnot Shock
                // Rist: Dynamik Realer Gase, Kap. 9.1.1

                // Eq. ( 9.15 )
                aU2 = aU1 * ( tK * tMa1 * tMa1 + aA2 / aA1 -
                   std::sqrt( std::pow( tMa1 * tMa1 - 1.0, 2 )
                   + ( aA2 / aA1 - 1.0 ) *
                     ( 1.0 + 2.0 * tK * tMa1 * tMa1 + aA2 / aA1 ) ) ) /
                     ( ( tK + 1.0 ) * tMa1 * tMa1 );


                real tHt = tCp * aT1 + 0.5 * aU1 * aU1;

                aT2 = ( tHt - 0.5 * aU2 * aU2 ) / tCp;
                aP2 = ( tMass * ( aU1 - aU2 ) + aP1 * aA2 ) / aA2;


                while ( tError > 1e-6 )
                {
                    // get values
                    aT2 -= tOmega * tF( 0 );
                    aP2 -= tOmega * tF( 1 );
                    aU2 -= tOmega * tF( 2 );

                    // calculate inverse density
                    tV2 = this->v( aT2, aP2 );

                    // Jacobian Terms

                    // d(mass)/dT
                    tJ( 0, 0 ) = -aA2 * aU2 * this->alpha( aT2, aP2 ) / ( tV2 * tMass );

                    // d(momentum)/dT
                    tJ( 1, 0 ) = tJ( 0, 0 ) * aU2 * tMass / tMomentum;

                    // d(energy)/dT
                    tJ( 2, 0 ) = this->cp( aT2, aP2 ) / tEnergy;

                    // d(mass)/dP
                    tJ( 0, 1 ) = -this->kappa( aT2, aP2 ) * aA2 * aU2 / ( tV2 * tMass );

                    // d(momentum)/dP
                    tJ( 1, 1 ) = ( aA2 - aU2 * tJ( 0, 1 ) * tMass ) / tMomentum;

                    //  d(energy)/dP
                    tJ( 2, 1 ) = mEoS->hdep( aT2, aP2 ) / tEnergy;

                    // d(mass)/du
                    tJ( 0, 2 ) = aA2 / ( tV2 * tMass );

                    // d(momentum)/du
                    tJ( 1, 2 ) = 2.0 * aU2 * tJ( 0, 2 ) * tMass / tMomentum;

                    // d(energy)/du
                    tJ( 2, 2 ) = aU2 / tEnergy;


                    // RHS
                    tF( 0 ) = ( aU2 * aA2 / tV2 - tMass ) / tMass;
                    tF( 1 ) = ( aU2 * aU2 * aA2 / tV2 + aP2 * aA2 - tMomentum ) / tMomentum;
                    tF( 2 ) = ( this->h( aT2, aP2 ) + 0.5 * aU2 * aU2 - tEnergy ) / tEnergy;

                    tError = norm( tF );

                    // solve system
                    gesv( tJ, tF, tPivot );

                    // increment counter
                    ++tCount;

                    BELFEM_ERROR( tCount < 100,
                                 "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                                 aT1,
                                 aP1,
                                 aU1,
                                 aA1,
                                 aA2 );

                }

            }
            else // supersonic
            {
                tError = 1.0;
                real tT = -1;
                tCount = 0;

                while( std::abs( aT2 - tT ) > 0.001 )
                {
                    tV2 = aA2 * aU2 / tMass;

                    aP2 = mEoS->p( aT2, tV2 );

                    //  Böckh, Saumweber: Fluidmechanik: Einführendes Lehrbuch

                    // Eq. (7.2)
                    aU2 = aU1 + ( aP1 - aP2 ) * aA2 / tMass;


                    tT = aT2;

                    aT2 = this->T_from_h( tEnergy - 0.5 * aU2 * aU2, aP2 );

                    ++tCount;
                    BELFEM_ERROR( tCount < 100,
                                 "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                                 aT1,
                                 aP1,
                                 aU1,
                                 aA1,
                                 aA2 );

                }
                tCount = 0;
                tT = -1.0;

                while( std::abs( aT2 - tT ) > 0.0001 )
                {
                    tT = aT2;

                    aT2 = this->isen_T( aT1, aP1, aP2 );

                    real tV2 = aA2 * aU2 / tMass;
                    aP2 = this->p( aT2, tV2 );

                    aU2 =  std::sqrt( 2.0 * ( tEnergy - this->h( aT2, aP2 ) ) );
                    ++tCount;
                    BELFEM_ERROR( tCount < 100,
                                 "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                                 aT1,
                                 aP1,
                                 aU1,
                                 aA1,
                                 aA2 );
                }

                while( tError > 1.0e-6 )
                {
                    // get values
                    aT2 -= tOmega * tF( 0 );
                    aP2 -= tOmega * tF( 1 );
                    aU2 -= tOmega * tF( 2 );

                    // calculate inverse density
                    tV2 = this->v( aT2, aP2 );

                    // Jacobian Terms

                    // d(mass)/dT
                    tJ( 0, 0 ) = -aA2 * aU2 * this->alpha( aT2, aP2 ) / ( tV2 * tMass );

                    // d(entropy)/dT
                    tJ( 1, 0 ) = this->dsdT( aT2, aP2 ) / tEntropy;

                    // d(energy)/dT
                    tJ( 2, 0 ) = this->cp( aT2, aP2 ) / tEnergy;

                    // d(mass)/dP
                    tJ( 0, 1 ) = -this->kappa( aT2, aP2 ) * aA2 * aU2 / ( tV2 * tMass );

                    // d(entropy)/dP
                    tJ( 1, 1 ) = this->dsdp( aT2, aP2 ) / tEntropy;

                    //  d(energy)/dP
                    tJ( 2, 1 ) = mEoS->hdep( aT2, aP2 ) / tEnergy;

                    // d(mass)/du
                    tJ( 0, 2 ) = aA2 / ( tV2 * tMass );

                    // d(entropy)/du
                    tJ( 1, 2 ) = 0.0;

                    // d(energy)/du
                    tJ( 2, 2 ) = aU2 / tEnergy;


                    // RHS
                    tF( 0 ) = ( aU2 * aA2 / tV2 - tMass ) / tMass;
                    tF( 1 ) = ( this->s( aT2, aP2 ) - tEntropy ) / tEntropy;
                    tF( 2 ) = ( this->h( aT2, aP2 ) + 0.5 * aU2 * aU2 - tEnergy ) / tEnergy;

                    tError = norm( tF );

                    // solve system
                    gesv( tJ, tF, tPivot );

                    // increment counter
                    ++tCount;

                    BELFEM_ERROR( tCount < 100,
                                 "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                                 aT1,
                                 aP1,
                                 aU1,
                                 aA1,
                                 aA2 );

                }

                // check mach number
                BELFEM_ERROR( aU2 / this->c( aT2, aP2 ) > 1.0,
                             "Could not find supersonic solition for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                             aT1,
                             aP1,
                             aU1,
                             aA1,
                             aA2 );

            }
        }
    }

//------------------------------------------------------------------------------

    real
    Gas::prandtl_meyer_angle( const real & aT, const real & aP, const real & aU )
    {

        real k = this->gamma( aT, aP );
        real Ma = aU / this->c( aT, aP );

        return std::sqrt( ( k + 1. ) / ( k - 1. ) )
               * std::atan( std::sqrt( (  k - 1. ) / ( k + 1. ) *
               ( Ma * Ma - 1. ) ) )
               - std::atan( std::sqrt(  Ma * Ma - 1. ) );

    }

//------------------------------------------------------------------------------
    real
    Gas::prandtl_meyer(
                    const real & aT1,
                    const real & aP1,
                    const real & aU1,
                    const real & aAlpha,
                          real & aT2,
                          real & aP2,
                          real & aU2 )
    {
        // minumum allowable temperature
        const real Tmin = gastables::gTmin ;

        // minimal allowabpe pressure
        const real pmin = gastables::gPmin ;

        real R = this->R( aT1, aP1 );

        const real ht = this->h( aT1, aP1 ) + 0.5 * aU1 * aU1;
        const real s = this->s( aT1, aP1 );

        const real aNu1 = this->prandtl_meyer_angle( aT1, aP1, aU1 );
        const real aNu2 = aAlpha + aNu1;

        // solution vector
        Vector<real> F(3,0.0);

        // Jacobian
        Matrix<real> J( 3,3  );

        // Help Vector
        Vector< int > P( 3 );

        // initial guess
        aT2 = aT1;
        aP2 = aP1;
        aU2 = aU1;

        // abbreviations
        real & T = aT2;
        real & p = aP2;
        real & u = aU2;

        real Ma = aU1 / this->c( aT1, aP1 );

        real Ma_old = Ma;

        real L;
        real G;
        real k = this->gamma( T, p );
        real dkdT;
        real cp = this->cp( T, p );

        real tError = 1.0;
        uint tCount = 0;

        real f = 1.0;
        real df;

        // initial guess  using perfect gas assumption
        G = std::sqrt(( k + 1. ) / ( k - 1. ));

        while( std::abs(f) > 0.001 )
        {
            L = std::sqrt( Ma * Ma - 1 );
            f = G * std::atan( L / G ) - std::atan( L ) - aNu2;
            df = ( G * G / ( L * L + G * G ) - 1. / ( 1. + L * L ))* Ma / L;

            Ma -= f/df;

            BELFEM_ERROR( tCount < 100,
                         "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, alpha=%12.3f",
                         aT1,
                         aP1,
                         aU1,
                         aAlpha );

            ++tCount;
        }

        BELFEM_ASSERT( std::abs( Ma ) > 1.0 ,"Fail in Prandly Meyer Expansion" );

        // total temperature
        real Tt = aT1 * ( 1.0 + 0.5*(k-1)*Ma_old * Ma_old ) ;

        T = Tt / ( 1.0 + 0.5 * ( k-1 ) * Ma * Ma );
        p = aP1 * std::pow( aT2/aT1, k/(k-1) );

        // check for temperature limits
        if( ! ( T > Tmin + 25.0 ) )
        {
            T = Tmin + 25.0  ;
            p = std::pow( T/aT1, k/(k-1) );
            Ma = std::sqrt( ( Tt / T - 1.0 ) * 2.0 / ( k-1.0 ) );
        }

        // reset counter
        tCount = 0;
        const real tOmega0 = 0.6;

        real tOmega = tOmega0 ;
        real tOmega1 ;

        while( tError > 1.0e-9 )
        {
            dkdT = -R * this->dcpdT( T, p ) / std::pow( cp - R, 2 );

            L = std::sqrt( Ma * Ma - 1 );
            G = std::sqrt(( k + 1 ) / ( k - 1. ));


            J( 0, 0 ) = ( L * G / ( L * L + G * G ) - std::atan( L / G ))
                        * dkdT / ( G * ( k - 1. ) * ( k - 1. ));

            J( 1, 0 ) = cp + 0.5 *  Ma * Ma * R * ( k + T * dkdT );

            J( 2, 0 ) = this->dsdT( T, p );


            J( 0, 1 ) = 0.0;

            J( 1, 1 ) = 0.0;

            J( 2, 1 ) =  this->dsdp( T, p );

            J( 0, 2 ) = ( G * G / ( L * L + G * G ) - 1. / ( 1. + L * L ))
                        * Ma / L;

            J( 1, 2 ) = Ma * k * R * T;

            J( 2, 2 ) = 0.0;

            F( 0 ) = this->prandtl_meyer_angle( T, p, u ) - aNu2;
            F( 1 ) = this->h( T, p ) + 0.5 * u * u - ht;
            F( 2 ) = this->s( T, p ) - s;

            tError = std::sqrt(
                      F( 0 ) * F( 0 ) / ( aNu2* aNu2 )
                    + F( 1 ) * F( 1 ) / ( ht * ht )
                    + F( 2 ) * F( 2 ) / ( s * s ) );

            gesv( J, F, P );

            // limit relaxation factor against lower limits of T, p and ma
            tOmega1 = std::abs( 0.9 * ( T - Tmin ) / F( 0 ) );
            tOmega1 = std::min(std::abs( 0.9 * ( p - pmin ) / F( 1 ) ), tOmega1 );
            tOmega1 = std::min(std::abs( 0.9 * ( Ma - 1.0 ) / F( 2 ) ), tOmega1 );


            if( tOmega < 1e-6 )
            {
                break ;
            }

            tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;



            if( tOmega < 1e-6 )
            {
                break ;
            }

            T  -= tOmega * F( 0 );
            p  -= tOmega * F( 1 );
            Ma -= tOmega * F( 2 );

            cp = this->cp( T, p );
            k = this->gamma( T, p );
            u = Ma * this->c( T, p );

            // update value of R ( has only an effect for TableGas )
            R = this->R( T, p );

            ++tCount;

            BELFEM_ERROR( tCount < 100,
                         "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, alpha=%12.3f",
                         aT1,
                         aP1,
                         aU1,
                         aAlpha );
        }

        return aNu2;
    }

//------------------------------------------------------------------------------

    void
    Gas::shock(  const real & aT1, const real & aP1, const real & aU1,
                       real & aT2,       real & aP2,       real & aU2 )
    {
        // relaxation factor
        real tOmega0 = 0.9;
        real tOmega1 ;
        real tOmega ;

        // mach number before shock
        real tMa1 = aU1 / this->c( aT1, aP1 );

        if ( tMa1 < 1.001 )
        {
            aT2 = aT1;
            aP2 = aP1;
            aU2 = aU1;
        }
        else
        {
            // total enthalpy
            real tHt = this->h( aT1, aP1 ) + 0.5 * aU1 * aU1;

            real tGamma = this->gamma( aT1, aP1 );

            // state vector: rho2, u2, T2
            Vector< real > tX( 3 );
            real & tRho2 = tX( 0 );
            real & tU2 = tX( 1 );
            real & tT2 = tX( 2 );

            // total temperature
            real tTt = aT1 * ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa1 * tMa1 );

            // initial guess for Ma2
            real tMa2 = std::sqrt(
                    ( tTt / aT1 ) / ( tGamma * tMa1 * tMa1 - 0.5 * ( tGamma - 1.0 )));



            // initial guess for temperature
            tT2 = tTt / ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa2 * tMa2 );

            // check for temperature limit
            if( tT2 > gastables::gTmax - 1000.0 )
            {
                tT2 = gastables::gTmax - 1000.0 ;
            }

            // initial guess for pressure
            aP2 = std::pow( tT2 / aT1, tGamma / ( tGamma - 1.0 ));

            // initial guess for density
            tRho2 = this->rho( tT2, aP2 );

            // initial guess for velocity
            tU2 = tMa2 * std::sqrt( this->R( tT2, aP2 ) * tGamma * tT2 );

            // mass
            real tM = this->rho( aT1, aP1 ) * aU1;

            // momentum
            real tI = aP1 + tM * aU1;

            // Jacobian
            Matrix< real > tJ( 3, 3 );

            // right hand side
            Vector< real > tY( 3 );

            real tResiduum = BELFEM_REAL_MAX;

            // counter for the loop
            uint tCount = 0;

            // for solver
            Vector< int > tPivot( 3 );

            real tTmax = gastables::gTmax - 1.0 ;

            while ( tResiduum > 1.0e-9 )
            {
                // compute Jacobian
                tJ( 0, 0 ) = tU2;
                //tJ( 1, 0 ) = 1.0 / ( this->kappa( tT2, aP2 ) * tRho2 ) + tU2 * tU2 ;
                tJ( 1, 0 ) = ( this->R( tT2, aP2 ) * tT2 + tU2 * tU2 );

                tJ( 2, 0 ) = 0.0;

                tJ( 0, 1 ) = tRho2;
                tJ( 1, 1 ) = 2.0 * tRho2 * tU2;
                tJ( 2, 1 ) = tU2;

                tJ( 0, 2 ) = 0.0;

                //tJ( 1, 2 ) = this->beta( tT2, aP2 ) ;
                tJ( 1, 2 ) = this->R( tT2, aP2 ) * tRho2;
                tJ( 2, 2 ) = this->cp( tT2, aP2 );

                // compute right hand side
                aP2 = this->p( tT2, 1.0 / tRho2 );

                tY( 0 ) = ( tRho2 * tU2 - tM );
                tY( 1 ) = ( aP2 + tRho2 * tU2 * tU2 - tI );
                tY( 2 ) = ( this->h( tT2, aP2 ) + 0.5 * tU2 * tU2 - tHt );

                // compute residuum
                tResiduum = std::sqrt(
                        std::pow( tY( 0 ) / tM, 2 )
                        + std::pow( tY( 1 ) / tI, 2 )
                        + std::pow( tY( 2 ) / tHt, 2 ));

                // solve system
                gesv( tJ, tY, tPivot );

                // limit relaxation
                tOmega1 = std::abs( 0.9 * ( tX( 2 ) - tTmax ) / tY( 2 ) );
                tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;

                // detect failure
                if ( tOmega < 1e-6 )
                {
                    // limit temperature
                    aT2 = tTmax ;

                    this->total( aT1, aP1, aU1, aT2, aP2 ) ;

                    // pressure
                    aP2 = this->p( tX( 2 ), 1.0 / tX( 0 ));

                    aU2 = tX( 1 );

                    return ;
                }
                // perform Newton Step
                tX -= tOmega * tY;

                BELFEM_ERROR( tCount++ < 1000,
                             "Infinite loop in Gas::shock" );
            }

            // postprocess output values
            aT2 = tX( 2 );
            aP2 = this->p( tX( 2 ), 1.0 / tX( 0 ));
            aU2 = tX( 1 );
        }
 }

//------------------------------------------------------------------------------

    void
    Gas::shock( const real & aT1, const real & aP1, const real & aU1, const real & aAlpha,
           real & aT2, real & aP2, real & aU2, real & aBeta )
    {
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // Step 1: indentify min and max possible beta-angle
        // - - - - - - - - - - - - - - - - - - - - - - - - -

        // compute mach number
        real tMa1 = aU1 / this->c( aT1, aP1 );

        BELFEM_ERROR( tMa1 > 0.0,
            "Mach number must be > 0 ( is %f )",
                     ( float ) tMa1 );

        // mach angle
        real tBetaA = std::asin( 1.0 / tMa1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // Step 1: indentify min and max possible beta-angle
        // - - - - - - - - - - - - - - - - - - - - - - - - -

        // critical angle ( correlation for air )
        real tBetaB = ( ( ( -0.1189 / tMa1 + 1.1617 ) / tMa1 - 0.7708 ) / ( tMa1 * tMa1 ) + 1.1832 ) ;

        real tDeltaB = ( tBetaB - tBetaA ) / 20.0;

        real tFA = this->shock_beta( aT1, aP1, aU1, aAlpha, aT2, aP2, aU2, tBetaA );
        real tFB = tFA;

        tBetaB = tBetaA;

        while( tFB > 0.0 )
        {
            // shift F
            tFA = tFB;

            // shift Beta
            tBetaA = tBetaB;

            // increment beta
            tBetaB += tDeltaB;

            tFB = this->shock_beta_simple( aT1, aP1, aU1, aAlpha, aT2, aP2, aU2, tBetaB );
        }

        // one more step for safety
        tBetaB += tDeltaB;


        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // Step 2: initial iteration using bisection
        // - - - - - - - - - - - - - - - - - - - - - - - - -

        real tF;

        // loop counter
        uint tCount = 0;
        real tBeta0 = tBetaA;
        aBeta  = tBetaB;

        while ( std::abs( tBeta0 - aBeta ) > 1.0e-3 )
        {
            // shift beta
            tBeta0 = aBeta;

            // new beta
            aBeta = 0.5 * ( tBetaA + tBetaB );


            // call beta function
            tF = this->shock_beta( aT1, aP1, aU1, aAlpha, aT2, aP2, aU2, aBeta );

            // test result
            if( tFA * tF > 0 )
            {
                tBetaA = aBeta;
                tFA = tF;
            }
            else
            {
                tBetaB = aBeta;
            }

            // increment counter
            BELFEM_ERROR( tCount++ < 1000, "too many iterations" );
        }


    }


//------------------------------------------------------------------------------
    real
    Gas::v( const uint & aIndex, const real & aT, const real & aP )
    {
        return mEoS->v( aIndex, aT, aP );
    }
//------------------------------------------------------------------------------

    real
    Gas::h( const uint & aIndex, const real & aT, const real & aP )
    {
        return mComponents( aIndex )->h( aT )
               / mComponents( aIndex )->data()->M()
               + mEoS->hdep( aIndex, aT, aP );
    }

//------------------------------------------------------------------------------

    real
    Gas::cp( const uint & aIndex, const real & aT, const real & aP )
    {
        return mComponents( aIndex )->heat_spline()->deval( aT )
            / mComponents( aIndex )->data()->M()
            + mEoS->cpdep( aIndex, aT, aP );
    }

    real
    Gas::dcpdT( const uint & aIndex, const real & aT, const real & aP )
    {
        // ideal gas contribution
        real tdCpdT0 = mComponents( aIndex )->heat_spline()->ddeval( aT )
                    /  mComponents( aIndex )->data()->M() ;

        // departure Hi
        real tT1 = std::min( 1.001 * aT, gastables::gTmax );
        real tdCpdepdT1 = mEoS->cpdep( aIndex, tT1, aP ) ;

        // departure Low
        real tT0 = std::max( 0.999 * aT, gastables::gTmin );
        real tdCpdepdT0 = mEoS->cpdep( aIndex, tT0, aP ) ;

        return tdCpdT0 + ( tdCpdepdT1 - tdCpdepdT0 ) / ( tT1 - tT0 );
    }

//------------------------------------------------------------------------------

    void
    Gas::print()
    {
        std::fprintf( stdout, "             Molar    Mass\n" );
        for( uint k=0; k<mNumberOfComponents; ++k )
        {
            std::fprintf( stdout, "%3d %8s %6.4e %6.4e \n",
                    k,
                    mComponents( k )->label().c_str(),
                    mMolarFractions( k ),
                    mMassFractions( k ) );
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::link_to_idgas_property_functions()
    {
        mFunctionCp = &Gas::idgas_cp;
        mFunctiondCpdT = &Gas::idgas_dcpdT;
        mFunctionCv = &Gas::idgas_cv;
        mFunctionGamma = &Gas::idgas_gamma;
        mFunctionC = &Gas::idgas_c;
        mFunctionH = &Gas::idgas_h;
        mFunctionS = &Gas::idgas_s;
        mFunctionDHDP = &Gas::dhdp_idgas ;
        mFunctionDSDT = &Gas::idgas_dsdT;
        mFunctionDSDP = &Gas::idgas_dsdp;
        mFunctionMU = &Gas::idgas_mu;
        mFunctionLAMBDA = &Gas::idgas_lambda;
    }

//------------------------------------------------------------------------------

    void
    Gas::link_to_realgas_property_functions()
    {
        mFunctionCp = &Gas::realgas_cp;
        mFunctiondCpdT = &Gas::realgas_dcpdT;
        mFunctionH = &Gas::realgas_h;
        mFunctionCv = &Gas::realgas_cv;
        mFunctionGamma = &Gas::realgas_gamma;
        mFunctionC = &Gas::realgas_c;
        mFunctionS = &Gas::realgas_s;
        mFunctionDHDP = &Gas::dhdp_eos_departure ;
        mFunctionDSDT = &Gas::realgas_dsdT;
        mFunctionDSDP = &Gas::realgas_dsdp;
        mFunctionMU = &Gas::realgas_mu;
        mFunctionLAMBDA = &Gas::realgas_lambda;
    }

//------------------------------------------------------------------------------

    // use the property functions from the equation of state
    void
    Gas::link_to_helmholtz_property_functions()
    {
        mFunctionCp = &Gas::helmholtz_cp;
        mFunctiondCpdT = &Gas::helmholtz_dcpdT;
        mFunctionH = &Gas::helmholtz_h;
        mFunctionCv = &Gas::helmholtz_cv;
        mFunctionGamma = &Gas::helmholtz_gamma;
        mFunctionC = &Gas::helmholtz_c;
        mFunctionS = &Gas::helmholtz_s;
        mFunctionDHDP = &Gas::dhdp_differential_quotient ;
        mFunctionDSDT = &Gas::helmholtz_dsdT;
        mFunctionDSDP = &Gas::helmholtz_dsdp;
        mFunctionMU = &Gas::helmholtz_mu;
        mFunctionLAMBDA = &Gas::helmholtz_lambda;
    }

//------------------------------------------------------------------------------

    gastables::GasData *
    Gas::data( const index_t aIndex )
    {
        BELFEM_ASSERT( aIndex < mNumberOfComponents,
                      "Requested Index is out of bounds ( %lu vs %lu )",
                      ( long unsigned int ) aIndex,
                      ( long unsigned int ) mNumberOfComponents );

        return mComponents( aIndex )->data();
    }

//------------------------------------------------------------------------------

    void
    Gas::check_thermo_exists()
    {
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            BELFEM_ERROR( mComponents( k )->has_thermo(),
                         "No thermodynamic data for %s found in database.",
                         mComponents( k )->label().c_str());
        }
    }

//------------------------------------------------------------------------------

    // special subroutine needed for shock
    real
    Gas::shock_beta_simple( const real & aT1,
                const real & aP1,
                const real & aU1,
                const real & aAlpha,
                real & aT2,
                real & aP2,
                real & aU2,
                const real & aBeta )
    {
        real tV  = aU1 * std::cos( aBeta );
        real tU1 = aU1 * std::sin( aBeta );
        real tU2;
        real tK = 1.4;

        // compute oblique shock for perfect gas
        real tMa1 = aU1 / std::sqrt( tK * this->R( aT1, aP1 ) * aT1 );
        real tTt  = aT1 * ( 1.0 + 0.5 * ( tK - 1.0 ) * tMa1 * tMa1 );

        real tMa2 = std::sqrt( tTt / ( aT1 * ( tK * tMa1 * tMa1 - 0.5 * ( tK - 1.0 ) ) ) );

        aT2 = tTt / ( 1.0 + 0.5 * ( tK - 1.0 ) * tMa2 * tMa2 );
        aP2 = aP1 * std::pow( ( aT2 / aT1 ), tK / ( tK - 1.0 ) );
        tU2 = tMa2 * std::sqrt( tK * this->R( aT1, aP1 ) * aT2 );

        // compute the velocity
        aU2 = std::sqrt( tU2 * tU2 + tV * tV );

        // return the function
        return tU2 / tU1 - std::tan( aBeta - aAlpha ) / std::tan( aBeta );
    }

//------------------------------------------------------------------------------

    // special subroutine needed for shock
    real
    Gas::shock_beta( const real & aT1,
                     const real & aP1,
                     const real & aU1,
                     const real & aAlpha,
                     real & aT2,
                     real & aP2,
                     real & aU2,
                     const real & aBeta )
    {
        real tV  = aU1 * std::cos( aBeta );
        real tU1 = aU1 * std::sin( aBeta );
        real tU2;

        // compute oblique shock
        if( tU1 > this->c( aT1, aP1 ) )
        {
            this->shock( aT1, aP1, tU1, aT2, aP2, tU2 );
        }
        else
        {
            aT2 = aT1;
            aP2 = aP1;
            tU2 = tU1;
        }

        // compute the velocity
        aU2 = std::sqrt( tU2 * tU2 + tV * tV );

        // return the function
        return tU2 / tU1 - std::tan( aBeta - aAlpha ) / std::tan( aBeta );
    }

//------------------------------------------------------------------------------

    // enthalpy derivative to pressure ( needed for total temperature )
    real
    Gas::dhdp( const real & aT, const real & aP )
    {
        mStatevals.update_Tp( aT, aP );

        if ( !mStatevals.test( BELFEM_STATEVAL_DHDP ) )
        {
            mStatevals.set( BELFEM_STATEVAL_DHDP,
                            ( this->*mFunctionDHDP )( aT, aP ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DHDP );
    }

//------------------------------------------------------------------------------

    void
    Gas::update_mixture_entropy()
    {
        mMixtureEntropy = 0.0 ;

        for( real tX : mMolarFractions )
        {
            if ( tX > 1e-9 )
            {
                mMixtureEntropy -= tX * std::log(tX );
            }
        }
    }

//------------------------------------------------------------------------------

}