/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

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
#include "cl_GM_EoS_Nitrogen.hpp"
#include "cl_GM_HelmholtzTransport.hpp"
#include "cl_GM_HelmholtzTransport_Methane.hpp"
#include "cl_GM_HelmholtzTransport_LemmonJacobsen.hpp"
#include "cl_GM_HelmholtzTransport_Hydrogen.hpp"
#include "cl_Timer.hpp"
#include "intpoints.hpp"

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
            else if ( aLabel == "N2" )
            {
                mHelmholzModel = HelmholtzModel::Nitrogen ;
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
            case( HelmholtzModel::Nitrogen ) :
            {
                tLabel = "N2" ;
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

            // create the help matrix for remixing
            tFactory.create_helpmatrix( mHelpMatrix );

            // reserve memory for heat spline. remix_heat() fills it through
            // matrix_data(), so the mode has to be declared here: the extra row
            // is mixed from the components, which all carry entropy tables.
            mHeatSpline.matrix_data().set_size( 5, gastables::gNumberOfSplinePoints );
            mHeatSpline.set_extra_mode( spline::ExtraMode::Entropy );

            // allocate memory for working vectors
            mWorkMu.set_size( gastables::gNumberOfSplinePoints );
            mWorkLambda.set_size( gastables::gNumberOfSplinePoints );

            // scratch for the duct solvers
            mFlowJacobian.set_size( 3, 3, BELFEM_QUIET_NAN );
            mFlowResidual.set_size( 3, BELFEM_QUIET_NAN );
            mFlowPivot.set_size( 3, 0 );

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
    Gas::element_to_molecule( const string & aElement ) const
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
    Gas::reference_element( const string & aElement ) const
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
            /* liquid, not crystal: bromine is liquid at the 298.15 K standard
             * state, and CEA anchors the assigned enthalpies of Br species to
             * Br2(L). The Br2(cr) record also ends below 298.15 K, so the
             * crystal would anchor the formation table on an extrapolation. */
            return string( "Br2(L)" );
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
        // allocate matrix; the sentinel must be larger than any possible pair
        // index ( up to N*(N-1)/2, which can exceed N )
        mViscosityInteractionTable.set_size(
                mNumberOfComponents,
                mNumberOfComponents,
                BELFEM_UINT_MAX );


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

        // the composition changes, so every cached state value is stale
        mStatevals.reset();

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

    real
    Gas::R( const real T, const real p ) const
    {
        return mR;
    }


//------------------------------------------------------------------------------

    real
    Gas::M( const real T, const real p ) const
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


        // the composition changes, so every cached state value is stale
        mStatevals.reset();

        // copy molar fractions
        mMolarFractions = aMolarFractions;

        // make molar fractions partition of unity
        mMolarFractions /= sum( mMolarFractions );

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

                        /* Muzny viscosity and Assael thermal conductivity.
                         * the viscosity correlation is for normal hydrogen and
                         * serves all three spin isomers, the conductivity
                         * tables are isomer aware, see the class header */
                        mTransport = new gasmodels::HelmholtzTransport_Hydrogen(
                                *this, mHelmholzModel );


                        break ;
                    }
                    case( HelmholtzModel::Oxygen ) :
                    {
                        mEoS = new gasmodels::EoS_Oxygen( *this );

                        mTransport = new gasmodels::HelmholtzTransport_LemmonJacobsen(
                                *this, HelmholtzModel::Oxygen );

                        break ;
                    }
                    case( HelmholtzModel::Methane ) :
                    {
                        mEoS = new gasmodels::EoS_Methane( *this );
                        mTransport = new gasmodels::HelmholtzTransport_Methane( *this );

                        break ;
                    }
                    case( HelmholtzModel::Nitrogen ) :
                    {
                        mEoS = new gasmodels::EoS_Nitrogen( *this );

                        mTransport = new gasmodels::HelmholtzTransport_LemmonJacobsen(
                                *this, HelmholtzModel::Nitrogen );

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

        // scale unit to J/kg
        tData /= mM;

        // additional term for idgas_s
        this->update_mixture_entropy() ;
    }
//------------------------------------------------------------------------------

    void
    Gas::remix_critical_point()
    {
        mEoS->eval_critical_point( mTcrit, mPcrit, mVcrit );

        // stiel thodos parameter
        mGamma = std::pow( mTcrit, 1.0 / 6.0 )
                 * std::pow( mM * 1000, 0.5 )               // scale: kg->g
                 * std::pow( mPcrit / 1.01325e5, -2.0 / 3.0 ) // scale: Pa->atm
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
    Gas::idgas_cp( const real T, const real p ) const
    {
        return mHeatSpline.deval( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dcpdT( const real T, const real p ) const
    {
        return mHeatSpline.ddeval( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_cv( const real T, const real p ) const
    {
        return this->cp( T, p ) - this->R( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_gamma( const real T, const real p ) const
    {
        real tCp = this->cp( T, p );
        return tCp / ( tCp - this->R( T, p ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_c( const real T, const real p ) const
    {
        return std::sqrt( this->idgas_gamma( T, p ) * this->R( T, p ) * T );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_h( const real T, const real p ) const
    {
        return mHeatSpline.eval( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_s( const real T, const real p ) const
    {
        return ( std::log( gastables::gPref / p ) + mMixtureEntropy )
            * this->R( T, p )
        + mHeatSpline.entropy( T, this->spline_col( T ) ) ;
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dsdT( const real T, const real p ) const
    {
        return mHeatSpline.dentropy( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_dsdp( const real T, const real p ) const
    {
        return -this->R( T, p ) / p;
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_mu( const real T, const real p ) const
    {
        return mViscositySpline.eval( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::idgas_lambda( const real T, const real p ) const
    {
        return mConductivitySpline.eval( T, this->spline_col( T ) );
    }

//------------------------------------------------------------------------------

    /**
     * REFERENCE PRESSURE CONVENTION OF THE REAL GAS FUNCTIONS
     *
     * The heat spline holds the CEA standard state: the IDEAL gas at
     * gPref = 1 bar. The ideal gas functions evaluate the spline directly and do
     * not use departure functions, so idgas_cp is the spline derivative. The real
     * gas functions below add the departure at ( T, p ) and subtract the departure
     * at ( T, gPref ).
     *
     * The subtraction enforces continuity between models. At 1 bar, the two
     * departure terms cancel, the result is the spline value, and SRK and PR reduce
     * exactly to the ideal gas model at the reference pressure. The departure
     * splines are rebuilt whenever the gas model changes, so the convention holds
     * for each cubic model. Helmholtz does not use this path because it takes
     * caloric properties from the equation of state.
     *
     * The subtraction also defines what the spline represents. If the spline held
     * the REAL gas value at 1 bar, this assembly would be exact:
     *
     *     h( T, p ) = h_real( T, gPref ) + hdep( T, p ) - hdep( T, gPref )
     *
     * The spline instead holds the ideal gas value at 1 bar. The result therefore
     * differs from the real gas enthalpy by the departure at 1 bar, which varies
     * with temperature. Dropping the subtraction gives the other self consistent
     * choice:
     *
     *     h( T, p ) = h_ideal( T ) + hdep( T, p )
     *
     * That expression gives the true real gas property, but it introduces a step
     * when switching from the ideal gas model to a cubic model at 1 bar. The step is
     * the real departure from ideality, not an artifact. This path cannot provide
     * both model continuity and the exact real gas property relative to the ideal
     * spline.
     */
    real
    Gas::realgas_cp( const real T, const real p ) const
    {
        return mHeatSpline.deval( T, this->spline_col( T ) )
               + mEoS->cpdep( T, p )
               - mEoS->cpdep0( T );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dcpdT( const real T, const real p ) const
    {
        real tCPDEP2 = mEoS->cpdep( T+1.0, p )  - mEoS->cpdep0( T+1.0 );
        real tCPDEP1 = mEoS->cpdep( T-1.0, p )  - mEoS->cpdep0( T-1.0 );

        return mHeatSpline.ddeval( T, this->spline_col( T ) ) + 0.5 * ( tCPDEP2 - tCPDEP1 );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_cv( const real T, const real p ) const
    {
        // ( 2.11 )
        return this->cp( T, p ) - p * this->v( T, p ) * T
                                    * this->alpha( T, p ) * this->beta( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_gamma( const real T, const real p ) const
    {
        // ( 2.20 )
        return this->cp( T, p ) * this->beta( T, p ) /
               ( this->cv( T, p ) * this->alpha( T, p ));
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_c( const real T, const real p ) const
    {
        // c^2 = gamma_pv * p * v with gamma_pv from ( 2.20 ):
        // the cp factor was missing here ( result had units of Kelvin )
        return std::sqrt( p * this->v( T, p ) * this->cp( T, p )
                          * this->beta( T, p ) /
                          ( this->alpha( T, p ) * this->cv( T, p )));
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_h( const real T, const real p ) const
    {
        return mHeatSpline.eval( T, this->spline_col( T ) )
               + mEoS->hdep( T, p )
               - mEoS->hdep0( T );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_s( const real T, const real p ) const
    {
        // Entropy includes an ideal gas pressure term; enthalpy and heat capacity do not.
        // The ideal gas change is:
        //
        //     s2^0 - s1^0 = int cp^0/T dT - R ln( p2/p1 )
        //
        // Therefore mR*log( gPref/p ) belongs to the IDEAL GAS contribution,
        // referenced to gPref = 1 bar. It is not a departure term and is not sdep.
        // See the reference pressure note above Gas::realgas_cp.
        // The pressure and mixture terms are identical to the ones in idgas_s,
        // so the cubic models reduce to the ideal gas at gPref.
        return ( std::log( gastables::gPref / p ) + mMixtureEntropy )
               * this->R( T, p )
               + mHeatSpline.entropy( T, this->spline_col( T ) )
               + mEoS->sdep( T, p )
               - mEoS->sdep0( T );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dsdT( const real T, const real p ) const
    {
        // This is the temperature derivative of realgas_s. The ideal gas pressure
        // term in entropy is constant in T, so it is absent here.
        return mHeatSpline.dentropy( T, this->spline_col( T ) )
               + mEoS->dsdepdT( T, p )
               - mEoS->dsdepdT0( T );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_dsdp( const real T, const real p ) const
    {
        return -mR / p + mEoS->dsdepdp( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_mu( const real T, const real p ) const
    {
        real mu = mViscositySpline.eval( T, this->spline_col( T ) );

        mu += this->mu_dep( mu, T, p );

        BELFEM_ASSERT( mu > 0.0 && mu < 1.0,
                      "Error in mu" );

        return mu;
    }

//------------------------------------------------------------------------------

    real
    Gas::realgas_lambda( const real T, const real p ) const
    {

        real lambda = mConductivitySpline.eval( T, this->spline_col( T ) )
                       + this->lambda_dep( T, p );

        BELFEM_ASSERT( lambda > 0.0 && lambda < 1.0,
                      "Error in lambda" );

        return lambda;
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_cp( const real T, const real p ) const
    {
        return mEoS->cp( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dcpdT( const real T, const real p ) const
    {
        real tT1 = T * 0.9999 ;
        real tT2 = T * 1.0001 ;

        return ( mEoS->cp( tT2, p ) - mEoS->cp( tT1, p ) ) / ( tT2 - tT1 );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_cv( const real T, const real p ) const
    {
        return mEoS->cv( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_gamma( const real T, const real p ) const
    {
        return this->cp( T, p ) * this->beta( T, p ) /
           ( this->cv( T, p ) * this->alpha( T, p ) );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_c( const real T, const real p ) const
    {
        return mEoS->w( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_h( const real T, const real p ) const
    {
        return mEoS->h( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_s( const real T, const real p ) const
    {
        return mEoS->s( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dsdT( const real T, const real p ) const
    {
        return mEoS->dsdT( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_dsdp( const real T, const real p ) const
    {
        return mEoS->dsdp( T, p ) ;
    }

//------------------------------------------------------------------------------

    real
    Gas::helmholtz_mu( const real T, const real p ) const
    {
        return mTransport->mu( T, p );
    }
//------------------------------------------------------------------------------

    real
    Gas::helmholtz_lambda( const real T, const real p ) const
    {
        return mTransport->lambda( T, p );
    }

//------------------------------------------------------------------------------
    void
    Gas::evaluate_viscosity_interaction( const real T ) const
    {
        // - - - - - - - - - - - - - - - - - - - -
        // Step 1: Update Viscosity of components
        // - - - - - - - - - - - - - - - - - - - -
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            const gastables::RefGas * tComponent = mComponents( k );

            mWorkVector( k ) = tComponent->mu( T );
        }

        uint k;

        // - - - - - - - - - - - - - - - - - - - - - - -
        // Step 2: Calculate interaction parameter Phi
        // - - - - - - - - - - - - - - - - - - - - - - -
        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            // read-only handle: this pass only evaluates the components
            const gastables::RefGas * tComponent_i = mComponents( i );

            // evaluate Interaction
            if ( tComponent_i->has_viscosity() )
            {
                const double & tM_i = tComponent_i->data()->M();
                const double & tMu_i = mWorkVector( i );

                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    if ( i != j )
                    {
                        const gastables::RefGas * tComponent_j = mComponents( j );

                        if ( tComponent_j->has_viscosity() )
                        {
                            const double & tM_j = tComponent_j->data()->M();
                            const double & tMu_j = mWorkVector( j );

                            k = mViscosityInteractionTable( i, j );
                            // test if interaction parameter exists
                            if ( k != BELFEM_UINT_MAX )
                            {
                                // NASA RP-1311 ( 5.7 )
                                mWorkMatrix( i, j ) = tMu_i /
                                                      mViscosityInteractionRefgas( k )->mu( T )
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
        mWorkTemperature = T;
    }

//------------------------------------------------------------------------------

    void
    Gas::evaluate_conductivity_interaction( const real T ) const
    {
        if ( T != mWorkTemperature )
        {
            this->evaluate_viscosity_interaction( T );
        }

        // - - - - - - - - - - - - - - - - - - - -
        // Step 1: Updatem Conductiviies
        // - - - - - - - - - - - - - - - - - - - -
        for ( uint k = 0; k < mNumberOfComponents; ++k )
        {
            const gastables::RefGas * tComponent = mComponents( k );

            mWorkVector2( k ) = tComponent->lambda( T );
        }

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            // read-only handle: this pass only evaluates the components
            const gastables::RefGas * tComponent_i = mComponents( i );

            if ( tComponent_i->has_conductivity() && tComponent_i->has_viscosity())
            {
                //const double & tMu_i = mWorkVector( i );
                const double & tM_i = tComponent_i->data()->M();

                for ( uint j = 0; j < mNumberOfComponents; ++j )
                {
                    if ( i != j )
                    {
                        const gastables::RefGas * tComponent_j = mComponents( j );

                        if ( tComponent_j->has_conductivity() && tComponent_j->has_viscosity())
                        {
                            //const double & tMu_j = mWorkVector( j );
                            const double & tM_j = tComponent_j->data()->M();
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
    Gas::cea_mu( const real T ) const
    {

        this->evaluate_viscosity_interaction( T );

        const Vector<real> & tX = mMolarFractions;

        // calculate value
        real mu = 0.0;

        // temporary vector for mixing
        mWorkVector2 = mWorkMatrix * tX;

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            if ( tX( i ) > BELFEM_EPSILON_X && mComponents( i )->has_viscosity() )
            {

                mu += tX( i ) * mWorkVector( i ) /
                       ( tX( i ) + mWorkVector2( i ));

            }

        }

        return mu;
    }

//------------------------------------------------------------------------------

    real
    Gas::cea_lambda( const real T ) const
    {

        this->evaluate_conductivity_interaction( T );

        const Vector<real> & tX = mMolarFractions;

        // calculate value
        real lambda = 0.0;

        // temporary vector for mixing
        mWorkVector = mWorkMatrix * tX;

        for ( uint i = 0; i < mNumberOfComponents; ++i )
        {
            if ( tX( i ) > BELFEM_EPSILON_X && mComponents( i )->has_conductivity())
            {
                lambda += tX( i ) * mWorkVector2( i ) /
                           ( tX( i ) + mWorkVector( i ));
            }

        }

        return lambda;
    }

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Thermodynamic States
//------------------------------------------------------------------------------

    real
    Gas::p( const real T, const real v ) const
    {
        mStatevals.update_Tv( T, v );

        // check if values are up to date
        if ( !mStatevals.test( BELFEM_STATEVAL_P ))
        {
            mStatevals.set(
                    BELFEM_STATEVAL_P,
                    mEoS->p( T, v ));
        }

        return mStatevals.get( BELFEM_STATEVAL_P );
    }

//------------------------------------------------------------------------------

    real
    Gas::v( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_V ))
        {
            mStatevals.set( BELFEM_STATEVAL_V,
                            mEoS->v( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_V );
    }

//------------------------------------------------------------------------------

    real
    Gas::rho( const real T, const real p ) const
    {
        return 1.0 / this->v( T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::T( const real p, const real v ) const
    {
        mStatevals.update_pv( p, v );

        if ( !mStatevals.test( BELFEM_STATEVAL_T ))
        {
            mStatevals.set( BELFEM_STATEVAL_T,
                            mEoS->T( p, v ));
        }

        return mStatevals.get( BELFEM_STATEVAL_T );
    }

//------------------------------------------------------------------------------

    real
    Gas::cp( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_CP ))
        {
            mStatevals.set( BELFEM_STATEVAL_CP,
                            ( this->*mFunctionCp )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_CP );
    }

//------------------------------------------------------------------------------

    // dissociation enthalpy ( only for tablegas at this time )
    real
    Gas::hd( const real T, const real p ) const
    {
        return 0.0 ;
    }

//------------------------------------------------------------------------------

    real
    Gas::dcpdT( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_DCPDT ))
        {
            mStatevals.set( BELFEM_STATEVAL_DCPDT,
                            ( this->*mFunctiondCpdT )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DCPDT );
    }

//------------------------------------------------------------------------------

    real
    Gas::cv( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_CV ))
        {
            mStatevals.set( BELFEM_STATEVAL_CV,
                            ( this->*mFunctionCv )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_CV );
    }

//------------------------------------------------------------------------------

    real
    Gas::gamma( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_GAMMA ))
        {
            mStatevals.set( BELFEM_STATEVAL_GAMMA,
                            ( this->*mFunctionGamma )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_GAMMA );
    }

//------------------------------------------------------------------------------

    real
    Gas::c( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_C ))
        {
            mStatevals.set( BELFEM_STATEVAL_C,
                            ( this->*mFunctionC )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_C );
    }

//------------------------------------------------------------------------------

    real
    Gas::u( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_U ))
        {
            mStatevals.set( BELFEM_STATEVAL_U,
                            ( this->*mFunctionH )( T, p )
                            - p * this->v( T, p ) );
        }

        return mStatevals.get( BELFEM_STATEVAL_U );
    }

//------------------------------------------------------------------------------

    real
    Gas::h( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_H ))
        {
            mStatevals.set( BELFEM_STATEVAL_H,
                            ( this->*mFunctionH )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_H );
    }

//------------------------------------------------------------------------------

    real
    Gas::s( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_S ))
        {
            mStatevals.set( BELFEM_STATEVAL_S,
                            ( this->*mFunctionS )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_S );
    }

//------------------------------------------------------------------------------

    real
    Gas::dsdT( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_DSDT ))
        {
            mStatevals.set( BELFEM_STATEVAL_DSDT,
                            ( this->*mFunctionDSDT )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DSDT );
    }

//------------------------------------------------------------------------------

    real
    Gas::dsdp( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_DSDP ))
        {
            mStatevals.set( BELFEM_STATEVAL_DSDP,
                            ( this->*mFunctionDSDP )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_DSDP );
    }

//------------------------------------------------------------------------------

    real
    Gas::mu( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_MU ))
        {
            mStatevals.set( BELFEM_STATEVAL_MU,
                            ( this->*mFunctionMU )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_MU );
    }

//------------------------------------------------------------------------------

    real
    Gas::lambda( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_LAMBDA ))
        {
            mStatevals.set( BELFEM_STATEVAL_LAMBDA,
                            ( this->*mFunctionLAMBDA )( T, p ));
        }

        return mStatevals.get( BELFEM_STATEVAL_LAMBDA );
    }

//------------------------------------------------------------------------------

    real
    Gas::Pr( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_PR ))
        {
            mStatevals.set( BELFEM_STATEVAL_PR,
                            this->cp( T, p ) * this->mu( T, p ) /
                            this->lambda( T, p ));
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
    Gas::Hf( const real T, Vector<real> & aHf ) const
    {
        // reset vector
        aHf.fill( 0.0 );

        // add species enthalpies to vector
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aHf( tCount++ ) +=  tSpecie->H( T ) -  tSpecie->H_ref() + tSpecie->data()->Hf();
        }

        // calculate element enthalpies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->H( T ) - tElement->H_ref();
        }

        aHf.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::Gibbs( const real T, Vector<real> & aGibbs ) const
    {
        // write species entropies into vector
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) = -tSpecie->S( T );
        }

        // calculate element entropies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->S( T );
        }

        // add entropies to vector
        aGibbs.vector_data() += mFormationTable * mFormationWork;

        // multiply vector with temperatures
        aGibbs *= T;

        // add species enthalpies to vector
        tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) +=  tSpecie->H( T ) -  tSpecie->H_ref() + tSpecie->data()->Hf();
        }

        // calculate element enthalpies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->H( T ) - tElement->H_ref();
        }

        aGibbs.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::dGibbsdT( const real T, Vector< real > & aGibbs ) const
    {

        // calculate entropy term for species
        uint tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) = -( tSpecie->S( T ) + T * tSpecie->dSdT( T ) );
        }

        // calculate element entropies
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->S( T ) + T * tElement->dSdT( T );
        }

        aGibbs.vector_data() += mFormationTable * mFormationWork;

        // add specific heats to vector
        tCount = 0;
        for ( gastables::RefGas * tSpecie: mComponents )
        {
            aGibbs( tCount++ ) += tSpecie->Cp( T );
        }

        // calculate specific heats
        tCount = 0;
        for ( gastables::RefGas * tElement: mElements )
        {
            mFormationWork( tCount++ ) = tElement->Cp( T );
        }

        aGibbs.vector_data() -= mFormationTable * mFormationWork;
    }

//------------------------------------------------------------------------------

    void
    Gas::compute_equilibrium( const real T, const real p, Vector< real > & aX )
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
            const uint tNumElements = mElements.size();

            if ( mPivotRAND.length() != tNumElements + 1 )
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
            this->Gibbs( T, tMu0 );

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
                    tMu( k ) = tMu0( k ) + constant::Rm * T * std::log( p / gastables::gPref * aX( k ));
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
                        tRHS( i ) += tA( k, i ) * aX( k ) * tMu( k ) / ( constant::Rm * T );
                    }
                }
                tRHS( tNumElements ) = dot( aX, tMu ) / ( constant::Rm * T );

                // solve system
                gesv( mWorkMatrixRAND, mWorkVectorRAND0, mPivotRAND );

                // compute change of Mols
                for ( uint k = 0; k < mNumberOfComponents; ++k )
                {
                    tDeltaX( k ) = tU - tMu( k ) / ( constant::Rm * T );
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
    Gas::remix_to_equilibrium( const real T, const real p,
                               bool aRemixHeat,
                               bool aRemixTransport )
    {
        if( mNumberOfComponents > 1 )
        {
            this->compute_equilibrium( T, p, mMolarFractions );
            this->remix( mMolarFractions, aRemixHeat, aRemixTransport );
        }
    }

//------------------------------------------------------------------------------

    real
    Gas::mu_dep( const real & mu, const real T, const real p ) const
    {
        // cutoff value
        real tMuMax = 10.0 * mu;

        // reduced temperature
        real tTr = T / mTcrit;
        real tPr = p / mPcrit;

        // help magnitude ( see Eq. 88 )
        real tX = ( 0.807 * std::pow( tTr, 0.618 )
                    - 0.357 * std::exp( -0.449 * tTr )
                    + 0.34 * std::exp( -4.058 * tTr )
                    + 0.018 ) * 1e-7;

        // correction factor
        real tFid = mu * mXi / tX;

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
                                                           ( mXi * mu ), -3 );
                return std::min( 1.0e-7 * tZ2 * tFp / mXi - mu, tMuMax );
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

            return std::min( mu * tZ2 * tFp - mu, tMuMax );

        }
    }

//----------------------------------------------------------------------------

    real
    Gas::lambda_dep( const real T, const real p ) const
    {
        real tX = mVcrit / mEoS->v( T, p );

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
    Gas::T_from_h( const real & h, const real p ) const
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
        real tHref = this->h( tTref, p );

        real tCp = this->cp( tTref, p );

        // initial guess
        real T = ( h - tHref ) / tCp + tTref;

        uint tCount = 0;

        real tT = 1;

        real tOmega = 0.95;

        T = std::max( std::min( T, 6000.0 ), 200.0 );

        while( std::abs( tT - T ) > BELFEM_EPSILON_T  )
        {
            // shift t
            tT = T;

            // newton step
            T -= tOmega * ( this->h( T, p ) - h ) / this->cp( T, p );

            // increment counter
            if( tCount++ == 100 )
            {
                // fallback: bisection over the full table range
                // ( h is monotone in T since cp > 0 )
                real tT0 = gastables::gDeltaT;
                real tT1 = gastables::gTmax;

                real tF0 = this->h( tT0, p ) - h ;

                // without a sign change the bisection would silently
                // converge to an endpoint
                BELFEM_ERROR( tF0 * ( this->h( tT1, p ) - h ) <= 0.0,
                        "T_from_h: enthalpy h=%12.3f at p=%12.3f is outside the table range",
                        h, p );

                while( std::abs( tT0 - tT1 ) > BELFEM_EPSILON_T )
                {
                    T = 0.5 * ( tT0 + tT1 );

                    real tF = this->h( T, p ) - h;

                    if( tF0 * tF >= 0 )
                    {
                        tT0 = T;
                        tF0 = tF;
                    }
                    else
                    {
                        tT1 = T;
                    }

                    BELFEM_ERROR( tCount++ < 1000,
                                 "T_from_h did not converge for h=%12.3f, p=%12.3f",
                                 h, p );

                }

                break;
            }


           BELFEM_ERROR( tCount < 1000,
                    "T_from_h did not converge for h=%12.3f, p=%12.3f",
                    h, p );
        }

        return T;
    }

//------------------------------------------------------------------------------

    real
    Gas::isen_T( const real T0, const real p0, const real p1 ) const
    {
        // guess value for new temperature
        real T1 = T0 * std::pow( p1 / p0, this->R( T0, p0 ) / this->cp( T0, p0 ));

        // entropy at this state
        real tS = this->s( T0, p0 );

        real tT1 = 0.0;

        uint tCount = 0;

        real tDeltaT ;

        while ( std::abs( tT1 - T1 ) > BELFEM_EPSILON_T )
        {
            tT1 = T1;
            tDeltaT = ( this->s( T1, p1 ) - tS ) / this->dsdT( T1, p1 ) ;

            if( T1 - tDeltaT < 15.0 )
            {
                tDeltaT = T1 - 15.0 ;
            }
            T1 -= 0.9 * tDeltaT ;

            ++tCount;


            BELFEM_ERROR( tCount < 1000,
                         "Too many iterations for isen_T ( T0=%f, p0=%f, p1=%f)",
                         ( float ) T0, ( float ) p0, ( float ) p1 );
        }

        return T1;
    }


// -----------------------------------------------------------------------------

    real
    Gas::isen_p( const real T0, const real p0, const real T1 ) const
    {
        // guess value for new pressure
        real p1 = p0 * std::pow( T1 / T0, this->idgas_cp( T0, p0 ) / this->R( T0, p0 ) );

        // entropy at this state
        real tS = this->s( T0, p0 );

        real tP1 = 0.0;

        uint tCount = 0;

        while ( std::abs( tP1 - p1 ) > BELFEM_EPSILON_P )
        {
            tP1 = p1;
            p1 -= ( this->s( T1, p1 ) - tS ) / this->dsdp( T1, p1 );

            ++tCount;

            BELFEM_ERROR( tCount < 1000,
                         "Too many iterations for isen_p ( T0=%f, p0=%f, T1=%f)",
                         ( float ) T0, ( float ) p0, ( float ) T1 );
        }

        return p1;
    }

// -----------------------------------------------------------------------------

    void
    Gas::total( const real T, const real p, const real & u,
                real & aTt, real & aPt ) const
    {
        // maximum temperature
        real tTmax = gastables::gTmax - BELFEM_EPSILON_T ;

        // Initial guesses
        real tCp = this->idgas_cp( T, p );

        // gas constant
        real tR = this->R( T, p );

        // guess ratio of specific heats
        real tGamma = tCp / ( tCp - tR );

        // guess speed of sound
        real tC = std::sqrt( tGamma * tR * T );

        // guess mach number
        real tMa = u / tC;

        // solution vector
        Vector<real> tX( 2 );

        real & tTt = tX( 0 );
        real & tPt = tX( 1 );

        // RHS
        Vector<real> tF( 2, 1.0 );

        // Jacobian
        Matrix<real> tJ( 2, 2 );

        // Pivot
        Vector< int_t > tPivot( 2 );

        // guess total temperature
        tTt = std::min( T * ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa * tMa ), tTmax -100.0 );

        // guess total preassure
        tPt = p * std::pow( tTt / T, tCp / tR );

        // calculate entropy
        real tS = this->s( T, p );

        // calculate enthalpy
        real tH = this->h( T, p ) + 0.5 * u * u;

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

    real
    Gas::area_mach( const real Ma, const real k ) const
    {
        return std::pow( ( 1.0 + 0.5 * ( k - 1.0 ) * Ma * Ma )
                         * 2.0 / ( k + 1.0 ),
                         0.5 * ( k + 1.0 ) / ( k - 1.0 ) ) / Ma ;
    }

//------------------------------------------------------------------------------

    void
    Gas::area_guess(
            const real   T1,
            const real   p1,
            const real   Ma1,
            const real & A1,
            const real & A2,
            const bool   aSupersonic,
                  real & T2,
                  real & p2,
                  real & u2 ) const
    {
        // the guess is a frozen gamma ideal gas, the Newton that follows does
        // the real gas work
        real k = this->gamma( T1, p1 );

        // outlet area, referenced to the sonic area of the inlet
        real tTarget = this->area_mach( Ma1, k ) * A2 / A1 ;

        // bracket the root. A/A* is monotonic on either side of Ma = 1, falling
        // below and rising above, so the branch fixes which end moves
        real tLow  = aSupersonic ?  1.0 : 1e-6 ;
        real tHigh = aSupersonic ? 50.0 :  1.0 ;

        // the open end of the bracket has to reach past the target, otherwise
        // the bisection pins itself to it and returns a guess that is not one
        BELFEM_ERROR( this->area_mach( aSupersonic ? tHigh : tLow, k ) >= tTarget,
                     "Area ratio A2/A1=%12.3f is outside the bracket of the area relation",
                     A2 / A1 );

        real Ma2 = 0.5 * ( tLow + tHigh );

        uint tCount = 0 ;

        while ( tHigh - tLow > BELFEM_EPSILON * tHigh )
        {
            Ma2 = 0.5 * ( tLow + tHigh );

            if ( ( this->area_mach( Ma2, k ) > tTarget ) == aSupersonic )
            {
                tHigh = Ma2 ;
            }
            else
            {
                tLow = Ma2 ;
            }

            ++tCount ;

            BELFEM_ERROR( tCount < 100,
                         "Bisection of the area relation did not close for A2/A1=%12.3f",
                         A2 / A1 );
        }

        // total temperature of the frozen gamma gas
        real tTt = T1 * ( 1.0 + 0.5 * ( k - 1.0 ) * Ma1 * Ma1 );

        T2 = tTt / ( 1.0 + 0.5 * ( k - 1.0 ) * Ma2 * Ma2 );
        p2 = p1 * std::pow( T2 / T1, k / ( k - 1.0 ) );
        u2 = Ma2 * std::sqrt( k * this->R( T1, p1 ) * T2 );
    }

//------------------------------------------------------------------------------

    void
    Gas::isentropic_duct(
            const real & aMass,
            const real & aEntropy,
            const real & aEnergy,
            const real & A2,
                  real & T2,
                  real & p2,
                  real & u2 ) const
    {
        const real tOmega = 0.9 ;

        Matrix< real > & tJ = mFlowJacobian ;
        Vector< real > & tF = mFlowResidual ;

        uint tCount = 0 ;

        while ( true )
        {
            // inverse density
            real tV2 = this->v( T2, p2 );

            // RHS. the entropy row is scaled with R and not with the entropy
            // itself, so that the tolerance does not depend on where the
            // entropy scale is anchored
            tF( 0 ) = ( u2 * A2 / tV2 - aMass ) / aMass ;
            tF( 1 ) = ( this->s( T2, p2 ) - aEntropy ) / mR ;
            tF( 2 ) = ( this->h( T2, p2 ) + 0.5 * u2 * u2 - aEnergy ) / aEnergy ;

            if ( norm( tF ) < 1.0e-10 )
            {
                break ;
            }

            // Jacobian Terms

            // d(mass)/dT
            tJ( 0, 0 ) = -A2 * u2 * this->alpha( T2, p2 ) / ( tV2 * aMass );

            // d(entropy)/dT
            tJ( 1, 0 ) = this->dsdT( T2, p2 ) / mR ;

            // d(energy)/dT
            tJ( 2, 0 ) = this->cp( T2, p2 ) / aEnergy ;

            // d(mass)/dP
            // d( 1/v )/dp = +kappa/v, since kappa = -( 1/v ) ( dv/dp )_T
            tJ( 0, 1 ) = this->kappa( T2, p2 ) * A2 * u2 / ( tV2 * aMass );

            // d(entropy)/dP
            tJ( 1, 1 ) = this->dsdp( T2, p2 ) / mR ;

            //  d(energy)/dP
            tJ( 2, 1 ) = this->dhdp( T2, p2 ) / aEnergy ;

            // d(mass)/du
            tJ( 0, 2 ) = A2 / ( tV2 * aMass );

            // d(entropy)/du
            tJ( 1, 2 ) = 0.0 ;

            // d(energy)/du
            tJ( 2, 2 ) = u2 / aEnergy ;

            // solve system
            gesv( tJ, tF, mFlowPivot );

            T2 -= tOmega * tF( 0 );
            p2 -= tOmega * tF( 1 );
            u2 -= tOmega * tF( 2 );

            ++tCount ;

            BELFEM_ERROR( tCount < 100,
                         "Too many iterations for A2=%12.3f, T2=%12.3f, p2=%12.3f, u2=%12.3f",
                         A2,
                         T2,
                         p2,
                         u2 );
        }
    }

//------------------------------------------------------------------------------
    void
    Gas::expand(
            const real & A1,
            const real T1,
            const real p1,
            const real & u1,
            const real & A2,
                  real & T2,
                  real & p2,
                  real & u2 ) const
    {

        T2 = T1;
        p2 = p1;
        u2 = u1;

        if( std::abs( A1 - A2 ) < BELFEM_EPSILON )
        {
            return;
        }
        else
        {
            BELFEM_ERROR( A2 > A1,
                         "expand needs A2 >= A1, use compress for a narrowing duct ( A1=%12.3f, A2=%12.3f )",
                         A1,
                         A2 );

            real tOmega = 0.9;

            // mass
            real tMass = this->rho( T1, p1 ) * u1 * A1;

            // momentum
            real tMomentum = tMass * u1 + p1 * A2; // <-- A2, not A1 !

            real tEntropy = this->s( T1, p1 );

            real tEnergy = this->h( T1, p1 ) + 0.5 * u1 * u1;

            // simplified constans for ideal gas


            real tMa1 = u1 / this->c( T1, p1 );

            // RHS
            Vector<real> & tF = mFlowResidual;
            tF.fill( 0.0 );

            // Jacobian
            Matrix<real> & tJ = mFlowJacobian;

            // Pivot
            Vector< int_t > & tPivot = mFlowPivot;

            real tV2;

            uint tCount = 0;

            real tError = 1.0;

            // initial guess
            if( tMa1 < 1.0 )
            {
                real tK = this->gamma( T1, p1 );
                real tCp = this->cp( T1, p1 );

                // calculate a Borda-Carnot Shock
                // Rist: Dynamik Realer Gase, Kap. 9.1.1

                // Eq. ( 9.15 )
                u2 = u1 * ( tK * tMa1 * tMa1 + A2 / A1 -
                   std::sqrt( std::pow( tMa1 * tMa1 - 1.0, 2 )
                   + ( A2 / A1 - 1.0 ) *
                     ( 1.0 + 2.0 * tK * tMa1 * tMa1 + A2 / A1 ) ) ) /
                     ( ( tK + 1.0 ) * tMa1 * tMa1 );


                real tHt = tCp * T1 + 0.5 * u1 * u1;

                T2 = ( tHt - 0.5 * u2 * u2 ) / tCp;
                p2 = ( tMass * ( u1 - u2 ) + p1 * A2 ) / A2;


                while ( tError > 1e-6 )
                {
                    // get values
                    T2 -= tOmega * tF( 0 );
                    p2 -= tOmega * tF( 1 );
                    u2 -= tOmega * tF( 2 );

                    // calculate inverse density
                    tV2 = this->v( T2, p2 );

                    // Jacobian Terms

                    // d(mass)/dT
                    tJ( 0, 0 ) = -A2 * u2 * this->alpha( T2, p2 ) / ( tV2 * tMass );

                    // d(momentum)/dT
                    tJ( 1, 0 ) = tJ( 0, 0 ) * u2 * tMass / tMomentum;

                    // d(energy)/dT
                    tJ( 2, 0 ) = this->cp( T2, p2 ) / tEnergy;

                    // d(mass)/dP
                    // d( 1/v )/dp = +kappa/v, since kappa = -( 1/v ) ( dv/dp )_T
                    tJ( 0, 1 ) = this->kappa( T2, p2 ) * A2 * u2 / ( tV2 * tMass );

                    // d(momentum)/dP
                    tJ( 1, 1 ) = ( A2 + u2 * tJ( 0, 1 ) * tMass ) / tMomentum;

                    //  d(energy)/dP
                    tJ( 2, 1 ) = this->dhdp( T2, p2 ) / tEnergy;

                    // d(mass)/du
                    tJ( 0, 2 ) = A2 / ( tV2 * tMass );

                    // d(momentum)/du
                    tJ( 1, 2 ) = 2.0 * u2 * tJ( 0, 2 ) * tMass / tMomentum;

                    // d(energy)/du
                    tJ( 2, 2 ) = u2 / tEnergy;


                    // RHS
                    tF( 0 ) = ( u2 * A2 / tV2 - tMass ) / tMass;
                    tF( 1 ) = ( u2 * u2 * A2 / tV2 + p2 * A2 - tMomentum ) / tMomentum;
                    tF( 2 ) = ( this->h( T2, p2 ) + 0.5 * u2 * u2 - tEnergy ) / tEnergy;

                    tError = norm( tF );

                    // solve system
                    gesv( tJ, tF, tPivot );

                    // increment counter
                    ++tCount;

                    BELFEM_ERROR( tCount < 100,
                                 "Too many iterations for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                                 T1,
                                 p1,
                                 u1,
                                 A1,
                                 A2 );

                }

            }
            else // supersonic
            {
                // the flow expands around the corner instead of shocking, so
                // the step face does work on it and momentum is not conserved
                // in the control volume. mass, entropy and energy close the
                // system, which is the problem a contraction poses as well
                this->area_guess( T1, p1, tMa1, A1, A2, true, T2, p2, u2 );

                this->isentropic_duct( tMass, tEntropy, tEnergy, A2, T2, p2, u2 );

                // check mach number
                BELFEM_ERROR( u2 / this->c( T2, p2 ) > 1.0,
                             "Could not find supersonic solution for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                             T1,
                             p1,
                             u1,
                             A1,
                             A2 );

            }
        }
    }

//------------------------------------------------------------------------------

    void
    Gas::compress(
            const real & A1,
            const real T1,
            const real p1,
            const real & u1,
            const real & A2,
                  real & T2,
                  real & p2,
                  real & u2 ) const
    {
        T2 = T1;
        p2 = p1;
        u2 = u1;

        if( std::abs( A1 - A2 ) < BELFEM_EPSILON )
        {
            return;
        }
        else
        {
            BELFEM_ERROR( A2 < A1,
                         "compress needs A2 <= A1, use expand for a widening duct ( A1=%12.3f, A2=%12.3f )",
                         A1,
                         A2 );

            // mass
            real tMass = this->rho( T1, p1 ) * u1 * A1;

            real tEntropy = this->s( T1, p1 );

            real tEnergy = this->h( T1, p1 ) + 0.5 * u1 * u1;

            real tMa1 = u1 / this->c( T1, p1 );

            bool tSupersonic = tMa1 > 1.0 ;

            // a contraction accelerates a subsonic flow and decelerates a
            // supersonic one. it does not separate, so there is no Borda-Carnot
            // loss and the flow stays isentropic on either branch, only the
            // root of the area relation differs
            real tSonicArea = A1 / this->area_mach( tMa1, this->gamma( T1, p1 ) );

            BELFEM_ERROR( A2 >= tSonicArea,
                         "The duct is choked, A2=%12.6f m^2 is below the sonic area of %12.6f m^2",
                         A2,
                         tSonicArea );

            this->area_guess( T1, p1, tMa1, A1, A2, tSupersonic, T2, p2, u2 );

            this->isentropic_duct( tMass, tEntropy, tEnergy, A2, T2, p2, u2 );

            // the contraction must not switch branches
            BELFEM_ERROR( ( u2 / this->c( T2, p2 ) > 1.0 ) == tSupersonic,
                         "The flow changed regime inside the contraction for T1=%12.3f, p1=%12.3f, u1=%12.3f, A1=%12.3f, A2=%12.3f",
                         T1,
                         p1,
                         u1,
                         A1,
                         A2 );
        }
    }

//------------------------------------------------------------------------------

    namespace prandtlmeyer
    {
        /**
         * characteristic integral for the Prandtl-Meyer turn of a thermally
         * perfect ideal gas, the working object behind Gas::prandtl_meyer.
         *
         * the simple wave relation d(nu) = sqrt( Ma^2 - 1 ) * dV / V and the
         * adiabatic energy conservation dh + V dV = 0 combine to
         *
         *     nu( T ) = int_T^T1  cp * sqrt( Ma^2 - 1 ) / V^2  dT'
         *
         * along the isentrope through the upstream state, with
         * V^2 = 2 * ( ht - h ) and Ma^2 = V^2 / ( gamma * R * T ). the
         * integral carries the stagnation enthalpy of the upstream state, so
         * unlike the calorically perfect closed form, nu is not a state
         * function of the Mach number alone. one object serves one upstream
         * state.
         *
         * the caloric data of an ideal gas mixture is a piecewise cubic
         * spline on a uniform temperature grid. the quadrature panels align
         * with the spline knots, so a fixed Gauss rule per panel integrates
         * to round-off, and the panel integrals accumulated from T1 downward
         * are cached over the Newton iterations.
         */
        class Wave
        {
            //! gas that provides the caloric model, must be an ideal gas
            const Gas & mGas ;

            //! Gauss rule on [ -1, 1 ], owned by the gas object
            const Vector< real > & mX ;
            const Vector< real > & mW ;

            //! upstream temperature
            const real mT1 ;

            //! pressure handed to the caloric calls ( dummy for an ideal gas )
            const real mP1 ;

            //! specific gas constant
            const real mR ;

            //! stagnation enthalpy of the upstream state
            const real mHt ;

            //! knot spacing of the caloric spline
            const real mDeltaT ;

            //! largest spline knot at or below T1
            real mTknot ;

            //! stagnation temperature, root of h( T ) = ht
            real mT0 ;

            //! sonic temperature, root of 2 * ( ht - h ) = gamma * R * T
            real mTsonic ;

            //! cached integrals from knot mTknot - j * mDeltaT up to T1
            Cell< real > mNu ;

            //! number of valid entries in mNu
            index_t mNumCached = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Wave(   const Gas & aGas,
                    const Vector< real > & aGaussPoints,
                    const Vector< real > & aGaussWeights,
                    const real T1,
                    const real p1,
                    const real u1 ) :
                mGas( aGas ),
                mX( aGaussPoints ),
                mW( aGaussWeights ),
                mT1( T1 ),
                mP1( p1 ),
                mR( aGas.R( T1, p1 ) ),
                mHt( aGas.h( T1, p1 ) + 0.5 * u1 * u1 ),
                mDeltaT( aGas.heat_spline().delta_x() )
            {
                BELFEM_ERROR( mGas.is_idgas(),
                        "prandtl_meyer needs an ideal gas model, the characteristic integral only closes if h is a function of T alone" );

                BELFEM_ERROR( T1 > gastables::gTmin,
                        "upstream temperature T1=%12.3f K undercuts the caloric table limit of %12.3f K",
                        T1,
                        gastables::gTmin );

                BELFEM_ERROR( u1 > mGas.c( T1, p1 ),
                        "upstream flow must be supersonic ( T1=%12.3f K, p1=%12.3f Pa, u1=%12.3f m/s )",
                        T1,
                        p1,
                        u1 );

                // largest spline knot at or below T1
                const real Tmin = mGas.heat_spline().x_min() ;
                mTknot = Tmin + std::floor( ( T1 - Tmin ) / mDeltaT ) * mDeltaT ;

                this->find_stagnation();
                this->find_sonic();

                // knot cache from T1 down to the table limit
                index_t tCapacity = ( index_t )
                        ( ( mTknot - gastables::gTmin ) / mDeltaT ) + 2 ;

                mNu.set_size( tCapacity, BELFEM_QUIET_NAN );

                // top piece from the anchor knot to T1
                mNu( 0 ) = this->panel( mTknot, mT1 );
                mNumCached = 1 ;
            }

//------------------------------------------------------------------------------

            ~Wave() = default ;

//------------------------------------------------------------------------------

            /**
             * integrand cp * sqrt( Ma^2 - 1 ) / V^2 , the negative slope
             * d(nu)/dT. vanishes at the sonic point, which makes the
             * quadrature benign there while the inverse map T( nu ) is
             * singular
             */
            real
            f( const real T ) const
            {
                const real cp = mGas.cp( T, mP1 );
                const real k  = cp / ( cp - mR );

                const real V2  = 2.0 * ( mHt - mGas.h( T, mP1 ) );
                const real Ma2 = V2 / ( k * mR * T );

                // Ma2 - 1 vanishes at the sonic point, clip round-off only
                return cp * std::sqrt( std::max( Ma2 - 1.0, 0.0 ) ) / V2 ;
            }

//------------------------------------------------------------------------------

            /**
             * turning angle in rad from the upstream state to temperature T
             * on the same isentrope, positive for T < T1
             */
            real
            nu( const real T )
            {
                // above the anchor knot: the compression side and the top
                // piece of the expansion side, both at most a few panels long
                if ( T >= mTknot )
                {
                    return T <= mT1 ?
                            this->integrate( T, mT1 ) :
                           -this->integrate( mT1, T );
                }

                // expansion below the anchor knot: smallest knot at or above T
                const index_t j = ( index_t )
                        std::floor( ( mTknot - T ) / mDeltaT );

                BELFEM_ASSERT( j + 1 < mNu.size(),
                        "T=%12.3f K undercuts the knot cache", T );

                // grow the cache knot by knot: the step from entry j - 1 to
                // entry j adds the interval [ Kj, Kj + dT ] just below the
                // already covered range
                while ( mNumCached <= j )
                {
                    const real Ka = mTknot - ( real ) mNumCached * mDeltaT ;

                    mNu( mNumCached ) = mNu( mNumCached - 1 )
                            + this->panel( Ka, Ka + mDeltaT );

                    ++mNumCached ;
                }

                return mNu( j )
                        + this->panel( T, mTknot - ( real ) j * mDeltaT );
            }

//------------------------------------------------------------------------------

            /**
             * temperature after turning by alpha, solved from
             * nu( T ) = alpha with a Newton iteration that a bisection
             * bracket safeguards against the sonic singularity of the
             * inverse map. Tguess enters as the initial iterate
             */
            real
            solve_T( const real alpha, const real Tguess )
            {
                if ( std::abs( alpha ) < 1e-14 )
                {
                    return mT1 ;
                }

                // bracket [ Ta, Tb ] with F( Ta ) > 0 and F( Tb ) < 0,
                // where F = nu - alpha decreases in T
                real Ta ;
                real Tb ;

                if ( alpha > 0.0 )
                {
                    // expansion towards lower temperature
                    Ta = gastables::gTmin ;
                    Tb = mT1 ;

                    BELFEM_ERROR( this->nu( Ta ) > alpha,
                            "expansion by alpha=%12.6f rad leaves the caloric table, T2 would undercut %12.3f K ( T1=%12.3f K )",
                            alpha,
                            gastables::gTmin,
                            mT1 );
                }
                else
                {
                    // smooth compression towards the sonic point
                    Ta = mT1 ;
                    Tb = mTsonic ;

                    BELFEM_ERROR( this->nu( Tb ) < alpha,
                            "compression by alpha=%12.6f rad reaches the sonic point, the isentropic limit is %12.6f rad ( T1=%12.3f K )",
                            alpha,
                            this->nu( Tb ),
                            mT1 );
                }

                real T = std::min( std::max( Tguess, Ta ), Tb );

                uint tCount = 0 ;

                while ( true )
                {
                    const real F = this->nu( T ) - alpha ;

                    if ( std::abs( F ) < 1e-12 )
                    {
                        break ;
                    }

                    // maintain the bracket
                    if ( F > 0.0 )
                    {
                        Ta = T ;
                    }
                    else
                    {
                        Tb = T ;
                    }

                    if ( Tb - Ta < 1e-12 * mT1 )
                    {
                        break ;
                    }

                    // Newton step, dF/dT = -f. near the sonic point f
                    // vanishes and the step degenerates, then the bisection
                    // takes over
                    const real df = this->f( T );

                    real Tn = df > BELFEM_EPSILON ?
                            T + F / df : 0.5 * ( Ta + Tb );

                    if ( ! ( Tn > Ta && Tn < Tb ) )
                    {
                        Tn = 0.5 * ( Ta + Tb );
                    }

                    T = Tn ;

                    BELFEM_ERROR( tCount++ < 200,
                            "no convergence for alpha=%12.6f rad ( T1=%12.3f K )",
                            alpha,
                            mT1 );
                }

                // the compression branch models smooth isentropic
                // compression only and must stay clear of the sonic point
                if ( alpha < 0.0 )
                {
                    const real cp = mGas.cp( T, mP1 );
                    const real k  = cp / ( cp - mR );
                    const real Ma2 = 2.0 * ( mHt - mGas.h( T, mP1 ) )
                            / ( k * mR * T );

                    BELFEM_ERROR( Ma2 > ( 1.0 + 1e-6 ),
                            "compression by alpha=%12.6f rad ends at Ma=%12.8f, too close to the sonic point ( T1=%12.3f K )",
                            alpha,
                            std::sqrt( std::max( Ma2, 0.0 ) ),
                            mT1 );
                }

                return T ;
            }

//------------------------------------------------------------------------------

            real
            T_sonic() const
            {
                return mTsonic ;
            }

//------------------------------------------------------------------------------

            real
            ht() const
            {
                return mHt ;
            }

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * Gauss integral of f over one panel, Ta and Tb must lie in
             * the same spline interval
             */
            real
            panel( const real Ta, const real Tb ) const
            {
                const real c0 = 0.5 * ( Tb + Ta );
                const real c1 = 0.5 * ( Tb - Ta );

                const uint n = mW.length() ;

                real tSum = 0.0 ;
                for ( uint k = 0; k < n; ++k )
                {
                    tSum += mW( k ) * this->f( c0 + c1 * mX( k ) );
                }

                return c1 * tSum ;
            }

//------------------------------------------------------------------------------

            /**
             * knot aligned integral of f over [ Ta, Tb ], Ta <= Tb
             */
            real
            integrate( const real Ta, const real Tb ) const
            {
                BELFEM_ASSERT( Ta <= Tb, "integrate needs Ta <= Tb" );

                // first knot strictly above Ta on the grid anchored at mTknot
                real tLo = Ta ;
                real tHi = mTknot
                        + ( std::floor( ( Ta - mTknot ) / mDeltaT ) + 1.0 )
                        * mDeltaT ;

                real tSum = 0.0 ;
                while ( tHi < Tb )
                {
                    tSum += this->panel( tLo, tHi );
                    tLo   = tHi ;
                    tHi  += mDeltaT ;
                }

                tSum += this->panel( tLo, Tb );

                return tSum ;
            }

//------------------------------------------------------------------------------

            /**
             * stagnation temperature from h( T0 ) = ht
             */
            void
            find_stagnation()
            {
                // h grows monotonically in T, Newton from below
                real T = mT1 ;

                uint tCount = 0 ;

                while ( true )
                {
                    const real dT = ( mHt - mGas.h( T, mP1 ) )
                            / mGas.cp( T, mP1 );

                    T += dT ;

                    // the sonic solve needs caloric data up to T0, and the
                    // spline extrapolates silently past its table limit
                    BELFEM_ERROR( T < gastables::gTmax,
                            "stagnation temperature exceeds the caloric table limit of %12.3f K ( T1=%12.3f K )",
                            gastables::gTmax,
                            mT1 );

                    if ( std::abs( dT ) < 1e-10 * T )
                    {
                        break ;
                    }

                    BELFEM_ERROR( tCount++ < 100,
                            "no convergence for the stagnation temperature ( T1=%12.3f K )",
                            mT1 );
                }

                mT0 = T ;
            }

//------------------------------------------------------------------------------

            /**
             * sonic temperature from 2 * ( ht - h ) = gamma * R * T,
             * bracketed by [ T1, T0 ]
             */
            void
            find_sonic()
            {
                // g = V^2 - a^2 falls monotonically from g( T1 ) > 0 in
                // the supersonic state to g( T0 ) = -a^2 < 0 at stagnation
                real Ta = mT1 ;
                real Tb = mT0 ;

                real T = 0.5 * ( Ta + Tb );

                uint tCount = 0 ;

                while ( true )
                {
                    const real cp = mGas.cp( T, mP1 );
                    const real k  = cp / ( cp - mR );
                    const real a2 = k * mR * T ;

                    const real g = 2.0 * ( mHt - mGas.h( T, mP1 ) ) - a2 ;

                    if ( std::abs( g ) < 1e-12 * a2 || Tb - Ta < 1e-12 * mT1 )
                    {
                        break ;
                    }

                    if ( g > 0.0 )
                    {
                        Ta = T ;
                    }
                    else
                    {
                        Tb = T ;
                    }

                    const real dkdT = -mR * mGas.dcpdT( T, mP1 )
                            / std::pow( cp - mR, 2 );

                    const real dg = -2.0 * cp - mR * ( k + T * dkdT );

                    real Tn = T - g / dg ;

                    if ( ! ( Tn > Ta && Tn < Tb ) )
                    {
                        Tn = 0.5 * ( Ta + Tb );
                    }

                    T = Tn ;

                    BELFEM_ERROR( tCount++ < 100,
                            "no convergence for the sonic temperature ( T1=%12.3f K )",
                            mT1 );
                }

                mTsonic = T ;
            }

//------------------------------------------------------------------------------
        };
    } /* namespace prandtlmeyer */

//------------------------------------------------------------------------------

    real
    Gas::prandtl_meyer_angle( const real T, const real p, const real & u ) const
    {

        real k = this->gamma( T, p );
        real Ma = u / this->c( T, p );

        return std::sqrt( ( k + 1. ) / ( k - 1. ) )
               * std::atan( std::sqrt( (  k - 1. ) / ( k + 1. ) *
               ( Ma * Ma - 1. ) ) )
               - std::atan( std::sqrt(  Ma * Ma - 1. ) );

    }

//------------------------------------------------------------------------------

    real
    Gas::prandtl_meyer(
                    const real T1,
                    const real p1,
                    const real & u1,
                    const real & alpha,
                          real & T2,
                          real & p2,
                          real & u2 ) const
    {
        // Gauss rule for the knot panels, computed on the first call.
        // fourteen points; on one spline interval cp is a quadratic, so far
        // fewer would already be exact and the surplus only buys round-off
        // margin
        if( mGaussWeights.length() == 0 )
        {


            int_t n = 14 ;
            mGaussPoints.set_size( n );
            mGaussWeights.set_size( n );

            intpoints_gauss( &n, mGaussWeights.data(), mGaussPoints.data() );
        }

        // characteristic integral along the isentrope of the upstream state,
        // errors out unless the gas model is ideal and the upstream flow is
        // supersonic
        prandtlmeyer::Wave tWave( *this, mGaussPoints, mGaussWeights,
                                  T1, p1, u1 );

        // initial guess from the calorically perfect closed form. on the
        // compression branch the perfect gas Newton would run into the
        // sonic singularity, the midpoint of the bracket serves instead
        real Tguess ;

        if( alpha < 0.0 )
        {
            Tguess = 0.5 * ( T1 + tWave.T_sonic() );
        }
        else
        {
            real k = this->gamma( T1, p1 );
            real Ma = u1 / this->c( T1, p1 );

            const real Ma1 = Ma;

            const real aNu2 = alpha + this->prandtl_meyer_angle( T1, p1, u1 );

            const real G = std::sqrt(( k + 1. ) / ( k - 1. ));

            real L;
            real f = 1.0;
            real df;

            uint tCount = 0;

            // best effort: the exact solver below only needs a starting
            // point, so bail out if the Newton overshoots past the sonic
            // singularity ( possible for a barely supersonic upstream )
            while( std::abs(f) > 0.001 && tCount++ < 50 )
            {
                L = std::sqrt( Ma * Ma - 1 );
                f = G * std::atan( L / G ) - std::atan( L ) - aNu2;
                df = ( G * G / ( L * L + G * G ) - 1. / ( 1. + L * L ))* Ma / L;

                Ma -= f/df;

                if( ! ( Ma > 1.0 + 1e-6 ) )   // also catches NaN
                {
                    Ma = Ma1;
                    break;
                }
            }

            // perfect gas temperature at the estimated Mach number
            const real Tt = T1 * ( 1.0 + 0.5 * ( k - 1. ) * Ma1 * Ma1 );

            Tguess = Tt / ( 1.0 + 0.5 * ( k - 1. ) * Ma * Ma );
        }

        // exact turn
        T2 = tWave.solve_T( alpha, Tguess );

        // pressure decouples and follows from ds = 0
        p2 = p1 * std::exp( ( this->s( T2, p1 ) - this->s( T1, p1 ) )
                / this->R( T1, p1 ) );

        // velocity from energy conservation
        u2 = std::sqrt( 2.0 * ( tWave.ht() - this->h( T2, p1 ) ) );

        return u2 / this->c( T2, p2 );
    }

//------------------------------------------------------------------------------

    void
    Gas::shock(  const real T1, const real p1, const real & u1,
                       real & T2,       real & p2,       real & u2 ) const
    {
        // relaxation factor
        real tOmega0 = 0.9;
        real tOmega1 ;
        real tOmega ;

        // mach number before shock
        real tMa1 = u1 / this->c( T1, p1 );

        if ( tMa1 < 1.001 )
        {
            T2 = T1;
            p2 = p1;
            u2 = u1;
        }
        else
        {
            // total enthalpy
            real tHt = this->h( T1, p1 ) + 0.5 * u1 * u1;

            real tGamma = this->gamma( T1, p1 );

            // state vector: rho2, u2, T2
            Vector< real > tX( 3 );
            real & tRho2 = tX( 0 );
            real & tU2 = tX( 1 );
            real & tT2 = tX( 2 );

            // total temperature
            real tTt = T1 * ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa1 * tMa1 );

            // initial guess for Ma2
            real tMa2 = std::sqrt(
                    ( tTt / T1 ) / ( tGamma * tMa1 * tMa1 - 0.5 * ( tGamma - 1.0 )));



            // initial guess for temperature
            tT2 = tTt / ( 1.0 + 0.5 * ( tGamma - 1.0 ) * tMa2 * tMa2 );

            // check for temperature limit
            if( tT2 > gastables::gTmax - 1000.0 )
            {
                tT2 = gastables::gTmax - 1000.0 ;
            }

            // initial guess for pressure
            p2 = p1 * std::pow( tT2 / T1, tGamma / ( tGamma - 1.0 ));

            // initial guess for density
            tRho2 = this->rho( tT2, p2 );

            // initial guess for velocity
            tU2 = tMa2 * std::sqrt( this->R( tT2, p2 ) * tGamma * tT2 );

            // mass
            real tM = this->rho( T1, p1 ) * u1;

            // momentum
            real tI = p1 + tM * u1;

            // Jacobian
            Matrix< real > tJ( 3, 3 );

            // right hand side
            Vector< real > tY( 3 );

            real tResiduum = BELFEM_REAL_MAX;

            // counter for the loop
            uint tCount = 0;

            // for solver
            Vector< int_t > tPivot( 3 );

            real tTmax = gastables::gTmax - 1.0 ;

            while ( tResiduum > 1.0e-9 )
            {
                // compute Jacobian
                tJ( 0, 0 ) = tU2;
                //tJ( 1, 0 ) = 1.0 / ( this->kappa( tT2, p2 ) * tRho2 ) + tU2 * tU2 ;
                tJ( 1, 0 ) = ( this->R( tT2, p2 ) * tT2 + tU2 * tU2 );

                tJ( 2, 0 ) = 0.0;

                tJ( 0, 1 ) = tRho2;
                tJ( 1, 1 ) = 2.0 * tRho2 * tU2;
                tJ( 2, 1 ) = tU2;

                tJ( 0, 2 ) = 0.0;

                //tJ( 1, 2 ) = this->beta( tT2, p2 ) ;
                tJ( 1, 2 ) = this->R( tT2, p2 ) * tRho2;
                tJ( 2, 2 ) = this->cp( tT2, p2 );

                // compute right hand side
                p2 = this->p( tT2, 1.0 / tRho2 );

                tY( 0 ) = ( tRho2 * tU2 - tM );
                tY( 1 ) = ( p2 + tRho2 * tU2 * tU2 - tI );
                tY( 2 ) = ( this->h( tT2, p2 ) + 0.5 * tU2 * tU2 - tHt );

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
                    this->total( T1, p1, u1, T2, p2 ) ;

                    // pressure
                    p2 = this->p( tX( 2 ), 1.0 / tX( 0 ));

                    u2 = tX( 1 );

                    return ;
                }
                // perform Newton Step
                tX -= tOmega * tY;

                BELFEM_ERROR( tCount++ < 1000,
                             "Infinite loop in Gas::shock" );
            }

            // postprocess output values
            T2 = tX( 2 );
            p2 = this->p( tX( 2 ), 1.0 / tX( 0 ));
            u2 = tX( 1 );
        }
 }

//------------------------------------------------------------------------------

    void
    Gas::shock( const real T1, const real p1, const real & u1, const real & alpha,
           real & T2, real & p2, real & u2, real & beta ) const
    {
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // Step 1: indentify min and max possible beta-angle
        // - - - - - - - - - - - - - - - - - - - - - - - - -

        // compute mach number
        real tMa1 = u1 / this->c( T1, p1 );

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

        real tFA = this->shock_beta( T1, p1, u1, alpha, T2, p2, u2, tBetaA );
        real tFB = tFA;

        tBetaB = tBetaA;

        uint tCountB = 0;

        while( tFB > 0.0 )
        {
            // shift F
            tFA = tFB;

            // shift Beta
            tBetaA = tBetaB;

            // increment beta
            tBetaB += tDeltaB;

            tFB = this->shock_beta_simple( T1, p1, u1, alpha, T2, p2, u2, tBetaB );

            BELFEM_ERROR( tCountB++ < 100,
                    "No sign change found while bracketing the shock angle" );
        }

        // one more step for safety
        tBetaB += tDeltaB;


        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // Step 2: initial iteration using bisection
        // - - - - - - - - - - - - - - - - - - - - - - - - -


        // loop counter
        uint tCount = 0;
        beta  = tBetaB;

        real tF =  this->shock_beta( T1, p1, u1, alpha, T2, p2, u2, beta );

        while ( std::abs( tF ) > 1.0e-3 )
        {

            // new beta
            beta = 0.5 * ( tBetaA + tBetaB );


            // call beta function
            tF = this->shock_beta( T1, p1, u1, alpha, T2, p2, u2, beta );

            // test result
            if( tFA * tF > 0 )
            {
                tBetaA = beta;
                tFA = tF;
            }
            else
            {
                tBetaB = beta;
            }

            // increment counter
            BELFEM_ERROR( tCount++ < 1000, "too many iterations" );
        }


    }


//------------------------------------------------------------------------------
    real
    Gas::v( const uint aIndex, const real T, const real p ) const
    {
        return mEoS->v( aIndex, T, p );
    }
//------------------------------------------------------------------------------

    real
    Gas::h( const uint aIndex, const real T, const real p ) const
    {
        // read-only handle: the component is only evaluated here
        const gastables::RefGas * tComponent = mComponents( aIndex );

        return tComponent->h( T )
               / tComponent->data()->M()
               + mEoS->hdep( aIndex, T, p );
    }

//------------------------------------------------------------------------------

    real
    Gas::cp( const uint aIndex, const real T, const real p ) const
    {
        // read-only handle: the component is only evaluated here
        const gastables::RefGas * tComponent = mComponents( aIndex );

        return tComponent->heat_spline()->deval( T )
            / tComponent->data()->M()
            + mEoS->cpdep( aIndex, T, p );
    }

    real
    Gas::dcpdT( const uint aIndex, const real T, const real p ) const
    {
        // read-only handle: the component is only evaluated here
        const gastables::RefGas * tComponent = mComponents( aIndex );

        // ideal gas contribution
        real tdCpdT0 = tComponent->heat_spline()->ddeval( T )
                    /  tComponent->data()->M() ;

        // departure Hi
        real tT1 = std::min( 1.001 * T, gastables::gTmax );
        real tdCpdepdT1 = mEoS->cpdep( aIndex, tT1, p ) ;

        // departure Low
        real tT0 = std::max( 0.999 * T, gastables::gTmin );
        real tdCpdepdT0 = mEoS->cpdep( aIndex, tT0, p ) ;

        return tdCpdT0 + ( tdCpdepdT1 - tdCpdepdT0 ) / ( tT1 - tT0 );
    }

//------------------------------------------------------------------------------

    void
    Gas::print() const
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    const gastables::GasData *
    Gas::data( const index_t aIndex ) const
    {
        BELFEM_ASSERT( aIndex < mNumberOfComponents,
                      "Requested Index is out of bounds ( %lu vs %lu )",
                      ( long unsigned int ) aIndex,
                      ( long unsigned int ) mNumberOfComponents );

        // read-only handle on the component, so the const overload of
        // RefGas::data() is the one that answers
        const gastables::RefGas * tComponent = mComponents( aIndex );

        return tComponent->data();
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
    Gas::shock_beta_simple( const real T1,
                const real p1,
                const real & u1,
                const real & alpha,
                real & T2,
                real & p2,
                real & u2,
                const real & beta ) const
    {
        real tV  = u1 * std::cos( beta );
        real tU1 = u1 * std::sin( beta );
        real tU2;
        real tK = 1.4;

        // compute oblique shock for perfect gas
        real tMa1 = u1 / std::sqrt( tK * this->R( T1, p1 ) * T1 );
        real tTt  = T1 * ( 1.0 + 0.5 * ( tK - 1.0 ) * tMa1 * tMa1 );

        real tMa2 = std::sqrt( tTt / ( T1 * ( tK * tMa1 * tMa1 - 0.5 * ( tK - 1.0 ) ) ) );

        T2 = tTt / ( 1.0 + 0.5 * ( tK - 1.0 ) * tMa2 * tMa2 );
        p2 = p1 * std::pow( ( T2 / T1 ), tK / ( tK - 1.0 ) );
        tU2 = tMa2 * std::sqrt( tK * this->R( T1, p1 ) * T2 );

        // compute the velocity
        u2 = std::sqrt( tU2 * tU2 + tV * tV );

        // return the function
        return tU2 / tU1 - std::tan( beta - alpha ) / std::tan( beta );
    }

//------------------------------------------------------------------------------

    // special subroutine needed for shock
    real
    Gas::shock_beta( const real T1,
                     const real p1,
                     const real & u1,
                     const real & alpha,
                     real & T2,
                     real & p2,
                     real & u2,
                     const real & beta ) const
    {
        real tV  = u1 * std::cos( beta );
        real tU1 = u1 * std::sin( beta );
        real tU2;

        // compute oblique shock
        if( tU1 > this->c( T1, p1 ) )
        {
            this->shock( T1, p1, tU1, T2, p2, tU2 );
        }
        else
        {
            T2 = T1;
            p2 = p1;
            tU2 = tU1;
        }

        // compute the velocity
        u2 = std::sqrt( tU2 * tU2 + tV * tV );

        // return the function
        return tU2 / tU1 - std::tan( beta - alpha ) / std::tan( beta );
    }

//------------------------------------------------------------------------------

    // enthalpy derivative to pressure ( needed for total temperature )
    real
    Gas::dhdp( const real T, const real p ) const
    {
        mStatevals.update_Tp( T, p );

        if ( !mStatevals.test( BELFEM_STATEVAL_DHDP ) )
        {
            mStatevals.set( BELFEM_STATEVAL_DHDP,
                            ( this->*mFunctionDHDP )( T, p ));
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