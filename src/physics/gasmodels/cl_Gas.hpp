//
// Created by Christian Messe on 01.09.19.
//

#ifndef BELFEM_CL_GAS_HPP
#define BELFEM_CL_GAS_HPP

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Spline.hpp"

#include "cl_GM_Statevals.hpp"
#include "cl_GT_GasData.hpp"
#include "cl_GM_EoS.hpp"
#include "en_GM_GasModel.hpp"
#include "en_Helmholtz.hpp"

namespace belfem
{
    namespace gastables
    {
        class RefGas;

        class RefGasFactory;
    }

//------------------------------------------------------------------------------

    namespace gasmodels
    {
        class Statevals;
        class EoS;
        class HelmholtzTransport ;
    }

//------------------------------------------------------------------------------
     /**
       * \brief The gas class that provides the fluid model
       */
    class Gas
    {
    protected:
        //const proc_t mMasterRank = 0;
        const proc_t mMyRank = comm_rank();

        //! size of components vector
        uint mNumberOfComponents;

        //! container for state variables
        gasmodels::Statevals mStatevals;

        //! x or chi
        Vector<real> mMolarFractions ;

        //! y or zeta
        Vector<real> mMassFractions;

        //! molar fractions at initialization
        Vector<real> mMolarFractions0 ;

        //! molar masses of components
        Vector<real> mMolarMasses;

        //! Molar Mass in kg/Mol
        const real & mM = mStatevals.get( BELFEM_STATEVAL_M ) ;

        //! Gas constant in J/(kg*K)
        const real & mR = mStatevals.get( BELFEM_STATEVAL_R ) ;

        //! Components of the mixgas
        Cell<gastables::RefGas *> mComponents;

        //! Reference gases for formation enthalpy
        //! ( may be redundant to mComponents )
        Cell<gastables::RefGas *> mElements;

        Cell<string> mElementNames;

        //! additional reference gases that are not present in
        //! the components list
        Cell<gastables::RefGas *> mExtra;

        //! help matrix for remixing
        SpMatrix mHelpMatrix;

        //! temperature steps for spline
        Vector<real> mTemperatureSteps;
        Spline mHeatSpline;
        Spline mViscositySpline;
        Spline mConductivitySpline;

        //! gasmodel type
        GasModel mGasModel = GasModel::UNDEFINED ;

        //! helmholz type, if used
        HelmholtzModel mHelmholzModel = HelmholtzModel::UNDEFINED ;

        //! Interaction polynomials for viscosity
        Cell<gastables::RefGas *> mViscosityInteractionRefgas;

        //! table telling if interaction parameter exists
        Matrix<uint> mViscosityInteractionTable;

        // Work matrices for viscosity and conductivity calculation of mixtures
        Matrix<real> mWorkMatrix;
        Vector<real> mWorkVector;
        Vector<real> mWorkVector2;
        Vector<real> mWorkMu;
        Vector<real> mWorkLambda;

        // work matrix and vectors for equilibrium
        Vector<real> mWorkVectorRAND0 ;
        Vector<real> mWorkVectorRAND1 ;
        Vector<real> mWorkVectorRAND2 ;
        Matrix<real> mWorkMatrixRAND ;
        Vector<int>  mPivotRAND ;

        real mWorkTemperature;

        // Work matrices for gibbs
        Matrix<real> mFormationTable;
        Vector<real> mFormationWork;

        // critical temperature
        real mTcrit = BELFEM_QUIET_NAN;

        // critical pressure
        real mPcrit = BELFEM_QUIET_NAN;

        // critical volume
        real mVcrit = BELFEM_QUIET_NAN;

        // constant for Stiel Thodos Equation
        real mGamma = 0.0;

        // constant for Lucas Equation
        real mXi    = 0.0;

        //! liquid flag. Currently not in use
        bool mLiquidFlag = false;

        //! the equation of state
        gasmodels::EoS    * mEoS;

        //! special class, only needed if this is a helmholtz eos
        gasmodels::HelmholtzTransport * mTransport = nullptr ;

        //! pointer to cp function
        real
        ( Gas::*mFunctionCp )
                ( const real & aT, const real & aP );

        //! pointer to dcpdT function
        real
        ( Gas::*mFunctiondCpdT )
                ( const real & aT, const real & aP );

        //! pointer to cv function
        real
        ( Gas::*mFunctionCv )
                ( const real & aT, const real & aP );

        //! pointer to gamma function
        real
        ( Gas::*mFunctionGamma )
                ( const real & aT, const real & aP );

        //! pointer to c function
        real
        ( Gas::*mFunctionC )
                ( const real & aT, const real & aP );

        //! pointer to h function
        real
        ( Gas::*mFunctionH )
                ( const real & aT, const real & aP );

        //! pointer to dhdp function
        real
        ( Gas::*mFunctionDHDP )
                ( const real & aT, const real & aP );

        //! pointer to s function
        real
        ( Gas::*mFunctionS )
                ( const real & aT, const real & aP );

        //! pointer to dsdT function
        real
        ( Gas::*mFunctionDSDT )
                ( const real & aT, const real & aP );

        //! pointer to dsdp function
        real
        ( Gas::*mFunctionDSDP )
                ( const real & aT, const real & aP );

        //! pointer to Mu function
        real
        ( Gas::*mFunctionMU )
                ( const real & aT, const real & aP );

        //! pointer to Lambda function
        real
        ( Gas::*mFunctionLAMBDA )
                ( const real & aT, const real & aP );

        //! reference enthalpy, needed for real gas s
        real mSref = 0.0;

        // support term for entropy
        real mMixtureEntropy = BELFEM_QUIET_NAN ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * default constructor, creates air as idgas
         */
        Gas();

//------------------------------------------------------------------------------

        /**
         * creates a pure gas
         */
         Gas( const string & aLabel, const GasModel aGasModel=GasModel::IDGAS );

//------------------------------------------------------------------------------

        /**
         * create a cryogenic gas
         */
         Gas( const HelmholtzModel aHelmholtzModel );

//------------------------------------------------------------------------------

        /**
         * advanced constructor
         * @param aSpecies          : List of Species names
         * @param aMolarFractions   : Molar fractions of species
         * @param aGasModel         : gas model ( default: IDGAS )
         */
        Gas(
                const Cell<string> & aSpecies,
                const Vector<real> & aMolarFractions,
                const GasModel       aGasModel=GasModel::IDGAS );

//------------------------------------------------------------------------------

        virtual ~Gas();

//------------------------------------------------------------------------------

        const uint &
        number_of_components() const;

//------------------------------------------------------------------------------

        virtual void
        remix( const Vector<real> & aMolarFractions,
               bool aRemixHeat=true,
               bool aRemixTransport=true );

//------------------------------------------------------------------------------

         virtual void
         remix_mass( const Vector<real> & aMassFractions,
                     bool aRemixHeat=true,
                     bool aRemixTransport=true );

//------------------------------------------------------------------------------

        // reset to mixture at initialization
        virtual void
        reset_mixture();

//------------------------------------------------------------------------------

        virtual const real &
        M( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        virtual const real &
        R( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        /**
         * expose statevals object
         */
        gasmodels::Statevals &
        statevals();

//------------------------------------------------------------------------------

        /**
         * expose the data object of a refgas
         */
         gastables::GasData *
         data( const index_t aIndex );
//------------------------------------------------------------------------------

         /**
          * expose state equation of a refgas
          */
         gasmodels::EoS *
         eos();

//------------------------------------------------------------------------------

        bool
        is_idgas() const;

//------------------------------------------------------------------------------

        /**
         * expose the molar fractions
         */
         const Vector< real > &
         molar_fractions() const;

//------------------------------------------------------------------------------

         /**
           * expose the mass fractions
           */
         const Vector< real > &
         mass_fractions() const;

//------------------------------------------------------------------------------

         /**
          * return one single milar fraction
          */
         const real &
         molar_fraction( const uint aIndex ) const;

//------------------------------------------------------------------------------

         /**
           * return one single mass fraction
           */
         const real &
         mass_fraction( const uint aIndex ) const;

//------------------------------------------------------------------------------

         /**
           * expose the element container
           */
         Cell< gastables::RefGas * > &
         elements();

//------------------------------------------------------------------------------

         /**
           * expose the component container
           */
         Cell< gastables::RefGas * > &
         components();
//------------------------------------------------------------------------------

         /**
           * expose one component
           */
         gastables::RefGas * &
         component( const index_t aIndex );

//------------------------------------------------------------------------------

         /**
           * expose the formation table, telling which component is
           * built from which element
           */
         const Matrix< real > &
         formation_table() const ;

//------------------------------------------------------------------------------

        /**
         * test if the gas is liquid. Currently always returns false.
         */
         bool
         is_liquid() const;

//------------------------------------------------------------------------------
        /**
         * set the liquid flag. Set e.g. by the Helmholtz EoS
         */
         void
         set_liquid_flag( const bool aFlag );

//------------------------------------------------------------------------------
// Thermodynamic States
//------------------------------------------------------------------------------

        real
        p( const real & aT, const real & aV );

//------------------------------------------------------------------------------

        real
        v( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        rho( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        T( const real & aP, const real & aV );

//------------------------------------------------------------------------------
// Caloric Properties
//------------------------------------------------------------------------------

        virtual real
        cp( const real & aT, const real & aP );

        virtual real
        cv( const real & aT, const real & aP );

        virtual real
        gamma( const real & aT, const real & aP );

        virtual real
        c( const real & aT, const real & aP );

        virtual real
        u( const real & aT, const real & aP );

        virtual real
        h( const real & aT, const real & aP );

        virtual real
        s( const real & aT, const real & aP );

        virtual real
        dsdT( const real & aT, const real & aP );

        virtual real
        dsdp( const real & aT, const real & aP );

        virtual real
        dcpdT( const real & aT, const real & aP );

        // dissociation enthalpy ( only for tablegas at this time )
        virtual real
        hd( const real & aT, const real & aP );

//------------------------------------------------------------------------------
// Transport Properties
//------------------------------------------------------------------------------

        /**
         * dynamic viscisity in Pa*s
         */
        real
        mu( const real & aT, const real & aP );

        /**
          * thermal conductivity in W/(m*K)
          */
        real
        lambda( const real & aT, const real & aP );

        /**
         * Prandtl Number
         */
        real
        Pr( const real & aT, const real & aP );

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

         /**
          * thermal expansion coefficient
          *
          * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
          */
        real
        alpha( const real & aT, const real & aP );

//------------------------------------------------------------------------------
         /**
          * isochoric stress coefficient
          *
          * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
          */
        real
        beta( const real & aT, const real & aP );

//------------------------------------------------------------------------------

         /**
          * isothermal compressibility coefficient
          *
          * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
          */
        real
        kappa( const real & aT, const real & aP );

//------------------------------------------------------------------------------
// CHEMISTRY
// -----------------------------------------------------------------------------

        /**
         * return the molar Gibbs potential at reference pressure
         */
         void
         Gibbs( const real & aT, Vector< real > & aGibbs );

         /**
         * return the molar formation enthalpy for each component
         */
         void
         Hf( const real & aT, Vector< real > & aHf );

         /**
          * return the temperature derivative of  gibs potentia at reference pressure
          */
         void
         dGibbsdT( const real & aT, Vector< real > & aGibbs );

         /**
          * remix to equilibrium
          */
         void
         remix_to_equilibrium(
                 const real & aT,
                 const real & aP,
                 const bool aRemixHeat=true,
                 const bool aRemixTransport=true );

         void
         compute_equilibrium( const real & aT, const real & aP, Vector< real > & aX );

//------------------------------------------------------------------------------
// State relevant methods
// -----------------------------------------------------------------------------

        real
        T_from_h( const real & aH, const real & aP );

// -----------------------------------------------------------------------------
        /**
         * get an isotropic temperature
         */
         real
         isen_T( const real & aT0, const real & aP0, const real & aP1 );

// -----------------------------------------------------------------------------

         /**
           * get an isotropic pressure
           */
         real
         isen_p( const real & aT0, const real & aP0, const real & aT1 );

// -----------------------------------------------------------------------------

        /**
         * calculate the total state
         */
         void
         total( const real & aT, const real & aP, const real & aU,
                      real & aTt, real & aPt );

// -----------------------------------------------------------------------------

        /**
         * isentropic expansion of a gas
         */
         void
         expand( const real & aA1,
                 const real & aT1,
                 const real & aP1,
                 const real & aU1,
                 const real & aA2,
                       real & aT2,
                       real & aP2,
                       real & aU2 );

// -----------------------------------------------------------------------------

        real
        prandtl_meyer_angle(
                const real & aT,
                const real & aP,
                const real & aU );

// -----------------------------------------------------------------------------

        /**
         * prandtl Meyer expansion, only for ideal gas
         * @param aIndex
         * @param aT
         * @param aP
         * @return
         */
         real
         prandtl_meyer(  const real & aT1,
                         const real & aP1,
                         const real & aU1,
                         const real & aAlpha,
                               real & aT2,
                               real & aP2,
                               real & aU2 );

// -----------------------------------------------------------------------------

        /**
         * perpendicular shock
         *
         */
         void
         shock(  const real & aT1, const real & aP1, const real & aU1,
                       real & aT2,       real & aP2,       real & aU2 );

// -----------------------------------------------------------------------------

        /**
         * oblique shock
         */
         void
         shock( const real & aT1, const real & aP1, const real & aU1, const real & aAlpha,
                real & aT2, real & aP2, real & aU2, real & aBeta );

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

         real
         v( const uint & aIndex, const real & aT, const real & aP );

         // specific enthalpy of gas, normalized to 298.15 K
         real
         h( const uint & aIndex, const real & aT, const real & aP );

         real
         cp( const uint & aIndex, const real & aT, const real & aP );

         real
         dcpdT( const uint & aIndex, const real & aT, const real & aP );

//------------------------------------------------------------------------------
// Print composition
//------------------------------------------------------------------------------

        void
        print();

//------------------------------------------------------------------------------
// Special access
//------------------------------------------------------------------------------

         /**
          * expose heat spline
          * @return
          */
         Spline &
         heat_spline() ;

         /**
          * expose viscosity spline
          * @return
          */
         Spline &
         viscosity_spline() ;

         /**
         * expose conductivity spline
         * @return
         */
         Spline &
         conductivity_spline() ;

         /**
          * returns what gas model is used
          */
         const GasModel &
         gas_model() const;

         /**
          * returns what helmhilzu model is used, of any
          */
         const HelmholtzModel &
         helmholtz_model() const;

//------------------------------------------------------------------------------

         // enthalpy derivative to pressure ( needed for total temperature )
         virtual real
         dhdp( const real & aT, const real & aP );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        create_eos( const GasModel & aGasModel );

//------------------------------------------------------------------------------

        void
        initialize(
                const Cell<string> & aSpecies,
                const Vector<real> & aMolarFractions,
                const GasModel       aGasModel );

//------------------------------------------------------------------------------

        void
        check_thermo_exists();

//------------------------------------------------------------------------------

        void
        remix_R( const Vector<real> & aMolarFractions );

//------------------------------------------------------------------------------

        void
        remix_heat();

//------------------------------------------------------------------------------

        void
        remix_critical_point();

//------------------------------------------------------------------------------

        void
        remix_transport();

//------------------------------------------------------------------------------

        void
        create_reference_gases(
                gastables::RefGasFactory & aFactory,
                const Cell<string> & aLables );

//------------------------------------------------------------------------------

        string
        element_to_molecule( const string & aElement );

//------------------------------------------------------------------------------

        string
        reference_element( const string & aElement );

//------------------------------------------------------------------------------

        void
        create_viscosity_table( gastables::RefGasFactory & aFactory );

//------------------------------------------------------------------------------

        void
        create_mass_properties( const Vector<real> & aMolarFractions );

//------------------------------------------------------------------------------

        void
        evaluate_viscosity_interaction( const real & aT );

//------------------------------------------------------------------------------

        void
        evaluate_conductivity_interaction( const real & aT );

//------------------------------------------------------------------------------

        real
        idgas_cp( const real & aT, const real & aP );

        real
        idgas_dcpdT( const real & aT, const real & aP );

        real
        idgas_cv( const real & aT, const real & aP );

        real
        idgas_gamma( const real & aT, const real & aP );

        real
        idgas_c( const real & aT, const real & aP );

        real
        idgas_h( const real & aT, const real & aP );

        real
        idgas_s( const real & aT, const real & aP );

        real
        idgas_dsdT( const real & aT, const real & aP );

        real
        idgas_dsdp( const real & aT, const real & aP );

        real
        idgas_mu( const real & aT, const real & aP );

        real
        idgas_lambda( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        realgas_cp( const real & aT, const real & aP );

        real
        realgas_dcpdT( const real & aT, const real & aP );

        real
        realgas_cv( const real & aT, const real & aP );

        real
        realgas_gamma( const real & aT, const real & aP );

        real
        realgas_c( const real & aT, const real & aP );

        real
        realgas_h( const real & aT, const real & aP );

        real
        realgas_s( const real & aT, const real & aP );

        real
        realgas_dsdT( const real & aT, const real & aP );

        real
        realgas_dsdp( const real & aT, const real & aP );

        real
        realgas_mu( const real & aT, const real & aP );

        real
        realgas_lambda( const real & aT, const real & aP );

//------------------------------------------------------------------------------

         real
         helmholtz_cp( const real & aT, const real & aP );

         real
         helmholtz_dcpdT( const real & aT, const real & aP );

         real
         helmholtz_cv( const real & aT, const real & aP );

         real
         helmholtz_gamma( const real & aT, const real & aP );

         real
         helmholtz_c( const real & aT, const real & aP );

         real
         helmholtz_h( const real & aT, const real & aP );

         real
         helmholtz_s( const real & aT, const real & aP );

         real
         helmholtz_dsdT( const real & aT, const real & aP );

         real
         helmholtz_dsdp( const real & aT, const real & aP );

         real
         helmholtz_mu( const real & aT, const real & aP );

         real
         helmholtz_lambda( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        // from cea, Eq. 5.3
        real
        cea_mu( const real & aT );

//------------------------------------------------------------------------------

        // from cea, Eq. 5.4
        real
        cea_lambda( const real & aT );

//------------------------------------------------------------------------------

        real
        lambda_dep( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        real
        mu_dep( const real & aMu, const real & aT, const real & aP  );

//------------------------------------------------------------------------------

        // link all caloric and transport functions to splines
        void
        link_to_idgas_property_functions();

//------------------------------------------------------------------------------

        // link all caloric and transport functions to splines
        // plus departure functions of underlying gas
        void
        link_to_realgas_property_functions();

//------------------------------------------------------------------------------

         // use the property functions from the equation of state
         void
         link_to_helmholtz_property_functions();

//------------------------------------------------------------------------------

        // create the table needed for formation enthalpy
        void
        create_formation_table();

//------------------------------------------------------------------------------

        // special subroutine needed for shock
        real
        shock_beta_simple( const real & aT1,
                    const real & aP1,
                    const real & aU1,
                    const real & aAlpha,
                          real & aT2,
                          real & aP2,
                          real & aU2,
                    const real & aBeta );

//------------------------------------------------------------------------------

         // special subroutine needed for shock
         real
         shock_beta( const real & aT1,
                     const real & aP1,
                     const real & aU1,
                     const real & aAlpha,
                     real & aT2,
                     real & aP2,
                     real & aU2,
                     const real & aBeta );

//------------------------------------------------------------------------------

         // just return zero for idgas
        real
        dhdp_idgas( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        // workaround if eos does that via departure function
        real
        dhdp_eos_departure( const real & aT, const real & aP );

//------------------------------------------------------------------------------

        // workaround if eos dies not support it
        real
        dhdp_differential_quotient( const real & aT, const real & aP );

//------------------------------------------------------------------------------

         void
         update_mixture_entropy();

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------

    inline bool
    Gas::is_idgas() const
    {
        return mGasModel == GasModel::IDGAS ;
    }

//------------------------------------------------------------------------------

    inline const Vector< real > &
    Gas::molar_fractions() const
    {
        return mMolarFractions;
    }

//------------------------------------------------------------------------------

    inline const Vector< real > &
    Gas::mass_fractions() const
    {
        return mMassFractions;
    }

//------------------------------------------------------------------------------


    inline const real &
    Gas::molar_fraction( const uint aIndex ) const
    {
        return mMolarFractions( aIndex );
    }

//------------------------------------------------------------------------------

    inline const real &
    Gas::mass_fraction( const uint aIndex ) const
    {
        return mMassFractions( aIndex );
    }

//------------------------------------------------------------------------------

    inline Cell< gastables::RefGas * > &
    Gas::elements()
    {
        return mElements ;
    }

//------------------------------------------------------------------------------

    inline Cell< gastables::RefGas * > &
    Gas::components()
    {
        return mComponents ;
    }

//------------------------------------------------------------------------------

    inline gastables::RefGas *&
    Gas::component( const index_t aIndex )
    {
        return mComponents( aIndex );
    }

//------------------------------------------------------------------------------

    inline const Matrix< real > &
    Gas::formation_table() const
    {
        return mFormationTable ;
    }

//------------------------------------------------------------------------------

    inline gasmodels::EoS *
    Gas::eos()
    {
        return mEoS;
    }

//------------------------------------------------------------------------------

    inline const uint &
    Gas::number_of_components() const
    {
        return mNumberOfComponents;
    }

//------------------------------------------------------------------------------

    inline gasmodels::Statevals &
    Gas::statevals()
    {
        return mStatevals;
    }

//------------------------------------------------------------------------------

    inline bool
    Gas::is_liquid() const
    {
        return mLiquidFlag;
    }

//------------------------------------------------------------------------------

    inline void
    Gas::set_liquid_flag( const bool aFlag )
    {
        mLiquidFlag = aFlag ;
    }

//------------------------------------------------------------------------------

    inline real
    Gas::alpha( const real & aT, const real & aP )
    {
        return mEoS->alpha( aT, aP );
    }

//------------------------------------------------------------------------------

    inline real
    Gas::beta( const real & aT, const real & aP )
    {
        return mEoS->beta( aT, aP );
    }

//------------------------------------------------------------------------------

    inline real
    Gas::kappa( const real & aT, const real & aP )
    {
        return mEoS->kappa( aT, aP );
    }

//------------------------------------------------------------------------------

    inline Spline &
    Gas::heat_spline()
    {
        return mHeatSpline ;
    }

//------------------------------------------------------------------------------

    inline Spline &
    Gas::viscosity_spline()
    {
        return mViscositySpline ;
    }

//------------------------------------------------------------------------------

    inline Spline &
    Gas::conductivity_spline()
    {
        return mConductivitySpline;
    }

//------------------------------------------------------------------------------

    inline const GasModel &
    Gas::gas_model() const
    {
        return mGasModel ;
    }

//------------------------------------------------------------------------------

    inline const HelmholtzModel &
    Gas::helmholtz_model() const
    {
        return mHelmholzModel ;
    }

//------------------------------------------------------------------------------

    inline real
    Gas::dhdp_idgas( const real & aT, const real & aP )
    {
        return 0 ;
    }

//------------------------------------------------------------------------------

    // workaround if eos does that via departure function
    inline real
    Gas::dhdp_eos_departure( const real & aT, const real & aP )
    {
        return mEoS->dhdepdp( aT, aP );
    }

//------------------------------------------------------------------------------

    // workaround if eos dies not support it
    inline real
    Gas::dhdp_differential_quotient( const real & aT, const real & aP )
    {
        return ( mEoS->h( aT, 1.01 * aP ) - mEoS->h( aT, 0.99 * aP ) ) /
                ( 0.02 * aP ) ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_GAS_HPP
