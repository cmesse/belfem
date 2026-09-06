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
     *
     * @ingroup grp_physics_gasmodels
     * @see @ref physics_gasmodels_gasmodels_usage_guide
       */
    /**
     * @par Const correctness
     *
     * The accessors that evaluate a property of a fixed mixture -- cp( T, p ),
     * h( T, p ), mu( T, p ), the flow routines, ... -- are const. They do not
     * change the identity of the gas: after the call it holds the same species
     * in the same proportions, and a second call with the same arguments
     * returns the same number. What they do write is the memoization cache and
     * the preallocated scratch, and those members carry `mutable` for exactly
     * that reason. The methods that change what the gas *is* -- remix(),
     * remix_mass(), the equilibrium solvers, the spline builders --
     * are the ones that stay non-const, so the signature alone separates the
     * two groups.
     *
     * Only genuine scratch is `mutable`. A work vector belonging to a
     * composition-changing path (mWorkMu, the RAND equilibrium block,
     * mFormationTable) deliberately is not, so the compiler keeps proving that
     * those paths stay out of the const ones.
     *
     * @warning `const` here means logically const, not thread safe. Because
     * the const evaluators write the shared cache, two threads must not call
     * them on the same Gas object even through a const reference. This matches
     * the framework-wide policy: BELFEM parallelises with MPI and is not
     * internally thread safe ( see doc/coding_philosophy.md ).
     */
    class Gas
    {
    protected:

        //! size of components vector
        uint mNumberOfComponents;

        //! container for state variables
        //! logically const memo of the state last asked for, hence mutable
        mutable gasmodels::Statevals mStatevals;

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


        Spline mHeatSpline;
        Spline mViscositySpline;
        Spline mConductivitySpline;

        mutable real    mLastSplineT   = BELFEM_REAL_MAX ;
        mutable index_t mLastSplineCol = BELFEM_UINT_MAX ;

        //! gasmodel type
        GasModel mGasModel = GasModel::UNDEFINED ;

        //! helmholz type, if used
        HelmholtzModel mHelmholzModel = HelmholtzModel::UNDEFINED ;

        //! Interaction polynomials for viscosity
        Cell<gastables::RefGas *> mViscosityInteractionRefgas;

        //! table telling if interaction parameter exists
        Matrix<uint> mViscosityInteractionTable;

        // Work matrices for viscosity and conductivity calculation of
        // mixtures. mutable: written by the const cea_mu / cea_lambda
        mutable Matrix<real> mWorkMatrix;
        mutable Vector<real> mWorkVector;
        mutable Vector<real> mWorkVector2;
        Vector<real> mWorkMu;
        Vector<real> mWorkLambda;

        // work matrix and vectors for equilibrium. NOT mutable: the RAND
        // solver changes the composition, so it is not a const path
        Vector<real> mWorkVectorRAND0 ;
        Vector<real> mWorkVectorRAND1 ;
        Vector<real> mWorkVectorRAND2 ;
        Matrix<real> mWorkMatrixRAND ;
        Vector< int_t >  mPivotRAND ;

        mutable real mWorkTemperature;

        // Work matrices for gibbs
        Matrix<real> mFormationTable;
        mutable Vector<real> mFormationWork;

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

        //! liquid flag, written by Helmholtz::v() and read by EoS_Cubic::v()
        bool mLiquidFlag = false;

        //! the equation of state
        gasmodels::EoS    * mEoS = nullptr;

        //! special class, only needed if this is a helmholtz eos
        gasmodels::HelmholtzTransport * mTransport = nullptr ;

        //! pointer to cp function
        real
        ( Gas::*mFunctionCp )
             ( const real T, const real p ) const ;

        //! pointer to dcpdT function
        real
        ( Gas::*mFunctiondCpdT )
                ( const real T, const real p ) const ;

        //! pointer to cv function
        real
        ( Gas::*mFunctionCv )
                ( const real T, const real p ) const ;

        //! pointer to gamma function
        real
        ( Gas::*mFunctionGamma )
                ( const real T, const real p ) const ;

        //! pointer to c function
        real
        ( Gas::*mFunctionC )
                ( const real T, const real p ) const ;

        //! pointer to h function
        real
        ( Gas::*mFunctionH )
                ( const real T, const real p ) const ;

        //! pointer to dhdp function
        real
        ( Gas::*mFunctionDHDP )
                ( const real T, const real p ) const ;

        //! pointer to s function
        real
        ( Gas::*mFunctionS )
                ( const real T, const real p ) const ;

        //! pointer to dsdT function
        real
        ( Gas::*mFunctionDSDT )
                ( const real T, const real p ) const ;

        //! pointer to dsdp function
        real
        ( Gas::*mFunctionDSDP )
                ( const real T, const real p ) const ;

        //! pointer to Mu function
        real
        ( Gas::*mFunctionMU )
                ( const real T, const real p ) const ;

        //! pointer to Lambda function
        real
        ( Gas::*mFunctionLAMBDA )
                ( const real T, const real p ) const ;

        // support term for entropy
        real mMixtureEntropy = BELFEM_QUIET_NAN ;

        //! scratch for the duct solvers, sized in initialize()
        mutable Matrix< real >  mFlowJacobian ;
        mutable Vector< real >  mFlowResidual ;
        mutable Vector< int_t > mFlowPivot ;

        //! Gauss rule for prandtl_meyer, allocated on first call
        mutable Vector< real >  mGaussPoints ;
        mutable Vector< real >  mGaussWeights ;

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

        // owns the component RefGas objects and the EoS
        Gas( const Gas & ) = delete;
        Gas & operator=( const Gas & ) = delete;

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

        /**
         * @name Mixture molar mass and specific gas constant
         *
         * Returned by value. Both used to hand out a reference into the state
         * cache, which made the result a live view that followed a later
         * remix; a caller wanting that had to know it, and a caller not
         * wanting it had no way to tell. A real is a register return, so the
         * copy costs nothing and the accessor now means what it says.
         *
         * @note Both currently ignore T and p; the parameters are there for a
         * derived model whose mixture dissociates, where M really does depend
         * on the state.
         * @{
         */
        virtual real
        M( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        virtual real
        R( const real T, const real p ) const ;

        /** @} */

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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * expose the data object of a refgas ( const version )
         */
         const gastables::GasData *
         data( const index_t aIndex ) const;

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
         * test if the last evaluated state is liquid. Set by the Helmholtz
         * EoS; stays false for ideal and cubic gases.
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
        p( const real T, const real v ) const ;

//------------------------------------------------------------------------------

        real
        v( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        real
        rho( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        real
        T( const real p, const real v ) const ;

//------------------------------------------------------------------------------
// Caloric Properties
//------------------------------------------------------------------------------

        virtual real
        cp( const real T, const real p ) const ;

        virtual real
        cv( const real T, const real p ) const ;

        virtual real
        gamma( const real T, const real p ) const ;

        virtual real
        c( const real T, const real p ) const ;

        virtual real
        u( const real T, const real p ) const ;

        virtual real
        h( const real T, const real p ) const ;

        virtual real
        s( const real T, const real p ) const ;

        virtual real
        dsdT( const real T, const real p ) const ;

        virtual real
        dsdp( const real T, const real p ) const ;

        virtual real
        dcpdT( const real T, const real p ) const ;

        // dissociation enthalpy ( only for tablegas at this time )
        virtual real
        hd( const real T, const real p ) const ;

//------------------------------------------------------------------------------
// Transport Properties
//------------------------------------------------------------------------------

        /**
         * dynamic viscosity in Pa*s
         */
        real
        mu( const real T, const real p ) const ;

        /**
          * thermal conductivity in W/(m*K)
          */
        real
        lambda( const real T, const real p ) const ;

        /**
         * Prandtl Number
         */
        real
        Pr( const real T, const real p ) const ;

//------------------------------------------------------------------------------
// Thermodynamic Coefficients
//------------------------------------------------------------------------------

         /**
          * thermal expansion coefficient
          *
          * \f$ \alpha = \frac{1}{v} \left( \frac{\partial v}{\partial T}\right)_p \f$
          */
        real
        alpha( const real T, const real p ) const ;

//------------------------------------------------------------------------------
         /**
          * isochoric stress coefficient
          *
          * \f$ \beta = \frac{1}{p} \left( \frac{\partial p}{\partial T}\right)_v \f$
          */
        real
        beta( const real T, const real p ) const ;

//------------------------------------------------------------------------------

         /**
          * isothermal compressibility coefficient
          *
          * \f$ \kappa = -\frac{1}{v} \left( \frac{\partial v}{\partial p}\right)_T \f$
          */
        real
        kappa( const real T, const real p ) const ;

//------------------------------------------------------------------------------
// CHEMISTRY
// -----------------------------------------------------------------------------

        /**
         * return the molar Gibbs potential at reference pressure
         */
         void
         Gibbs( const real T, Vector< real > & aGibbs ) const ;

         /**
         * return the molar formation enthalpy for each component
         */
         void
         Hf( const real T, Vector< real > & aHf ) const ;

         /**
          * return the temperature derivative of  gibs potentia at reference pressure
          */
         void
         dGibbsdT( const real T, Vector< real > & aGibbs ) const ;

         /**
          * remix to equilibrium
          */
         void
         remix_to_equilibrium(
                 const real T,
                 const real p,
                 const bool aRemixHeat=true,
                 const bool aRemixTransport=true );

         void
         compute_equilibrium( const real T, const real p, Vector< real > & aX );

//------------------------------------------------------------------------------
// State relevant methods
// -----------------------------------------------------------------------------

        real
        T_from_h( const real & h, const real p ) const ;

// -----------------------------------------------------------------------------
        /**
         * get an isentropic temperature
         */
         real
         isen_T( const real T0, const real p0, const real p1 ) const ;

// -----------------------------------------------------------------------------

         /**
           * get an isentropic pressure
           */
         real
         isen_p( const real T0, const real p0, const real T1 ) const ;

// -----------------------------------------------------------------------------

        /**
         * calculate the total state
         */
         void
         total( const real T, const real p, const real & u,
                      real & aTt, real & aPt ) const ;

// -----------------------------------------------------------------------------

        /**
         * expansion of a gas into a widening duct, A2 >= A1
         *
         * a subsonic flow separates at the step and loses total pressure over a
         * Borda-Carnot shock, so mass, momentum and energy are solved. a
         * supersonic flow expands around the corner instead, the step face does
         * work on it and momentum is not a control volume invariant, so mass,
         * entropy and energy are solved
         */
         void
         expand( const real & A1,
                 const real T1,
                 const real p1,
                 const real & u1,
                 const real & A2,
                       real & T2,
                       real & p2,
                       real & u2 ) const ;

// -----------------------------------------------------------------------------

        /**
         * compression of a gas into a narrowing duct, A2 <= A1
         *
         * a contraction does not separate, so the flow stays isentropic on
         * either branch and only the root of the area relation differs. errors
         * out if the duct is choked, that is if A2 undercuts the sonic area
         */
         void
         compress( const real & A1,
                   const real T1,
                   const real p1,
                   const real & u1,
                   const real & A2,
                         real & T2,
                         real & p2,
                         real & u2 ) const ;

// -----------------------------------------------------------------------------

        /**
         * Prandtl-Meyer turn of a supersonic stream around a corner, for a
         * thermally perfect ideal gas.
         *
         * integrates the exact simple wave relation
         * d(nu) = sqrt( Ma^2 - 1 ) * dV / V along the isentrope of the
         * upstream state, then recovers the pressure from ds = 0. both
         * laws of thermodynamics hold
         * identically, because isentropy is built into the characteristic
         * derivation instead of being enforced as a constraint afterwards.
         *
         * a negative angle models a smooth isentropic compression and is
         * admissible while the flow stays supersonic. coalescing
         * characteristics, that is an embedded shock on a concave wall, are
         * not detected by any isentropic method.
         *
         * @param T1     upstream temperature in K
         * @param p1     upstream pressure in Pa
         * @param u1     upstream velocity in m/s, must be supersonic
         * @param alpha  turning angle in rad, positive for an expansion
         * @param T2     downstream temperature in K, written
         * @param p2     downstream pressure in Pa, written
         * @param u2     downstream velocity in m/s, written
         * @return       the Mach number reached downstream
         */
         real
         prandtl_meyer(  const real T1,
                         const real p1,
                         const real & u1,
                         const real & alpha,
                               real & T2,
                               real & p2,
                               real & u2 ) const ;

// -----------------------------------------------------------------------------

        /**
         * perpendicular shock
         *
         */
         void
         shock(  const real T1, const real p1, const real & u1,
                       real & T2,       real & p2,       real & u2 ) const ;

// -----------------------------------------------------------------------------

        /**
         * oblique shock
         */
         void
         shock( const real T1, const real p1, const real & u1, const real & alpha,
                real & T2, real & p2, real & u2, real & beta ) const ;

//------------------------------------------------------------------------------
// Component Volume and Departure Functions
//------------------------------------------------------------------------------

         /**
          * @name Per-component properties
          *
          * Properties of one mixture component, not of the mixture as a whole.
          *
          * @warning At the reference pressure, the caloric properties h, cp and
          * dcpdT return the ideal gas value, not the component's real gas property.
          *
          * The indexed departure functions in the cubic equation of state include the
          * reference pressure subtraction, so these functions use the same convention
          * as the mixture level real gas functions. See the note on Gas::realgas_cp in
          * cl_Gas.cpp. The two levels agree and must change together if BELFEM changes
          * this convention.
          *
          * The finite rate combustion solver uses these functions. It takes its
          * chemistry from the ideal gas Gibbs energy of the thermo tables. Holding
          * component enthalpies at the ideal gas value at the reference pressure keeps
          * the energy equation on the same footing as the equilibrium constants. For
          * combustion at moderate pressure, the ideal gas assumption is a reasonable
          * engineering approximation.
          *
          * @param aIndex component index
          * @param T     temperature in K
          * @param p     pressure in Pa
          * @{
          */

         real
         v( const uint aIndex, const real T, const real p ) const ;

         //! specific enthalpy of the component on the RefGas scale
         //! ( H( 0 K ) = formation enthalpy at 298.15 K, see cl_GT_RefGas.hpp )
         //! plus the departure
         real
         h( const uint aIndex, const real T, const real p ) const ;

         real
         cp( const uint aIndex, const real T, const real p ) const ;

         real
         dcpdT( const uint aIndex, const real T, const real p ) const ;

         /** @} */

//------------------------------------------------------------------------------
// Print composition
//------------------------------------------------------------------------------

        void
        print() const ;

//------------------------------------------------------------------------------
// Special access
//------------------------------------------------------------------------------

         /**
          * expose heat spline
          * @return
          */
         Spline &
         heat_spline() ;

         const Spline &
         heat_spline() const ;

         /**
          * expose viscosity spline
          * @return
          */
         Spline &
         viscosity_spline() ;

         const Spline &
         viscosity_spline() const ;

         /**
         * expose conductivity spline
         * @return
         */
         Spline &
         conductivity_spline() ;

         const Spline &
         conductivity_spline() const ;

         /**
          * returns what gas model is used
          */
         const GasModel &
         gas_model() const;

         /**
          * which Helmholtz model this gas uses, if any
          */
         const HelmholtzModel &
         helmholtz_model() const;

//------------------------------------------------------------------------------

         // enthalpy derivative to pressure ( needed for total temperature )
         virtual real
         dhdp( const real T, const real p ) const ;

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
        element_to_molecule( const string & aElement ) const;

//------------------------------------------------------------------------------

        string
        reference_element( const string & aElement ) const ;

//------------------------------------------------------------------------------

        void
        create_viscosity_table( gastables::RefGasFactory & aFactory );

//------------------------------------------------------------------------------

        void
        create_mass_properties( const Vector<real> & aMolarFractions );

//------------------------------------------------------------------------------

        void
        evaluate_viscosity_interaction( const real T ) const ;

//------------------------------------------------------------------------------

        void
        evaluate_conductivity_interaction( const real T ) const ;

//------------------------------------------------------------------------------

        real
        idgas_cp( const real T, const real p ) const ;

        real
        idgas_dcpdT( const real T, const real p ) const ;

        real
        idgas_cv( const real T, const real p ) const ;

        real
        idgas_gamma( const real T, const real p ) const ;

        real
        idgas_c( const real T, const real p ) const ;

        real
        idgas_h( const real T, const real p ) const ;

        real
        idgas_s( const real T, const real p ) const ;

        real
        idgas_dsdT( const real T, const real p ) const ;

        real
        idgas_dsdp( const real T, const real p ) const ;

        real
        idgas_mu( const real T, const real p ) const ;

        real
        idgas_lambda( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        real
        realgas_cp( const real T, const real p ) const ;

        real
        realgas_dcpdT( const real T, const real p ) const ;

        real
        realgas_cv( const real T, const real p ) const ;

        real
        realgas_gamma( const real T, const real p ) const ;

        real
        realgas_c( const real T, const real p ) const ;

        real
        realgas_h( const real T, const real p ) const ;

        real
        realgas_s( const real T, const real p ) const ;

        real
        realgas_dsdT( const real T, const real p ) const ;

        real
        realgas_dsdp( const real T, const real p ) const ;

        real
        realgas_mu( const real T, const real p ) const ;

        real
        realgas_lambda( const real T, const real p ) const ;

//------------------------------------------------------------------------------

         real
         helmholtz_cp( const real T, const real p ) const ;

         real
         helmholtz_dcpdT( const real T, const real p ) const ;

         real
         helmholtz_cv( const real T, const real p ) const ;

         real
         helmholtz_gamma( const real T, const real p ) const ;

         real
         helmholtz_c( const real T, const real p ) const ;

         real
         helmholtz_h( const real T, const real p ) const ;

         real
         helmholtz_s( const real T, const real p ) const ;

         real
         helmholtz_dsdT( const real T, const real p ) const ;

         real
         helmholtz_dsdp( const real T, const real p ) const ;

         real
         helmholtz_mu( const real T, const real p ) const ;

         real
         helmholtz_lambda( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        // from cea, Eq. 5.3
        real
        cea_mu( const real T ) const ;

//------------------------------------------------------------------------------

        // from cea, Eq. 5.4
        real
        cea_lambda( const real T ) const ;

//------------------------------------------------------------------------------

        real
        lambda_dep( const real T, const real p  ) const ;

//------------------------------------------------------------------------------

        real
        mu_dep( const real & mu, const real T, const real p  ) const ;

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

        /**
         * closed form Prandtl-Meyer angle of a calorically perfect gas,
         * evaluated with the local gamma( T, p ). for a thermally perfect
         * gas this is not a state function and must not be differenced
         * between two states, it only serves the initial guess of
         * prandtl_meyer
         */
        real
        prandtl_meyer_angle(
                const real T,
                const real p,
                const real & u ) const ;

//------------------------------------------------------------------------------

        // special subroutine needed for shock
        real
        shock_beta_simple( const real T1,
                    const real p1,
                    const real & u1,
                    const real & alpha,
                          real & T2,
                          real & p2,
                          real & u2,
                    const real & beta ) const ;

//------------------------------------------------------------------------------

         // special subroutine needed for shock
         real
         shock_beta( const real T1,
                     const real p1,
                     const real & u1,
                     const real & alpha,
                     real & T2,
                     real & p2,
                     real & u2,
                     const real & beta ) const ;

//------------------------------------------------------------------------------

         // just return zero for idgas
        real
        dhdp_idgas( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        // workaround if eos does that via departure function
        real
        dhdp_eos_departure( const real T, const real p ) const ;

//------------------------------------------------------------------------------

        // workaround if eos dies not support it
        real
        dhdp_differential_quotient( const real T, const real p ) const ;

//------------------------------------------------------------------------------

         void
         update_mixture_entropy();

//------------------------------------------------------------------------------

         /**
          * isentropic area relation A/A* of a frozen gamma ideal gas
          */
         real
         area_mach( const real Ma, const real k ) const ;

//------------------------------------------------------------------------------

         /**
          * starting values for a duct, bisected out of the isentropic area
          * relation with a frozen gamma. the branch is selected by the caller,
          * since both roots of the relation are physical
          */
         void
         area_guess( const real   T1,
                     const real   p1,
                     const real   Ma1,
                     const real & A1,
                     const real & A2,
                     const bool   aSupersonic,
                           real & T2,
                           real & p2,
                           real & u2 ) const ;

//------------------------------------------------------------------------------

         /**
          * solves mass, entropy and energy for a duct in which the flow stays
          * isentropic. T2, p2 and u2 enter as the initial guess
          */
         void
         isentropic_duct( const real & aMass,
                          const real & aEntropy,
                          const real & aEnergy,
                          const real & A2,
                                real & T2,
                                real & p2,
                                real & u2 ) const ;

//------------------------------------------------------------------------------

        index_t
        spline_col( const real T ) const ;

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
    Gas::alpha( const real T, const real p ) const
    {
        return mEoS->alpha( T, p );
    }

//------------------------------------------------------------------------------

    inline real
    Gas::beta( const real T, const real p ) const
    {
        return mEoS->beta( T, p );
    }

//------------------------------------------------------------------------------

    inline real
    Gas::kappa( const real T, const real p ) const
    {
        return mEoS->kappa( T, p );
    }

//------------------------------------------------------------------------------

    inline Spline &
    Gas::heat_spline()
    {
        return mHeatSpline ;
    }

//------------------------------------------------------------------------------

    inline const Spline &
    Gas::heat_spline() const
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

    inline const Spline &
    Gas::viscosity_spline() const
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

    inline const Spline &
    Gas::conductivity_spline() const
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
    Gas::dhdp_idgas( const real T, const real p ) const
    {
        return 0 ;
    }

//------------------------------------------------------------------------------

    // workaround if eos does that via departure function
    inline real
    Gas::dhdp_eos_departure( const real T, const real p ) const
    {
        return mEoS->dhdepdp( T, p );
    }

//------------------------------------------------------------------------------

    // workaround if eos does not support it
    inline real
    Gas::dhdp_differential_quotient( const real T, const real p ) const
    {
        return ( mEoS->h( T, 1.01 * p ) - mEoS->h( T, 0.99 * p ) ) /
                ( 0.02 * p ) ;
    }

//------------------------------------------------------------------------------

    /**
     * Cached column lookup for the heat, viscosity and conductivity splines.
     *
     * The lookup uses the heat spline column and reuses it for viscosity and
     * conductivity. This is valid because all three splines use the same grid; see
     * the constructor initializer list. The assertion below guards that invariant,
     * because a mismatched grid would return values from the wrong interval.
     *
     * The cache cannot go stale. The column depends only on temperature and grid,
     * and remix_heat rewrites the spline coefficients without changing the grid.
     */
    inline index_t
    Gas::spline_col( const real T ) const
    {
        BELFEM_ASSERT(    mViscositySpline.x_min() == mHeatSpline.x_min()
                       && mViscositySpline.x_max() == mHeatSpline.x_max()
                       && mConductivitySpline.x_min() == mHeatSpline.x_min()
                       && mConductivitySpline.x_max() == mHeatSpline.x_max(),
            "Gas splines do not share a grid, the cached column cannot be reused" );

        if ( T != mLastSplineT )
        {
            mLastSplineT = T ;
            mLastSplineCol = mHeatSpline.find_col( T );
        }
        return mLastSplineCol ;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_GAS_HPP
