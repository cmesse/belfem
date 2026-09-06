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

#ifndef BELFEM_CL_MATERIAL_HPP
#define BELFEM_CL_MATERIAL_HPP
#include "globals.hpp"

#include <cmath>

#include "typedefs.hpp"
#include "constants.hpp"
#include "assert.hpp"
#include "cl_Bitset.hpp"
#include "cl_Cell.hpp"

// note: NEVER include cl_Vector.hpp or cl_Matrix.hpp here, or any class that uses it.
//       Doing so would break the API for the user defined materials.

namespace belfem
{
    /**
     * @brief Upper bound for the thermal expansion split temperature [K]
     *
     * The split is 0.618 * theta_Debye, capped here. The Debye factor keeps the
     * anchor above the steep part of cp, where dln(C)/dT would be negative; the
     * cap keeps it inside the range where expansion measurements still carry
     * signal, and limits how far a single anchor has to extrapolate. For a
     * high-theta material the cap binds, for a low-theta one the Debye factor does.
     *
     * The value must not undercut the temperature at which dln(C)/dT changes
     * sign, which sits at roughly 0.43 - 0.51 * theta for the fitted metals -
     * 0.618 * theta clears it by construction, but a flat cap breaks that
     * relationship once theta exceeds about 425 K ( aluminum, chromium ). The
     * ice point clears the sign change for every material in the roster, and
     * the choice of anchor inside the admissible window is uncritical: moving
     * it shifts the cryogenic alpha by less than the scatter of the underlying
     * expansion data ( under 5 % for copper, aluminum, iron and nickel ).
     */
    constexpr real gTAlphaSwitchMax = 273.15 ;
}

//#include "cl_Mesh.hpp"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#pragma clang diagnostic ignored "-Wformat"
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"

#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#elif defined(BELFEM_INTEL)
    // Intel compiler diagnostics
#pragma warning(push)
    // disable format warnings
#pragma warning(disable: 1011)
    // disable unused variable warnings
#pragma warning(disable: 177)
    // keep your original one
#pragma warning(disable: 1595)

#endif



// Forward declaration (full definition only needed in implementation files)

#include "cl_JcFunction.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    namespace material
    {
        class BhCurve ;
        class SplineLookupTable ;
    }

//------------------------------------------------------------------------------

    /**
     * @brief Material classification based on physical behavior
     *
     * Material types determine which properties are available and their dependencies:
     * - UserDefined: Custom material with user-specified properties
     * - Ferro: Ferromagnetic materials where permeability depends on H (T dependence planned)
     * - HTS: High-temperature superconductors where rho and lambda depend on j, T, n×B and n·B
     * - PureMetal: Noble metals where rho and lambda depend on j, T, and B
     * - LookupAlloy: Metallic alloys where rho and lambda depend only on T (tabulated)
     * - CompositeAlloy: Composite materials with orthotropic rho and lambda that depend on T
     * - NonMetal: Non-metallic materials (e.g. Magnesia)
     */
    enum class MaterialType
    {
        UserDefined,
        Ferro,
        HTS,
        PureMetal,
        LookupAlloy,
        CompositeAlloy,
        NonMetal
    };

    /**
     * @brief Dependencies that material properties can have
     */
    enum class MaterialDependency
    {
        T         = 0,    // temperature [K]
        normB     = 1,    // magnitude of in-plane magnetic field [T] (HTS, for lambda and rho)
        angleBxJ  = 2,    // angle between B and j [rad]
        angleNxB  = 3,    // angle between surface normal and B [rad]
        normH     = 4,    // magnetic field strength [A/m] (for mu in ferromagnetic materials)
        normJ     = 5,
        Jc        = 6,
        rho       = 7,
        UNDEFINED = 8
    };

    /**
     * @brief Bitset for tracking material property dependencies
     */
    typedef Bitset< static_cast< size_t > ( MaterialDependency::UNDEFINED ) > MaterialDependencyBitset ;

    /**
     * @brief All material properties that can be defined
     *
     * Properties are organized by category:
     * - Physical: density, M (molar mass), q (atoms per molecule)
     * - Mechanical: E, nu, alpha, Rp02
     * - Thermal: cp, lambda, debye, gamma, beta
     * - Electric: rho, rho_i, rho_0, kohler_trans, kohler_long, RRR
     * - Magnetic: mu
     * - Superconductor: jc, n, ec, T_crit, layer_thickness
     * - Reference values: ref_density, T_ref_density, rho_i_ref, T_ref_rho_i, debye0K
     * - Gas: R (specific gas constant), Gamma (impurity parameter)
     * - Limits: T_max
     */
    enum class MaterialProperty
    {
        density             =  0,    // density [kg/m³]
        E                   =  1,    // Young's modulus [Pa]
        nu                  =  2,    // Poisson's ratio [-]
        cp                  =  3,    // specific heat capacity [J/(kg·K)]
        lambda              =  4,    // thermal conductivity [W/(m·K)]
        mu                  =  5,    // magnetic permeability [H/m], defined as ∂B/∂H
        rho                 =  6,    // electric resistivity [Ω·m]
        alpha               =  7,    // thermal expansion coefficient [1/K], defined as (1/l)·∂l/∂T, not (1/l)·Δl/ΔT!
        Rp02                =  8,    // yield stress [Pa]
        debye               =  9,    // Debye temperature [K]
        rho_i               = 10,    // inner resistivity of noble metal [Ω·m]
        kohler_trans        = 11,    // Kohler parameter for transverse magnetoresistance [-]
        kohler_long         = 12,    // Kohler parameter for longitudinal magnetoresistance [-]
        T_crit              = 13,    // critical temperature [K]
        M                   = 14,    // molar mass [kg/mol]
        Gamma               = 15,    // impurity parameter [-]
        R                   = 16,    // specific gas constant [J/(kg·K)]
        T_max               = 17,    // maximum temperature [K]
        ref_density         = 18,    // reference density [kg/m³]
        T_ref_density       = 19,    // temperature at reference density [K]
        gamma               = 20,    // linear Debye parameter [J/(kg·K²)], cv = γ·T + β·T³ = ∂cp/∂T at T=0
        beta                = 21,    // cubic Debye parameter [J/(kg·K⁴)], cv = γ·T + β·T³
        q                   = 22,    // number of atoms per molecule [-] (default: 1)
        debye0K             = 23,    // Debye temperature at 0 K [K]
        rho_i_ref           = 24,    // reference inner resistivity [Ω·m]
        T_ref_rho_i         = 25,    // temperature at reference inner resistivity [K]
        A_bloch_gruen       = 26,    // A-parameter for Bloch-Grüneisen law [-]
        n_bloch_gruen       = 27,    // n-parameter for Bloch-Grüneisen law [-]
        RRR                 = 28,    // residual resistivity ratio [-]
        rho_0               = 29,    // residual resistivity [Ω·m]
        layer_thickness     = 30,    // characteristic length, e.g., layer thickness [m]
        ec                  = 31,    // critical electric field [V/m]
        jc                  = 32,    // critical current density [A/m²]
        n                   = 33,    // exponent for power law [-]
        Tcurie              = 34,    // Curie temperature [K]
        A_electron_magnon   = 35,    // electron-magnon scattering parameter [ Ω · m / K² ]
        A_spin_disorder     = 36,    // spin disorder amplitude [ Ω · m ]
        grueneisen          = 37,    // Grüneisen parameter
        density_correction  = 38,    // density scaling factor, to correct solder thickness [-]
        rho0_pure           = 39,    // resistivity from TPCR for pure metals
        UNDEFINED           = 40
    };

    /**
     * @brief Number of non-constant material properties
     *
     * Does not include jc and n because they are handled differently
     */
    constexpr size_t gNumNonConstantMaterialProperties = static_cast< size_t >( MaterialProperty::T_crit ) ;

    /**
     * @brief Total number of material properties
     */
    constexpr size_t gNumMaterialProperties = static_cast< size_t >( MaterialProperty::UNDEFINED ) ;

    /**
     * @brief Convert material property enum to string
     */
    string
    to_string( const MaterialProperty aProperty ) ;

    /**
     * @brief Get the SI unit for a material property
     */
    unit
    get_unit( const MaterialProperty aProperty ) ;

    /**
     * @brief E-J law used for a superconductor's resistivity
     *
     * PowerLaw  — power-law channel parallel to the normal-state channel
     *             ( Duron et al. 2004 ); the historical default.
     * Piecewise — three-regime law: raw power law, Bézier flux-flow blend,
     *             normal state ( Rhyner 1993 exponent, log-log blend ).
     * Riva      — same parallel model as PowerLaw, but total over the full
     *             jc/n table range ( Riva 2021, EPFL thesis 8754 ): guarded
     *             against dead defects, under-/overflow and n → 1.
     */
    enum class ResistivityLaw
    {
        PowerLaw  = 0,
        Piecewise = 1,
        Riva      = 2
    };

    class Material ;

    typedef real ( MatFunc1 )( const Material * , const real ) ; // T
    typedef real ( MatFunc2 )( const Material *, const real, const real ) ;
    typedef real ( MatFunc3 )( const Material *, const real, const real, const real ) ;
    typedef real ( DefectFunc )( const real, const real, const real, const real ) ;
    typedef real ( HeatFunc )( const real, const real, const real, const real ) ;

    /**
     * @brief Base class for all materials in BELFEM
     *
     * The Material class provides a unified interface for accessing physical properties
     * of materials. Properties can be constant, temperature-dependent, or depend on
     * multiple variables (field, current, etc.) depending on the material type.
     *
     * Key features:
     * - Automatic property function dispatching based on dependencies
     * - Support for both constant and variable properties
     * - Spline-based interpolation for tabulated data
     * - Custom functions for complex material behavior
     * - Ownership management for B-H curves and Jc functions
     *
     * Usage:
     * @code
     * MaterialFactory factory;
     * Material* cu = factory.create_material("copper");
     * cu->set_RRR(100);
     * real rho = cu->rho(77.0);  // Resistivity at 77K
     * @endcode
     *
     * IMPORTANT NOTES:
     * - For HTS materials, use rho_powerlaw() instead of rho()
     * - For thermal M-matrix, use ref_density() or density(gTroom), not density(T)
     * - Material owns and deletes assigned B-H curves and Jc/n functions
     * - Variable names in material functions are exempt from BELFEM naming rules
     *
     * @ingroup grp_physics_materials
     * @see @ref physics_materials_materials_usage_guide
     */
    class Material
    {
        const proc_t mCommRank ;

        const MaterialType mType ;

        const bool mIsIsotropic ;  // Flag for isotropic properties (default: true)

        string mLabel ;   // Material label (e.g., "Copper", "YBCO")

        string mNumber = "" ;  // Material number/identifier if it exists

        // bitwise flags for this vertex
        // 1: used by kernel
        // 7: reserved for factory
        uint8_t mFlags = 0;

        Bitset< gNumMaterialProperties > mHaveProperty ;  // Flags which properties are implemented

        Cell< real > mConstantProperties ;  // Storage for constant property values

        Cell< MaterialDependencyBitset * > mPropertyDependencies ;  // Tracks dependencies for each property


        //Piecewise resistivity parameters
        real mNff = 3.0 ;
        real mD = 2.5 ;

        //! numerical floor on the intrinsic power-law resistivity. ZERO by
        real mRhoMin = 0.0 ;

        //! temperature below which alpha is taken from the Grueneisen branch,
        //! min( 0.618 * theta_Debye, gTAlphaSwitchMax ). Set by
        //! set_alpha_switch_temperature(), NaN until then; 0 K means "no
        //! cryogenic branch" ( the plain curve is used everywhere ).
        real mTAlphaSwitch = BELFEM_QUIET_NAN ;

        //! which E-J law the resistivity dispatch uses ( replaces the old
        //! mUsePiecewise bool; use_piecewise() tests == Piecewise ONLY, so
        //! Riva never falls into the Bézier branch by accident )
        ResistivityLaw mResistivityLaw = ResistivityLaw::PowerLaw ;

        // Function pointers for property evaluation (dispatched based on dependencies)
        real ( Material::*mFunctionDensity ) ( const real T )  const = nullptr ;

        real ( Material::*mFunctionCp ) ( const real T )  const = nullptr ;
        real ( Material::*mFunctionLambda ) ( const real T ) const = nullptr ;

        real ( Material::*mFunctiondCpdT ) ( const real T )  const = nullptr ;
        real ( Material::*mFunctiondLambdadT ) ( const real T ) const = nullptr ;

        real ( Material::*mFunctiond2CpdT2 ) ( const real T )  const = nullptr ;

        real ( Material::*mFunctionRho ) ( const real T ) const = nullptr ;
        real ( Material::*mFunctionRhoI ) ( const real T ) const = nullptr ;
        real ( Material::*mFunctiondRhodT ) ( const real T ) const = nullptr ;

        real ( Material::*mFunctionE ) ( const real T ) const = nullptr ;
        real ( Material::*mFunctionNu ) ( const real T ) const = nullptr ;

        real ( Material::*mFunctionAlpha ) ( const real T ) const = nullptr ;

        real ( Material::*mFunctionRp02 ) ( const real T ) const = nullptr ;
        real ( Material::*mFunctionDebye ) ( const real T ) const = nullptr ;

        DefectFunc * mDefectFunction  = nullptr ;
        HeatFunc   * mHeatFunction   = nullptr ;

        //! dlopen handle of the defect plugin that supplied mDefectFunction,
        //! null when no defect was read. Owned: closed by ~Material, which
        //! must outlive every mDefectFunction call
        void * mDefectHandle = nullptr ;

        //! dlopen handle of the heating plugin, same ownership as mDefectHandle
        void * mHeatHandle = nullptr ;

        // we don't want all materials to access private functions
        // but for this one, it's OK
        friend class material::SplineLookupTable ;

    protected:

        /**
         * O1 "full-signature policy" helpers: the assembly path always passes
         * the full ( T, normB, angleNxB ) set; the material consumes what its
         * Jc / n functions depend on and ignores the rest. JcFunction
         * subclasses override either eval( B, angle ) or eval( B, angle, T ) —
         * the respective other default throws — so we route on the declared
         * T-dependence. Falls back to the plain constants when no function is
         * attached.
         */
        real
        jc_eval( const real T, const real normB, const real angleNxB ) const ;

        //! d(jc)/d|B| routed like jc_eval: JcFunction::deval_dB when a
        //! function is attached, 0 for the constants fallback ( the exact
        //! derivative of a constant ) — so constant-jc decks are bit-identical
        real
        djc_eval_dB( const real T, const real normB, const real angleNxB ) const ;

        //! d(n)/d|B|, same routing as djc_eval_dB
        real
        dn_eval_dB( const real T, const real normB, const real angleNxB ) const ;

        //! d(jc)/dT, T-leg twin of djc_eval_dB: JcFunction::deval_dT when a
        //! function is attached, 0 for the constants fallback. For a
        //! T-dependent UserDefined without a deval_dT override the base-class
        //! 0 is a conservative fallback, not an exact derivative
        real
        djc_eval_dT( const real T, const real normB, const real angleNxB ) const ;

        //! d(n)/dT, same routing as djc_eval_dT
        real
        dn_eval_dT( const real T, const real normB, const real angleNxB ) const ;

        real
        n_eval( const real T, const real normB, const real angleNxB ) const ;

        //! setup gate for the riva law: refuses a PROVABLY bad n source
        //! ( table bound or constant at or below 1 ). Called from
        //! set_resistivity_law and set_n_function, so the outcome does not
        //! depend on which of the two is called first.
        void
        check_riva_n_source() const ;

        real ( Material::*mFunctionRhoKohler ) ( const real T, const real normB, const real angle ) const = nullptr ;
        real ( Material::*mFunctiondRhoKohlerdT ) ( const real T, const real normB, const real angle ) const = nullptr ;
        real ( Material::*mFunctiondRhoKohlerdB ) ( const real T, const real normB, const real angle ) const = nullptr ;
        real ( Material::*mFunctiondRhoKohlerdbeta ) ( const real T, const real normB, const real angle ) const = nullptr ;

        real ( Material::*mFunctionH ) ( const real B ) const = nullptr ;
        real ( Material::*mFunctionMu ) ( const real H, const real T ) const = nullptr ;
        void ( Material::*mFunctionDMuDH )( double, double&, double& ) const = nullptr ;

        // Function objects (material owns and deletes these)

        const material::JcFunction * mJcFunction = nullptr ;
        const material::JcFunction * mNFunction  = nullptr ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * @brief Constructor
         * @param aType Material type (Ferro, HTS, PureMetal, etc.)
         * @param aIsIsotropic True if material properties are isotropic (default)
         */
        Material( const MaterialType aType, const bool aIsIsotropic=true ) ;

        /**
         * @brief Destructor - deletes owned B-H curves and Jc/n functions
         */
        virtual ~Material() ;

        // Non-copyable and non-movable (Rule of Five)
        Material( const Material & ) = delete ;
        Material & operator=( const Material & ) = delete ;
        Material( Material && ) = delete ;
        Material & operator=( Material && ) = delete ;

//------------------------------------------------------------------------------
// Property Checks
//------------------------------------------------------------------------------

        /**
         * @brief Check if a material property is available
         * @param aProperty The property to check
         * @return True if the property is implemented for this material
         */
        bool
        have( const MaterialProperty aProperty ) const ;

        /**
         * @brief Check if the material has a defect
         * @return True if a defect function was defined
         */
        bool
        have_defect() const ;

        /**
         * @brief Check if the material has a heating function
         * @return True if a heating function was defined
         */
        bool
        have_heating() const ;


        /**
         * @brief Check if we use piecewise instead of power-law
         * @return True if a we use piecewise resistivity
         */
        bool
        use_piecewise() const ;

        /**
         * @brief Check if a property depends on a specific variable
         * @param aProperty The material property
         * @param aDependency The dependency to check (T, B, angle, etc.)
         * @return True if the property depends on the specified variable
         */
        bool
        depends( const MaterialProperty aProperty, const MaterialDependency aDependency ) const ;

        /**
         * @brief Get all dependencies for a property
         * @param aProperty The material property
         * @return Pointer to bitset containing all dependencies
         */
        const MaterialDependencyBitset *
        dependencies( const MaterialProperty aProperty ) const ;

        /**
         * @brief Check if material is isotropic
         * @return True if all properties are isotropic
         */
        bool
        is_isotropic() const ;

        /**
         * @brief Get the material type
         * @return Material type (Ferro, HTS, PureMetal, etc.)
         */
        MaterialType
        type() const ;

//------------------------------------------------------------------------------
// Parameters and Tools
//------------------------------------------------------------------------------

        /**
         * @brief Get material label
         * @return Material label string (e.g., "Copper", "YBCO")
         */
        const string &
        label() const ;

        /**
         * @brief Get material number/identifier
         * @return Material number if it exists, empty string otherwise
         */
        const string &
        number() const ;

        /**
         * @brief Activate a B-H curve as the permeability source (ferromagnets)
         * @param aCurve Concrete B-H curve; the material takes ownership.
         *
         * Kept free of the concrete curve type: the caller (MaterialFactory)
         * builds the linalg-backed BhSplineCurve, so base Material never depends
         * on Spline/linalg. Stores the curve via set_bh_curve() and routes
         * mu/H/dmudH through it.
         */
        void
        load_bh_curve( const material::BhCurve * aCurve );

        /**
         * @brief Set flag (multi-purpose flag used by Kernel)
         */
        void
        flag( const uint8_t aIndex = 0 );

        /**
         * @brief Clear flag
         */
        void
        unflag( const uint8_t aIndex = 0 );

        /**
         * @brief Check if material is flagged
         * @return True if flag is set
         */
        bool
        is_flagged( const uint8_t aIndex = 0 ) const ;

        /**
         * @brief Check if a property is constant (temperature-independent)
         * @param aProperty The property to check
         * @return True if the property is constant
         */
        bool
        is_constant( const MaterialProperty aProperty ) const ;

        /**
         * @brief Set the residual resistivity ratio (for noble metals)
         * @param RRR Residual resistivity ratio (ρ(273K)/ρ(0K))
         */
        virtual void
        set_RRR( const real RRR ) ;

//------------------------------------------------------------------------------
// Mass and Weight
//------------------------------------------------------------------------------

        /**
         * @brief Density as a function of temperature
         * @param T Temperature [K] (default: room temperature)
         * @return Density [kg/m³]
         *
         * IMPORTANT: For thermal M-matrix, use ref_density() or density(gTroom),
         * not density(T), because BELFEM computes on the undeformed mesh.
         */
        virtual real
        density( const real T=gTroom ) const ;

        /**
         * @brief Reference density at reference temperature
         * @return Reference density [kg/m³]
         *
         * Use this for thermal M-matrix calculations.
         */
        real
        ref_density() const ;

        /**
         * @brief Molar mass
         * @return Molar mass [kg/mol]
         */
        real
        M() const ;

//------------------------------------------------------------------------------
// Thermal Properties
//------------------------------------------------------------------------------

        /**
         * @brief Specific heat capacity
         * @param T Temperature [K] (default: room temperature)
         * @return Specific heat capacity [J/(kg·K)]
         */
        virtual real
        cp( const real T=gTroom ) const ;

        virtual real
        dcpdT( const real T=gTroom ) const ;

        virtual real
        d2cpdT2( const real T=gTroom ) const ;

        /**
         * @brief Thermal conductivity (isotropic)
         * @param T Temperature [K] (default: room temperature)
         * @return Thermal conductivity [W/(m·K)]
         */
        virtual real
        lambda( const real T=gTroom ) const ;

        virtual real
        dlambdadT( const real T=gTroom ) const ;


        /**
         * @brief Thermal conductivity for noble metals with magnetoresistance
         * @param T Temperature [K]
         * @param B Magnetic field magnitude [T]
         * @param beta Angle between field and current [rad]
         * @return Thermal conductivity [W/(m·K)]
         *
         * NOTE: metals evaluate Kohler's rule directly; when a metal was
         * constructed with tables enabled, set_RRR() builds <label>_RRR<n>.hdf5
         * and the field-dependent evaluation reads it instead. Alloys require
         * the table.
         */
        virtual real
        lambda( const real T, const real B, const real beta ) const ;

        virtual real
        dlambdadT( const real T, const real B, const real beta ) const ;

        virtual real
        dlambdadB( const real T, const real B, const real beta ) const ;

        virtual real
        dlambdadbeta( const real T, const real B, const real beta ) const ;

        /**
         * @brief Thermal conductivity for HTS materials
         * @param T Temperature [K]
         * @param B_par Parallel magnetic field component [T]
         * @param B_perp Perpendicular magnetic field component [T]
         * @param J Current density magnitude [A/m²]
         * @return Thermal conductivity [W/(m·K)]
         */
        virtual real
        lambda( const real T, const real B_par, const real B_perp, const real J ) const ;

        virtual real
        dlambdadT( const real T, const real B_par, const real B_perp, const real J ) const ;


//------------------------------------------------------------------------------
// Electric Properties
//------------------------------------------------------------------------------

        /**
         * @brief Electrical resistivity (isotropic)
         * @param T Temperature [K]
         * @return Electrical resistivity [Ω·m]
         *
         * For HTS materials, use rho_powerlaw() instead.
         */
        virtual real
        rho( const real T ) const ;

        virtual real
        drhodT( const real T ) const ;

        /**
         * @brief Electrical resistivity for noble metals with magnetoresistance
         * @param T Temperature [K]
         * @param B Magnetic field magnitude [T]
         * @param beta Angle between field and current [rad]
         * @return Electrical resistivity [Ω·m]
         *
         * NOTE: metals evaluate Kohler's rule directly; when a metal was
         * constructed with tables enabled, set_RRR() builds <label>_RRR<n>.hdf5
         * and the field-dependent evaluation reads it instead. Alloys require
         * the table.
         */
        virtual real
        rho( const real T, const real B, const real beta ) const ;

        virtual real
        drhodT( const real T, const real B, const real beta ) const ;

        virtual real
        drhodB( const real T, const real B, const real beta ) const ;

        virtual real
        drhodbeta( const real T, const real B, const real beta ) const ;

        /**
         * @brief Power law resistivity for HTS (constant jc and n)
         *
         * Parameters and return value are documented on the definition in
         * powerlaws.hpp, together with the derivation and its references.
         *
         * Computes ρ = ( 1/ρ_n + 1/ρ_PL )^-1 with ρ_PL = (Ec/jc)·(J/jc)^(n-1)
         * for constant jc and n ( parallel combination with the normal-state
         * resistivity ρ_n, see powerlaws.hpp ).
         * Use this for HTS materials instead of rho(T).
         */
        real
        rho_powerlaw( const real normJ ) const ;

        real
        rho_powerlaw( const real normJ, const real x, const real y, const real z, const real t ) const ;

        /**
         * @brief Power law resistivity for HTS (field-dependent jc, constant n)
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @return Electrical resistivity [Ω·m]
         *
         * Computes ρ = ( 1/ρ_n + 1/ρ_PL )^-1 with ρ_PL = (Ec/jc)·(J/jc)^(n-1)
         * where jc = jc(B, angle).
         */
        real
        rho_powerlaw( const real normJ, real normB, const real angleNxB  ) const ;

        real
        rho_powerlaw( const real normJ, real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Power law resistivity for HTS (full temperature and field dependence)
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @param T Temperature [K]
         * @return Electrical resistivity [Ω·m]
         *
         * Computes ρ = ( 1/ρ_n + 1/ρ_PL )^-1 with ρ_PL = (Ec/jc)·(J/jc)^(n-1)
         * where jc = jc(B, angle, T) and n = n(B, angle, T).
         */
        real
        rho_powerlaw( const real normJ, const real T, real normB, const real angleNxB  ) const ;

        real
        rho_powerlaw( const real normJ, const real T, real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        //Temperature dependent only cases (custom)
        real
        rho_powerlaw( const real normJ, const real T  ) const ;

        real
        rho_powerlaw( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Piecewise resistivity for HTS (constant jc and n)
         *
         * Parameters and return value are documented on the definition in
         * powerlaws.hpp, together with the derivation and its references.
         *
         * Use this for HTS materials instead of rho(T).
         */
        real
        rho_piecewise( const real normJ ) const ;

        real
        rho_piecewise( const real normJ, const real x, const real y, const real z, const real t ) const ;

        /**
         * @brief Piecewise resistivity for HTS (field-dependent jc, constant n)
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @return Electrical resistivity [Ω·m]
         *
         */
        real
        rho_piecewise( const real normJ, real normB, const real angleNxB  ) const ;

        real
        rho_piecewise( const real normJ, real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Piecewise resistivity for HTS (full temperature and field dependence)
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @param T Temperature [K]
         * @return Electrical resistivity [Ω·m]
         *
         */
        real
        rho_piecewise( const real normJ, const real T, real normB, const real angleNxB  ) const ;

        real
        rho_piecewise( const real normJ, const real T, real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        //Temperature dependent only cases (custom)
        real
        rho_piecewise( const real normJ, const real T  ) const ;

        real
        rho_piecewise( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Riva-law resistivity: the superconducting power-law channel
         *        in parallel with the normal-state channel
         * @return ρ = ρ_PL·ρ_n / ( ρ_PL + ρ_n ) [Ω·m]
         *
         * The parallel combination of Duron et al. 2004, used by Riva 2021
         * (EPFL thesis 8754, Eq. 5.4). Same model core as rho_powerlaw, but
         * total over the full table range: jc_eff ≤ 0 or nonfinite (dead
         * defect, underflow) falls back to the fully normal branch, the
         * power-law channel is evaluated in log10 space with an overflow
         * early-out, and n arrives pre-floored at 1 from n_eval (ohmic
         * limit ρ_PL = ec/jc, J-independent). Only the assembly-path
         * signatures (J,T,|B|,θ) ± defect exist; jc and n arrive through
         * the usual dependency routing ( table or constants ).
         */
        real
        rho_riva( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        //! defect overload: jc_eff = D(x,y,z,t)·jc throughout
        real
        rho_riva( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        //! dρ/d|J| of rho_riva: w²·dρ_PL/dJ with w = ρ_n/(ρ_PL+ρ_n)
        real
        drho_riva_dJ( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        real
        drho_riva_dJ( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        //! dρ/d|B| of rho_riva: w²·dρ_PL/dB ( dρ_n/dB = 0 on this path )
        real
        drho_riva_dB( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        real
        drho_riva_dB( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        //! dρ/dT of rho_riva: w²·dρ_PL/dT + (1−w)²·dρ_n/dT
        real
        drho_riva_dT( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        real
        drho_riva_dT( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        real
        n( const real normB, const real angleNxB, const real T  ) const ;

        real
        n( const real normB, const real angleNxB ) const ;

        real
        jc( const real normB, const real angleNxB, const real T  ) const ;

        real
        jc( const real normB, const real angleNxB  ) const ;

        //Case of jc with user-defined defect
        real
        jc( const real normB, const real angleNxB, const real T, const real x, const real y, const real z, const real t  ) const ;

        real
        jc( const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of power-law resistivity with respect to current density magnitude
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @param T Temperature [K]
         * @return Derivative dρ/d||J|| [Ω·m³/A]
         *
         * For Newton-Raphson: computes dρ/d||J|| where ρ = ρ(||J||, B, θ, T)
         */
        real
        drho_powerlaw_dJ( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        // Overloads matching rho_powerlaw signatures
        real
        drho_powerlaw_dJ( const real normJ ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real x, const real y, const real z, const real t ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real T  ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real normB, const real angleNxB  ) const ;

        real
        drho_powerlaw_dJ( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of piecewise resistivity with respect to current density magnitude
         * @param normJ Current density magnitude [A/m²]
         * @param normB Magnetic field magnitude [T]
         * @param angleNxB Angle between surface normal and B [rad]
         * @param T Temperature [K]
         * @return Derivative dρ/d||J|| [Ω·m³/A]
         *
         * For Newton-Raphson: computes dρ/d||J|| for the piecewise model
         */
        real
        drho_piecewise_dJ( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of power-law resistivity with respect to |B|
         *        at fixed J, T, θ ( jc = jc(T,|B|,θ), n = n(T,|B|,θ) )
         * @return dρ/d|B| [Ω·m/T]
         *
         * Chain: dρ/d|B| = (r/p)² · [ ∂p/∂jc·djc/d|B| + ∂p/∂n·dn/d|B| ]
         * with p the power-law resistivity, r the parallel combination with
         * ρ_n, ∂p/∂jc = −n·p/jc and ∂p/∂n = p·ln(J/jc). Parallel factor and
         * floor convention mirror drho_powerlaw_dJ exactly ( floored p in
         * the factor, unfloored law differentiated ).
         */
        real
        drho_powerlaw_dB( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        //! defect overload: jc_eff = D(x)·jc throughout, djc_eff = D·djc
        real
        drho_powerlaw_dB( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of power-law resistivity with respect to T
         *        at fixed J, |B|, θ ( T-leg: jc(T), n(T) AND ρ_n(T)
         *        all move — the quench-feedback tangent )
         * @return dρ/dT [Ω·m/K]
         *
         * Chain: dρ/dT = dp₀/dT / (1+ρ_PL/ρ_n)² + dρ_n/dT · (ρ_PL/(ρ_n+ρ_PL))²
         * with dp₀/dT = p₀·( −(n/jc)·djc/dT + ln(J/jc)·dn/dT ). Floor
         * convention mirrors drho_powerlaw_dB ( floored ρ_PL in the parallel
         * factors, unfloored law differentiated ). No (djc==0 && dn==0)
         * early-out: the ρ_n term is the only correct T-dependence of a
         * constant-jc material and must survive.
         */
        real
        drho_powerlaw_dT( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        //! defect overload: jc_eff = D(x)·jc throughout, djc_eff = D·djc
        real
        drho_powerlaw_dT( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of piecewise resistivity with respect to |B|
         *        at fixed J, T, θ
         * @return dρ/d|B| [Ω·m/T]
         *
         * Exact per regime where the regime's law permits, staged where it
         * does not ( 2026-08-13 audit, C5 accepted by both voices ):
         * - power-law regime ( J ≤ j1 ): dp/d|B| WITHOUT the parallel
         *   factor, because rho_piecewise returns the RAW power law there —
         *   consistency is with the residual, not with drho_powerlaw_dB.
         * - normal regime ( J > j3 ) and T > T_crit: exactly 0 ( ρ = ρ_n is
         *   independent of jc, n ).
         * - Bézier blend ( j1 < J ≤ j3 ): 0 for now — the boundaries j1..j3
         *   and the control points also move with jc, and that derivative is
         *   deferred; a zero here is today's global behaviour, locally.
         *   The tangent therefore has a jump at j1.
         */
        real
        drho_piecewise_dB( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        //! defect overload: jc_eff = D(x)·jc throughout, djc_eff = D·djc
        real
        drho_piecewise_dB( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Derivative of piecewise resistivity with respect to T
         *        at fixed J, |B|, θ ( T-leg )
         * @return dρ/dT [Ω·m/K]
         *
         * Follows rho_piecewise's OWN residual branch by branch:
         * - T > T_crit and normal regime ( J > j3 ): exactly dρ_n/dT.
         * - power-law regime ( J ≤ j1 ): raw-law derivative, NO parallel
         *   factor and NO ρ_n term — consistency is with the residual.
         * - Bézier blend ( j1 < J ≤ j3 ): BOTH parts — the frozen-knot
         *   partial and the knot motion dtParam/dT through j1(T)..j3(T);
         *   the knot motion dominates. Degrades to the frozen part alone at
         *   the degenerate n−1 == mNff transition. Jumps at j1 and T_crit
         *   are inherited from the residual itself.
         */
        real
        drho_piecewise_dT( const real normJ, const real T, const real normB, const real angleNxB  ) const ;

        //! defect overload: jc_eff = D(x)·jc throughout, djc_eff = D·djc
        real
        drho_piecewise_dT( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        // Overloads matching rho_piecewise signatures
        real
        drho_piecewise_dJ( const real normJ ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real x, const real y, const real z, const real t ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real T  ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real normB, const real angleNxB  ) const ;

        real
        drho_piecewise_dJ( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const ;

        /**
         * @brief Intrinsic (phonon) electrical resistivity of a metal, ρ = ρ_i(T) + ρ_0
         * @param T Temperature [K]
         * @return Intrinsic electrical resistivity [Ω·m]
         */
        virtual real
        rho_i( const real T ) const ;

//------------------------------------------------------------------------------
// Magnetic Properties
//------------------------------------------------------------------------------

        /**
         * @brief Magnetic field strength from flux density
         * @param B Magnetic flux density [T]
         * @return Magnetic field strength H [A/m]
         */
        virtual real
        H( const real B ) const ;

        /**
         * @brief Magnetic permeability
         * @param H Magnetic field strength [A/m]
         * @param T Temperature [K] (used by user-defined / polynomial mu(H,T); ignored by the constant and B-H-curve bindings)
         * @return Magnetic permeability [H/m]
         */
        virtual real
        mu( const real H, const real T=BELFEM_QUIET_NAN ) const ;

        /**
         * @brief Magnetic permeability and its derivative
         * @param H Magnetic field strength [A/m]
         * @param mu Output: magnetic permeability [H/m]
         * @param dmudH Output: derivative dμ/dH [H/m²]
         */
        virtual void
        dmudH( const real H, real & mu, real & dmudH ) const ;

//------------------------------------------------------------------------------
// Mechanical Properties
//------------------------------------------------------------------------------

        /**
         * @brief Young's modulus
         * @param T Temperature [K] (default: room temperature)
         * @return Young's modulus [Pa]
         */
        virtual real
        E( real T=gTroom ) const ;

        /**
         * @brief Poisson's ratio
         * @param T Temperature [K] (default: room temperature)
         * @return Poisson's ratio [-]
         */
        virtual real
        nu( real T=gTroom ) const ;

        /**
         * @brief Shear modulus
         * @param T Temperature [K] (default: room temperature)
         * @return Shear modulus G = E/(2(1+ν)) [Pa]
         */
        virtual real
        G( real T=gTroom ) const ;

        /**
         * @brief Bulk modulus
         * @param T Temperature [K] (default: room temperature)
         * @return Bulk modulus K = E/(3(1-2ν)) [Pa]
         */
        virtual real
        K( real T=gTroom ) const ;

        /**
         * @brief Thermal expansion coefficient
         * @param T Temperature [K] (default: room temperature)
         * @return Thermal expansion coefficient [1/K]
         *
         * Defined as α = (1/l)·∂l/∂T, NOT (1/l)·Δl/ΔT!
         */
        virtual real
        alpha( real T=gTroom ) const ;

        /**
         * @brief Relative length after thermal expansion
         * @param T Temperature [K]
         * @return Relative length l/l₀ [-]
         */
        virtual real
        l( real T ) const ;

        /**
         * @brief Yield stress (0.2% offset)
         * @param T Temperature [K] (default: room temperature)
         * @return Yield stress Rp0.2 [Pa]
         */
        virtual real
        Rp02( real T=gTroom ) const ;

//------------------------------------------------------------------------------
// Other Properties
//------------------------------------------------------------------------------

        /**
         * @brief Debye temperature
         * @param T Temperature [K]
         * @return Debye temperature [K]
         */
        virtual real
        debye( const real T ) const ;

        /**
         * @brief Assign a B-H curve to this material
         * @param aCurve Pointer to B-H curve (material takes ownership)
         *
         * IMPORTANT: Material takes ownership and will delete the curve.
         *
         * Virtual so that load_bh_curve() (defined on the base) dispatches to
         * the B-H-capable override in Metal; the base implementation errors.
         */
        virtual void
        set_bh_curve( const material::BhCurve * aCurve );

        /**
         * @brief Assign a critical current density function
         * @param aFunction Pointer to Jc function (material takes ownership)
         *
         * IMPORTANT: Material takes ownership and will delete the function.
         */
        void
        set_jc_function( const material::JcFunction * aFunction );

        /**
         * @brief Assign a power law exponent function
         * @param aFunction Pointer to n function (material takes ownership)
         *
         * IMPORTANT: Material takes ownership and will delete the function.
         */
        void
        set_n_function( const material::JcFunction * aFunction );

        void
        set_piecewise( const bool aUsePiecewise );

        void
        set_resistivity_law( const ResistivityLaw aLaw );

        ResistivityLaw
        resistivity_law() const ;

        bool
        use_riva() const ;
//------------------------------------------------------------------------------
// API for user defined functions
//------------------------------------------------------------------------------

        /**
         * @brief Set a user-defined function with one dependency (for UserDefinedMaterial)
         *
         * This method is primarily used by UserDefinedMaterial to assign custom property
         * evaluation functions. Most users should not call this directly; instead, use
         * MaterialFactory::create_material() with a user library.
         *
         * @param Property The material property to define
         * @param Dependency The dependency type (typically MaterialDependency::T)
         * @param Function Pointer to user function with signature: real(const Material*, real)
         *
         * @note This is a virtual method overridden by UserDefinedMaterial
         */
        virtual void
        set_user_defined_function(
            const MaterialProperty   Property,
            const MaterialDependency Dependency,
                  MatFunc1 * Function ) ;

        /**
         * @brief Set a user-defined function with two dependencies (for UserDefinedMaterial)
         *
         * This method is primarily used by UserDefinedMaterial to assign custom property
         * evaluation functions for two-parameter properties like mu(H,T) or jc(B,angle).
         *
         * @param Property The material property to define
         * @param Dependency1 First dependency
         * @param Dependency2 Second dependency
         * @param Function Pointer to user function with signature: real(const Material*, real, real)
         *
         * @note This is a virtual method overridden by UserDefinedMaterial
         */
        virtual void
        set_user_defined_function(
            const MaterialProperty   Property,
            const MaterialDependency Dependency1,
            const MaterialDependency Dependency2,
                  MatFunc2 * Function ) ;

        /**
         * @brief Set a user-defined function with three dependencies (for UserDefinedMaterial)
         *
         * This method is primarily used by UserDefinedMaterial to assign custom property
         * evaluation functions for three-parameter properties like rho(T,B,angle).
         *
         * @param Property The material property to define
         * @param Dependency1 First dependency
         * @param Dependency2 Second dependency
         * @param Dependency3 Third dependency
         * @param Function Pointer to user function with signature: real(const Material*, real, real, real)
         *
         * @note This is a virtual method overridden by UserDefinedMaterial
         */
        virtual void
        set_user_defined_function(
            const MaterialProperty   Property,
            const MaterialDependency Dependency1,
            const MaterialDependency Dependency2,
            const MaterialDependency Dependency3,
                  MatFunc3 * Function ) ;

        /**
         * @brief Set a user-defined defect function with x, y, z, t dependencies
         *
         *
         * @param Function Pointer to user function with signature: real(real, real, real, real)
         *
         */
        void
        set_user_defined_defect( DefectFunc * Function ) ;

        /**
         * @brief Read a defect from the library
         *
         * @param aLibraryPath path to the shared library holding the defect,
         *                     resolved through material::data_file(): the run
         *                     directory first, then $BELFEM_DATA/material, and
         *                     otherwise handed to dlopen unchanged
         * @param aLabel       name of the defect to read
         */
        void
        read_defect( const string & aLibraryPath,
                    const string & aLabel ) ;

        /**
         * @brief Set a user-defined volumetric heat load [ W/m³ ] as a
         *        function of x, y, z [ m ] and t [ s ]
         *
         * @param Function Pointer to user function with signature: real(real, real, real, real)
         */
        void
        set_user_defined_heating( HeatFunc * Function ) ;

        /**
         * @brief Read a heating function from a plugin library
         *
         * @param aLibraryPath path to the shared library holding the heating
         *                     function, resolved through material::data_file():
         *                     the run directory first, then $BELFEM_DATA/material,
         *                     and otherwise handed to dlopen unchanged
         * @param aLabel       name of the heating function; the library must
         *                     export <label>_init( Material * )
         */
        void
        read_heating( const string & aLibraryPath,
                    const string & aLabel ) ;


        /**
         * @brief Evaluate a polynomial for a given property (for UserDefinedMaterial)
         *
         * Internal method used by UserDefinedMaterial to evaluate polynomial-based
         * property functions. Uses Horner's method for efficient evaluation.
         *
         * @param Property The material property
         * @param T Temperature [K]
         * @return Property value at temperature T
         *
         * @note This is a virtual method overridden by UserDefinedMaterial
         */
        virtual real
        evaluate_polynomial(
            const MaterialProperty  Property,
            const real T ) const ;

        virtual real
        evaluate_derivative_of_polynomial(
            const MaterialProperty  Property,
            const real T ) const ;


        /**
         * @brief Set a polynomial function for a property (for UserDefinedMaterial)
         *
         * Defines a property as a polynomial in temperature: f(T) = c₀T^n + c₁T^(n-1) + ... + c_n
         * This is a convenience function for user-defined materials.
         *
         * IMPORTANT: Coefficients are in DESCENDING order (MATLAB style), highest degree first.
         *
         * @param Property The material property to define
         * @param Coefficients Polynomial coefficients [c₀, c₁, c₂, ...] in DESCENDING order (highest degree first)
         *
         * @note This is a virtual method overridden by UserDefinedMaterial
         */
        virtual void
        set_user_defined_polynomial(
            const MaterialProperty Property,
            const std::vector< real > & Coefficients ) ;

        virtual void
        set_user_defined_polynomial(
            const MaterialProperty Property,
            const Cell< real > & Coefficients ) ;

        /**
         * PROTECTED INTERFACE FOR DERIVED MATERIALS
         *
         * This section contains utility functions for derived material classes
         * to set up their properties. The Material class uses a strategy pattern
         * where properties can be:
         * 1. Constant - stored directly in mConstantProperties
         * 2. Spline-interpolated - using tabulated data
         * 3. Custom - using specialized functions in derived classes
         *
         * When a property is set up, the appropriate function pointer is assigned
         * (e.g., mFunctionRho points to rho_const, rho_spline, or rho_custom).
         */

        /**
         * @brief Set material label (for derived classes)
         */
        void
        set_label( const string & aLabel ) ;

        /**
         * @brief Set material number/identifier (for derived classes)
         */
        void
        set_number( const string & aNumber ) ;

        /**
         * @brief Mark a property as available
         */
        void
        set_have( const MaterialProperty aProperty, const bool aHave=true ) ;

        /**
         * @brief Reset all dependencies for a property
         */
        void
        reset_dependencies( const MaterialProperty aProperty );

        /**
         * @brief Add a dependency to a property
         */
        void
        set_dependency( const MaterialProperty aProperty, const MaterialDependency aDependency ) ;

        /**
         * @brief Define a property as constant
         */
        void
        set_constant( const MaterialProperty aProperty, const real aValue ) ;

        /**
         * @brief Mark property as using custom evaluation function
         */
        void
        set_custom( const MaterialProperty aProperty ) ;

        /**
         * @brief Get constant property value
         */
        real
        constant_property( const MaterialProperty aProperty ) const ;

        /**
         * @brief Evaluate property using spline interpolation
         */
        virtual real
        spline_property( const MaterialProperty aProperty, const real aX ) const ;

        virtual real
        dspline_property( const MaterialProperty aProperty, const real aX ) const ;

        virtual real
        ddspline_property( const MaterialProperty aProperty, const real aX ) const ;

        virtual void
        set_table_flags( const bool aFlag );

        virtual real
        jc_custom( const real T ) const ;

        virtual real
        n_custom( const real T ) const ;

        /**
         * @brief Artificial volumetric heat load [ W/m³ ] from the heating plugin
         */
        real
        volumetric_heatload( const real x, const real y, const real z, const real time ) const ;

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

        //! unfloored n as routed from table / custom / constant — n_eval
        //! clamps this at 1, and the dn_eval_* derivatives return 0 while
        //! the clamp binds
        real
        n_eval_raw( const real T, const real normB, const real angleNxB ) const ;

        //! shared riva kernel: power-law channel resistivity in log10 space.
        //! Returns false when the channel has overflowed past the normal
        //! state ( caller takes the fully-normal branch: ρ = ρn )
        bool
        riva_rho_pl( const real normJ, const real jc, const real n, const real ec, real & rhoPL ) const ;

//------------------------------------------------------------------------------
        /**
         * PROPERTY EVALUATION FUNCTIONS
         *
         * These functions are called via function pointers for property evaluation.
         * Each property has up to four variants:
         * - *_const:  Returns constant value
         * - *_spline: Interpolates from spline
         * - *_custom: Calls derived class's custom implementation
         * - *_table:  Uses lookup table (for field-dependent properties)
         *
         * The appropriate function is selected during material initialization
         * based on how the property is defined.
         */


        /**
        * @brief Create spline for a property with optional boundary derivatives
        */
        void
        create_spline(  const MaterialProperty aProperty,
            const real adYdX0=BELFEM_QUIET_NAN,
            const real adXdX1=BELFEM_QUIET_NAN ) ;


        /**
         * @brief Temperature below which alpha is taken from the Grueneisen branch
         *
         * min( 0.618 * theta_Debye, gTAlphaSwitchMax ). The Debye factor is a
         * convention, not a derivation - it sits comfortably above the
         * temperature at which dln(C)/dT changes sign, which is where the
         * fitted expansion curve stops being consistent with the heat capacity.
         */
        real
        alpha_switch_temperature() const ;

        real
        density_const( const real T ) const ;

        virtual real
        density_custom( const real T ) const ;

        /**
         * @brief Shared binding for derivative channels that are identically zero
         *
         * Bound wherever a property is held constant, so that its first and
         * higher derivatives vanish. One function rather than one clone per
         * property and per derivative order.
         */
        real
        return_zero( const real T ) const ;

        real
        cp_const( const real T ) const ;

        real
        dcpdT_finite_difference( const real T ) const ;

        real
        d2cpdT2_finite_difference( const real T ) const ;

        real
        cp_spline( const real T ) const ;

        real
        dcpdT_spline( const real T ) const ;

        real
        d2cpdT2_spline( const real T ) const ;


        virtual real
        cp_custom( const real T ) const ;

        virtual real
        dcpdT_custom( const real T ) const ;

        virtual real
        d2cpdT2_custom( const real T ) const ;

        real
        lambda_const( const real T ) const  ;

        real
        dlambdadT_finite_difference( const real T ) const ;

        real
        lambda_spline( const real T ) const ;

        real
        dlambdadT_spline( const real T ) const ;

        virtual real
        lambda_custom( const real T ) const  ;

        virtual real
        dlambdadT_custom( const real T ) const  ;

        real
        drhodT_finite_difference( const real T ) const ;

        real
        rho_const( const real T ) const ;

        virtual real
        rho_custom( const real T ) const ;

        virtual real
        drhodT_custom( const real T ) const ;

        virtual real
        rho_spline( const real T ) const ;

        virtual real
        drhodT_spline( const real T ) const ;

        virtual real
        rho_kohler( const real T, const real normB, const real angleJxB  ) const ;

        virtual real
        drhodT_kohler( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        drhodB_kohler( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        drhodbeta_kohler( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        rho_table( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        drhodT_table( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        drhodB_table( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        drhodbeta_table( const real T, const real normB, const real angleJxB ) const ;

        virtual real
        lambda_custom( const real T, const real normB, const real angle ) const ;

        virtual real
        lambda_table( const real T, const real normB, const real angle ) const ;

        real
        H_const( const real B ) const ;

        virtual real
        H_bhcurve( const real B ) const ;

        real
        mu_const( const real H, const real T ) const ;

        virtual real
        mu_bhcurve( const real H, const real T ) const ;

        virtual real
        mu_custom( const real H, const real T ) const ;

        void
        dmudH_const( const real H, real & mu, real & dmudH ) const ;

        virtual void
        dmudH_bhcurve( const real H, real & mu, real & dmudH ) const ;

        real
        E_const( const real T ) const ;

        real
        E_spline( const real T) const ;

        virtual real
        E_custom( const real T ) const ;

        virtual real
        dEdT_custom( const real T ) const ;

        real
        nu_const( const real T ) const ;

        real
        nu_spline( const real T ) const ;

        virtual real
        nu_custom( const real T ) const ;

        real
        alpha_const( const real T ) const ;

        real
        alpha_spline( const real T ) const ;

        virtual real
        alpha_custom( const real T ) const ;

        real
        Rp02_const( const real T ) const ;

        real
        Rp02_spline( const real T ) const ;

        virtual real
        Rp02_custom( const real T ) const ;

        real
        debye_const( const real T ) const ;

        real
        debye_spline( const real T ) const ;

        virtual real
        debye_custom( const real T ) const ;

        virtual real
        rho_i_custom( const real T ) const ;

        real
        rho_i_spline( const real T ) const ;

        virtual void
        reset_spline( const MaterialProperty aProperty ) ;

        virtual void
        create_spline( real (Material::*aFunction)(const real aT) const,
            const MaterialProperty aProperty,
            const uint aStartBC,
            const uint aEndBC,
            const real adYdX0=BELFEM_QUIET_NAN,
            const real adYdX1=BELFEM_QUIET_NAN ) ;

    };

    inline const string &
    Material::label() const
    {
        return mLabel ;
    }

    inline const string &
    Material::number() const
    {
        return mNumber ;
    }

    inline MaterialType
    Material::type() const
    {
        return mType ;
    }


    inline void
    Material::flag( const uint8_t aIndex )
    {
        BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                                  ( unsigned int ) aIndex );
        mFlags |= ( 1 << aIndex );
    }


    inline void
    Material::unflag( const uint8_t aIndex )
    {
        BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                                   ( unsigned int ) aIndex );
        mFlags &= ~( 1 << aIndex );
    }

    inline bool
    Material::is_flagged( const uint8_t aIndex ) const
    {
        BELFEM_ASSERT( aIndex < 8, "Flag index %u out of bounds (must be < 8)",
                                ( unsigned int ) aIndex );
        return mFlags & ( 1 << aIndex );
    }

    inline bool
    Material::have( const MaterialProperty aProperty ) const
    {
        return mHaveProperty.test( static_cast< size_t >( aProperty ) ) ;
    }

    inline bool
    Material::have_defect() const
    {
        return mDefectFunction != nullptr ;
    }

    inline bool
    Material::have_heating() const
    {
        return mHeatFunction != nullptr ;
    }

    inline bool
    Material::use_piecewise() const
    {
        return mResistivityLaw == ResistivityLaw::Piecewise ;
    }

    inline bool
    Material::use_riva() const
    {
        return mResistivityLaw == ResistivityLaw::Riva ;
    }

    inline ResistivityLaw
    Material::resistivity_law() const
    {
        return mResistivityLaw ;
    }

    inline bool
    Material::depends( const MaterialProperty aProperty, const MaterialDependency aDependency ) const
    {
        return mPropertyDependencies( static_cast< size_t >( aProperty ) )->test( static_cast< uint >( aDependency ));
    }

    inline const MaterialDependencyBitset *
    Material::dependencies( const MaterialProperty aProperty ) const
    {
        return mPropertyDependencies( static_cast< size_t >( aProperty ) ) ;
    }

    inline bool
    Material::is_isotropic() const
    {
        return mIsIsotropic ;
    }

    inline bool
    Material::is_constant( const MaterialProperty aProperty ) const
    {
        return ! std::isnan( mConstantProperties( static_cast< index_t >( aProperty ) ) ) ;
    }

    inline real
    Material::constant_property( const MaterialProperty aProperty ) const
    {
        BELFEM_ASSERT( this->is_constant( aProperty ), "Property is not constant" ) ;
        return mConstantProperties( static_cast< size_t >( aProperty ) ) ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::density( real T ) const
    {
        BELFEM_ASSERT( this->have( MaterialProperty::density ),
            "Material '%s' does not have density assigned", this->label().c_str() ) ;

        return ( this->*mFunctionDensity )( T );
    }

    inline real
    Material::density_const( real T ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::density ),
           "Wrong function call for Material::density");
        return this->constant_property( MaterialProperty::density );
    }

    inline real
    Material::density_custom( real T ) const
    {
        BELFEM_ASSERT( this->have( MaterialProperty::ref_density ),
            "Material %s does not have a reference density set", mLabel.c_str() );
        BELFEM_ASSERT( this->have( MaterialProperty::T_ref_density ),
         "Material %s does not have a temperature for the reference density set", mLabel.c_str() );

        if ( T == mConstantProperties( static_cast< index_t >( MaterialProperty::T_ref_density ) ) )
        {
            return mConstantProperties( static_cast< index_t >( MaterialProperty::ref_density ) ) ;
        }

        return mConstantProperties( static_cast< index_t >( MaterialProperty::ref_density ) ) / std::pow( this->l( T ), 3.0 ) ;
    }

    inline real
    Material::ref_density() const
    {
        BELFEM_ASSERT( this->have( MaterialProperty::ref_density ),
            "Material %s does not have a reference density set", mLabel.c_str() );
        return this->constant_property( MaterialProperty::ref_density );
    }

    inline real Material::M() const
    {
        return this->constant_property( MaterialProperty::M );
    }

    inline real
    Material::cp( real T ) const
    {
       BELFEM_ASSERT( this->have( MaterialProperty::cp ), "Material does not have specific heat capacity" ) ;

       return ( this->*mFunctionCp )( T );
    }


    inline real
    Material::cp_const( real T ) const
    {
        return this->constant_property( MaterialProperty::cp );
    }

    inline real Material::dcpdT( const real T ) const
    {
        BELFEM_ASSERT( this->have( MaterialProperty::cp ), "Material does not have specific heat capacity" ) ;

        return (this->*mFunctiondCpdT )( T ) ;
    }

    inline real Material::d2cpdT2( const real T ) const
    {
        BELFEM_ASSERT( this->have( MaterialProperty::cp ), "Material does not have specific heat capacity" ) ;

        return (this->*mFunctiond2CpdT2 )( T ) ;
    }

    inline real
    Material::return_zero( const real T ) const
    {
        return 0. ;
    }

    inline real
    Material::dcpdT_finite_difference( const real T ) const
    {
        return 0.5 * (  this->cp( T + gFinDiffDeltaT )
                      - this->cp( T - gFinDiffDeltaT ) ) /gFinDiffDeltaT ;
    }

    inline real
    Material::d2cpdT2_finite_difference( const real T ) const
    {
        return (  this->cp( T + gFinDiffDeltaT ) - 2. * this->cp( T )
                      + this->cp( T - gFinDiffDeltaT ) ) / ( gFinDiffDeltaT * gFinDiffDeltaT );
    }

    inline real
    Material::cp_spline( real T ) const
    {
       return this->spline_property( MaterialProperty::cp, T );
    }

    inline real
    Material::dcpdT_spline( real T ) const
    {
        return this->dspline_property( MaterialProperty::cp, T );
    }

    inline real
    Material::d2cpdT2_spline( real T ) const
    {
        return this->ddspline_property( MaterialProperty::cp, T );
    }


//------------------------------------------------------------------------------

    inline real
    Material::lambda( const real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
           "Material does not have thermal conductivity" ) ;
       return (this->*mFunctionLambda )( T ) ;
    }

    inline real
    Material::lambda_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::lambda );
    }

    inline real
    Material::dlambdadT( const real T ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
            "Material does not have thermal conductivity" ) ;

        return (this->*mFunctiondLambdadT )( T ) ;
    }

    inline real
    Material::dlambdadT_finite_difference( const real T ) const
    {
        return 0.5 * (  this->lambda( T + gFinDiffDeltaT )
                      - this->lambda( T - gFinDiffDeltaT ) ) /gFinDiffDeltaT ;
    }

    inline real
    Material::lambda_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::lambda, T );
    }

    inline real
    Material::dlambdadT_spline( const real T ) const
    {
        return this->dspline_property( MaterialProperty::lambda, T );
    }

    inline real
    Material::lambda( const real T, const real B, const real beta ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
           "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

       BELFEM_ASSERT( mType == MaterialType::PureMetal,
           "wrong function call for Material::lambda in material %s", mLabel.c_str() ) ;

       return this->lambda( T );
    }

    inline real
    Material::dlambdadT( const real T, const real B, const real beta ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
            "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::PureMetal,
            "wrong function call for Material::dlambdadT in material %s", mLabel.c_str() ) ;

        return this->dlambdadT( T );
    }

    inline real
    Material::dlambdadB( const real T, const real B, const real beta ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
            "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::PureMetal,
            "wrong function call for Material::dlambdadB in material %s", mLabel.c_str() ) ;

        return 0.;
    }

    inline real
    Material::dlambdadbeta( const real T, const real B, const real beta ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
            "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::PureMetal,
            "wrong function call for Material::dlambdadbeta in material %s", mLabel.c_str() ) ;

        return 0.;
    }


//------------------------------------------------------------------------------

    inline real
    Material::lambda( const real T, const real B_par, const real B_perp, const real J ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
           "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::HTS, "wrong function call for Material::lambda" ) ;

        return this->lambda( T );
    }

    inline real
    Material::dlambdadT( const real T, const real B_par, const real B_perp, const real J )const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::lambda ) ,
            "Material %s does not have thermal conductivity", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::HTS,
            "wrong function call for Material::dlambdadT in material %s", mLabel.c_str() ) ;

        return this->dlambdadT( T );
    }

//------------------------------------------------------------------------------

    inline real
    Material::rho( const real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::rho ) ,
           "Material %s does not have electric resistivity", mLabel.c_str() ) ;

       return ( this->*mFunctionRho )( T ) ;
    }


    inline real
    Material::drhodT( const real T ) const
    {
        BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->have( MaterialProperty::rho ) ,
            "Material %s does not have electric resistivity", mLabel.c_str() ) ;

        return ( this->*mFunctiondRhodT )( T ) ;
    }

    inline real
    Material::drhodT_finite_difference( const real T ) const
    {
        return 0.5 * (  this->rho( T + gFinDiffDeltaT )
                      - this->rho( T - gFinDiffDeltaT ) ) / gFinDiffDeltaT ;
    }

    inline real
    Material::rho_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::rho );
    }

    inline real
    Material::drhodT_custom( const real T ) const
    {
        return this->drhodT_finite_difference( T );
    }

    inline real
    Material::rho_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::rho, T );
    }

    inline real
    Material::drhodT_spline( const real T ) const
    {
        return this->dspline_property( MaterialProperty::rho, T );
    }

    inline real
    Material::rho( const real T, const real B, const real beta ) const
    {

        BELFEM_ASSERT( mFunctionRhoKohler != nullptr, "Material %s does not have a Kohler function", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::PureMetal,
           "wrong function call for Material::rho in %s", mLabel.c_str() ) ;

       return this->rho( T );
    }

    inline real
    Material::drhodT( const real T, const real B, const real beta ) const
    {

        BELFEM_ASSERT( mFunctionRhoKohler != nullptr, "Material %s does not have a Kohler function", mLabel.c_str() ) ;

        BELFEM_ASSERT( mType == MaterialType::PureMetal,
           "wrong function call for Material::drhodT in %s", mLabel.c_str() ) ;

        return this->drhodT( T );
    }

    // base fallback: a metal without field-dependent resistivity has no
    // B or beta sensitivity — no Kohler function required here, the
    // Metal override handles the field-dependent case
    inline real
    Material::drhodB( const real T, const real B, const real beta ) const
    {
        BELFEM_ASSERT( mType == MaterialType::PureMetal,
           "wrong function call for Material::drhodB in %s", mLabel.c_str() ) ;

        return 0. ;
    }

    inline real
    Material::drhodbeta( const real T, const real B, const real beta ) const
    {
        BELFEM_ASSERT( mType == MaterialType::PureMetal,
           "wrong function call for Material::drhodbeta in %s", mLabel.c_str() ) ;

        return 0. ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::H( const real B ) const
    {
       return (this->*mFunctionH ) ( B );
    }

    inline real
    Material::H_const( const real B ) const
    {
        return B * constant::nu0 ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::mu( real H, const real T  ) const
    {
        return ( this->*mFunctionMu )( H, T );
    }

    inline real
    Material::mu_const( real H, const real T ) const
    {
       return this->constant_property( MaterialProperty::mu );
    }

//------------------------------------------------------------------------------

    inline void
    Material::dmudH( const real H, real & mu, real & dmudH ) const
    {
       ( this->*mFunctionDMuDH )( H, mu, dmudH );
    }

    inline void
    Material::dmudH_const( const real H, real & mu, real & dmudH ) const
    {
       mu = this->constant_property( MaterialProperty::mu );
       dmudH = 0.0 ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::E( real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::E ) ,
           "Material %s does not have Youngs Modulus", mLabel.c_str() ) ;
       return ( this->*mFunctionE )( T );
    }

    inline real
    Material::E_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::E );
    }

    inline real
    Material::E_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::E, T );
    }

//------------------------------------------------------------------------------

    inline real
    Material::nu( real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(),
           "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::nu ) ,
           "Material does not have Poisson Ratio" ) ;

        return ( this->*mFunctionNu )( T );
    }

    inline real
    Material::nu_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::nu );
    }

    inline real
    Material::nu_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::nu, T );
    }

//------------------------------------------------------------------------------

    inline real
    Material::G( real T ) const
    {
        return this->E( T ) / ( 2.0 * ( 1.0 + this->nu( T ) ) ) ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::K( real T ) const
    {
        return this->E( T ) / ( 3.0 * ( 1.0 - 2.0 * this->nu( T ) ) ) ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::alpha( real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::alpha ) ,
           "Material %s does not have thermal expansion", mLabel.c_str() ) ;

        return ( this->*mFunctionAlpha )( T );
    }

    inline real
    Material::alpha_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::alpha );
    }

    inline real
    Material::alpha_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::alpha, T );
    }

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------

    inline real
    Material::Rp02( real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::Rp02 ) ,
           "Material %s does not have yield stress", mLabel.c_str() ) ;

       return ( this->*mFunctionRp02 )( T );
    }

    inline real
    Material::Rp02_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::Rp02 );
    }

    inline real
    Material::Rp02_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::Rp02, T );
    }

//------------------------------------------------------------------------------

    inline real
    Material::debye( real T ) const
    {
       BELFEM_ASSERT( this->is_isotropic(), "Material %s is not isotropic", mLabel.c_str() ) ;
       BELFEM_ASSERT( this->have( MaterialProperty::debye ) ,
           "Material %s does not have Debye temperature", mLabel.c_str() ) ;

       return ( this->*mFunctionDebye )( T );
    }

    inline real
    Material::debye_const( const real T ) const
    {
       return this->constant_property( MaterialProperty::debye );
    }

    inline real
    Material::debye_spline( const real T ) const
    {
       return this->spline_property( MaterialProperty::debye, T );
    }

    inline real
    Material::rho_i( const real T ) const
    {
        return ( this->*mFunctionRhoI )( T );
    }

    inline real
    Material::rho_i_spline( const real T ) const
    {
        return this->spline_property( MaterialProperty::rho_i, T );
    }

//------------------------------------------------------------------------------


    inline real
    Material::alpha_switch_temperature() const
    {
        BELFEM_ASSERT( ! std::isnan( mTAlphaSwitch ),
            "%s: alpha split temperature is not set, create_low_temperature_alpha() "
            "must run before alpha is evaluated", mLabel.c_str() );
        return mTAlphaSwitch ;
    }

//------------------------------------------------------------------------------

    inline real
    Material::volumetric_heatload( const real x, const real y, const real z, const real time ) const
    {
        BELFEM_ASSERT( mHeatFunction != nullptr , "No heat function set for material %s", mLabel.c_str() ) ;
        return this->mHeatFunction( x, y, z, time ) ;
    }

}

#include "powerlaws.hpp"

// The pop MUST mirror the push at the top of this header condition for
// condition. It used to test BELFEM_CLANG / BELFEM_GCC / BELFEM_INTEL, which
// are BELFEM build defines, while the push tests compiler builtins. A
// translation unit that includes this header without those defines -- every
// user-material plugin built from UserMaterialTemplate.cmake -- therefore
// pushed and never popped, silently disabling -Wformat, -Wunused-variable and
// -Wunused-parameter for the rest of the file.
#if defined(__clang__)
#pragma clang diagnostic pop

#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop

#elif defined(BELFEM_INTEL)
#pragma warning(pop)

#endif

#endif //BELFEM_CL_MATERIAL_HPP