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

#ifndef CL_MATERIALFACTORY_HPP
#define CL_MATERIALFACTORY_HPP

#include "cl_Material.hpp"
#include "cl_BhCurve.hpp"
#include "cl_JcFunction.hpp"
#include "cl_Input_Section.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Factory class for creating material objects and associated property functions
     *
     * The MaterialFactory provides a centralized interface for creating:
     * - Material objects (Copper, Silver, YBCO, etc.)
     * - B-H curves for ferromagnetic materials
     * - Critical current density (Jc) functions for superconductors
     * - Power law exponent (n) functions for superconductors
     *
     * IMPORTANT USAGE NOTES:
     *
     * 1. OWNERSHIP: When you assign a B-H curve, Jc function, or n function to a
     *    material using set_bh_curve(), set_jc_function(), or set_n_function(),
     *    the material takes OWNERSHIP and will delete the object in its destructor.
     *    DO NOT manually delete these objects after assignment.
     *
     * 2. MANUAL ASSIGNMENT with the create_*() methods: they only create the
     *    object; call the material's load_bh_curve() / set_jc_function() /
     *    set_n_function() yourself. The input-section constructor does this
     *    wiring for you.
     *
     * 3. CONSTANT PROPERTIES: If Jc or n are constant values, no function object
     *    should be created. Instead, the constant value is stored directly in the
     *    material class. Only create functions when properties vary with field,
     *    angle, or temperature.
     *
     * Example usage:
     * @code
     * MaterialFactory factory;
     *
     * // Create a copper material
     * Material* copper = factory.create_material("copper");
     * copper->set_RRR(2000);
     *
     * // Create HTS with Jc function
     * Material* ybco = factory.create_material("YBCO");
     * material::JcFunction* jc = factory.create_jc_function(1e9, 5.0, 0.5, 2.0);
     * ybco->set_jc_function(jc);  // YBCO now owns jc, will delete it
     *
     * // Create ferromagnetic material with B-H curve
     * Material* ferro = factory.create_material("SomeFerro");
     * material::BhCurve* bh = factory.create_bh_curve("bhdata.txt", "Ferro");
     * ferro->load_bh_curve(bh);  // Ferro now owns bh and routes mu/H/dmudH through it
     *                            // ( set_bh_curve() alone only stores the pointer )
     * @endcode
     *
     * @ingroup grp_physics_materials
     * @see @ref physics_materials_materials_usage_guide
     */
    class MaterialFactory
    {

        Map < string ,Material * > mMaterialsMap ; //the material map, with labels as keys
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * @brief Default constructor
         */
        MaterialFactory() = default;

        /**
         * @brief Constructor with input section
         */
        MaterialFactory( const input::Section * aSection );

        /**
         * @brief Default destructor
         */
        ~MaterialFactory() = default;

//------------------------------------------------------------------------------

        /**
         * @brief Create a material by label
         *
         * Creates a material object from a predefined material database.
         * The label is case-insensitive.
         *
         * Currently supported labels ( aliases in parentheses ):
         * - aluminum (aluminium, al), chromium (cr), iron (ferro, fe),
         *   nickel (ni), copper (cu), silver (ag), indium (in), tin (sn),
         *   lead (pb)                    : PureMetal
         * - hastelloyc276 (hastelloy)    : Hastelloy C-276 nickel alloy (LookupAlloy)
         * - ybco                         : YBCO high-temperature superconductor (HTS)
         * - magnesia (mgo, buffer)       : NonMetal
         * - solder compositions of the form <Element><percentage>..., e.g.
         *   Sn40Pb60                     : class Alloy ( typed PureMetal )
         * See print_material_list().
         *
         * @param aLabel Material label (case-insensitive)
         * @param aRRR Residual resistivity ratio; NaN keeps the material default
         * @param aBuildTables Build the property lookup tables on construction
         * @return Pointer to created material (caller owns the material)
         * @throws Error if material label is not recognized
         *
         * Example:
         * @code
         * MaterialFactory factory;
         * Material* cu = factory.create_material("copper");
         * Material* ag = factory.create_material("SILVER");  // case-insensitive
         * @endcode
         */
        Material *
        create_material( const string & aLabel,
                         const real aRRR = BELFEM_QUIET_NAN,
                         const bool aBuildTables = true );

//------------------------------------------------------------------------------

        void
        print_material_list( std::ostream & aStream );

//------------------------------------------------------------------------------

        /**
         * @brief Create a user-defined material from shared library
         *
         * Loads a user-defined material from an external shared library (.so, .dylib, .dll).
         * The library must contain an initialization function with the signature:
         * @code
         * extern "C" void \<aLabel\>_init(Material* mat);
         * @endcode
         *
         * Within the init function, users define material properties using:
         * - mat->set_constant() - For constant values
         * - mat->set_user_defined_function() - For custom functions
         * - mat->set_user_defined_polynomial() - For polynomial functions
         * - mat->set_jc_function() / mat->set_n_function() - For superconductor properties
         *
         * @param aLibraryPath Path to shared library file
         * @param aLabel Material label (must match the init function name: \<aLabel\>_init)
         * @return Pointer to created user-defined material (caller owns the material)
         * @throws Error if library cannot be loaded or init function not found
         *
         * Example user library (myalloy.cpp):
         * @code
         * #include "cl_Material.hpp"
         * using namespace belfem;
         *
         * real my_rho(const Material* mat, real T) {
         *     return 1.7e-8 * (1.0 + 0.004 * (T - 293.0));
         * }
         *
         * extern "C" void MyAlloy_init(Material* mat) {
         *     mat->set_constant(MaterialProperty::E, 200e9);
         *     mat->set_constant(MaterialProperty::nu, 0.3);
         *     mat->set_user_defined_function(MaterialProperty::rho,
         *                                     MaterialDependency::T,
         *                                     &my_rho);
         *     // cp(T) = 0.12*T + 385  (descending order: T¹, T⁰)
         *     std::vector<real> cp_coeffs = {0.12, 385.0};
         *     mat->set_user_defined_polynomial(MaterialProperty::cp, cp_coeffs);
         * }
         * @endcode
         *
         * Compile: g++ -shared -fPIC myalloy.cpp -o libmyalloy.so -I/path/to/belfem/include
         *
         * The library path is resolved through material::data_file(): the run
         * directory first, then $BELFEM_DATA/material. A name that matches
         * nothing is still passed to dlopen, which searches $LD_LIBRARY_PATH.
         *
         * Usage:
         * @code
         * MaterialFactory factory;
         * Material* mat = factory.create_material("libmyalloy.so", "MyAlloy");
         * @endcode
         */
        Material *
        create_material( const string & aLibraryPath, const string & aLabel );

//------------------------------------------------------------------------------

        Map <string, Material *> &
        materials();

//------------------------------------------------------------------------------

        /**
         * @brief Create a B-H curve from file for ferromagnetic materials
         *
         * Creates a B-H curve object by loading data from a file. This curve
         * defines the relationship between magnetic flux density (B) and
         * magnetic field strength (H) for ferromagnetic materials.
         *
         * IMPORTANT: This function only CREATES the B-H curve. You must manually
         * assign it to a material using material->load_bh_curve(). Once assigned,
         * the material takes ownership and will delete the curve.
         *
         * @param aPath  Path to the B-H curve data file, resolved through
         *               material::data_file(): the run directory first, then
         *               $BELFEM_DATA/material
         * @param aLabel Identifier label for this B-H curve
         * @return Pointer to created B-H curve (ownership transfers to material upon assignment)
         *
         * Example:
         * @code
         * MaterialFactory factory;
         * Material* ferro = factory.create_material("SomeFerro");
         *
         * // Create B-H curve
         * material::BhCurve* bh = factory.create_bh_curve("iron_bh.dat", "Iron");
         *
         * // Assign to material (material now owns bh and uses it for mu/H)
         * ferro->load_bh_curve(bh);
         * @endcode
         */
        material::BhCurve *
        create_bh_curve( const string & aPath, const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * @brief Create a Jc function using the modified Kim analytical model
         *
         * Creates a critical current density function for superconductors using
         * the modified Kim model. This is an analytical model suitable for
         * quick calculations when experimental data is not available.
         *
         * The modified Kim model expresses Jc as a function of magnetic field.
         *
         * IMPORTANT: This function only CREATES the Jc function. You must manually
         * assign it to a material using material->set_jc_function() or
         * material->set_n_function(). Once assigned, the material takes ownership
         * and will delete the function.
         *
         * @param aJc0 Critical current density at zero field [A/m²]
         * @param aB   Characteristic field scale B0 [T]
         * @param aBc  Anisotropy parameter k² (dimensionless)
         * @param aK   Field-dependence exponent α (dimensionless)
         * @return Pointer to created Jc function (ownership transfers to material upon assignment)
         *
         * Example:
         * @code
         * MaterialFactory factory;
         * Material* ybco = factory.create_material("YBCO");
         *
         * // Create Jc function using modified Kim model
         * real jc0 = 1e9;    // A/m²
         * real B0 = 5.0;     // T
         * real k2 = 0.5;     // anisotropy k², dimensionless
         * real alpha = 2.0;  // exponent
         * material::JcFunction* jc = factory.create_jc_function(jc0, B0, k2, alpha);
         *
         * // Assign to material (material now owns jc)
         * ybco->set_jc_function(jc);
         * @endcode
         */
        material::JcFunction *
        create_jc_function( const real aJc0, const real aB, const real aBc, const real aK );

//------------------------------------------------------------------------------

        /**
         * @brief Create a Jc (or n) function from database file
         *
         * Creates a critical current density function (or power law exponent function)
         * by loading experimental data from a file. This provides more accurate
         * results than analytical models when experimental data is available.
         *
         * The function can represent:
         * - Jc(B, angle, T) : Critical current density
         * - n(B, angle, T)  : Power law exponent
         *
         * IMPORTANT: This function only CREATES the function object. You must manually
         * assign it to a material using material->set_jc_function() or
         * material->set_n_function(). Once assigned, the material takes ownership
         * and will delete the function.
         *
         * NOTE: If Jc or n are constant, do NOT create a function. Instead, store
         * the constant value directly in the material using set_constant().
         *
         * @param aPath  Path to the database file containing Jc or n data,
         *               resolved through material::data_file(): the run
         *               directory first, then $BELFEM_DATA/material
         * @param aLabel Identifier label for this function
         * @return Pointer to created function (ownership transfers to material upon assignment)
         *
         * Example:
         * @code
         * MaterialFactory factory;
         * Material* ybco = factory.create_material("YBCO");
         *
         * // Create Jc function from database
         * material::JcFunction* jc = factory.create_jc_function("ybco.hdf5", "jc");
         * ybco->set_jc_function(jc);  // YBCO now owns jc
         *
         * // Create n function from database
         * material::JcFunction* n = factory.create_jc_function("ybco.hdf5", "n");
         * ybco->set_n_function(n);    // YBCO now owns n
         * @endcode
         */
        material::JcFunction *
        create_jc_function( const string & aPath, const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * refuse any key or subsection the resolved material shape does not
         * read
         *
         * A material section selects exactly one shape -- a b-h curve, a
         * builtin, or a plugin -- and each shape consumes a different set of
         * keys. Everything else used to be dropped without a word, so a
         * misspelled key, a key belonging to another shape, or a constant
         * that lost to a file all behaved as if the deck had never mentioned
         * them.
         *
         * This is an ALLOW-LIST rather than a record of what happened to be
         * read, because the two are not the same predicate: `RRR` on a ybco
         * IS read here and then dropped by the constructor, so only an
         * explicit per-shape contract catches it.
         *
         * @param aSection   the material section, or one of its subsections
         * @param aLabel     material name, for the error message
         * @param aShape     shape name, for the error message
         * @param aKeys      keys this shape reads
         * @param aSections  subsection types this shape enters
         */
        void
        check_unused_input(
                const input::Section * aSection,
                const string         & aLabel,
                const string         & aShape,
                const Cell< string > & aKeys,
                const Cell< string > & aSections );

//------------------------------------------------------------------------------
    };

    inline Map <string, Material *> &
    MaterialFactory::materials()
    {
        return mMaterialsMap ;
    }

//------------------------------------------------------------------------------
}
#endif //CL_MATERIALFACTORY_HPP
