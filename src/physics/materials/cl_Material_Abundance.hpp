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

#ifndef BELFEM_CL_MATERIAL_ABUNDANCE_HPP
#define BELFEM_CL_MATERIAL_ABUNDANCE_HPP

/**
 * @file cl_Material_Abundance.hpp
 * @brief Natural isotope abundance database for computing alloy properties
 *
 * This file provides the Abundance class for computing effective molar masses
 * and mass-difference impurity parameters (Γ) for alloys and composites.
 *
 * ## Physical Background
 *
 * The impurity parameter Γ quantifies phonon scattering from mass differences
 * between isotopes and alloying elements. It appears in the mass-difference
 * scattering rate τ_M^(-1) in thermal conductivity models (Callaway model).
 *
 * For an alloy with composition {X_i} (molar fractions):
 *
 *   Γ = Σ_i Σ_j X_i · A_j^i · [(M_j^i - M̄_i) / M̄]²
 *
 * where:
 * - X_i: Molar fraction of element i
 * - A_j^i: Natural abundance of isotope j of element i [%]
 * - M_j^i: Atomic mass of isotope j of element i [u]
 * - M̄_i: Average atomic mass of element i [u]
 * - M̄: Average molar mass of the alloy [u]
 *
 * ## Database Source
 *
 * All isotope masses and abundances are taken from the CRC Handbook of
 * Chemistry and Physics. The database covers elements H (Z=1) through Bi (Z=83).
 *
 * ## Usage
 *
 * @code
 * // Example: Hastelloy C-276 (Ni-Mo-Cr alloy)
 * Abundance db;
 * Cell<string> elements = {"Ni", "Mo", "Cr"};
 * Vector<real> mass_fractions = {0.0, 0.16, 0.155};  // fractions; the 0 marks the balance element (Ni)
 * real M, Gamma;
 * db.compute_molar_mass_and_impurity_from_masses(
 *     elements, mass_fractions, M, Gamma
 * );
 * // M ≈ 0.0593 kg/mol, Γ used in Callaway conductivity model
 * @endcode
 *
 * Reference: Zou & Balandin, J. Appl. Phys. 89, 2932 (2001) for τ_M formulation
 */

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace material
    {
        /**
         * @brief Database entry for an element's isotopes
         *
         * Stores parallel arrays:
         * - first: Isotope masses [atomic mass units]
         * - second: Natural abundances [%]
         */
        typedef std::pair< Vector< real >, Vector< real  > > AbundanceEntry ;

        /**
         * @brief Natural isotope abundance database for alloy property calculations
         *
         * Provides isotope mass and abundance data for computing effective
         * molar masses and mass-difference impurity parameters (Γ) used in
         * phonon scattering models.
         *
         * The database contains isotope data for all stable elements from
         * H (Z=1) through Bi (Z=83), sourced from the CRC Handbook of
         * Chemistry and Physics.
         *
         * ## Applications
         *
         * 1. **Molar Mass Calculation**: M̄ = Σ_i X_i · M̄_i
         *    - Used for density calculations
         *    - Required for Callaway thermal conductivity
         *
         * 2. **Impurity Parameter (Γ)**: Quantifies mass-difference scattering
         *    - Appears in τ_M^(-1) ∝ ω⁴ · Γ (Callaway model)
         *    - Accounts for both isotope and alloy scattering
         *
         * ## Usage Patterns
         *
         * For single elements:
         * @code
         * Abundance db;
         * real M_Cu = db.compute_molar_mass("Cu");  // 0.06355 kg/mol
         * @endcode
         *
         * For alloys with mass fractions:
         * @code
         * Cell<string> elements = {"Ni", "Cr", "Fe"};
         * Vector<real> wt_percent = {70.0, 20.0, 10.0};  // Can sum to 100
         * real M, Gamma;
         * db.compute_molar_mass_and_impurity_from_masses(
         *     elements, wt_percent, M, Gamma
         * );
         * @endcode
         *
         * For composites with volume fractions:
         * @code
         * Cell<string> elements = {"Y", "Ba", "Cu", "O"};
         * Vector<real> vol_fractions = {1.0, 2.0, 3.0, 7.0};  // YBCO stoichiometry ( not normalized )
         * real M, Gamma;
         * db.compute_molar_mass_and_impurity_from_volumes(
         *     elements, vol_fractions, M, Gamma
         * );
         * @endcode
         */
        class Abundance
        {
            /**
             * @brief Internal database mapping element symbols to isotope data
             *
             * Key: Element symbol (e.g., "Cu", "Ni", "O")
             * Value: AbundanceEntry with masses [u] and abundances [%]
             */
            Map< string, AbundanceEntry > mDatabase ;

        public:

            /**
             * @brief Constructor - initializes isotope database
             *
             * Loads natural isotope masses and abundances for all elements
             * from H (Z=1) through Bi (Z=83) from CRC Handbook data.
             */
            Abundance();

            /**
             * @brief Destructor
             */
            ~Abundance() = default;

            /**
             * @brief Compute molar mass and impurity parameter from volume fractions
             * @param[in] aElements Element symbols (e.g., {"Ni", "Cr", "Fe"})
             * @param[in] aVolumeFractions Stoichiometric counts per formula unit
             *            ( NOT normalized: the returned molar mass is Σ X_i·M̄_i,
             *            i.e. the formula-unit mass when the counts are integers )
             * @param[out] aMolarMass Average molar mass [kg/mol]
             * @param[out] aImpurity Mass-difference impurity parameter Γ [-]
             *
             * Computes effective properties for composites where volume fractions
             * are known (e.g., from stoichiometry or microstructure).
             *
             * Algorithm:
             * 1. Compute M̄_i for each element from isotope abundances
             * 2. Compute alloy molar mass: M̄ = Σ_i X_i · M̄_i
             * 3. Compute Γ = Σ_i Σ_j X_i · A_j^i · [(M_j^i - M̄_i) / M̄]²
             *
             * Example (YBCO: Y₁Ba₂Cu₃O₇):
             * @code
             * Cell<string> elem = {"Y", "Ba", "Cu", "O"};
             * Vector<real> vol = {1.0, 2.0, 3.0, 7.0};  // Stoichiometric ratios
             * real M, Gamma;
             * db.compute_molar_mass_and_impurity_from_volumes(elem, vol, M, Gamma);
             * // M ≈ 0.6658 kg/mol (for YBa₂Cu₃O₇)
             * @endcode
             *
             * Used by: the YBCO constructor ( cl_Material_YBCO.cpp )
             */
            void
            compute_molar_mass_and_impurity_from_volumes(
                const Cell< string > & aElements,
                const Vector< real > & aVolumeFractions,
                  real & aMolarMass,
                  real & aImpurity ) const;

            /**
             * @brief Compute molar mass and impurity parameter from mass fractions
             * @param[in] aElements Element symbols (e.g., {"Ni", "Mo", "Cr"})
             * @param[in] aMassFractions Mass fractions [wt%] or [fraction]
             *            Can sum to 100 (wt%), 1.0 (fraction), or be unbalanced
             * @param[out] aMolarMass Average molar mass [kg/mol]
             * @param[out] aImpurity Mass-difference impurity parameter Γ [-]
             *
             * Computes effective properties for alloys where weight percentages
             * are known from composition specifications.
             *
             * Input handling:
             * - If Σ(wt%) ≈ 100: Treats as percentages, normalizes to fractions
             * - If Σ(frac) ≈ 1.0: Treats as fractions
             * - If unbalanced: Exactly one element must be zero (balance element)
             *
             * Algorithm:
             * 1. Convert mass fractions Y_i to molar fractions: X_i = Y_i / M̄_i
             * 2. Normalize: X_i ← X_i / Σ(X_i)
             * 3. Compute M̄ = Σ_i X_i · M̄_i
             * 4. Compute Γ = Σ_i Σ_j X_i · A_j^i · [(M_j^i - M̄_i) / M̄]²
             *
             * Example (Hastelloy C-276):
             * @code
             * Cell<string> elem = {"Ni", "Mo", "Cr", "Fe", "W"};
             * Vector<real> wt = {0.0, 0.16, 0.155, 0.055, 0.04};  // fractions, Ni = balance
             * // percent lists must sum to exactly 100; an unbalanced list needs exactly one 0
             * real M, Gamma;
             * db.compute_molar_mass_and_impurity_from_masses(elem, wt, M, Gamma);
             * @endcode
             *
             * Used by: HastelloyC276::set_constants()
             */
            void
            compute_molar_mass_and_impurity_from_masses(
               const Cell< string > & aElements,
               const Vector< real > & aMassFractions,
                 real & aMolarMass,
                 real & aImpurity ) const;

            /**
             * @brief Get molar mass of a single element
             * @param aElement Element symbol (e.g., "Cu", "Ni", "O")
             * @return Molar mass [kg/mol]
             *
             * Returns the average molar mass weighted by natural isotope
             * abundances: M̄ = Σ_j A_j · M_j
             *
             * Example:
             * @code
             * Abundance db;
             * real M_Cu = db.compute_molar_mass("Cu");
             * // Returns 0.063546 kg/mol (natural Cu: 69.17% ⁶³Cu, 30.83% ⁶⁵Cu)
             * @endcode
             *
             * Used by: Material constructors for pure elements
             */
            real
            compute_molar_mass(
                const string & aElement ) const;

        };
    }
}
#endif //BELFEM_CL_MATERIAL_ABUNDANCE_HPP