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

#ifndef BELFEM_DEBYE_HPP
#define BELFEM_DEBYE_HPP

/**
 * @file debye.hpp
 * @brief Fortran interface for Debye integrals and Callaway thermal conductivity
 *
 * This header provides C++ access to high-precision Fortran implementations
 * of computationally intensive functions used in material property calculations.
 *
 * Functions:
 * - debye_table: Computes Debye integral J_n(Z) for electrical resistivity
 * - callaway_conductivity: Computes phonon thermal conductivity
 *
 * Implementation language: Fortran (for numerical stability and performance)
 */

#ifdef __cplusplus
extern"C"
{
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * @brief Compute Debye integral lookup table
     *
     * Computes J_n(Z) = ∫₀^Z x^n · exp(x) / (exp(x) - 1)² dx
     *
     * Used for Bloch-Grüneisen electrical resistivity model:
     * ρᵢ(T) = C · J_n(θ/T) / (θ/T)^n
     *
     * The integral is computed using adaptive Gauss-Legendre quadrature
     * with variable step sizes to capture both the sharp peak and the
     * exponential tail. The peak location is determined analytically
     * via Newton-Raphson (z_peak ≈ n for large n).
     *
     * @param[in] n Exponent of the integrand ( any real ≥ 3; the metals use
     *              3, 4, 5, Iron uses 4.5 )
     *              Peak locations: J₃ at 2.576, J₄ at 3.830, J₅ at 4.928
     * @param[in] m Number of integration points per interval
     *              Recommended: 8-16 for good accuracy
     * @param[in] p Number of points to peak value (set to ~20)
     *              Controls fine mesh density near peak
     * @param[in] q Total number of points in lookup table
     *              Recommended: n=3: 330, n=4: 237, n=5: 196
     * @param[out] z Z-values (abscissas) for lookup table [dimensionless]
     * @param[out] y J_n(z) integrated values (ordinates) [dimensionless]
     *
     * Grid construction:
     * - Constant step size dz = z_peak / (p - 1)
     * - [0, z_peak]: p points with uniform spacing
     * - [z_peak, ∞]: q-p points with uniform spacing (captures tail)
     *
     * The integration uses quadruple precision internally for numerical
     * stability, then converts to double precision for output.
     *
     * Used by: Metal::rho_i_custom() ( Bloch–Grüneisen intrinsic resistivity )
     *          via the J₃/J₄/J₅/Jₙ splines built in Metal::create_J_spline()
     */
    void
    debye_table(
        const double * n,
        const int * m,
        const int * p,
        const int * q,
        double * z,
        double * y
        );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    /**
     * @brief Evaluate phonon thermal conductivity using Callaway model
     *
     * Computes λ_ph(T) using the Callaway model with multiple scattering
     * mechanisms combined via Matthiessen's rule:
     * 1/τ_total = 1/τ_U + 1/τ_M + 1/τ_B + 1/τ_ph-e
     *
     * Scattering mechanisms:
     * 1. Umklapp (phonon-phonon): τ_U⁻¹ ∝ ω² · T · exp(-θ/(b·(1+d·T/θ)·T))
     * 2. Mass difference impurity: τ_M⁻¹ ∝ ω⁴ · Γ
     * 3. Boundary scattering: τ_B⁻¹ = v_g / L₀
     * 4. Phonon-electron scattering: τ_ph-e⁻¹ ∝ ω (acoustic or optical)
     * 5. Superconducting gap reduction (d-wave for YBCO)
     *
     * The integral computed is:
     * λ_ph = (kB/(2π²v)) · (kB·T/ℏ)³ · ∫₀^(θ/T) τ(z)·z^n·exp(z)/[exp(z)-1]² dz
     *
     * Reference: J. Callaway, Phys. Rev. 113, 1046 (1959)
     * Scattering: Zou & Balandin, J. Appl. Phys. 89, 2932 (2001)
     *
     * @param[in] params Array of 20 parameters (see below)
     * @param[out] result Phonon thermal conductivity λ_ph [W/(m·K)]
     * @param[out] status Error status (0 = success, 1 = error)
     *
     * Parameter array layout (C++ zero-based / Fortran one-based):
     *
     * Temperature-dependent material properties:
     * - params[0]  / params(1):  T       - Temperature [K]
     * - params[1]  / params(2):  θ       - Debye temperature [K]
     * - params[2]  / params(3):  v_g     - Group velocity (phonon velocity) [m/s]
     * - params[3]  / params(4):  ρ       - Density [kg/m³], used for τ_ph-e
     * - params[4]  / params(5):  G       - Shear modulus [Pa], used for τ_U
     * - params[5]  / params(6):  γ       - Grüneisen parameter [-], used for τ_U
     *
     * Fixed material properties (composition-based):
     * - params[6]  / params(7):  M       - Molar mass [kg/mol]
     * - params[7]  / params(8):  Γ       - Impurity parameter [-], used for τ_M
     * - params[8]  / params(9):  T_c     - Critical temperature [K] (0 if not superconducting)
     * - params[9]  / params(10): L₀      - Thickness of REBCO layer [m], used for τ_B
     *
     * Model parameters:
     * - params[10] / params(11): n       - Exponent for the function (usually n = 4)
     * - params[11] / params(12): b       - Umklapp parameter [-], typically 1.5-3.0 (ignored if 0)
     * - params[12] / params(13): d       - Umklapp temperature correction [-]
     * - params[13] / params(14): ε       - Deformation potential [eV], used for τ_ph-e
     * - params[14] / params(15): m_e,eff - Effective electron mass [factor × m_e], used for τ_ph-e
     * - params[15] / params(16): n_e     - Conduction electron concentration [1/m³]
     * - params[16] / params(17): Δ₀/(k_B·T_c) - D-wave gap ratio [-] (0 if no SC reduction)
     * - params[17] / params(18): λ_pe    - Acoustic electron-phonon coupling [-]
     * - params[18] / params(19): λ_opt   - Optical electron-phonon coupling [-]
     * - params[19] / params(20): ω_opt   - Optical band frequency [1/cm]
     *
     * IMPORTANT NOTES:
     * - Cannot use both λ_opt and λ_pe (returns status=1 if both nonzero)
     * - Acoustic model: τ_ph-e⁻¹ based on deformation potential ε
     * - Optical model: τ_ph-e⁻¹ based on Raman frequency ω_opt
     * - If T < T_c and Δ₀ > 0 and an e-ph model is active: applies d-wave SC gap reduction
     * Typical values for YBCO:
     * - θ ≈ 275 K, v_g ≈ 3000 m/s, L₀ ≈ 1 μm, n = 4
     * - b ≈ 2.0, T_c ≈ 92 K, Δ₀/(k_B·T_c) ≈ 2.2, λ_opt ≈ 0.1
     *
     * Used by: YBCO::lambda_custom() (override of Material_Metal::lambda_custom);
     *          also YBCO::test_callaway()
     */
    void
    callaway_conductivity (
        const double * params,
              double * result,
              int    * status );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef __cplusplus
}
#endif

#endif //BELFEM_DEBYE_HPP
