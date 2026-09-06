/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/**
 * \file powerlaws.hpp
 * \brief Inline E–J resistivity models for HTS materials and their Jacobians.
 *
 * This header provides three resistivity models that an HTS-aware
 * Material may evaluate, together with their analytic derivatives for
 * use in the nonlinear (Newton–Raphson) solver:
 *
 *  - Material::rho_powerlaw():   standard E–J power law parallel-combined
 *                                with the normal-state resistivity;
 *  - Material::rho_piecewise():  three-regime extension that smoothly
 *                                transitions from the power law into the
 *                                flux-flow / normal regime through a
 *                                quadratic Bezier curve in log–log space;
 *  - Material::rho_riva():       the same parallel model as rho_powerlaw(),
 *                                made total over the full jc/n table range
 *                                ( Riva 2021 ); selected via ResistivityLaw.
 *
 * Each model is offered as a family of overloads that differ only in how
 * the critical current density \f$J_c\f$ and exponent \f$n\f$ are obtained
 * (plain constants, function of \f$(B,\angle)\f$, function of
 * \f$(B,\angle,T)\f$, or user-supplied temperature-only callbacks) and in
 * whether the spatial defect modulation \f$d(x,y,z,t)\f$ is applied to
 * \f$J_c\f$. The math is identical across overloads; only the source of
 * \f$J_c\f$ and \f$n\f$ changes.
 *
 * \section powerlaws_refs References
 *  - J. Rhyner, "Magnetic properties and AC-losses of superconductors with
 *    power-law current-voltage characteristics," Physica C 212 (1993)
 *    pp. 292–300. — original derivation of the E–J power-law model used
 *    throughout the HTS community.
 *  - C. P. Plummer and J. E. Evetts, "Dependence of the shape of the
 *    resistive transition on composite inhomogeneity in multifilamentary
 *    wires," IEEE Trans. Magn. 23 (2) (1987) pp. 1179–1182. — phenomenology
 *    of the n-value characterizing the transition sharpness.
 *  - C. Messe et al., "BELFEM: a finite-element framework for HTS magnet
 *    quench analysis," Supercond. Sci. Technol. 36 (2023) 114001. — BELFEM
 *    material database (Sec. 2.6) and the convergence strategy that drives
 *    the residual below \f$10^{-11}\f$ (Sec. 2.7) for which these
 *    derivatives are required.
 */

#ifndef BELFEM_POWERLAWS_HPP
#define BELFEM_POWERLAWS_HPP

#include "cl_Material.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
// Power-law resistivity model
//------------------------------------------------------------------------------

    inline real
    Material::jc_eval( const real T, const real normB, const real angleNxB ) const
    {
        return ( mJcFunction != nullptr )
            ? ( mJcFunction->depends_on( material::JcParameter::T )
                ? mJcFunction->eval( normB, angleNxB, T )
                : mJcFunction->eval( normB, angleNxB ) )
            : this->constant_property( MaterialProperty::jc ) ;
    }

    inline real
    Material::n_eval_raw( const real T, const real normB, const real angleNxB ) const
    {
        return ( mNFunction != nullptr )
            ? ( mNFunction->depends_on( material::JcParameter::T )
                ? mNFunction->eval( normB, angleNxB, T )
                : mNFunction->eval( normB, angleNxB ) )
            : this->constant_property( MaterialProperty::n ) ;
    }

    inline real
    Material::n_eval( const real T, const real normB, const real angleNxB ) const
    {
        // ohmic floor ( 2026-08-27 ): measured tables soften through n = 1
        // near T_crit; below that the raw power law is sub-ohmic and its
        // J -> 0 limit flips. At n = 1 the tape is a plain resistor
        // ec/jc. dn_eval_dB / dn_eval_dT return 0 while the floor binds,
        // so the tangents differentiate the SAME clamped law
        return std::max( this->n_eval_raw( T, normB, angleNxB ), 1.0 ) ;
    }

    // derivative routing mirrors jc_eval / n_eval: with no function attached
    // the value is a constant and its derivative is EXACTLY zero — which also
    // makes constant-jc decks bit-identical, since the Newton consumer
    // early-outs on zero ( audited 2026-08-13 )

    inline real
    Material::djc_eval_dB( const real T, const real normB, const real angleNxB ) const
    {
        return ( mJcFunction != nullptr )
            ? mJcFunction->deval_dB( normB, angleNxB, T )
            : 0.0 ;
    }

    inline real
    Material::dn_eval_dB( const real T, const real normB, const real angleNxB ) const
    {
        if ( mNFunction == nullptr ) return 0.0 ;

        // n_eval clamps at 1; the clamped law is locally constant there
        if ( ! ( this->n_eval_raw( T, normB, angleNxB ) > 1.0 ) ) return 0.0 ;

        return mNFunction->deval_dB( normB, angleNxB, T ) ;
    }

    // T-leg twins ( audited 2026-08-13 ): same null-check
    // routing. Exact zero for constants, ModifiedKim and 2-arg functions;
    // for a 3-arg ( T-dependent ) UserDefined WITHOUT a deval_dT override
    // the base-class 0.0 is a conservative fallback, not an exact derivative

    inline real
    Material::djc_eval_dT( const real T, const real normB, const real angleNxB ) const
    {
        return ( mJcFunction != nullptr )
            ? mJcFunction->deval_dT( normB, angleNxB, T )
            : 0.0 ;
    }

    inline real
    Material::dn_eval_dT( const real T, const real normB, const real angleNxB ) const
    {
        if ( mNFunction == nullptr ) return 0.0 ;

        // n_eval clamps at 1; the clamped law is locally constant there
        if ( ! ( this->n_eval_raw( T, normB, angleNxB ) > 1.0 ) ) return 0.0 ;

        return mNFunction->deval_dT( normB, angleNxB, T ) ;
    }

    /**
     * \brief Effective resistivity from the E–J power law (Rhyner 1993),
     *        parallel-combined with the normal-state resistivity.
     *
     * The intrinsic HTS power law is
     * \f[
     *     \rho_{PL}(J) \;=\; \frac{E_c}{J_c}\left(\frac{|J|}{J_c}\right)^{n-1},
     * \f]
     * where \f$E_c\f$ is the critical-field criterion (typically
     * \f$10^{-4}\f$ V/m), \f$J_c\f$ is the critical current density and
     * \f$n\f$ is the resistive-transition exponent
     * (\f$n \rightarrow 1\f$ → linear, \f$n \rightarrow \infty\f$ → ideal
     * critical-state).
     *
     * To bound \f$\rho\f$ from above as \f$J \gg J_c\f$ and from below as
     * \f$J \rightarrow 0\f$, the effective resistivity is the parallel
     * combination of \f$\rho_{PL}\f$ with the normal-state resistivity
     * \f$\rho_n(T)\f$:
     * \f[
     *     \rho_{eff}(J) \;=\;
     *       \left(\frac{1}{\rho_n} + \frac{1}{\rho_{PL}}\right)^{-1}.
     * \f]
     * \f$\rho_{PL}\f$ is floored at mRhoMin, which is ZERO by default
     * (2026-08-10): a positive floor desynchronizes the value from
     * drho_powerlaw_dJ, which differentiates the unfloored law — with the
     * old 1e-16 the Newton tangent was wrong for all
     * \f$J \lesssim 0.87\,J_c\f$ on typical tape constants. Anyone raising
     * mRhoMin again must also zero the derivative while the floor binds.
     *
     * In this overload \f$J_c\f$ and \f$n\f$ are read as plain constants
     * (configured via Material::set_constant) and \f$\rho_n\f$ is evaluated
     * at the global bulk temperature gTbulk.
     *
     * \param normJ  current-density magnitude \f$|J|\f$ [A/m²]
     * \return       effective electric resistivity [Ω·m]
     *
     * \see Rhyner, Physica C 212 (1993) 292–300.
     * \see Messe et al., Supercond. Sci. Technol. 36 (2023) 114001, Sec. 2.6.
     */
    inline real
    Material::rho_powerlaw( const real normJ  ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Constant \f$J_c\f$ / \f$n\f$ with spatially-dependent defect modulation
     * \f$J_c \rightarrow J_c \cdot d(x,y,z,t)\f$.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real x, const real y, const real z, const real t  ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc )*((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Field-dependent critical current: \f$J_c = J_c(|B|, \angle(n,B))\f$
     * supplied via mJcFunction. \f$n\f$ remains a constant.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real normB, const real angleNxB  ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;


        BELFEM_ASSERT(
            mJcFunction->depends_on( material::JcParameter::normB ) &&
            mJcFunction->depends_on( material::JcParameter::angleNxB ) &&
            ! mJcFunction->depends_on( material::JcParameter::T ), "wrong powerlaw for material %s", mLabel.c_str() );


        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ with spatial defect modulation.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;


        BELFEM_ASSERT(
            mJcFunction->depends_on( material::JcParameter::normB ) &&
            mJcFunction->depends_on( material::JcParameter::angleNxB ) &&
            ! mJcFunction->depends_on( material::JcParameter::T ), "wrong powerlaw for material %s", mLabel.c_str() );


        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB )*((this->mDefectFunction) (x,y,z,t)) ; ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Fully field- and temperature-dependent variant:
     * \f$J_c = J_c(|B|, \angle, T)\f$ and \f$n = n(|B|, \angle, T)\f$,
     * both supplied via mJcFunction / mNFunction. \f$\rho_n\f$ is evaluated
     * at the local temperature \p T (not gTbulk).
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real T, const real normB, const real angleNxB  ) const
    {
        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant with spatial defect
     * modulation applied to \f$J_c\f$.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const
    {
        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Temperature-only variant. \f$J_c(T)\f$ and \f$n(T)\f$ are obtained
     * from the user-supplied callbacks Material::jc_custom() and
     * Material::n_custom(), allowing materials with bespoke analytic
     * \f$T\f$-laws (e.g. Kim-style fits) without going through a
     * JcFunction.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real T  ) const
    {

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) ;
        real  n = this->n_custom( T ) ;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

    /**
     * \overload
     * Temperature-only custom-callback variant with spatial defect modulation.
     */
    inline real
    Material::rho_powerlaw( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const
    {
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_custom( T ) ;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;
    }

//------------------------------------------------------------------------------
// Three-regime piecewise resistivity model
//------------------------------------------------------------------------------

    /**
     * \brief Three-regime resistivity (power law / flux-flow / normal) with a
     *        smooth Bezier transition in log–log space.
     *
     * The pure power law diverges as \f$J \rightarrow \infty\f$, which is
     * unphysical above the critical regime: an HTS conductor saturates at
     * the resistivity of its normal-state matrix once the superconductor is
     * fully driven. This model splits the \f$J\f$-axis into three regimes
     * connected smoothly in \f$(\log_{10} J,\,\log_{10}\rho)\f$ space:
     *
     *  - **Power-law regime** \f$(J \le J_1)\f$: standard
     *    \f$\rho_{PL}(J) = (E_c/J_c)\,(J/J_c)^{n-1}\f$.
     *  - **Flux-flow transition** \f$(J_1 < J \le J_3)\f$: a quadratic
     *    Bezier curve in \f$(\log_{10} J,\,\log_{10}\rho)\f$ connecting
     *    \f$(\log J_1, \log \rho_1)\f$ to \f$(\log J_3, \log \rho_n)\f$
     *    with the middle control point at \f$(\log J_2, \log \rho_n)\f$.
     *  - **Normal regime** \f$(J > J_3)\f$: \f$\rho = \rho_n(T)\f$.
     *
     * The transition currents are
     * \f[
     *     J_1 \;=\; J_c \cdot 10^{D/n}, \qquad
     *     J_2 \;=\; J_1 \cdot \left(\rho_n / \rho_1\right)^{1/(n-1)}, \qquad
     *     J_3 \;=\; J_1 \cdot \left(\rho_n / \rho_1\right)^{1/N_{ff}},
     * \f]
     * where the knee is placed a fixed \f$D = 2.5\f$ decades above
     * \f$J_c\f$ ( hardcoded; the member mD is not consulted ) and mNff
     * (default 3) controls the slope of the flux-flow regime. Above the assumed critical
     * temperature \f$T_{crit}\f$, only \f$\rho_n\f$ is returned.
     *
     * The Bezier parameter \f$t \in [0,1]\f$ for a queried \f$J\f$ is the
     * positive root of \f$a\,t^2 - 2\,b\,t - c = 0\f$ with
     * \f$a = \log J_1 - 2\log J_2 + \log J_3\f$,
     * \f$b = \log J_1 - \log J_2\f$, and
     * \f$c = \log J - \log J_1\f$.
     * The actual resistivity is then
     * \f[
     *   \log_{10}\rho_{FF}(J) = (1-t)^2 \log\rho_1
     *                         + 2(1-t)\,t \log\rho_n
     *                         + t^2 \log\rho_n .
     * \f]
     *
     * This three-regime extension is BELFEM-specific; it complements the
     * pure Rhyner power law (Messe et al. 2023, Sec. 2.6) so that the same
     * material can be evaluated above the critical regime in quench
     * simulations without numerical blow-up.
     *
     * \param normJ  current-density magnitude \f$|J|\f$ [A/m²]
     * \return       effective electric resistivity [Ω·m]
     *
     * \pre  \f$J_c > 0\f$ and \f$n > 1\f$ (checked in debug mode).
     * \pre  The Bezier coefficient \f$|a|\f$ must be non-degenerate
     *       (checked in debug mode for the flux-flow branch).
     */
    inline real
    Material::rho_piecewise( const real normJ  ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) ) // (assumed critical temperature)
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Constant \f$J_c\f$ / \f$n\f$ with spatially-dependent defect modulation
     * \f$J_c \rightarrow J_c \cdot d(x,y,z,t)\f$.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real x, const real y, const real z, const real t  ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc )*((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ via mJcFunction; \f$n\f$ remains constant.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real normB, const real angleNxB  ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;


        BELFEM_ASSERT(
            mJcFunction->depends_on( material::JcParameter::normB ) &&
            mJcFunction->depends_on( material::JcParameter::angleNxB ) &&
            ! mJcFunction->depends_on( material::JcParameter::T ), "wrong powerlaw for material %s", mLabel.c_str() );

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }


        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ with spatial defect modulation.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;


        BELFEM_ASSERT(
            mJcFunction->depends_on( material::JcParameter::normB ) &&
            mJcFunction->depends_on( material::JcParameter::angleNxB ) &&
            ! mJcFunction->depends_on( material::JcParameter::T ), "wrong powerlaw for material %s", mLabel.c_str() );

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }


        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB )*((this->mDefectFunction) (x,y,z,t)) ; ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant: both \f$J_c\f$ and
     * \f$n\f$ come from JcFunction objects, and \f$\rho_n\f$ is evaluated at
     * the local temperature \p T.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real T, const real normB, const real angleNxB  ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin);

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant with spatial defect
     * modulation applied to \f$J_c\f$.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t  ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) ) // (assumed critical temperature)
        {
            return rhon ;
        }

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Temperature-only variant: \f$J_c(T)\f$ and \f$n(T)\f$ from the
     * user-supplied callbacks Material::jc_custom() / Material::n_custom().
     */
    inline real
    Material::rho_piecewise( const real normJ, const real T  ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) ;
        real  n = this->n_custom( T ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

    /**
     * \overload
     * Temperature-only custom-callback variant with spatial defect modulation.
     */
    inline real
    Material::rho_piecewise( const real normJ, const real T, const real x, const real y, const real z, const real t  ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) ) // (assumed critical temperature)
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_custom( T ) ;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        //Power law resistivity
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        //Power-law regime
        real j1 = jc * std::pow(10.0,2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            return rhoPL ;
        }

        //Normal regime
        real j3 = j1 * std::pow((rhon/rho1),1.0/(mNff)) ;
        real rho3 = rhon ;

        if ( normJ > j3 )
        {
            return rhon ;
        }

        //Flux-flow regime
        real j2 = j1 * std::pow((rhon/rho1),1.0/(n-1.0)) ;
        real rho2 = rhon ;

        //Bezier control points
        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;

        real a = logj1 - 2.0*logj2 + std::log10(j3) ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b+std::pow(b*b+a*c,0.5))/a ;

        real rhoFF = std::pow(10.0,std::log10(rho1)*(1.0-tParam)*(1.0-tParam) + std::log10(rho2)*2.0*(1.0-tParam)*tParam + std::log10(rho3)*tParam*tParam );

        return rhoFF ;
    }

//------------------------------------------------------------------------------
// Function-backed Jc and n accessors
//------------------------------------------------------------------------------

    /**
     * \brief Direct evaluation of the n-value as a function of \f$|B|\f$,
     *        \f$\angle(n,B)\f$ and \f$T\f$.
     *
     * Used by post-processing and by overloads of rho_powerlaw() / rho_piecewise()
     * that take a JcFunction-supplied \f$n(B,\angle,T)\f$. For materials with
     * a constant \f$n\f$, use Material::constant_property(MaterialProperty::n)
     * instead.
     *
     * \pre mNFunction must have been installed (checked in debug mode).
     */
    inline real
    Material::n( const real normB, const real angleNxB, const real T  ) const
    {
        BELFEM_ASSERT( mNFunction != nullptr,
            "Material %s does not have an N function", mLabel.c_str() ) ;

        return mNFunction->eval( normB, angleNxB, T ) ;
    }

    /**
     * \overload
     * Field-dependent \f$n(|B|, \angle(n,B))\f$ without temperature dependence.
     */
    inline real
    Material::n( const real normB, const real angleNxB  ) const
    {
        BELFEM_ASSERT( mNFunction != nullptr,
            "Material %s does not have an N function", mLabel.c_str() ) ;

        return mNFunction->eval( normB, angleNxB ) ;
    }

    /**
     * \brief Direct evaluation of the critical current density
     *        \f$J_c(|B|, \angle(n,B), T)\f$.
     *
     * Used by post-processing and by overloads of rho_powerlaw() /
     * rho_piecewise() that take a JcFunction-supplied
     * \f$J_c(B,\angle,T)\f$. For materials with a constant \f$J_c\f$,
     * use Material::constant_property(MaterialProperty::jc) instead.
     *
     * \pre mJcFunction must have been installed (checked in debug mode).
     *      Note: Material::have(jc) is true even when \f$J_c\f$ was set as a
     *      plain constant, so it is *not* sufficient to guard the function
     *      pointer — the explicit nullptr check is required.
     */
    inline real
    Material::jc( const real B, const real angleNxB, const real T ) const
    {
        // Constants-fallback: when jc was set as a plain constant via the
        // MaterialFactory constants path, mJcFunction is null but
        // constant_property(jc) is valid. See rho_powerlaw(...) for the
        // matching pattern in the assembly path.
        return ( mJcFunction != nullptr )
               ? mJcFunction->eval( B, angleNxB, T )
               : this->constant_property( MaterialProperty::jc ) ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle(n,B))\f$ without temperature dependence.
     */
    inline real
    Material::jc( const real B, const real angleNxB ) const
    {
        // Constants-fallback (see jc(B,angle,T)): have(jc) is true even for a
        // plain-constant jc, where mJcFunction is null.
        return ( mJcFunction != nullptr )
               ? mJcFunction->eval( B, angleNxB )
               : this->constant_property( MaterialProperty::jc ) ;
    }

    /**
     * \overload
     * Variant with spatial defect modulation
     * \f$J_c \rightarrow J_c \cdot d(x,y,z,t)\f$.
     */
    inline real
    Material::jc( const real B, const real angleNxB, const real T, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr,
            "Material %s does not have a Jc function", mLabel.c_str() ) ;

        return mJcFunction->eval( B, angleNxB, T )*((this->mDefectFunction) (x,y,z,t)) ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle(n,B))\f$ with spatial defect
     * modulation \f$J_c \rightarrow J_c \cdot d(x,y,z,t)\f$.
     */
    inline real
    Material::jc( const real B, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr,
            "Material %s does not have a Jc function", mLabel.c_str() ) ;

        return mJcFunction->eval( B, angleNxB )*((this->mDefectFunction) (x,y,z,t)) ;
    }


//------------------------------------------------------------------------------
// Power-law resistivity derivatives for Newton-Raphson
//------------------------------------------------------------------------------

    /**
     * \brief Analytic Jacobian \f$d\rho_{eff}/dJ\f$ of the power-law model.
     *
     * For the parallel-combined model
     * \f$\rho_{eff} = (1/\rho_n + 1/\rho_{PL})^{-1}\f$ the chain rule gives
     * \f[
     *     \frac{d\rho_{eff}}{dJ}
     *       \;=\; \frac{d\rho_{PL}/dJ}{(1 + \rho_{PL}/\rho_n)^2}, \qquad
     *     \frac{d\rho_{PL}}{dJ}
     *       \;=\; \frac{E_c}{J_c^{2}}\,(n-1)\,\left(\frac{|J|}{J_c}\right)^{n-2}.
     * \f]
     * Required by the nonlinear (Newton-Raphson) iteration in the HTS solver
     * to drive the residual below \f$10^{-11}\f$ as recommended in
     * Messe et al. 2023, Sec. 2.7. For \f$|J| < \f$BELFEM_EPSILON the
     * derivative is taken as zero to avoid the
     * \f$0^{n-2}\f$ singularity when \f$n < 2\f$.
     *
     * This overload assumes constant \f$J_c\f$ and \f$n\f$.
     *
     * \param normJ  current-density magnitude \f$|J|\f$ [A/m²]
     * \return       \f$d\rho_{eff}/dJ\f$ [Ω·m / (A/m²)]
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        // dρ_PL/dJ = (Ec/Jc^n)·(n-1)·J^(n-2)
        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        // dρ/dJ = (ρ²/ρ_PL²)·dρ_PL/dJ
        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Constant \f$J_c\f$ / \f$n\f$ with spatial defect modulation applied to \f$J_c\f$.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        // dρ_PL/dJ = (Ec/Jc^n)·(n-1)·J^(n-2)
        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        // dρ/dJ = (ρ²/ρ_PL²)·dρ_PL/dJ
        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Temperature-only variant: \f$J_c(T)\f$ and \f$n(T)\f$ from
     * Material::jc_custom() / Material::n_custom().
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real T ) const
    {
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) ;
        real  n = this->n_custom( T ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Temperature-only custom-callback variant with spatial defect modulation.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real T, const real x, const real y, const real z, const real t ) const
    {
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_custom( T ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ via mJcFunction; \f$n\f$ remains constant.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real normB, const real angleNxB ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( ! mJcFunction->depends_on( material::JcParameter::T ),
            "wrong powerlaw derivative for material %s", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ with spatial defect modulation;
     * \f$n\f$ remains constant.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( gTbulk ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant: both \f$J_c\f$ and
     * \f$n\f$ come from JcFunction objects.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant with spatial defect
     * modulation applied to \f$J_c\f$.
     */
    inline real
    Material::drho_powerlaw_dJ( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        real rhon = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );

        real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );

        return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
    }

//------------------------------------------------------------------------------
// Piecewise resistivity derivatives for Newton-Raphson
//------------------------------------------------------------------------------

    /**
     * \brief Branch-aware Jacobian \f$d\rho/dJ\f$ of the three-regime
     *        piecewise model.
     *
     * Branches mirror those of rho_piecewise():
     *
     *  - **Power-law regime** \f$(J \le J_1)\f$: the parallel-combined
     *    chain rule of drho_powerlaw_dJ() is applied here. NOTE:
     *    rho_piecewise() returns the RAW power law in this regime (no
     *    parallel combination), so unlike drho_piecewise_dB/dT this
     *    tangent is not consistent with the residual there.
     *  - **Normal regime** \f$(J > J_3)\f$: \f$\rho\f$ is constant, so the
     *    derivative is exactly zero.
     *  - **Flux-flow transition** \f$(J_1 < J \le J_3)\f$: differentiate the
     *    Bezier blend in log–log space. With
     *    \f$L_\rho \equiv \log_{10}\rho\f$,
     *    \f[
     *      \frac{d\rho}{dJ}
     *        \;=\; \rho\,\ln 10\,
     *              \frac{dL_\rho}{dt}\,
     *              \frac{dt}{dJ},
     *    \f]
     *    where the parametric derivatives follow from the implicit relation
     *    \f$a t^2 - 2bt - c(J) = 0\f$ defined in rho_piecewise(). Above
     *    \f$T_{crit}\f$ the derivative is taken as zero (only \f$\rho_n\f$
     *    is returned by rho_piecewise()), and for \f$|J| < \f$BELFEM_EPSILON
     *    the derivative is zero to avoid the singular \f$\log J\f$ term.
     *
     * Required by the nonlinear iteration in the HTS solver
     * (Messe et al. 2023, Sec. 2.7) so that Newton-Raphson can drive the
     * residual through the curved transition without losing quadratic
     * convergence.
     *
     * This overload assumes constant \f$J_c\f$ and \f$n\f$.
     *
     * \param normJ  current-density magnitude \f$|J|\f$ [A/m²]
     * \return       \f$d\rho/dJ\f$ [Ω·m / (A/m²)]
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        // Power-law regime thresholds
        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            // Power-law region: same derivative as power-law
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            // Normal regime: constant resistivity
            return 0.0 ;
        }

        // Flux-flow regime: derivative of Bezier curve in log-space
        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;
        real rho2 = rhon ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        // t parameter from quadratic formula
        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        // Derivative dt/dJ
        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        // Bezier curve: log(ρ) = (1-t)² log(rho1) + 2(1-t)t log(rho2) + t² log(rho3)
        // d(log(ρ))/dt = -2(1-t) log(rho1) + 2(1-2t) log(rho2) + 2t log(rho3)
        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rho2) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        // ρ = 10^(log(ρ))
        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        // dρ/dJ = ρ · ln(10) · d(log(ρ))/dt · dt/dJ
        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Constant \f$J_c\f$ / \f$n\f$ with spatial defect modulation applied to \f$J_c\f$.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( this->is_constant( MaterialProperty::jc ), "Material %s does not have a constant jc parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->constant_property( MaterialProperty::jc ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->constant_property( MaterialProperty::n ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        // Power-law regime thresholds
        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            // Power-law region: same derivative as power-law
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            // Normal regime: constant resistivity
            return 0.0 ;
        }

        // Flux-flow regime: derivative of Bezier curve in log-space
        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;
        real rho2 = rhon ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        // t parameter from quadratic formula
        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        // Derivative dt/dJ
        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        // Bezier curve: log(ρ) = (1-t)² log(rho1) + 2(1-t)t log(rho2) + t² log(rho3)
        // d(log(ρ))/dt = -2(1-t) log(rho1) + 2(1-2t) log(rho2) + 2t log(rho3)
        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rho2) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        // ρ = 10^(log(ρ))
        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        // dρ/dJ = ρ · ln(10) · d(log(ρ))/dt · dt/dJ
        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Temperature-only variant: \f$J_c(T)\f$ and \f$n(T)\f$ from
     * Material::jc_custom() / Material::n_custom().
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real T ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) ;
        real  n = this->n_custom( T ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Temperature-only custom-callback variant with spatial defect modulation.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real T, const real x, const real y, const real z, const real t ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_custom( T ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_custom( T ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ via mJcFunction; \f$n\f$ remains constant.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real normB, const real angleNxB ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;
        BELFEM_ASSERT( ! mJcFunction->depends_on( material::JcParameter::T ),
            "wrong powerlaw derivative for material %s", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Field-dependent \f$J_c(|B|, \angle)\f$ with spatial defect modulation;
     * \f$n\f$ remains constant.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        BELFEM_ASSERT( mJcFunction != nullptr, "Material %s does not have Jc function", mLabel.c_str() ) ;
        BELFEM_ASSERT( this->is_constant( MaterialProperty::n ),   "Material %s does not have a constant n parameter defined", mLabel.c_str() ) ;

        real rhon = this->rho( gTbulk ) ;

        if ( gTbulk > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = mJcFunction->eval( normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = mNFunction->eval( normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin );
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant: both \f$J_c\f$ and
     * \f$n\f$ come from JcFunction objects.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ERROR( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        cplx d = b*b + a*c ;
        BELFEM_ERROR( std::abs( d ) > BELFEM_EPSILON,
           "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
           ( double ) std::abs( d ), mLabel.c_str() ) ;

        cplx tParam = ( b + std::sqrt( d ) ) / a ;

        cplx dt_dc = 0.5 / std::sqrt( d ) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        cplx dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        cplx dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        cplx rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        cplx drho_dJ = rho * std::log(10.0) * dlogrho_dt * dt_dJ ;

        BELFEM_ERROR( std::imag( drho_dJ ) < BELFEM_EPSILON,
           "Piecewise power law derivative: imaginary part of drho_dJ is non-zero (|a| = %g) for material %s, %g + i * %g",
           ( double ) std::imag( drho_dJ ), mLabel.c_str(), std::real( drho_dJ ) , std::imag( drho_dJ ) ) ;

        return std::real( drho_dJ ) ;
    }

    /**
     * \overload
     * Fully \f$(|B|, \angle, T)\f$-dependent variant with spatial defect
     * modulation applied to \f$J_c\f$.
     */
    inline real
    Material::drho_piecewise_dJ( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ((this->mDefectFunction) (x,y,z,t)) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if (normJ < BELFEM_EPSILON) return 0.0;

        BELFEM_ASSERT( jc > 0.0, "Piecewise power law derivative requires jc > 0 (got jc = %g) for material %s", ( double ) jc, mLabel.c_str() ) ;
        BELFEM_ASSERT( n  > 1.0, "Piecewise power law derivative requires n > 1 (got n = %g) for material %s",  ( double ) n,  mLabel.c_str() ) ;

        real j1 = jc * std::pow(10.0, 2.5/n) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin) ;

        if (normJ <= j1)
        {
            real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;
            real drho_PL_dJ = ( ec / std::pow(jc, 2) ) * ( n - 1.0 ) * std::pow( normJ/jc, n - 2.0 );
            return drho_PL_dJ/std::pow((1+rhoPL/rhon),2.0) ;
        }

        real j3 = j1 * std::pow((rhon/rho1), 1.0/(mNff)) ;

        if ( normJ > j3 )
        {
            return 0.0 ;
        }

        real j2 = j1 * std::pow((rhon/rho1), 1.0/(n-1.0)) ;

        real logj1 = std::log10(j1) ;
        real logj2 = std::log10(j2) ;
        real logj3 = std::log10(j3) ;

        real a = logj1 - 2.0*logj2 + logj3 ;
        real b = logj1 - logj2 ;
        real c = std::log10(normJ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law derivative: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tParam = (b + std::pow(b*b + a*c, 0.5)) / a ;

        real dt_dc = 0.5 / std::pow(b*b + a*c, 0.5) ;
        real dc_dJ = 1.0 / (normJ * std::log(10.0)) ;
        real dt_dJ = dt_dc * dc_dJ ;

        real logrho1 = std::log10(rho1) ;
        real logrho2 = std::log10(rhon) ;
        real logrho3 = std::log10(rhon) ;

        real dlogrho_dt = -2.0*(1.0-tParam)*logrho1 + 2.0*(1.0-2.0*tParam)*logrho2 + 2.0*tParam*logrho3 ;

        real rho = std::pow(10.0, (1.0-tParam)*(1.0-tParam)*logrho1 + 2.0*(1.0-tParam)*tParam*logrho2 + tParam*tParam*logrho3) ;

        return rho * std::log(10.0) * dlogrho_dt * dt_dJ ;
    }

//------------------------------------------------------------------------------
// Field derivatives of the HTS laws ( audited 2026-08-13 )
//
// jc = jc(T,|B|,θ) and n = n(T,|B|,θ) through the JcFunction hooks; at fixed
// J, T, θ the unfloored power law p0 = (ec/jc)·(J/jc)^(n−1) has
//
//     ∂p0/∂jc = −n·p0/jc          ∂p0/∂n = p0·ln(J/jc)
//
// so  dp0/d|B| = p0·[ −(n/jc)·djc/d|B| + ln(J/jc)·dn/d|B| ].
//
// Conventions per the audit round ( tmp/ai_exchange/jc_derivative_plumbing.md,
// both voices ):
//  - drho_powerlaw_dB mirrors drho_powerlaw_dJ exactly: parallel factor
//    ( 1 + rhoPL/rhon )^-2 with the FLOORED rhoPL, unfloored law in the
//    numerator, J < eps guard.
//  - drho_piecewise_dB is consistent with rho_piecewise's OWN residual: the
//    power-law regime returns the RAW dp0 ( rho_piecewise returns raw rhoPL
//    there, no parallel combination ); the normal regime and T > T_crit are
//    exactly 0; the Bézier blend is 0 FOR NOW ( staged — the regime
//    boundaries also move with jc, that derivative is deferred, and the
//    tangent has a documented jump at j1 ).
//  - defect overloads modulate jc AND its derivative by the same D(x,y,z,t).
//------------------------------------------------------------------------------

    inline real
    Material::drho_powerlaw_dB( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        real rhon  = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        // unfloored law differentiated, cf. drho_powerlaw_dJ
        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        real dp0dB = p0 * ( -( n / jc ) * djcdB
                            + std::log( normJ / jc ) * dndB ) ;

        return dp0dB / std::pow( ( 1.0 + rhoPL / rhon ), 2.0 ) ;
    }

    inline real
    Material::drho_powerlaw_dB( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        real ec = this->constant_property( MaterialProperty::ec ) ;

        // defect modulation applies to jc AND its derivative ( product rule
        // with D independent of |B| ); n is not modulated, cf. rho_powerlaw
        real tD    = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc    = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n    = this->n_eval( T, normB, angleNxB ) ;

        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) * tD ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        real rhon  = this->rho( T ) ;
        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        real dp0dB = p0 * ( -( n / jc ) * djcdB
                            + std::log( normJ / jc ) * dndB ) ;

        return dp0dB / std::pow( ( 1.0 + rhoPL / rhon ), 2.0 ) ;
    }

    inline real
    Material::drho_piecewise_dB( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        // above T_crit rho_piecewise returns rhon, independent of jc and n
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        // outside the power-law regime: normal regime is exactly 0, the
        // Bézier blend is staged 0 ( see the block comment above )
        real j1 = jc * std::pow( 10.0, 2.5 / n ) ;
        if ( normJ > j1 ) return 0.0 ;

        // power-law regime: rho_piecewise returns the RAW power law here,
        // so the exact derivative carries NO parallel factor
        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        return p0 * ( -( n / jc ) * djcdB
                      + std::log( normJ / jc ) * dndB ) ;
    }

    inline real
    Material::drho_piecewise_dB( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;

        real tD    = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc    = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n    = this->n_eval( T, normB, angleNxB ) ;

        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) * tD ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        real j1 = jc * std::pow( 10.0, 2.5 / n ) ;
        if ( normJ > j1 ) return 0.0 ;

        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        return p0 * ( -( n / jc ) * djcdB
                      + std::log( normJ / jc ) * dndB ) ;
    }

//------------------------------------------------------------------------------
// Temperature derivatives of the HTS laws ( T-leg, audited 2026-08-13 )
//
// The quench-feedback tangent: jc(T), n(T) AND rho_n(T) all move. With
// a = rho_n, b = p0 = (ec/jc)·(J/jc)^(n−1) and c = parallel(a,b),
//
//     ∂c/∂a = b²/(a+b)² = ( c/a )²        ∂c/∂b = a²/(a+b)² = 1/(1+b/a)²
//
// dp0/dT = p0·( −(n/jc)·djc/dT + ln(J/jc)·dn/dT ), same chain as the |B|
// channel. These closed forms replace the retired b = a·c/(c−a)
// reconstruction in the old compute_drhodT_hts, which was sign-flipped
// ( correct: a·c/(a−c) ) and whose dadT weighting c²/(a−2c)² diverged at
// the flux-flow crossover rho_PL = rho_n. Conventions:
//  - floored rhoPL in the parallel factors, unfloored law differentiated
//    ( dB/dJ convention );
//  - NO ( djc==0 && dn==0 ) early-out: the rho_n term is the only correct
//    T-dependence of a constant-jc material and must survive;
//  - drho_piecewise_dT follows rho_piecewise's OWN residual branch for
//    branch; the Bézier blend carries BOTH parts ( 2026-08-14 ): the
//    frozen-knot partial ( rho_n control points + the explicit rho1
//    dependence ) AND the knot motion, i.e. dtParam/dT through j1(T),
//    j2(T), j3(T). The knot motion dominates — with it staged out the
//    frozen part alone reproduced only 4-33 % of dρ/dT across the blend,
//    which is why the thermal Newton stalled in flux-flow. It degrades to
//    the frozen part alone at the degenerate n−1 == mNff transition ( see
//    the guard at the discriminant ). Jumps at j1 and T_crit are
//    inherited from the residual itself;
//  - defect overloads modulate jc AND djc/dT by the same D(x,y,z,t).
//------------------------------------------------------------------------------

    inline real
    Material::drho_powerlaw_dT( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        real rhon  = this->rho( T ) ;
        real dadT  = this->drhodT( T ) ;

        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        // unfloored law differentiated, cf. drho_powerlaw_dB
        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        real dp0dT = p0 * ( -( n / jc ) * djcdT
                            + std::log( normJ / jc ) * dndT ) ;

        return dp0dT / std::pow( ( 1.0 + rhoPL / rhon ), 2.0 )
             + dadT * std::pow( rhoPL / ( rhon + rhoPL ), 2.0 ) ;
    }

    inline real
    Material::drho_powerlaw_dT( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        real ec = this->constant_property( MaterialProperty::ec ) ;

        // defect modulation applies to jc AND its derivative ( product rule
        // with D independent of T ); n is not modulated, cf. rho_powerlaw
        real tD    = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc    = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n    = this->n_eval( T, normB, angleNxB ) ;

        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) * tD ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        real rhon  = this->rho( T ) ;
        real dadT  = this->drhodT( T ) ;

        real rhoPL = std::max(( ec / jc ) * std::pow( normJ / jc, n - 1 ), mRhoMin ) ;

        real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

        real dp0dT = p0 * ( -( n / jc ) * djcdT
                            + std::log( normJ / jc ) * dndT ) ;

        return dp0dT / std::pow( ( 1.0 + rhoPL / rhon ), 2.0 )
             + dadT * std::pow( rhoPL / ( rhon + rhoPL ), 2.0 ) ;
    }

    inline real
    Material::drho_piecewise_dT( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        real dadT = this->drhodT( T ) ;

        // above T_crit rho_piecewise returns rhon(T) exactly
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return dadT ;
        }

        // dependency-routed evaluation with constants-fallback ( O1 policy )
        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        // power-law regime: rho_piecewise returns the RAW power law here,
        // so the exact derivative carries NO parallel factor and NO rho_n
        // term ( J < eps is inside this branch: p0 -> 0 for n > 1 )
        real j1 = jc * std::pow( 10.0, 2.5 / n ) ;
        if ( normJ <= j1 )
        {
            if ( normJ < BELFEM_EPSILON ) return 0.0 ;

            real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

            return p0 * ( -( n / jc ) * djcdT
                          + std::log( normJ / jc ) * dndT ) ;
        }

        real rhon = this->rho( T ) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin ) ;

        // normal regime: rho_piecewise returns rhon(T) exactly
        real j3 = j1 * std::pow( rhon / rho1, 1.0 / ( mNff ) ) ;
        if ( normJ > j3 )
        {
            return dadT ;
        }

        // Bézier blend: FULL derivative ( 2026-08-14 ). The residual's Bézier
        // machinery is reproduced verbatim ( cf. rho_piecewise ) to obtain
        // tParam and rhoFF, then differentiated in two parts — (i) the
        // control points rho1 and rho_n at frozen knots, and (ii) the knot
        // motion, i.e. dtParam/dT through j1(T), j2(T), j3(T). Part (ii)
        // dominates: it was staged out when the T-leg first landed, and a
        // finite-difference check then showed the frozen part alone
        // reproduces only 4-33 % of dρ/dT across the blend. With both parts
        // the derivative is exact to roundoff.
        real j2 = j1 * std::pow( rhon / rho1, 1.0 / ( n - 1.0 ) ) ;

        real logj1 = std::log10( j1 ) ;
        real logj2 = std::log10( j2 ) ;

        real a = logj1 - 2.0 * logj2 + std::log10( j3 ) ;
        real b = logj1 - logj2 ;
        real c = std::log10( normJ ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tDisc  = b * b + a * c ;
        real tSqrt  = std::pow( tDisc, 0.5 ) ;
        real tParam = ( b + tSqrt ) / a ;

        real rhoFF = std::pow( 10.0,
              std::log10( rho1 ) * ( 1.0 - tParam ) * ( 1.0 - tParam )
            + std::log10( rhon ) * 2.0 * ( 1.0 - tParam ) * tParam
            + std::log10( rhon ) * tParam * tParam ) ;

        // (i) control points move: rho1 = (ec/jc)·10^{2.5(n−1)/n} and rho_n(T)
        real dlnrho1dT = -djcdT / jc + 2.5 * std::log( 10.0 ) / ( n * n ) * dndT ;

        real tFrozen = rhoFF * ( ( 2.0 * tParam - tParam * tParam ) * dadT / rhon
                                 + ( 1.0 - tParam ) * ( 1.0 - tParam ) * dlnrho1dT ) ;

        // (ii) KNOT MOTION — the dominant part in this regime, and the reason
        // the frozen-knot form alone was not usable: j1, j2, j3 all slide as
        // jc(T) falls, so an element at fixed J moves DEEPER into the
        // transition and tParam itself carries a T-derivative. Measured
        // against a finite difference at 77 K REBCO constants, the frozen
        // terms above reproduce only 4-33 % of dρ/dT across the blend; with
        // this term the agreement is exact to roundoff.
        //
        // ln j1 = ln jc + 2.5·ln10/n, hence dln(j1)/dT = −dln(rho1)/dT
        // exactly; j3 and j2 follow from their definitions, j2 additionally
        // through the explicit 1/(n−1) exponent.
        real dlnj1dT   = -dlnrho1dT ;
        real dlnrhondT = dadT / rhon ;
        real dlnj3dT   = dlnj1dT + ( dlnrhondT - dlnrho1dT ) / mNff ;
        real dlnj2dT   = dlnj1dT + ( dlnrhondT - dlnrho1dT ) / ( n - 1.0 )
                         - dndT / ( ( n - 1.0 ) * ( n - 1.0 ) )
                           * std::log( rhon / rho1 ) ;

        // a, b, c are base-10 logs of the knots: d/dT picks up 1/ln10
        real tInvLn10 = 1.0 / std::log( 10.0 ) ;
        real dA = ( dlnj1dT - 2.0 * dlnj2dT + dlnj3dT ) * tInvLn10 ;
        real dB = ( dlnj1dT - dlnj2dT ) * tInvLn10 ;
        real dC = -dlnj1dT * tInvLn10 ;

        // tSqrt is the Bézier root's discriminant, and it is NOT bounded away
        // from zero for every legal material. With R = log10( rho_n / rho1 ),
        // p = n−1 and q = mNff, at the upper knot J = j3 it reduces to
        //
        //     b² + ac = R² ( 1/p − 1/q )²
        //
        // which VANISHES when p == q — n = 4 at the default mNff = 3, legal
        // under the n > 1 precondition, and exactly where j2 == j3, i.e. the
        // Bézier control points coincide. The existing |a| > eps assert does
        // not cover this ( a = −R/q there, nonzero ).
        //
        // The singularity is REMOVABLE, not a pole: dt/dT carries a 1/tSqrt
        // factor, but the weight derivative it multiplies is ∝ ( 1 − t ), and
        // at p == q one has a == b, hence 1 − t = −tSqrt/a. The product is
        // finite ( checked numerically: it converges as J → j3 ). What the
        // guard avoids is therefore an indeterminate 0/0 in floating point,
        // not an infinite physical derivative.
        //
        // The fallback is exact at the endpoint rather than merely safe: at
        // tSqrt = 0 the root gives t = 1, where rhoFF = rho_n and tFrozen
        // reduces to dadT — the normal-branch derivative — so the tangent
        // stays continuous with the J > j3 branch. In the infinitesimal
        // neighborhood where the guard bites, a finite knot contribution is
        // dropped; that is a documented approximation for a degenerate
        // material ( n = 4 ), not a mechanism any REBCO deck reaches.
        // Guarded relative to |b|, the natural scale of the root
        // ( tSqrt = |b| at c = 0 ).
        // Tested on the DISCRIMINANT, not on its root, and with the
        // comparison negated: roundoff can push tDisc slightly negative at
        // the degenerate point, which makes tSqrt a NaN — and a NaN fails
        // every ordinary comparison, so a `tSqrt <= tol` guard would be
        // bypassed and the NaN would propagate into the tangent. `!( x > tol )`
        // takes the fallback for NaN as well.
        if ( ! ( tDisc > BELFEM_EPSILON * ( b * b + 1.0 ) ) )
        {
            return tFrozen ;
        }

        real dSdT = ( 2.0 * b * dB + dA * c + a * dC ) / ( 2.0 * tSqrt ) ;
        real dtdT = ( ( dB + dSdT ) - tParam * dA ) / a ;

        // d(rhoFF)/dt at frozen control points: the weight derivatives
        // −2(1−t) on log ρ1 and (2−2t) on log ρ_n collapse to
        // 2(1−t)·ln( ρ_n / ρ1 ) once the ln10 factors cancel
        return tFrozen
             + rhoFF * 2.0 * ( 1.0 - tParam ) * std::log( rhon / rho1 ) * dtdT ;
    }

    inline real
    Material::drho_piecewise_dT( const real normJ, const real T, const real normB, const real angleNxB, const real x, const real y, const real z, const real t ) const
    {
        real dadT = this->drhodT( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return dadT ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;

        // defect modulation applies to jc AND its derivative ( product rule
        // with D independent of T ); n is not modulated, cf. rho_piecewise
        real tD    = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc    = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n    = this->n_eval( T, normB, angleNxB ) ;

        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) * tD ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        real j1 = jc * std::pow( 10.0, 2.5 / n ) ;
        if ( normJ <= j1 )
        {
            if ( normJ < BELFEM_EPSILON ) return 0.0 ;

            real p0 = ( ec / jc ) * std::pow( normJ / jc, n - 1 ) ;

            return p0 * ( -( n / jc ) * djcdT
                          + std::log( normJ / jc ) * dndT ) ;
        }

        real rhon = this->rho( T ) ;
        real rho1 = std::max(( ec / jc ) * std::pow( j1 / jc, n - 1 ), mRhoMin ) ;

        real j3 = j1 * std::pow( rhon / rho1, 1.0 / ( mNff ) ) ;
        if ( normJ > j3 )
        {
            return dadT ;
        }

        real j2 = j1 * std::pow( rhon / rho1, 1.0 / ( n - 1.0 ) ) ;

        real logj1 = std::log10( j1 ) ;
        real logj2 = std::log10( j2 ) ;

        real a = logj1 - 2.0 * logj2 + std::log10( j3 ) ;
        real b = logj1 - logj2 ;
        real c = std::log10( normJ ) - logj1 ;

        BELFEM_ASSERT( std::abs( a ) > BELFEM_EPSILON,
            "Piecewise power law: degenerate Bezier transition (|a| = %g) for material %s",
            ( double ) std::abs( a ), mLabel.c_str() ) ;

        real tDisc  = b * b + a * c ;
        real tSqrt  = std::pow( tDisc, 0.5 ) ;
        real tParam = ( b + tSqrt ) / a ;

        real rhoFF = std::pow( 10.0,
              std::log10( rho1 ) * ( 1.0 - tParam ) * ( 1.0 - tParam )
            + std::log10( rhon ) * 2.0 * ( 1.0 - tParam ) * tParam
            + std::log10( rhon ) * tParam * tParam ) ;

        // (i) control points move: rho1 = (ec/jc)·10^{2.5(n−1)/n} and rho_n(T)
        real dlnrho1dT = -djcdT / jc + 2.5 * std::log( 10.0 ) / ( n * n ) * dndT ;

        real tFrozen = rhoFF * ( ( 2.0 * tParam - tParam * tParam ) * dadT / rhon
                                 + ( 1.0 - tParam ) * ( 1.0 - tParam ) * dlnrho1dT ) ;

        // (ii) KNOT MOTION — the dominant part in this regime, and the reason
        // the frozen-knot form alone was not usable: j1, j2, j3 all slide as
        // jc(T) falls, so an element at fixed J moves DEEPER into the
        // transition and tParam itself carries a T-derivative. Measured
        // against a finite difference at 77 K REBCO constants, the frozen
        // terms above reproduce only 4-33 % of dρ/dT across the blend; with
        // this term the agreement is exact to roundoff.
        //
        // ln j1 = ln jc + 2.5·ln10/n, hence dln(j1)/dT = −dln(rho1)/dT
        // exactly; j3 and j2 follow from their definitions, j2 additionally
        // through the explicit 1/(n−1) exponent.
        real dlnj1dT   = -dlnrho1dT ;
        real dlnrhondT = dadT / rhon ;
        real dlnj3dT   = dlnj1dT + ( dlnrhondT - dlnrho1dT ) / mNff ;
        real dlnj2dT   = dlnj1dT + ( dlnrhondT - dlnrho1dT ) / ( n - 1.0 )
                         - dndT / ( ( n - 1.0 ) * ( n - 1.0 ) )
                           * std::log( rhon / rho1 ) ;

        // a, b, c are base-10 logs of the knots: d/dT picks up 1/ln10
        real tInvLn10 = 1.0 / std::log( 10.0 ) ;
        real dA = ( dlnj1dT - 2.0 * dlnj2dT + dlnj3dT ) * tInvLn10 ;
        real dB = ( dlnj1dT - dlnj2dT ) * tInvLn10 ;
        real dC = -dlnj1dT * tInvLn10 ;

        // tSqrt is the Bézier root's discriminant, and it is NOT bounded away
        // from zero for every legal material. With R = log10( rho_n / rho1 ),
        // p = n−1 and q = mNff, at the upper knot J = j3 it reduces to
        //
        //     b² + ac = R² ( 1/p − 1/q )²
        //
        // which VANISHES when p == q — n = 4 at the default mNff = 3, legal
        // under the n > 1 precondition, and exactly where j2 == j3, i.e. the
        // Bézier control points coincide. The existing |a| > eps assert does
        // not cover this ( a = −R/q there, nonzero ).
        //
        // The singularity is REMOVABLE, not a pole: dt/dT carries a 1/tSqrt
        // factor, but the weight derivative it multiplies is ∝ ( 1 − t ), and
        // at p == q one has a == b, hence 1 − t = −tSqrt/a. The product is
        // finite ( checked numerically: it converges as J → j3 ). What the
        // guard avoids is therefore an indeterminate 0/0 in floating point,
        // not an infinite physical derivative.
        //
        // The fallback is exact at the endpoint rather than merely safe: at
        // tSqrt = 0 the root gives t = 1, where rhoFF = rho_n and tFrozen
        // reduces to dadT — the normal-branch derivative — so the tangent
        // stays continuous with the J > j3 branch. In the infinitesimal
        // neighborhood where the guard bites, a finite knot contribution is
        // dropped; that is a documented approximation for a degenerate
        // material ( n = 4 ), not a mechanism any REBCO deck reaches.
        // Guarded relative to |b|, the natural scale of the root
        // ( tSqrt = |b| at c = 0 ).
        // Tested on the DISCRIMINANT, not on its root, and with the
        // comparison negated: roundoff can push tDisc slightly negative at
        // the degenerate point, which makes tSqrt a NaN — and a NaN fails
        // every ordinary comparison, so a `tSqrt <= tol` guard would be
        // bypassed and the NaN would propagate into the tangent. `!( x > tol )`
        // takes the fallback for NaN as well.
        if ( ! ( tDisc > BELFEM_EPSILON * ( b * b + 1.0 ) ) )
        {
            return tFrozen ;
        }

        real dSdT = ( 2.0 * b * dB + dA * c + a * dC ) / ( 2.0 * tSqrt ) ;
        real dtdT = ( ( dB + dSdT ) - tParam * dA ) / a ;

        // d(rhoFF)/dt at frozen control points: the weight derivatives
        // −2(1−t) on log ρ1 and (2−2t) on log ρ_n collapse to
        // 2(1−t)·ln( ρ_n / ρ1 ) once the ln10 factors cancel
        return tFrozen
             + rhoFF * 2.0 * ( 1.0 - tParam ) * std::log( rhon / rho1 ) * dtdT ;
    }

//------------------------------------------------------------------------------
// The riva law ( 2026-08-27 )
//
// The same parallel model as rho_powerlaw — the superconducting power-law
// channel in parallel with the normal-state channel ( Duron et al. 2004;
// Riva 2021, EPFL thesis 8754, Eq. 5.4 ) — but TOTAL over the full range a
// measured jc/n table can produce mid-iterate:
//  - jc_eff ≤ 0 or nonfinite ( dead defect D(x)=0, spline underflow ):
//    the superconducting channel is gone, the material is fully normal;
//  - the power-law channel is evaluated in log10 space; past a cap the
//    parallel combination is ρn to machine precision, so the residual
//    returns ρn and the tangents return the matching normal-branch values
//    instead of the raw inf/inf;
//  - n arrives from n_eval pre-floored at 1 ( ohmic limit ρPL = ec/jc,
//    J-independent; dn_eval_* are 0 while the floor binds );
//  - no n > 1 precondition anywhere. ( A PROVABLY bad n source -- a table
//    or constant whose bound sits at or below 1 -- is refused at setup in
//    set_resistivity_law; the runtime stays total for what remains. )
// Weight algebra as in the powerlaw T-leg: with w = ρn/(ρPL+ρn),
// ∂ρ/∂ρPL = w², ∂ρ/∂ρn = (1−w)².
// Deviation from Riva Eq. 5.1/5.4: BELFEM keeps its mRhoMin FLOOR semantics
// ( zero by default ) instead of Riva's additive 1e-17 regularization — an
// additive term would reopen the 2026-08-10 value/tangent desync.
//------------------------------------------------------------------------------

    inline bool
    Material::riva_rho_pl( const real normJ, const real jc, const real n, const real ec, real & rhoPL ) const
    {
        if ( n <= 1.0 )
        {
            // ohmic floor: J-independent resistor. A subnormal jc passes the
            // caller's positivity gate but overflows ec/jc, and a bad deck
            // can set ec <= 0 — an infinite or negative channel resistance
            // takes the fully-normal branch
            rhoPL = ec / jc ;
            return std::isfinite( rhoPL ) && rhoPL >= 0.0 ;
        }
        else if ( normJ < BELFEM_EPSILON )
        {
            rhoPL = mRhoMin ;
        }
        else
        {
            real lg = std::log10( ec / jc )
                    + ( n - 1.0 ) * std::log10( normJ / jc ) ;

            // past this cap 1/ρPL vanishes to machine precision against any
            // physical ρn — the caller takes the fully-normal branch
            // negated NaN-aware comparison: a NaN lg ( NaN n from the
            // table, negative ec, or inf·0 at J == jc with infinite n )
            // must also land here, and NaN fails every ordinary comparison
            if ( ! ( lg <= 250.0 ) )
            {
                return false ;
            }
            rhoPL = std::max( std::pow( 10.0, lg ), mRhoMin ) ;
        }
        return true ;
    }

    inline real
    Material::rho_riva( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) )
        {
            return rhon ;
        }

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return rhon ;
        }

        // parallel combination, branch-stable against a huge ρPL
        return rhoPL > rhon ? rhon  / ( 1.0 + rhon  / rhoPL )
                            : rhoPL / ( 1.0 + rhoPL / rhon  ) ;
    }

    inline real
    Material::drho_riva_dJ( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        // fully-normal fallbacks: ρ = ρn there, which is J-independent
        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return 0.0 ;
        if ( n <= 1.0 ) return 0.0 ;
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return 0.0 ;
        }

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;

        // grouped so each factor stays bounded: w·ρPL ≤ ρn even when ρPL
        // is huge, and w ≤ 1
        return ( w * rhoPL ) * ( w * ( n - 1.0 ) / normJ ) ;
    }

    inline real
    Material::drho_riva_dB( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return 0.0 ;

        // the J < eps early-out is only valid for n > 1 ( rhoPL -> mRhoMin );
        // at the ohmic floor rhoPL = ec/jc is J-independent but still
        // B-dependent through jc, so the tangent must follow the residual
        if ( n > 1.0 && normJ < BELFEM_EPSILON ) return 0.0 ;

        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return 0.0 ;
        }

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;

        // dρPL/dB = ρPL·( −(n/jc)·djc/dB + ln(J/jc)·dn/dB ); the ln term
        // vanishes with dndB = 0 while the n-floor binds
        // djc/jc grouped as a ratio ( bounded for every real jc source ) so
        // a tiny jc cannot overflow n/jc on its own; weights grouped so each
        // factor stays bounded ( w·ρPL ≤ ρn, w ≤ 1 )
        return ( w * rhoPL ) * ( w * ( -( n * ( djcdB / jc ) )
                                       + ( normJ < BELFEM_EPSILON ? 0.0
                                           : std::log( normJ ) - std::log( jc ) ) * dndB ) ) ;
    }

    inline real
    Material::drho_riva_dT( const real normJ, const real T, const real normB, const real angleNxB ) const
    {
        real dadT = this->drhodT( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return dadT ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        // fully-normal fallbacks: ρ = ρn there, tangent follows it
        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return dadT ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return dadT ;
        }

        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;
        real v    = rhoPL / ( rhoPL + rhon ) ;

        // dρPL/dT = ρPL·( −(n/jc)·djc/dT + ln(J/jc)·dn/dT ); the ln term is
        // suppressed below BELFEM_EPSILON where ρPL is 0 ( n > 1 ) or
        // J-independent ( n = 1, dndT = 0 under the floor )
        // djc/jc grouped as a ratio and the weights as (w·ρPL)·(w·dln) so
        // no factor can overflow on a tiny-but-finite jc
        real dlnp = ( normJ < BELFEM_EPSILON && n > 1.0 ) ? 0.0 :
            -( n * ( djcdT / jc ) )
                + ( normJ < BELFEM_EPSILON ? 0.0
                    : std::log( normJ ) - std::log( jc ) ) * dndT ;

        return ( w * rhoPL ) * ( w * dlnp ) + v * v * dadT ;
    }

    inline real
    Material::rho_riva( const real normJ, const real T, const real normB, const real angleNxB,
                        const real x, const real y, const real z, const real t ) const
    {
        real rhon = this->rho( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return rhon ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ( ( this->mDefectFunction )( x, y, z, t ) ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) )
        {
            return rhon ;
        }

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return rhon ;
        }

        return rhoPL > rhon ? rhon  / ( 1.0 + rhon  / rhoPL )
                            : rhoPL / ( 1.0 + rhoPL / rhon  ) ;
    }

    inline real
    Material::drho_riva_dJ( const real normJ, const real T, const real normB, const real angleNxB,
                            const real x, const real y, const real z, const real t ) const
    {
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * ( ( this->mDefectFunction )( x, y, z, t ) ) ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return 0.0 ;
        if ( n <= 1.0 ) return 0.0 ;
        if ( normJ < BELFEM_EPSILON ) return 0.0 ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return 0.0 ;
        }

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;

        // grouped so each factor stays bounded: w·ρPL ≤ ρn even when ρPL
        // is huge, and w ≤ 1
        return ( w * rhoPL ) * ( w * ( n - 1.0 ) / normJ ) ;
    }

    inline real
    Material::drho_riva_dB( const real normJ, const real T, const real normB, const real angleNxB,
                            const real x, const real y, const real z, const real t ) const
    {
        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return 0.0 ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real tD = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return 0.0 ;

        // J < eps early-out only for n > 1, cf. the non-defect overload
        if ( n > 1.0 && normJ < BELFEM_EPSILON ) return 0.0 ;

        // defect modulates jc AND its derivative by the same D
        real djcdB = this->djc_eval_dB( T, normB, angleNxB ) * tD ;
        real dndB  = this->dn_eval_dB( T, normB, angleNxB ) ;

        if ( djcdB == 0.0 && dndB == 0.0 ) return 0.0 ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return 0.0 ;
        }

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;

        // djc/jc grouped as a ratio ( bounded for every real jc source ) so
        // a tiny jc cannot overflow n/jc on its own; weights grouped so each
        // factor stays bounded ( w·ρPL ≤ ρn, w ≤ 1 )
        return ( w * rhoPL ) * ( w * ( -( n * ( djcdB / jc ) )
                                       + ( normJ < BELFEM_EPSILON ? 0.0
                                           : std::log( normJ ) - std::log( jc ) ) * dndB ) ) ;
    }

    inline real
    Material::drho_riva_dT( const real normJ, const real T, const real normB, const real angleNxB,
                            const real x, const real y, const real z, const real t ) const
    {
        real dadT = this->drhodT( T ) ;

        if ( T > this->constant_property( MaterialProperty::T_crit ) )
        {
            return dadT ;
        }

        real ec = this->constant_property( MaterialProperty::ec ) ;
        real tD = ( this->mDefectFunction )( x, y, z, t ) ;
        real jc = this->jc_eval( T, normB, angleNxB ) * tD ;
        real  n = this->n_eval( T, normB, angleNxB ) ;

        if ( ! ( jc > 0.0 && std::isfinite( jc ) ) ) return dadT ;

        real rhoPL ;
        if ( ! this->riva_rho_pl( normJ, jc, n, ec, rhoPL ) )
        {
            return dadT ;
        }

        // defect modulates jc AND its derivative by the same D
        real djcdT = this->djc_eval_dT( T, normB, angleNxB ) * tD ;
        real dndT  = this->dn_eval_dT( T, normB, angleNxB ) ;

        real rhon = this->rho( T ) ;
        real w    = rhon / ( rhoPL + rhon ) ;
        real v    = rhoPL / ( rhoPL + rhon ) ;

        // djc/jc grouped as a ratio and the weights as (w·ρPL)·(w·dln) so
        // no factor can overflow on a tiny-but-finite jc
        real dlnp = ( normJ < BELFEM_EPSILON && n > 1.0 ) ? 0.0 :
            -( n * ( djcdT / jc ) )
                + ( normJ < BELFEM_EPSILON ? 0.0
                    : std::log( normJ ) - std::log( jc ) ) * dndT ;

        return ( w * rhoPL ) * ( w * dlnp ) + v * v * dadT ;
    }

}
#endif //BELFEM_POWERLAWS_HPP
