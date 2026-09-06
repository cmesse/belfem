# Nitrogen Helmholtz EoS: Completion and Wiring

**Date:** 2026-08-06
**Purpose:** Everything still needed to finish `EoS_Nitrogen` and wire it into
`Gas`. Written as reconstruction insurance after an editor save failed to reach
disk — the three outstanding derivative bodies are given in full below, derived
from Span et al. 2000, so nothing depends on session context surviving.
**Module:** `src/physics/gasmodels`
**Status:** **THE PREMISE OF THIS FILE IS NO LONGER TRUE — corrected 2026-09-03 by the
`todo/` currentness sweep (Grok found it; every claim below re-checked against the tree).**
All nine functions are written, the file is in the build, and it is wired into `Gas`.
What remains is a *correctness* question, not a completion one — and it has never been
gated. See the correction block below before reading anything else in this file.

**2026-09-03 correction.** Each of the three "outstanding" bodies now exists and each named
syntax defect is gone: `cl_GM_EoS_Nitrogen.cpp` is **13691 bytes** (was 7050),
`compute_phir_dd` (`:246`), `compute_phir_tt` (`:304`) and `compute_phir_dt` (`:338`) all
carry real implementations, the stray `* *` is absent (zero matches in the file), and `kk`
is properly declared at `:232` and incremented in the `k=32..36` loop at `:233`. The file is
**in** `SOURCES` (`src/physics/gasmodels/CMakeLists.txt:14`) — the opposite of what the
paragraph below asserted — and `Gas` constructs it at `cl_Gas.cpp:845`, so it is reachable.
The completion work was evidently done during the gasmodels open-source migration and never
recorded here.

**Therefore the "reconstruction insurance" framing below is spent**: the corrected bodies no
longer exist only in this file, and the versions in this file are now the *older* text. Do
not paste them over the source. **What is genuinely owed** is what this plan never had: a
numerical gate on the nine bodies against Span et al. 2000 reference values. Nothing here
has been compiled or run by this sweep — reviewed, not verified.

*Superseded status text, kept for the record:* "OPEN — six of nine functions correct on
disk; three outstanding. Not wired; `EoS_Nitrogen` is unreachable from `Gas` today.
Re-verified unchanged 2026-08-09 (currentness sweep). Nothing has moved since 2026-08-06:
`cl_GM_EoS_Nitrogen.cpp` is still 7050 bytes and the three bodies still do not compile …
deliberately absent from the `SOURCES` list …"

---

## State on disk (verified 2026-08-06 19:27)

`cl_GM_EoS_Nitrogen.cpp`, 7128 bytes. Only one copy exists on the machine and
only one git worktree (`sideconnectors`) — a save intended for this file did
not land, so anything typed into an editor buffer is NOT recoverable from disk.

**Correct and verified:** `compute_phi0`, `compute_phi0_d`, `compute_phi0_t`,
`compute_phi0_tt`, `compute_phir`, `compute_phir_d`, `compute_phir_dt`.
`phir_dt` was rewritten and finite-difference verified (agrees with FD of
`phir_d` to ~1e-10; the previous version was wrong by 100-650 %).

**Outstanding — will not compile:** `compute_phir_dd` (`:150` stray `*`),
`compute_phir_t` (`:175` `uint k = 0` should declare `kk`; `kk` then undeclared
at `:176`/`:178`), `compute_phir_tt` (`:198` stray `*`, undeclared `kk`, and
`:194` missing the leading `mJ(k) *` factor).

---

## D1 — the three outstanding bodies

Conventions in this class: groups are `k=0..5` polynomial, `k=6..31`
`exp(-δ^l)`, `k=32..35` Gaussian with `F = exp( -φ(δ-1)² - β(τ-γ)² )` and
**positive** `mPhi`/`mBeta` (the minus sits inside the exponential).
`mDeltaPowL(k)` holds δ^l, filled by `update_f` for `k=6..31` only.
Sign-checked against the equivalent methane code, which uses the opposite
convention (negative stored coefficients) and was verified numerically.

```cpp
        real
        EoS_Nitrogen::compute_phir_dd()
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k )
                        * mI( k ) * ( mI( k ) - 1. );
            }

            /* ∂²/∂δ²( δ^i exp( -δ^l ) ) = δ^i exp( -δ^l )
             *     * ( ( i - l δ^l )( i - 1 - l δ^l ) - l² δ^l ) / δ² */
            for ( uint k=6; k<32; ++k )
            {
                const real tU = mI( k ) - mL( k ) * mDeltaPowL( k ) ;

                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( tU * ( tU - 1. )
                            - mL( k ) * mL( k ) * mDeltaPowL( k ) );
            }

            /* ∂²/∂δ²( δ^i F ) = δ^i F
             *     * ( ( i - 2φδ(δ-1) )² - i - 2φδ² ) / δ² */
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                const real tU = mI( k )
                        - 2. * mPhi( kk ) * mDelta * ( mDelta - 1. ) ;

                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( tU * tU - mI( k )
                            - 2. * mPhi( kk ) * mDelta * mDelta );
            }

            alpha_r /= mDelta * mDelta ;

            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_t()
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mJ( k );
            }

            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * mJ( k );
            }

            // ∂/∂τ( τ^j F ) = τ^j F * ( j - 2βτ(τ-γ) ) / τ
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( mJ( k )
                            - 2. * mBeta( kk ) * mTau * ( mTau - mGamma( kk ) ) );
            }

            alpha_r /= mTau ;

            return alpha_r;
        }

        real
        EoS_Nitrogen::compute_phir_tt()
        {
            this->update_f();

            real alpha_r = 0. ;

            for ( uint k=0; k<6; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k )
                        * mJ( k ) * ( mJ( k ) - 1. );
            }

            for ( uint k=6; k<32; ++k )
            {
                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * mJ( k ) * ( mJ( k ) - 1. );
            }

            /* ∂²/∂τ²( τ^j F ) = τ^j F
             *     * ( ( j - 2βτ(τ-γ) )² - j - 2βτ² ) / τ² */
            uint kk = 0 ;
            for ( uint k=32; k<36; ++k, ++kk )
            {
                const real tU = mJ( k )
                        - 2. * mBeta( kk ) * mTau * ( mTau - mGamma( kk ) ) ;

                alpha_r += mN( k ) * mDeltaPowI( k ) * mTauPowJ( k ) * mF( k )
                        * ( tU * tU - mJ( k )
                            - 2. * mBeta( kk ) * mTau * mTau );
            }

            alpha_r /= mTau * mTau ;

            return alpha_r;
        }
```

- [ ] Paste the three bodies (or Christian's own, if the buffer is recovered —
      then diff them against these, do not assume either is right).
- [ ] Re-run the FD gate below before wiring.

## D2 — verification gate (run before wiring)

`scratchpad/n2_check.py` already mirrors this EoS in Python and FD-checks
`phir_dt` against `phir_d`. Extend the same way for the three new ones:
`phir_dd` vs d(`phir_d`)/dδ, `phir_t` vs d(`phir`)/dτ, `phir_tt` vs
d(`phir_t`)/dτ. Test at least (δ,τ) = (0.8,1.1), (1.0,1.16), (1.5,1.3),
(1.0,1.0) **and** a state where the Gaussians are alive — the earlier
`phir_dt` defect was invisible at (2.2,0.9) because both Gaussian exponentials
underflow to ~3e-13 there. Expect agreement at the FD floor, ~1e-9.

- [ ] FD gate green for all six `phir_*`.

## R1 — wiring (after D1 and D2)

Verified sites, 2026-08-06:

- [ ] **Constructor.** Today the class has `EoS_Nitrogen()` taking no parent,
      unlike its three siblings. Needs
      `EoS_Nitrogen( Gas & aParent ) : Helmholtz( aParent, "N2" )`, keeping the
      existing const-member initializer list, then `init_tables()` and
      `set_reference_point()`. Destructor should call `delete_cubic_eos()` as
      `EoS_Oxygen` does.
- [ ] **`init_tables()` does not exist.** Needs, from Span et al. 2000 p. 1364
      (values verified against the PDF): `mTtriple = 63.151`,
      `mTcrit = 126.192`, `mPcrit = 3.3958e6`,
      `mRhocrit = mM * 11.1839e3`, `mVcrit = 1.0 / mRhocrit`,
      then `set_critical_point_in_data_object()`. Validity range is
      63.151-1000 K to 2200 MPa, so `mTmax = 1000.0` and `mPmax = 2200.0e6`
      (the latter is now enforced by `BELFEM_ERROR` in `Helmholtz::v`).
      M = 28.01348 g/mol comes from the gas tables, do not hard-code it.
- [ ] **Vapor ancillary.** `mNvap` / `mKvap` are unset; `pi_vap`/`psi_vap` in
      the base need them, and `init_Tvap_poly()` must be called. Span gives the
      ancillary equations in §4.3; `mVliq` (liquid volume initial guess) also
      needs fitting — see how `EoS_Oxygen` does it.
- [ ] **Enum.** Add `Nitrogen` to `HelmholtzModel`
      (`en_Helmholtz.hpp`) **before** `UNDEFINED`.
- [ ] **Dispatch.** Add a `case( HelmholtzModel::Nitrogen )` to the switch in
      `cl_Gas.cpp:790-820`, constructing `EoS_Nitrogen( *this )` plus a
      `HelmholtzTransport( *this )` placeholder (as Oxygen and Hydrogen do —
      nitrogen has no dedicated transport class).
- [ ] **Includes.** `cl_GM_EoS_Nitrogen.hpp` into `cl_Gas.cpp` (alongside the
      three at `:31-33`) and into `helmholtz.cpp` (`:27-29`) if the example
      should cover it.
- [ ] Wherever the model is selected by name from input, add the `N2` string.

## O — open questions

- [ ] **O1** Dead `mSwap` (`cl_GM_EoS_Nitrogen.hpp:42,76`) — allocated, never
      read or written. Both auditors flagged it. Remove, or use it in the three
      functions above.
- [ ] **O2** `update_delta_pow_i` uses the shape `x = δ; multiply n-1 times`
      guarded by `BELFEM_ASSERT(n>0)`, while methane/oxygen/hydrogen use the
      zero-safe `x = 1; multiply n times`. Decide whether nitrogen adopts the
      zero-safe form for uniformity.
- [ ] **O3** Once nitrogen lands, unify all four classes on the base-class
      `ipow` helper — currently `mDeltaPowD` uses a naive loop while
      `mDeltaPowC` uses `ipow`, which Grok flagged as two algorithms for one
      job. Deferred deliberately so as not to create a third variant.
- [ ] **O4** Nitrogen's τ exponents sit on a 1/8 grid, so the `hpow`
      half-integer trick does not apply directly; a `qpow` on `τ^(1/8)` would
      need its own accuracy argument (methane/oxygen measured 9 ulp on the
      1/2 grid; three nested roots will be worse).
- [ ] **O5** `compute_phi0_d` is declared without `override` because the
      Helmholtz base declares no such virtual — it is currently dead code.
      Either wire it or drop it.
