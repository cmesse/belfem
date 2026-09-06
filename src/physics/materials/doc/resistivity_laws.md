# HTS Resistivity Laws: powerlaw, piecewise, riva {#physics_materials_resistivity_laws}

**Date:** 2026-08-27
**Purpose:** Reference for the three E-J constitutive laws a superconductor can use in BELFEM —
what each one computes, how they differ, how the solver consumes them, and how to choose.
**Module:** `src/physics/materials` (laws), `src/fem/kernel` (dispatch)

---

## 1. Where the resistivity law sits

The h-φ solver never sees a critical current directly. The weak form uses the local resistivity
ρ(|J|, T, |B|, θ) of the superconducting layer, while the Newton tangent uses its three derivatives:
dρ/d|J|, dρ/d|B| and dρ/dT. The *law* is the rule that turns the material data — the critical field
criterion `ec`, the critical current density `jc` and the transition exponent `n` — into those four
numbers.

All three laws share the same inputs:

- `ec` — constant, default `1e-4` V/m (the 1 µV/cm criterion).
- `jc(T, |B|, θ)` and `n(T, |B|, θ)` — either constants from the deck or measured tables loaded
  from an HDF5 file (the `file` key). Both are evaluated through the same dependency routing
  (`Material::jc_eval` / `n_eval`, `powerlaws.hpp`). Table lookups are clamped to the table hull
  in all three coordinates.
- θ is the **field-to-tape-normal angle** (`bn_angle`, unfolded to [0, π] since 2026-08-16). The
  measured jc(θ) tables are asymmetric about π/2 and consume the angle as is.
- The intrinsic power-law channel shared by all three laws is Rhyner's E-J power law
  (Rhyner 1993):

      ρ_PL(J) = (ec / jc) · (|J| / jc)^(n−1)

The deck selects the law separately for each material:

```
materials
{
    ybco
    {
        builtin : ybco ;
        file : sp-ap.hdf5 ;
        resistivity type : riva ;    // powerlaw | piecewise | riva
    }
}
```

`MaterialFactory` parses the key into a `ResistivityLaw` enum. During setup, the FEM `Calculator`
binds the four evaluation channels (ρ, dρ/dJ, dρ/dB, dρ/dT) to the corresponding member functions
(`cl_FEM_Calculator.cpp`, four dispatch sites: thin-shell and bulk, each with and without a defect
function). The default is `powerlaw`.

**The ohmic floor (2026-08-27).** Measured n(T, |B|, θ) tables soften through n = 1 near T_crit.
That is real physics: the transition approaches ohmic behavior. Below n = 1, however, the raw
power law is sub-ohmic and its J → 0 limit flips. `Material::n_eval` therefore floors the exponent
at 1 globally for every law. While the floor binds, `dn_eval_dB` and `dn_eval_dT` return exactly
zero, so the tangents differentiate the same clamped law as the residual. At n = 1, the power-law
channel becomes a plain resistor, ρ_PL = ec/jc, independent of J.

---

## 2. powerlaw — the parallel model (default)

    ρ = ( 1/ρ_n + 1/ρ_PL )^(−1)

This law places the power-law channel in parallel with the normal-state channel ρ_n(T) (Duron et
al. 2004). This standard construction makes Rhyner's law usable "in an arbitrary current range."
Far below jc, the parallel combination is indistinguishable from ρ_PL: the relative difference is
ρ_PL/ρ_n, about 4·10⁻⁸ at 1.05 jc for 77 K REBCO constants. Past jc, it saturates smoothly to ρ_n.
With n ≈ 19, saturation is essentially complete by |J| ≈ 3 jc.

The derivative legs carry the parallel weights. With w = ρ_n/(ρ_PL + ρ_n),

    ∂ρ/∂ρ_PL = w² ,   ∂ρ/∂ρ_n = (1−w)²

Therefore, dρ/dT = w²·dρ_PL/dT + (1−w)²·dρ_n/dT. The dρ/dJ and dρ/dB legs scale the
corresponding ρ_PL derivatives by w² because ρ_n does not depend on J or, on this path, on B.

**Caveat — not total at the extremes.** The residual degrades gracefully when ρ_PL overflows
because the reciprocal form returns ρ_n. The derivative legs, however, evaluate `ρ_PL` and its
ratios directly and can produce NaN when (n−1)·log₁₀(J/jc) overflows the double range. This
requires an extreme combination: it was observed only at a table-corner artifact with n ≈ 86 at
|B| = 10 T together with |J| ≳ 2·10¹², far outside self-field tape simulations. This limitation is
the reason `riva` exists.

---

## 3. piecewise — three regimes with a flux-flow blend

The piecewise law divides the resistive transition into three explicit regimes:

1. **Power-law regime**, |J| ≤ j1 with j1 = jc·10^(2.5/n): the raw Rhyner law ρ_PL, *not*
   parallel-combined.
2. **Flux-flow blend**, j1 < |J| ≤ j3: a quadratic Bézier curve in log ρ – log J space connecting
   (j1, ρ1) to (j3, ρ_n), with the upper knot j3 = j1·(ρ_n/ρ1)^(1/n_ff) set by the flux-flow
   exponent n_ff (member `mNff`, default 3).
3. **Normal regime**, |J| > j3 or T > T_crit: ρ = ρ_n(T) exactly.

The transition is much wider than in the parallel model. With 77 K REBCO constants, the blend
extends from ≈ 1.35 jc to ≈ 76 jc, and the two laws differ by up to two orders of magnitude within
that window. Determining which shape is closer to reality requires measurement. Pulsed-current
data on REBCO tapes (Riva 2021, Ch. 3–5) show that the true overcritical resistivity rises *more
slowly* than the raw power law. Neither the hard parallel knee nor the Bézier stretch reproduces
that behavior exactly.

**Degenerate windows.** The Bézier construction assumes that the power-law slope remains steeper
than the flux-flow slope (n − 1 > n_ff). A measured table violates this assumption as n approaches
1 near T_crit. The machinery then degenerates in two distinct ways: the blend's discriminant
vanishes when n − 1 = n_ff (n = 4 at the default n_ff), and the middle knot
j2 = j1·(ρ_n/ρ1)^(1/(n−1)) overflows as n → 1⁺. The 2026-08-27 quench-onset abort of the
tapestack3d deck came from this second window: NaN passed through the dT leg's degeneracy guard
because the guard computed its fallback before testing. Hardening these windows in place is
planned but **not yet landed** (see `todo/fix_piecewise_degenerate_window.md`, step R8). Until it
lands, decks whose tables reach n ≤ 4 near T_crit should prefer `riva`. Piecewise also still asserts
n > 1 in debug builds, so the global n-floor causes an abort at exactly n = 1. This is by design
until R8 relaxes the assertion.

---

## 4. riva — the parallel model, made total

`riva` computes the *same* parallel model as `powerlaw`. For ordinary inputs, the two return
identical values. Unlike `powerlaw`, `riva` is hardened to stay finite for every input that a
measured table, a defect function, or a wild solver iterate can produce:

| situation | behavior |
|---|---|
| `jc_eff ≤ 0` or nonfinite (dead defect `D = 0`, spline underflow, NaN) | fully normal: ρ = ρ_n, dρ/dT = dρ_n/dT, dρ/dJ = dρ/dB = 0 |
| ρ_PL past the overflow cap (evaluated in log₁₀ space) | fully normal, residual and all tangents consistently |
| n floored at 1 (measured softening near T_crit) | ohmic closed form ρ_PL = ec/jc, J-independent — evaluated without `pow`, and the dB/dT tangents keep tracking jc(B, T) |
| NaN or infinite n, negative `ec` | land in the NaN-aware overflow guard → fully normal |
| `|J| < ε` with n > 1 | ρ_PL = 0 exactly; dρ/dJ = 0 |

The name credits Nicolò Riva's thesis (Riva 2021), which uses the Duron parallel construction in
its models of the overcritical current regime and analyzes the continuity of the relevant limits.
BELFEM makes two deliberate deviations from the thesis formulation, both documented in
`doc/input_file_reference.md`:

- BELFEM keeps its **floor** semantics for the minimum resistivity (`mRhoMin`, zero by default)
  instead of the thesis's additive 10⁻¹⁷ Ω·m regularization — an additive term is
  tangent-consistent, but a positive *floor* combined with an unfloored tangent is the exact
  value/derivative desynchronization that was removed from BELFEM on 2026-08-10, and mixing the
  two conventions in one code invites it back.
- The n-floor at 1 is BELFEM's totalization of measured sub-ohmic table values, not part of the
  published model.

`riva` is available through the assembly-path signatures (|J|, T, |B|, θ), with and without defect
modulation. Its jc and n values follow the same routing as those of the other laws, whether they
come from tables or constants.

---

## 5. Choosing a law

| | `powerlaw` | `piecewise` | `riva` |
|---|---|---|---|
| model core | Rhyner ∥ ρ_n (Duron) | Rhyner + Bézier flux-flow + ρ_n | Rhyner ∥ ρ_n (Duron) |
| transition width | knee at ρ_PL = ρ_n, saturated by ≈ 3 jc | stretched to ≈ 76 jc | as powerlaw |
| valid n range | n > 1 practical (no assert, no totality guarantee) | n > 1 asserted | runtime: any table output (floored at 1); setup refuses a stored table min <= 1 or a constant n not finite-and-greater-than-1 (2026-08-31) |
| behavior at table edges / dead defects | residual safe, tangents can NaN | NaN / abort windows until R8 | total |
| since | original | original | 2026-08-27 |

Practical guidance:

- **Sub-critical operation (|J| ≲ jc everywhere):** all three laws agree to within ρ_PL/ρ_n. The
  choice is irrelevant to the physics, and `powerlaw` (the default) is suitable.
- **Quench and current-sharing studies with measured jc/n tables**, where the solver visits
  T near T_crit and |J| past jc during an iteration: use `riva`. The tapestack3d abort occurred in
  this regime, where finite tangents help keep Newton alive through wild iterates.
- **Legacy piecewise decks** away from the degenerate windows reproduce their previous results
  unchanged. If a deck's table reaches n ≤ 4 near T_crit, switch to `riva` or wait for R8.
- The transition-zone difference between the parallel and piecewise shapes (up to ~200× in ρ)
  is a *model* difference, not a numerical one. If loss or quench-speed predictions matter in the
  overcritical window, validate the choice against measurements. Riva 2021 (Ch. 6) shows that
  power-law-shaped models quench faster than measured tapes.

---

## 6. Newton-tangent conventions

- Each law's tangent legs differentiate **that law's own residual, branch by branch**. Fallback
  branches also match: wherever the riva residual returns ρ_n, its tangents return dρ_n/dT and
  zeros.
- The β (angle) channel of the tangent is deliberately **not** bound for HTS materials. The HTS
  angle is `bn_angle` (field to tape normal), while the generic angle-tangent machinery
  differentiates `bj_angle` (field to current). Binding it would apply the wrong ∂β/∂q rows.
  (`cl_FEM_Calculator.cpp`, dispatch comment, 2026-08-13 audit.)
- While the n-floor binds, dn/dT = dn/dB = 0 by convention (one-sided at the kink). The floor is
  C0 but not C1 in raw n. The branch *switch* at T_crit is exact in every law, but continuity of
  the *value* across it is not automatic: it holds only when jc has become negligible by T_crit,
  and a table whose hull ends below T_crit freezes jc at the edge value, which can leave a jump.
  Check this per table.

## 7. Pitfalls

| pitfall | symptom | rule |
|---|---|---|
| raising `mRhoMin` above 0 | Newton stalls below ≈ 0.87 jc | the floor desynchronizes value and tangent; leave at 0 |
| `ec ≤ 0` in a deck | riva silently fully normal; other laws NaN | keep ec > 0; a loud setup check is a planned follow-up |
| piecewise + table reaching n ≤ 4 near T_crit | debug assert at n = 1, or NaN in the blend windows | use riva until R8 lands |
| assuming θ is folded at 90° | wrong jc lobe on measured tables | θ is unfolded [0, π]; the sideset sign selects the lobe |
| reading `Material::n()` for physics | unfloored table value | the floored exponent is `n_eval`; `n()` is the raw accessor |

## 8. References

- J. Rhyner, "Magnetic properties and AC-losses of superconductors with power-law current-voltage
  characteristics," Physica C 212 (1993) 292–300. — the E-J power law.
- C. P. Plummer, J. E. Evetts, "Dependence of the shape of the resistive transition on composite
  inhomogeneity in multifilamentary wires," IEEE Trans. Magn. 23 (2) (1987) 1179–1182. — the
  n-value phenomenology.
- J. Duron, F. Grilli, B. Dutoit, S. Stavrev, "Modeling the E–J relation of high-Tc
  superconductors in an arbitrary current range," Physica C 401 (2004) 231–235. — the parallel
  combination used by `powerlaw` and `riva`.
- N. Riva, "Quench behavior of high-temperature superconductor tapes for power applications"
  (EPFL thesis 8754, 2021). — measured overcritical resistivity 77–90 K (Ch. 3–5), the ρ_ηβ
  model and continuity analysis (§5.1), quench-speed impact of the law choice (Ch. 6).
- C. Messe et al. 2023, §2.6–2.7 — the BELFEM material database and the convergence strategy the
  tangents serve.

Full citations: `doc/literature_references.md`.
