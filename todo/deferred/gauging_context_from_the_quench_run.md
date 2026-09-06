# Context for the gauging / stabilization work: what the 2026-08-21/22 quench run measured

**Date:** 2026-08-22
**Purpose:** Hand the measured behaviour of the near-critical h-φ system to whoever
implements a gauge/stabilization term, so the design starts from evidence rather than
from the usual assumptions about ill-conditioning.
**Status:** evidence only, no design. Nothing here is a recommendation to gauge.

## Why this file exists

The tapestack3d coarse-mesh campaign spent 2026-08-21/22 chasing a timestep collapse
near I ≈ 1.1·Ic. Four candidate causes were excluded BY MEASUREMENT, and the one that
turned out to matter was not conditioning. If a gauging term is added, these are the
numbers it has to be judged against, and the one open question that decides whether it
helps at all.

## Measured: κ(A) of the magnetic system

`compute conditioning : true` in `linear magnetic` ( MUMPS supplies κ from its own error
analysis, free; do NOT put the key in a nonlinear block, it is not read there ).
97 completed steps, t = 6350.7 … 6386.3 ms:

| quantity | value |
|---|---|
| range | 1.1e17 … 7.4e19 |
| median, first third | 7.4e17 |
| median, middle third | 4.2e17 |
| median, last third ( into the collapse ) | 6.8e17 |
| median on steps converging in ≤ 6 iterates | 5.9e17 ( n = 85 ) |
| median on steps needing ≥ 12 iterates | 3.0e17 ( n = 5 ) |

**Two negative results, both load-bearing:**

1. **No trend into the failure.** κ is flat across the whole approach to the wall.
2. **κ ANTI-correlates with difficulty** — hard steps had roughly HALF the κ of easy
   ones ( ratio 0.51 ). The worst-conditioned step of the whole set ( 7.4e19 ) converged
   in five iterates to −113.7 dB.

So κ ~ 1e17–1e19 is the resting state of this formulation, not a pathology, and it does
not predict which steps are hard. **A gauging term that only improves κ would, on this
evidence, change nothing about solver behaviour here.** That is the central warning.

## Measured: what the collapse actually was

The residual plateaued and refused to fall below roughly −90 to −100 dB. Excluded:

| candidate | how it was excluded |
|---|---|
| timestep | plateau identical from Δt = 0.016 ms down to 1e-7 ms ( 5 orders ) |
| linear solver | STRUMPACK and MUMPS both, same plateau |
| linear precision | rel tol 1e-10, 1e-12 ( maxit-capped, undelivered ), and MUMPS exact backsolve — same plateau |
| conditioning | the table above |
| nonlinear tolerance | **this was it** |

Loosening the magnetic nonlinear tolerance 1e-11 → 1e-9:

| | 1e-11 arm | 1e-9 arm |
|---|---|---|
| reached | died at 6386.35 ms | past 6407.6 ms, running |
| rejections | 41 / 139 attempts | 3 / 139 |
| watchdog cuts | 24 | 3 |
| median iterates | ( collapsed ) | 4 |
| accepted magnetic final | — | **−90 … −96 dB, and it stops there** |

The accepted steps in the 1e-9 arm converge to a genuine plateau just past the criterion.
They are not being accepted early; they stop improving at ≈ −92 dB.

## THE question a gauging design should answer first

**Is that −92 dB plateau living in the gradient null space of the curl-curl operator?**

- If YES: gauging removes it, 1e-11 becomes reachable again, and the tolerance
  loosening is a workaround that gauging retires. This is the case worth building for.
- If NO: gauging improves κ ( a number already shown not to matter here ) and the
  plateau stays. The formulation would then be limited by something else, and the
  stabilization would be solving a problem this deck does not have.

**Cheap test, no new code:** take a converged state from the 1e-9 arm, form the residual,
and check whether it is concentrated on cotree edges / in the range of the discrete
gradient. A projection or a tree-cotree split of the residual vector answers it directly.
Do this before designing the term.

## Physical state where all of this happens

I ≈ 1846 A against a self-field-degraded stack Ic ≈ 1670 A ( **≈ 1.1·Ic** ), T_max ≈
+273 mK above the 77 K bath, dissipating volume frozen at 375 911 elements for many
saves. Every element sits on the steep flank of the E–J power law with local n ≈ 17–22.
Whatever the stabilization does must survive that regime, not just a benign one.

## Literature situation: the corpus will NOT help

Searched `literature/books/*.txt` and `literature/papers/fem/*.txt` for
`coulomb gauge`, `gauging`, `tree-cotree`, `tree cotree`:

- **One hit total:** `papers/fem/dular2021.txt:541-551`, and it is a 2-D special case —
  `a_δ = a ψ_n z` satisfies `div a_δ = 0` **automatically by construction** for a
  perpendicular edge function. Not a general 3-D scheme.
- `dular2021.txt:380-385` gauges the current vector potential **t** by choosing it along
  the tape normal ( `t = t n`, citing its ref [7] ) — again a construction, not a
  penalty or a tree.
- **Monk 2003 has no hits** under these terms, which is worth knowing before anyone
  goes looking there. No general gauging treatment in Bathe, Zienkiewicz, Belytschko,
  Hughes, Brenner, Boffi, Arnold either.

So the theoretical grounding has to come from outside the local library ( Bossavit and
the tree-cotree literature are the obvious starting points, unverified here ). Under the
project's literature-first rule that means the gauging design needs sources FETCHED and
cited, not paraphrased from training knowledge.

## BELFEM-specific constraints the design must respect

- **Cuts and cohomology.** The air region carries cohomology generators ( the cut around
  the stack, plus the free axial generator from periodicity ), and the current boundary
  condition consumes one generator per terminal group. A tree-cotree gauge interacts
  with those generators directly. See the `topology` block comments in the campaign deck.
- **Thin-shell hanging edges and the layer stack.** Enriched meshes carry twin edges that
  share both end nodes, so node-pair keys are NOT unique. Any tree construction keyed on
  edge endpoints will be wrong on this mesh.
- **A penalty term changes the discrete operator**, and this deck's jc / lobe / thin-shell
  apparatus is calibrated against the current one ( the `bn_angle` lobe selection is
  physics here, not a convention — see DR-69 ). Tree-cotree is exact by comparison
  ( it removes dofs rather than adding terms ), at the cost of the interactions above.
- **The residual is `‖Ax − b‖/‖b‖`** ( Messe et al. 2023, Eqs. 10–11 ) and is evaluated
  PRE-update on the Picard path. If a stabilization term changes what `A` is, the
  meaning of every residual number quoted above changes with it.

## Related records

`tmp/ai_exchange/jjc_bias.md` ( postproc/assembly jc consistency, three-AI audit ),
`todo/debt_register.md` DR-88 ( raw factor gives ~4 digits on this matrix ),
`src/fem/kernel/doc/nonlinear_controller_theory.md` ( residual definition, watchdog ),
`doc/parallel_execution.md`.
