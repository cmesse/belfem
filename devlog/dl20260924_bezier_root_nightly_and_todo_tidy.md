# One clamped Bézier root for the piecewise power law, G5 on the 2-D tape stacks, and five plans closed

**Date:** 2026-09-24
**Purpose:** Three threads of one session. (1) The `todo/` tidy-up Christian asked for while a peer
session held the materials cluster: five finished plans moved to `closed/`. (2) The nightly's one
failing test, relayed by the CI session: `rho_piecewise` and `drho_piecewise_dT` returned NaN at the
coincident knots; fixed with one shared helper behind plan and code juries. (3) The G5 gate owed by
the thin-shell recovery plan: the post-fix postprocessor against the pre-fix one on the two shipped
2-D (QUAD4TS) tape-stack decks.
**Module:** `src/physics/materials` (`powerlaws.hpp`, `cl_Material.hpp`), `tests/physics`; `todo/`;
`src/fem/kernel` (gate only, no edit)
**AIs involved:** Claude Fable 5.1 (plan, implementation, gates, commit). Plan jury and code jury:
Codex `gpt-5.6-terra`/high, Grok `grok-4.7`/high, `tmp/ai_exchange/review_piecewise_bezier_root.md`.
The CI session (a peer Claude on the runner side) supplied the failure diagnosis.
**Status:** COMPLETE — fix landed `aa468c3`; five plans closed; G5 passed on `tapestack2d_layered`, the racetrack deck raised O1 (§3–4)
**Verification:** executable — `test_physics` 151/151 in release and debug Armadillo scratch builds of the fix (red-then-green on the new sweep); G5 A/B: layered deck ≤ 2.1e-9 with a bitwise-zero floor over 150 saves (pass); racetrack deck diverges through the solve from 1.4 ms (§3, O1). Not run: `make check-fast` in the shared trees, the pipeline (Christian).
**Exchange:** `tmp/ai_exchange/review_piecewise_bezier_root.md` (plan and code rounds), `.diff`

---

## 1. The `todo/` tidy-up

The directory rule: a plan stays active only while code or tests are missing from the tree. Five
plans had landed and were carried by their unchecked boxes alone. Moved to `todo/closed/`, each
with a `CLOSED 2026-09-24` Status line naming its residue, and the index rewritten:

| plan | why closed | residue |
|---|---|---|
| `smith_inverse_grouped_terminals` | report to Gregory delivered 2026-09-19; branch stays unmerged by design | none |
| `test_gap_tier1_execution` | run closed 2026-09-22; its rulings are the rows of `test_gap_followups.md` | none |
| `citation_index_completion` | checker green, 0 unresolved | two rulings for Christian: gate or not (R6), an `unheld` marker (O2) |
| `thinshell_nodal_recovery` | landed `5ab4d41`, G1–G4 passed | G5, a QUAD4TS deck — run below |
| `thinshell_recovery_unit_test` | landed `5ab4d41`, 11 tests, mutations each red on their own | T5, the optional 2-D analogue, not built |

The Tier 1 run had no index entry; it has one now. Links between the moved files and the follow-ups
plan were repointed; devlog links were left as written. `materials_tangent_cluster.md` and both
registers were the peer's live files and were not touched. Still active: the Tier 2 run, the
test-gap backlog and follow-ups, the periodic-postprocessor plan, the mutation-suite follow-ups.

## 2. The nightly's NaN at the coincident knots

**Report.** Pipeline 41 (job 150, RELEASE/Armadillo) failed `HtsLaws.PiecewiseDegenerateBlendStaysFinite`
at `test_MaterialDerivatives.cpp:570` and `:572`: value and dT at J = j3 for n = 4 non-finite, the
dJ line fine. The CI session's diagnosis: with p = n − 1 equal to q = mNff the Bézier discriminant
b² + ac is exactly 0 at the upper knot, roundoff lands it below, `pow( d, 0.5 )` is NaN, and only the
value path is unguarded.

**Confirmed at HEAD, with one addition.** All eight `rho_piecewise` overloads took the root unguarded.
The dT overloads *have* a guard, but it sits after `tSqrt`, `tParam`, `rhoFF` and `tFrozen` are
computed from the root, and returns `tFrozen`, which already carries the NaN; its comment claimed the
negated comparison "takes the fallback for NaN as well", which protects the branch and not the value.
That is why both lines failed. The dJ overloads return the flux-flow limit before any root (yesterday's
`9373e44`); one of them then took a complex square root.

**Why not here.** A release Armadillo scratch build of the pristine tree passes the test on this
machine. With the test's constants `logj2 == log10( j3 )` bitwise, and `( logj1 − 2·logj2 ) + log10( j3 )`
is exact while the partial sum stays inside the binade of `logj1`, so `a == b`, `c == −b`, and the
discriminant is +0. It goes negative when that sum crosses 16, from jc ≈ 2.5e11 at n = 4, or when a
compiler contracts `b*b + a*c` into a fused multiply-add and the result is the rounding error of `b²`.
The runner does one of those; which one was not determined and the fix does not depend on it. A jc sweep
through the in-tree overloads on the pristine library: **63 of 488 points NaN**, value, value( J ) and
dT together at every one.

**Fix (`aa468c3`).** `Material::piecewise_bezier_sqrt( a, b, c )`, a protected static defined inline at
the top of `powerlaws.hpp`: `sqrt( d < 0 ? 0 : d )` with `BELFEM_ASSERT( !( d < −BELFEM_EPSILON ( b² + 1 ) ) )`,
the band the existing dJ and dT guards use. It serves all 26 root evaluations in the 18 functions: 8
value sites, 7 dJ sites with two real `pow` each, the eighth dJ site whose `cplx` intermediates and
`BELFEM_ERROR` on the imaginary part are gone (behind its guard the root was real), and both dT
copies, whose six-line comment tail became three true lines. The clamp is legal because the
discriminant is affine in c and both blend endpoints are squares, so it is non-negative on the whole
blend in exact arithmetic; a negative value is roundoff or a broken invariant, never a second root.

**Tests.** `PiecewiseDegenerateBlendFiniteAcrossJc`: n = 4, jc 1e9…1e15 (61) × T 20…90 K (8), value,
value( J ), dT and dJ finite at the knot and the value equal to `rho( T )` within 1e-9 — red on the
pristine library (63 points), green now, 24 ms. `PiecewiseBezierRootClampsRoundoff`: the clamp's
contract through a derived-type `using` seam (`d = ∓4.4e-16` exactly by Sterbenz, a plain root, a NaN
passing through, and the assert's throw under `BELFEM_ASSERTIONS_ACTIVE`), so the pipeline's verdict
no longer depends on which side of zero a platform rounds to.

**Juries.** Plan round: the helper snippet's linkage (a `static` at an out-of-class definition), the
complex eighth site the plan had counted as `pow`, the value-only endpoint assertion (the blend-side
dT at the knot is 0 for constant jc and n, the normal side is `drhodT`; the two differ, so only the
value is asserted), an unsupported "correctly rounded" claim withdrawn, and the debug assert. Code
round: no blocking defect; `std::max( d, 0.0 )` is not specified for NaN, so a ternary; the helper's
comment overstated the band; one material per temperature in the sweep (1.5 s → 24 ms); the index `m`
renamed. Accepted residue: a discriminant below the band aborts in debug on the value and dT paths and
is the silent flux-flow limit on the dJ paths, whose guard returns first — unreachable for legal knots,
and release agrees across all three.

**Gates run by Claude** (scratch worktrees, own build directories): `test_physics` 151/151 in release
and debug Armadillo, zero warnings under `-Wall -Werror -pedantic-errors`; Blaze release flags,
syntax-clean. **Owed to Christian:** `make check-fast` in both shared trees, and a manual pipeline on
devel before tonight's schedule. The CI session was sent the commit.

## 3. G5: the post-fix recovery on the 2-D tape stacks

Scratch worktree at HEAD `ad8064f`, its own Blaze/MKL release build (the `cmake-build-debug`
configure), serial. A = HEAD with `cl_FEM_Postprocessor.{hpp,cpp}` from before `5ab4d41`, B = HEAD,
B2 = HEAD again for the floor. Both shipped 2-D decks carry `thinshell : tape`, which the mesh types
as QUAD4TS, so both exercise the curvilinear shell mode; both have straight tapes, so G5 is an
agreement gate like G2 (tolerance 1e-3 of the field), not a curvature gate. Comparison per block and
nodal field, max over nodes of |A − B| / max|B|, max over the shared saves.

**`tapestack2d_layered`** (11 tapes, seven resolved layers, constant jc = 1e10, n = 20, 160 A ramp,
15 s, 150 saves): **PASS.** A → B at most 2.1e-9 (buffer B), 1.4e-9 in YBCO for B, H, J and J/Jc,
1.3e-10 in air; the B → B2 floor is **bitwise zero** on every field at every save, so the 1e-9 is
the postprocessor change itself: the old absolute-coordinate fit's round-off on a millimetre mesh.
The solve is untouched (identical step sequence, identical φ).

**`tapestack2d_racetrack`** (one 1 µm YBCO layer, `sp-ap.hdf5` table, copper bulk, 3040 A at 50 Hz,
cut to 5 ms, 25 saves; 15 shared at close, both B runs still going in the scratch): **not a
postprocessor comparison after 2.3 ms, and a finding before it.** A → B is 1e-10 (copper J) to
2e-8 (YBCO J/Jc) through 1.2 ms, then grows: 3e-7 at 1.4 ms, 1e-5 at 2.2 ms, 1e-3 at 2.4 ms, 4e-2
in copper J and 2e-2 in J/Jc at 3 ms, while B → B2 stays bitwise zero. The step sequences are
identical to step 38 (2.2 ms); at step 37 A took 21 Newton iterations and B 23, the controller's
next Δt differs (0.1055 vs 0.1044 ms), and from step 39 the two runs are different BDF trajectories
(78 steps vs 49 to 3.2 ms). So on this deck **the recovered nodal fields feed the solve**: with the
same Δt the residual differs at 1e-7 from 1.4 ms, which only the postprocessor can have caused, and
the discrete controller amplifies that into a step-size fork. The layered deck, whose jc is a
constant, shows no such coupling. The consumer was not traced in this session (the Maxwell IWGs'
`field_data` reads are all in `custom_postprocess`; the candidate is the table's jc( B, θ ) or the
thin shell's normal field through a nodal quantity). Which trajectory is "right" is not a question
G5 answers; B's fit is the better conditioned one on a metre-scaled mesh.

**Verdict:** G5 passed on the deck that isolates the postprocessor; the racetrack deck exposes a
dependency of the solve on the recovered fields that the rutherford A/B could not have seen (its
solver was not bitwise reproducible, U₂ 3e-5 run to run). Filed as an open question, not a defect,
until the consumer is named.

## 4. Open questions

- **O1** Which quantity of the h-φ solve reads a postprocessor-recovered nodal field on the
  racetrack deck, and is that by design? Evidence: identical Δt to 2.2 ms with 1e-7 divergence from
  1.4 ms, floor 0. A bisect over the deck's features (constant jc in place of `sp-ap.hdf5`; the
  copper bulk removed) would name it in two runs of 3 ms each.
- **O2** (from the closed citation plan) whether `check_citations.py` becomes a gate, and whether
  the index marks unheld works. Rulings for Christian.
- **O3** The racetrack B and B2 runs finish in the scratch after this close; their last ten saves
  are not in the table above and would not change the verdict.

## 5. Files

- `src/physics/materials/powerlaws.hpp`, `cl_Material.hpp`, `tests/physics/test_MaterialDerivatives.cpp` — `aa468c3`
- `todo/piecewise_bezier_root.md` (plan, closed), `todo/closed/{smith_inverse_grouped_terminals, test_gap_tier1_execution, citation_index_completion, thinshell_nodal_recovery, thinshell_recovery_unit_test}.md`, `todo/README.md`, `todo/test_gap_followups.md`
- `tmp/ai_exchange/review_piecewise_bezier_root.{md,diff}` — both rounds, pre-registrations frozen
- scratch: `g5/` (worktree, A/B runs, `g5compare.py`), `fix/` (worktree, release and debug builds), `disc/disc.cpp`
