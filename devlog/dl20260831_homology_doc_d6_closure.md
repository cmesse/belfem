# Devlog 2026-08-31 — Homology Documentation: D6 Closed, Plus a Second Defect of the Same Class

**Date:** 2026-08-31
**Purpose:** Close the D6 residue left open by the documentation repair campaign — the homology
usage guide still teaching a cleanup pipeline that no code implements — and record a second,
previously unrecorded defect of the same class found while doing it
**Module:** src/homology (documentation only)

---

## Authorization and scope

Christian authorized editing the **documentation of the homology module**, explicitly not the
source. Nothing under `src/homology/*.{cpp,hpp}` was touched; the cohomology core edit ban
(`doc/ai_collaboration_protocol.md` §7.1) was never approached. Files changed:

| File | Sites |
|---|---|
| `src/homology/doc/homology_usage_guide.md` | 7 (4 for D6, 3 for the debug-output defect) |
| `src/homology/doc/README.md` | 1 |
| `src/homology/doc/handoff_for_gregory_20260831.md` | updated, new item **d** |
| `todo/doc_currentness_fixes.md` | D6 ticked, one stale claim inside it corrected |
| `todo/handoff_20260831_open_defects.md` | F4 struck |

## D6 — the residue was four sites, not two

The prior session's note named `homology_usage_guide.md:335` and `:1329-1350`. Two more existed:
the **Internal Workflow** step 3, and — the one that matters — the troubleshooting entry at
`:2205-2208`, which instructed the reader to *tune* the nonexistent filter:

```
3. **Manual manifold filtering:**
   - Check `manifold_filter_3d()` parameters
   - Adjust cycle density threshold (<0.3 → <0.4)
```

That is the most actionable form the error can take: not a stale description a reader skims, but a
debugging procedure they would follow, fail at, and blame themselves for. All four now carry the
correction. The "Manifold Cleanup" section gained the same design-note blockquote the two sibling
documents already used, so the module speaks with one voice; the design text was **marked, not
deleted**, on the standing assumption that it may be the plan of record.

One claim inside the tracker was itself wrong and is corrected: D6 asserted the Tarjan phases had
"zero hits, `archive/` included". Tarjan **does** exist — `archive/graph/fn_Graph_tarjan.{hpp,cpp}`,
outside the build. That was the `ugrep --ignore-files` false negative (F5) reaching one more
artifact than F5's audit had found. The glossary row now says "archived prototype" rather than
implying it ships.

## Second defect, same class: debug outputs that no code emits

Found by checking an adjacent claim rather than by being told to. `README.md` and the guide
advertised three debug files:

- `thick_cut_*.vtk`, `thin_cut_*.vtk`, `manifold_*.vtk` — **no code writes any of these names.**
  `manifold` appears **zero times** across `src/homology/*.{cpp,hpp}`.

What actually exists: `CutProcessor::save_debug_meshes()` (public, `cl_CutProcessor.hpp:116`)
writes `cut_<index>.vtk`, one per cut (`cl_CutProcessor.cpp:139-147`). The guide additionally
called it on `CutFactory`, which has no such member — already noted for Gregory at
`handoff_for_gregory_20260831.md:39`, but the guide's own snippet had not been fixed.

A practical consequence worth Gregory's eye, recorded in his handoff rather than acted on: the
factory's only call to `save_debug_meshes()` is **commented out** (`cl_CutFactory.cpp:471`), so the
module's debug output is unreachable without a source edit and a rebuild. That is a source
question, so it was left alone.

## The completeness check

Reading was not trusted to have found everything. All 44 code identifiers asserted in backticks
across the six homology documents were extracted and diffed against the concatenated sources of
`src/` and `archive/`. Three were absent:

| Symbol | Verdict |
|---|---|
| `manifold_filter_3d` | real defect — fixed |
| `check_surface` | real defect — fixed |
| `rectify_to_unit_or_certify` | **false positive** — both sites already label it "Proposed", and `thin_cut_nonunit_rectification.md:343-350` carries a banner naming the shipped path |

No further missing symbols. Every search used `command grep`, not the shell wrapper.

## Posture

**Reviewed, not verified.** Nothing was compiled or run; these are documentation edits and no gate
applies to them. The content decision — whether the manifold-filtering design is revived, retired,
or kept as the plan of record — remains Gregory Giard's, and the handoff says so.

## Rulings received, and the two plans they produced

Christian ruled on both open items the same evening; neither was implemented, both are planned for
2026-09-01.

**F2 → wire it.** `todo/parmetis_ptscotch_wiring.md`. The consequence worth stating: F1 stops being
a latent curiosity the moment there is a caller, so it is **D1** of that plan and a precondition for
R3, not follow-up work. Confirmed while planning that `gNoOwner` is `INT_MAX` = 2147483647 exactly
(`typedefs.hpp:42,59`, `proc_t` is `int`), which is what makes the release-build increment an
in-range-looking out-of-bounds write rather than a reliable fault. Build gates
`BELFEM_PARMETIS`/`BELFEM_PTSCOTCH` already exist, so the wiring has its guards. The test gap is the
real cost: nothing in `make check` enters the path and the suite runs multi-rank only for
`tests/comm`.

**Q1 → a real gap.** `todo/ac_loss_postprocessing.md`. Christian's framing was the useful part: the
heat imposed on the structure is already computed, so summing it is about half an hour. Confirmed —
`rho * norm_j * norm_j`, already carrying `w(k)` and `dV(k)`, at `mt_thermal_h.cpp:49` (Picard) and
`:124` (Newton). The plan is therefore short, but it logs three ways the sum can be quietly wrong,
because a plausible loss number is worse than none:

1. **Which resistivity.** The thermal path is per-integration-point and exists only in coupled runs;
   the magnetic side keeps a per-element `element_rho` field (`mt_maxwell_h.hpp:27-37`) that works
   isothermally — which is how AC loss is usually measured. This decision sets the task's scope.
2. **The clamp.** `compute_rho` clamps into `[gRhoMin, gRhoMax]` = `[0, 1e10]`
   (`cl_FEM_Calculator.hpp:2755-2756`). Right for solver stability, wrong for a reported engineering
   quantity. Cheap mitigation: `mRhoClamped` (`:182`) already records where it bit, so the sum can
   report its own contamination instead of hiding it.
3. **Thin shells.** For stacked-tape and CORC decks the loss lives in the shell, and
   `src/fem/thermal/matrices/` has no shell kernel — so the headline number for exactly the
   geometries BELFEM targets could come out zero.

Both plans name a validation gate rather than ending at "it produced a number": a published strip
solution for the loss, a multi-rank test plus a benchmark for the ordering.


---

## Addendum 2026-09-01 — the AC-loss plan, and a method lesson worth more than it

The Q1 plan written above was overtaken within hours. `7323d528` on Christian's Mac
("saving heatlosses, currents and voltages to exodus", two commits ahead of this checkout and on no
branch here) already lands the summation and the MPI reduction. The field was then renamed
`heatloss` → `dotQ`, for exactly the reason the review surfaced: the old name did not say whether
the number was a power or an energy. It is a power. The plan was rewritten to its real remaining
scope — time integration, per-block decomposition, a strip-solution gate, the owed two-rank run —
rather than executed as written.

**The method lesson is the durable part, and it is not about resistivity.**

A peer session and I spent two rounds arguing whether the resistivity clamp could fire. I claimed it
bit the quench regime; wrong, a quenched HTS is ~1e-6…1e-2 Ω·m against a 1e10 cap. I retracted and
claimed a deep over-critical power law would exceed the cap, since `rho ∝ (J/Jc)^(n-1)`; also wrong.
The peer had earlier supplied the parallel combination as background without noticing it settled the
question, so both of us were wrong in the same way at the same time.

Both arguments were about `rhoPL`. The functions **return** `1/((1/rhon) + (1/rhoPL))` — a parallel
combination, which discards `rhoPL` precisely when it grows. `rho ≤ rhon` unconditionally, verified
here across all eight `rho_powerlaw` overloads and the branch-stable `rho_riva`. We had both
reasoned about what the power law *computes* and neither had read what the function *returns*, and
the return statement was four lines of arithmetic in a file neither had opened.

The author had already written the answer down, twice — in prose (*"past a cap the parallel
combination is ρn to machine precision"*) and in code, at the cap itself (*"past this cap 1/ρPL
vanishes to machine precision against any physical ρn"*, guarding a deliberately NaN-aware
`if ( ! ( lg <= 250.0 ) )` so that a corrupt jc/n table degrades to normal-state instead of
poisoning the assembly). **One grep for an author comment would have beaten both sides' reasoning.**

Recorded in the plan as two struck FALSE POSITIVEs with their reasons, so the dead end is marked
rather than rediscovered a third time, and saved as a memory. Both citations are now anchored by
greppable sentence rather than `file:line`, per the anchor convention.

Net effect: the clamp flag drops from a correctness fix to low-priority insurance against a plugin
material with `rhon > 1e10` — the one route left, and the one where nobody would otherwise find out.
