# Campaign: Release 1.0 (September lens)

**Status seed:** 2026-08-05 — every claim `[seeded — confirm]` until Christian's
correction pass. This page aggregates; the row-level truth is
`todo/debt_register.md` (blocking-1.0 column).
**Branch:** `main` + `sideconnectors` (merge unit TBD)

## What 1.0 means `[seeded — confirm]`

Open-source release of the framework (literature/ and nonfree/ stay out; gasmodels
migrates in via its own campaign). METHODOLOGY.md + cross-review + falsification
tooling are part of the public story.

## Blocking rows — refreshed 2026-08-10 against the register

**Still open (14).** Most are *run* gates, not code:

| kind | rows |
|---|---|
| needs a run | DR-02 (kernel-collapse §4.1 — Gate B done on helix; thermal `T_h` coverage owed, and its decks need repair first) · DR-06 (−1.9 K dip) · DR-08 (−58.74 dB freeze) · DR-09 (4-proc coupled verify) · DR-15 (greg2 smoke, Gregory) · DR-18 (side-connector smoke) · DR-22 (bfm save→reload) · DR-34 (STRUMPACK S4 rank sweep) · DR-38 (circuit restart) · DR-52 (greg3 A/B — code is committed, only the run remains) |
| needs code | DR-19 (side-connector FEM wiring — far advanced) · DR-21 (bfm checksum ignores processing options) · DR-23 (SPFA rectifier has zero checked-in tests) · DR-42 (LAPACK/spline test debt, Christian) |
| needs a ruling | *(none — DR-54 was the last one, closed 2026-08-10)* |

**Closed since the seed:** ~~DR-01~~ (stale decision) · ~~DR-03/04/05~~ (premises dissolved
or fixed) · ~~DR-10~~ (Anderson landed, `6b2a2b98`) · ~~DR-12~~ (orientation trust flag) ·
~~DR-20~~ · ~~DR-24~~ · ~~DR-25~~ · ~~DR-32~~ · ~~DR-44~~ · ~~DR-47~~ (HEX8TB factory case;
`Hex8TbUnitCirculation` re-enabled 2026-08-10) · ~~DR-48~~ · ~~DR-50~~ · ~~DR-51~~ ·
~~DR-55~~ · ~~DR-57/58~~ (first-run repair, `bc578b5e`) · ~~DR-54~~ (`examples/corc` replaced
with the working 6-tape model, verified from a clean dir). **Downgraded:** DR-13 (uint8 →
uint16 facet capacity). **Non-blocking but live:** DR-53 (fusing flags stay `false`) ·
DR-59 (latent landmines in the retired parallel rho builder).

**Commit state (2026-08-10):** working tree clean; the bundle that gated several of these —
Anderson residual fix, PID controller, material first-run repair — is in history
(`4f2c11cd`, `1d6ef305`, `bc578b5e`). No open row is waiting on a commit any more.

## Dormant campaigns folded in `[seeded — confirm]`

- **periodic** (branches `periodic`, `periodic_new` both unmerged — diff before trusting
  any periodic audit): validation debt = DR-23/24/25/27/28; twisted-helix + single-layer
  corc are the regression decks.
- **double-corc** (seam campaign complete 2026-07-16, uncommitted parts since landed):
  residual = D1 kernel-collapse handoff (`todo/handoff_double_corc_d1_session.md`),
  half-cut conductor watch, T8 port (= DR-23).

## Open decisions for the cut `[seeded — confirm]`

- Merge order: sideconnectors → main vs feature-gated; periodic branches' fate.
- ~~FVM module in or out of 1.0 (DR-41 assumes out).~~ **Ruled 2026-08-14: out —
  moved to `nonfree/fvm` (unfinishable before release).**
- Umbrella naming question (dl20260615_architecture_structure_review).
- Version at tag time: CMake says 0.9.0; bump to 1.0.0 (README/banner/doxygen inherit).
- CITATION.cff co-author list to be completed from the SUST DOI record.

## Release packaging state — 2026-08-14

Rulings: **`USE_TEST=ON` is the release default** (docs + bidirectional claim
probe updated) and **FVM is out** (see above). Landed in the working tree:
rewritten root `README.md` (front door), `CITATION.cff`, published
`doc/getting_started.md` user entry point, `PROJECT_NUMBER` from the CMake
version, the doc-sweep's false-claim fixes (Maxwell XML fiction, fixed-bug
banners, wrong CMake flags/class names/geo names), and a genuine 0-warning
doxygen state (four pre-existing CWD-resolved README-link warnings fixed).
Details: dl20260814_doc_sweep_fixes_and_release_packaging. Owed gates: real
`make doc` and `make check` on the reconfigured tree (Christian).

## Dated entries

todo/debt_register.md (the row-level lens) · dl20260615_architecture_structure_review ·
dl20260805_cross_review_tooling · dl20260809_todo_currentness_sweep ·
dl20260810_examples_and_rho_database_repair (first-run repair; DR-54 is the residue)
