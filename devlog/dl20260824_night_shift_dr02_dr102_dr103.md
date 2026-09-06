# Night shift: DR-02 resolved to deck physics; DR-102 found+fixed; DR-103 found

**Date:** 2026-08-24
**Purpose:** Overnight autonomous session on Christian's mandate ("you have codex and
grok, so let's see what you can figure out in option 3 over the night"). Also records
the post-midnight register closures and the DR-02 discriminating experiments that
preceded it.

## Register closures rolled across midnight (each on Christian's explicit ruling)

- **DR-11** closed won't-do (both halves symptom-free in production).
- **DR-13** FIXED — `Vertex::allocate_facet_container` uint→uint16_t clamp, full
  plan+audit round (both vendors: land as written), check-fast 9/9.
- **DR-16** closed as by-design (2D same-type propagation deliberately 3D-gated).
- **DR-17 / DR-42 / DR-49** signed off — fresh check-fast plus name-filtered runs:
  28/28 gesvd (4 LAPACK flavors), 3/3 spline semantic + 2/2 debug-assert (debug
  tree), 4/4 Penta6Ts.
- **DR-18** closed — tapestack3d's live campaign runs `edge coating : on` with 16
  `coating_*_copper` blocks in the saved mesh (real constructed side-connector
  geometry under coupled physics).
- **DR-19** first gate discharged by the same evidence; R7's five checks + P8/B7/D2
  stay open (precisely recorded, not overclaimed).
- **DR-22** closed — the current tapestack3d `.bfm` carries the coating blocks from
  its first save, closing the side-connector residual the 2026-08-15 note held out.
- **DR-23** T9 discharged — Christian's live 4-rank periodic corc (`out.txt`): cuts
  built clean, solve converging. Only the layer-tier forced-repair fixture remains.
- **DR-26** FIXED — dead Step-6c machinery deleted (`mPairVerdict`,
  `classify_periodic_pairs`, `CutPairVerdict`/`CutPairInfo`, plus Grok's addition
  `cut_pattern()`); comments corrected to cite the closed plan; check-fast 9/9.
- **DR-27** closed by scope (λ-guard half refuted; memory-watch half unactionable).
- **DR-06** run executed: purpose-built coupled corc scratch deck, 95 timesteps to
  the full 2 s target, T pinned 77.000000–77.000004 K, zero warnings — clean at
  iteration level and modest heating; closure is Christian's call (marathon-scale
  and quench-level heating not covered).

## DR-02: the tape_quench wall is deck physics, not a kernel defect

Three discriminating experiments (Christian's pick): E1 conservative-Picard =
floor-stuck at step 1 (settings not the cause); E2 defect-off = advanced to
t≈1.7 ms (the wall IS the defect); E3 database jc/n with defect = floor-stuck
(plugin fits not the cause). Root cause read from `lib/defect.cpp`: a hard
10⁴× jc suppression step replaced an earlier smooth Gaussian (still present,
commented out) — with n≈37 an effectively discontinuous resistivity wall, matching
the observed bit-identical Picard limit cycle. On Christian's ruling the smooth
Gaussian was restored: the deck then passed steps 1–5 at full Δt (never before past
step 1) and now stalls at step 6, t≈0.77 ms — exactly the measured current table's
first jump. A smooth+conservative-Picard combination was WORSE (step-1 stuck), so
the CustomMat Newton settings are correct and restored. **Remaining for Christian:
the step-6 wall at the current-table jump** (smooth the waveform, or accept the
envelope). Not copied to `examples/` — Christian's converge-first condition unmet.

## DR-102 (new, FIXED): order-2 gmsh meshes could never be loaded

DR-31's gate attempt exposed that no order-2 gmsh mesh loads at all —
`Face::compute_orientation` aborts during mesh setup. Investigation
(`tmp/ai_exchange/dr102_face_key_order2.md`, executable probe in
`scratchpad/dr31_order2/probe_tet10.cpp`):

- Ruled out: the gmsh→exo TET10 8/9 swap (probe: zero raw-gmsh elements in
  memory), the facet tables (consistent), the failing element (perfectly ordered).
- Root cause: `FaceFactory::face_key_3d` keyed faces by the 3 lowest-indexed of the
  FULL facet node set. At order 2, midside indices interleave with corner indices
  (27% of facet instances on the probe mesh), and two adjacent facets of one
  element sharing an edge (2 corners + its midside) can collide on the same key —
  the factory pairs non-coincident faces and the orientation corner-match rightly
  aborts (8,862 aborts in replication; Grok: silent wrong pairings also existed).
- Fix, pre-verified at data level by the probe (corner-only keys: 358,954 pairs,
  0 failures, 0 triple-owner keys, geometrically consistent): key from
  `get_corner_nodes_of_facet` — bit-identical at order 1 by construction. Landed
  with both audits clearing it, plus `face_key_2d` `set_size` hygiene (Codex's
  scratch-buffer catch, Grok concurring as latent). check-fast 9/9.
- Green gate: the orientation abort is gone; the mesh loads through edges, faces,
  and the cohomology solve — into DR-103.

## DR-103 (new, OPEN): cut-pipeline TET10 support has never executed

Directly behind DR-102: `CutProcessor::flip_node_bitsets_tet10` feeds midside
nodes carrying the `collect_nodes` `gNoIndex` sentinel into
`node_bitset()->set()` — debug assert; release would be an out-of-range write.
Mechanism and open design question (flag midsides? skip sentinels? TET10 cut
support half-built by design?) recorded in the register row. Deliberately not
fixed overnight — homology-layer semantics, Christian's call. No in-tree 3D deck
is simply-connected, so no cut-free order-2 bypass exists for DR-31's gate.

## DR-74 fixed (morning session, Christian's ruling)

Christian ruled segregated mode should write a memdump like the coupled path. The
audit round earned its keep again: both vendors independently found that the
one-liner alone would have created a silent-wrong-answer path — `load_memdump`
never restored the thermal clock (`mTime2`), so a segregated warm restart would
have re-integrated thermal from t=0, rotating the restored BDF history into
garbage. Landed both edits (`save_memdump` in the segregated save block +
`mTime2 = mTime` in the loader). Gates: check-fast 9/9, RED confirmed pre-fix,
GREEN 1 (216 MB dump written by a segregated run), GREEN 2 (warm restart banner,
full BDF5 order retained, DR-92 re-entry cap applied, exactly one thermal block
per restarted step — no catch-up avalanche). Exchange:
`tmp/ai_exchange/dr74_segregated_memdump.md`.

## Day session: DR-94, DR-100 probe, DR-103 — and the order-2 arc concludes

On Christian's "let's do DR-94 and then the other two" + "make it so":

- **DR-94 FIXED** — both input-contract artifacts corrected (the `restart` row now
  states the DR-92 cap and the `adapt timestep : false` permanence; the "raise the
  knob" advice warns instead of pointing at the sick configuration), Codex language
  sweep applied. Design call resolved as warning-not-ceiling (doc-tier row; an
  enforced ceiling would be a fresh feature row).
- **DR-100 mechanism class PINNED by the probe** — the per-request status decode
  (landed in `collect(Vector<T>&)`'s Waitall, marked for removal) shows rank 0's
  sizes collect failing with `MPI_ERR_TRUNCATE` from ranks 2 and 3 on exactly the
  size-exchange tags: a larger, stale message from an earlier exchange (prime
  suspect: the exodus save's field collects moments before) is matched into the
  1-element receives. Tag-collision class, real framework defect, precise target
  for the fix round.
- **DR-103 FIXED** (full round, both auditors: land, high confidence) — the
  one-guard fix: `check_midside` gains the `is_flagged` guard its sibling
  `check_edge` always had; skipping unflagged midsides is physically correct
  (Grok: the abort shape is the diagonal cut cases, and the guard implements the
  documented push-to-face rule; `relink_element` is a third consumer the placement
  covers). 2D/TRI6's identical latent hole fixed by the same guard (both auditors
  confirmed order-2 edge flagging was already complete). `ctest -L fast` 9/9 (one
  transient Error 8 under X2 contention, clean on rerun); order-2 inductor green
  through cuts/cohomology/thin-cut creation.
- **The order-2 excavation ends at design intent**: behind DR-103 sits
  `create_hanging_edges_and_facets`'s deliberate, always-active
  `BELFEM_ERROR( max_element_order() == 1, "Not implemented for higher order" )`
  (`cl_MaxwellFactory.cpp:1399`). Every in-tree 3D h-φ deck computes cohomologies,
  so order-2 h-φ terminates at an honest feature boundary — settling DR-31's scope
  ask (order-2 is not a 1.0 feature; its committed synch fix stays
  dormant-but-correct) and closing the DR-102→DR-103 chain exactly where it
  should: at a wall someone built on purpose.
- **DR-78's np≥2 gate PASSED** (autopilot): 2-rank warm restart from a 235 MB
  history dump resumed at full BDF5 with the banner, clean steps, zero errors.
  The standing safe-word is retired.
- **DR-89 X2**: the specified gate (50/100 ms dump) is unrunnable — that dump no
  longer exists, the same evidence decay that closed DR-76. The analog on the
  surviving t=7.1 s dump: restart mechanics PASS (full BDF5 resume, ~25 steps at
  3 iterates, Δt recovering — no startup lottery); at step 3378 the residual
  floors at ~−80 dB with Δt collapsing — NOT the DR-89 signature, consistent with
  the open DR-88 conditioning-floor class at deep quench, recorded there.
  Closure is Christian's call (accept the analog, or ride the next cold start).
- **DR-69 gate 1**: attempted and found infeasible by warm restart — swapping to
  constant jc/n under a state evolved with database jc/n is physically
  inconsistent; steps rejected, no frame in 90 min. Valid gate needs a cold
  control start (days) or a smaller purpose-built deck; gate remains open, the
  design lesson is in the row.
- **DR-06 CLOSED** (Christian's ruling on his live run's evidence): his coupled
  `cmake-build-debug/corc` campaign — 93 frames to t = 9.3 s at ~149 A with real
  Joule heating (max T 79.84 K) — shows minimum T across ALL frames of
  76.999965 K: a transient 35 µK undershoot, recovered to exactly 77.000000,
  five orders of magnitude below the −1.9 K symptom and quantitatively
  consistent with the fixed drift mechanism's µK-class residual at rtol 1e-10.
  Combined with the scratch run's iteration-level evidence, both suspects are
  cleared in their own regimes. The T-floor guard proposal lapses with the row.

## DR-19 residuals advanced (Christian: "What would you like to do for DR-19?")

The feasible slice, executed:

- **P8 doc sync DONE**: `side_connector_wall_element.md` header, §5, and §6 rewritten
  against the landed state (§5's recovery is implemented since the wall kernel;
  §6's stale 2026-07-31 WIP-stop status replaced with the real end-to-end
  implemented list). Language sweep applied. The rewrite surfaced a precisely
  stated gap: **wall Joule heating never reaches the thermal problem** (connector
  domain types absent from `src/fem/thermal/`, grep-verified) — now a named
  pending item instead of a buried clause.
- **B7 finding**: the reader-side version gate is half-built already — the loader
  prints writer provenance and warns on newer-file datasets-ignored
  (`cl_Mesh_BfmFile.cpp:291-323`). Only the strict `format_version` hard-error
  variant remains; strict-vs-warn is Christian's design call.
- **R7(a) presence half CONFIRMED**: a temporary probe on the sidecoating deck
  printed "2 coating walls dispatched ( 1 LeftCoating, 1 RightCoating )" — both
  signs constructed and dispatched. Probe stripped after the answer. The numeric
  handedness half (the O1-resolving check) stays open with the rest of R7.

Left deliberately for Christian + Prof. Sirous: R7(b) analytic r′ (needs the wall-term
theory first), R7(c/d/e) (dedicated decks/tests), D2, the coated variant.

**DR-19 then STRUCK on Christian's ruling** — the titular wiring is complete, verified,
doc-synced, and production-exercised. Surviving items stay tracked in the plan file's
live checkboxes; the wall-thermal-Joule gap is promoted to **DR-104** (missing physics:
the wall's r′ dissipation never reaches the thermal problem — grows exactly in the
scenarios the connectors exist to model; design pass needed before code).

## DR-23 closed: the forced-repair fixture

On Christian's relay of Codex's recommendation (one focused test, not another
production run), the last DR-23 gap landed as `Cohomology.AnnulusForcedRepair`
in `tests/homology/test_Cohomology.cpp`: perturb the cleaned annulus H¹ generator
by an exact node coboundary (same cohomology class), assert a genuine ±2 exists,
`clean()`, assert the unit-representative + winding-pairing contract. Full round:
both auditors cleared the pre-registered fixture with converging refinements
(drop the ring exclusion — the annulus has only ring nodes; closed-form sign
s = −g₀; index comparison; `original_edges` domain filter), all adopted. Proven
red AND green: the production log prints "rectified 1 non-unit coefficients"
where every prior fixture printed 0, and with the repair skipped the test fails
on "edge index 1 carries coefficient 2". Homology suite 9/9 (was 8), check-fast
9/9. DR-23 struck; DR-19 struck earlier the same session on Christian's ruling
(wiring complete; residuals live in the plan file; thermal-Joule gap promoted
to DR-104).

## DR-02 closed by ruling; the tape_quench four-run comparison; DR-105/DR-106

Christian ran the tape_quench deck four ways (`build/tape_quench/out*.txt`):
plain BDF1/mumps (step 22, t=1.38 ms, genuine floor stall, zero solver errors),
BDF1+gauge/mumps (14× MUMPS −9), BDF5 ungauged/mumps (41× −9 — **exonerating the
gauge**), and BDF1+gauge/strumpack (zero solver failures — the clean view of the
nonlinear map: 5–8 dB oscillation at Δt=0.1 ms, 2-iterate −75 dB convergence at
Δt≈1–2 µs). Findings: (1) the −9s are workspace exhaustion misread by the
controller as timestep failures — the code's own comment prescribes the
retry-on-(−9) fix, now **DR-106**; (2) the gauge neither causes the failures nor
enlarges the contraction radius — consistent with the closed campaign's
vacuous-on-simplex/TS finding for this mesh family; (3) the deck's µs-scale
contraction radius is invariant across every configuration — deck physics,
now **DR-105** with the four-run table as evidence.

**DR-02 then CLOSED BY RULING** (Christian, on Codex's recommendation: production
evidence over the no-longer-meaningful historical §4.1 matrix; Gate A stays dead).
The row was **compacted at strike time** — its four generations of amendments had
grown beyond human readability (Christian had to run Codex over it to process it);
the closure row carries the six-point evidence list, and the full history lives in
the register's git history and the session devlogs. Lesson adopted: long-runner
rows get dated compaction before they become machine-only.

## tape_quench tuned (DR-105): six-variant matrix, recipe applied

On Christian's "let's see if you can tweak the input.conf": a six-variant matrix
at equal wall budgets overturned the "µs-scale contraction radius" reading — it
was a cold-start artifact. Winner **T5 (start 2 µs, cap 20 µs): t = 0.954 ms in
12 min, riding the cap with near-zero rejections, straight through the 0.77 ms
wall**; T6's free 100 µs cap located the natural ceiling at ≈16–30 µs (capping
just under beats PID hunting, 4.8 vs 3.2 ms/hour). Large cold first steps are
pure waste (T2/T3: 100–160× worse); the ramp-BC control was inconclusive (it is
a stronger drive, not a gentler one). Recipe applied to
`build/tape_quench/input.conf`: initial 0.002 ms, maximum 0.02 ms, rest of
Christian's settings unchanged. Projection ~60 h serial for the 0.3 s target;
remaining levers are physics (defect/current-table softening). Key mechanism:
the radius is state-dependent — the deck that rejects a cold 50 µs first step
grows happily into 20 µs steps once BDF5 history exists.

## Standing state for the morning

- Uncommitted in the working tree: DR-13's clamp, DR-26's deletion, DR-101's
  PARDISO restore (from the previous session), DR-102's face-key fix, the
  tape_quench plugin/deck repairs.
- Open decisions queued: DR-06 closure, DR-31/DR-102/DR-103 order-2 scope call,
  DR-103's fix direction, tape_quench's current-table wall, DR-23's repair
  fixture, PARDISO-on for DR-83's last site (superseded — done), DR-100's
  comm_check probe.
- Session records: exchange threads `dr13_*`, `dr26_*`, `dr101_*`, `dr102_*`;
  probe + all order-2 artifacts in `scratchpad/dr31_order2/`.
