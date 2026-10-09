# Devlog 2026-10-09 — No automatic gauge pin in a φ component that a symmetry plane already fixes

**Date:** 2026-10-09
**Topic:** Gregory Giard's report, relayed by the Server session on Christian's word ("this needs fixing before 0.9.3"): with a symmetry plane in the deck and no bearing, the automatic gauge pin of 0.9.0 was still imposed at an interior node of the gauged component. Mathematics round, plan and code juries, the fix, a red-on-HEAD test, the shipped-deck gate
**AIs involved:** Claude Fable 5.1 (literature, plan, implementation, gates); Codex gpt-5.6-terra (mathematics, light pass), gpt-6-astra (plan and code audits); Grok grok-4.7 (all three rounds); Server session (brief, rulings relayed); DrClean session (00b7972 landed meanwhile, rebased onto)
**Claude Confidence:** high on the mechanism and the gate numbers (executed); medium on the attribution of the largest far-field changes (observed, not traced)
**Literature References:** Brenner & Scott 2008 §5.2 (pure Neumann: unique up to a constant, one scalar gauge, Eq. 5.2.3); Monk 2003 Thm. 3.37 (curl-free ⇒ gradient of a potential unique up to a constant, simply connected); Bathe 2016 §4.2.2 Eqs. 4.42–4.45 (a prescribed dof's row leaves the system; its residual is the reaction); Arsenault et al. 2021, Appendix (φ gauged by the Dirichlet data of the applied field); Dular et al. 1999 §VI Fig. 3 (the two symmetry surface kinds, h×n = 0 and b·n = 0)

## Summary

Christian challenged the premise ("I challenge the mathematical reasoning; where would the
offset come from?"). The mathematics round settled it, and the deck gate then showed the
effect. In a non-conducting φ component the problem is pure Neumann and φ is unique up to one
constant; an h×n = 0 symmetry plane imposes φ = 0 on every node of the plane and so uses up that
constant. The automatic pin, by construction the interior node farthest from any sideset,
then adds φ(x_p) = 0 as a second condition. Its "offset" is φ*(x_p) = −∫ h·dl from the plane to
the pin, the magnetomotive force between them, generically nonzero; forcing it to zero removes
the pin's row from the solved system and its residual becomes a reaction, a point source at
the pin. Both auditors confirmed; neither found a generic harmless case beyond φ*(x_p) = 0 by
coincidence. On costheta the gauged potential at the pin is −757 A and the reference run shows
|h| = 1.4e5 A/m in the pin's element patch where the fixed run has 100 A/m.

The fix: `MaxwellFactory::impose_autopins()` decides per pin, when the pins are imposed and on
the full mesh of rank 0, whether the pin's φ component holds a node whose φ is already fixed —
a non-duplicate node of an Air/Buffer/Ferro Symmetry sideset or of a `dirichlet` sideset — and
leaves such a component to its plane. Selection (`find_autopins`) and the `.bfm` cache are
untouched, so stored pins stay valid candidates and a `dirichlet` condition added under a
cached mesh is honored. Antisymmetry planes fix nothing and keep their pin; a deck `bearing`
still disables every automatic pin. One Default-level line per component names the pin or the
contact.

## Key Findings

- **Where the rule has to live** (plan jury, Codex astra/high + Grok 4.7/high, both "request
  changes"): the symmetry switch at the factory iterates `dofmgr()->sidesets()`, which exist
  only for sidesets with a locally owned facet, and the dof-level fixed state on rank 0 is
  partial until `synchronize_dirichlet_bcs()` runs inside the first `initialize()`. A plane
  owned entirely by workers is invisible to both. The contacts are therefore read from
  `mMesh->sidesets()` by domain type on rank 0, whose mesh is complete. Confirmed by the
  two-rank gate: gantry and costheta on `mpirun -np 2` print the skip line and reproduce the
  serial fixed run to 1.7e-9 (costheta terminals) and 4.5e-9 median (gantry's 913 voltages).
- **A duplicate's fix is inert** (`SideSet::impose_dirichlet` keeps a duplicate's alias fix on
  purpose; the Dirichlet loop skips duplicates), so only non-duplicate sideset nodes seed the
  contact set; a plane touching a component only through duplicates does not suppress its pin.
- **What a `gauge` is not** (code jury, Codex astra/xhigh + Grok 4.7/xhigh, two P1s): my first
  diff treated `gauge` nodes as contacts. A gauge fixes its dof only in `impose_bc`, after
  `shift_fields()` has frozen the free/fixed split, and names vertex ids, not node ids. The
  branch is removed; a gauge keeps HEAD's behavior. The late flip itself is a pre-existing
  hazard of the `gauge` type (no shipped deck uses it) — for Christian, not fixed here.
- **Gregory's workaround** (a bearing on a node of the plane) works because any deck bearing
  disables the whole autopin loop, not because a bearing on the plane repeats the plane's
  condition — my first message said the latter; both auditors corrected it.
- **Error scaling** is capacity scaling, not the universal O(h) I first wrote: in 3-D the
  monopole's far field is O(φ* h/r); in 2-D it decays only logarithmically, so on the 2-D decks
  (costheta, gantry) the change is wide, not local — as the gate shows.
- **Walk adjacency** mirrors `find_autopins` exactly: the flagged elements of the node, of its
  duplicates, of its periodic partner, and of the partner's duplicates, all through
  `original()`; visited nodes on flag plane 0, contacts on plane 1, φ blocks on element flag 0;
  all cleared before and after the loop, plane 0 per search through the queue.

## Changes Made

- `src/fem/maxwell/cl_MaxwellFactory.cpp/.hpp`: the autopin block of `create_magnetic_kernel()`
  becomes `impose_autopins()` (the hanging/regular fix logic moved unchanged; the pre-feature
  warning stays at the call site); new `find_fixed_phi_contact()` (plane-0 BFS, queue = visited
  `Cell`) and `push_phi_neighbors()`. `BELFEM_ERROR` on the rank precondition (setup code) and on
  a `dirichlet` sideset id missing from the mesh.
- `tests/fem/test_AutoPinGauge.cpp` (registered under `USE_MAXWELL`): a 2-D TRI3 strip with two
  air components separated by a copper column, sideset 24 at x = 0, a point element at node 11.
  Five cases: both pinned without a plane; symmetry plane → left pin free, right fixed, skip line
  with sideset 24; `dirichlet` likewise; the `.bfm` of the plain deck reused under the
  `dirichlet` deck ("Reusing processed mesh" asserted) → left pin free, and both pinned again
  under the plain deck; a deck bearing → no autopin fixed, node 11 fixed. The symmetry,
  dirichlet and cache cases were RED on the unfixed factory (built with the change stashed),
  all five green with it.
- `doc/input_file_reference.md` (bearing paragraph), `doc/input_schema.yaml` (`bearing_note`),
  `examples/README.md` §"Do not write a bearing", `CHANGELOG.md` [Unreleased] `### Fixed`
  naming gantry, costheta and racetrack_usermat as the decks whose results change.
- `todo/closed/autopin_symmetry_gauge.md` (plan, D1–D5 from the plan jury folded in),
  `todo/README.md`.

## The shipped-deck gate

Serial release build (MKL flavor, `-O2`), one regenerated mesh per deck shared by all runs, no
`.bfm`, simulation times cut to 1 s (gantry, costheta), 20 s (undulator2d) and 100 s
(racetrack_usermat, capped at 25 min wall time); the reference binary is devel 53b67ec, the
fixed one this branch. Pins read from the `.bfm` (identical by construction), fields from the
last Exodus save at the pin, in its patch (2–10 nearest-neighbor distances) and mesh-wide.

| deck | pin (contact) | terminal history ref → fix | at the pin, ref → fix | mesh-wide |
|---|---|---|---|---|
| costheta (2-D, air symmetry 330–332) | 5843 at (0.113, 0.093); left to node 3269 of sideset 331 | U_01 0.03264 → 0.03236 V (−0.9 %), worst 9.0e-3 over 6 rows | φ 0 → −757 A; \|h\| 6.3e3 → 95 A/m, patch max 1.42e5 → 102 A/m | median \|Δh\| 61 A/m (3.3 % rel); largest 1.8e5 A/m on the y = 0 line in air block 72, 9 cm from the pin |
| gantry (2-D, air + ferro symmetry) | 12891 at (0.239, 0.081); one component, left to node 1551 of sideset 8 | U_001 2.907e-4 → 3.196e-4 V (+10 %); largest column U_456 −7.979e-2 → −7.811e-2 (2.1 %); median 2.2 % over 912 voltages; currents identical | φ 0 → +28.2 A; \|h\| 120 → 88 A/m, patch max 1.2e3 → 90 A/m | median rel 0.28 %; largest 1.4e3 A/m (0.5 %) at a tape |
| undulator2d (2-D, no plane — control) | 14624, pinned in both | 1.06e-5 (Picard at 1e-6) | identical to 1.9e-10 rel | max 4.7e-5 rel |
| racetrack_usermat (3-D, air symmetry `sym1` = sideset 278) | 204629; left to node 204193 of sideset 278 | none comparable: both runs capped at 25 min after 14 steps to 50 s, the reference before its first save | — | base mesh identical (245542 nodes); prescribed dofs 1023 → 1022, free 1286964 → 1286965: exactly the released pin |

Two ranks (`mpirun -np 2`, fixed binary): costheta terminals to 1.7e-9 and fields to 2.7e-12
of the serial fix; gantry median 4.5e-9 over the voltages (one small column 2.7e-3), fields to
8.8e-6. The symmetry planes' ownership is whatever the partitioner gave — the production path.

Reading: the reference runs carry a point source at the pin whose near field is 10× (gantry)
to 1000× (costheta) the physical field there; in 2-D its footprint reaches the terminals at the
percent level. The largest far-field changes on costheta sit in air on the midplane next to
the iron — attributed to the nonlinear yoke redistributing a percent-level MMF change, not
traced (medium confidence). The control deck is unchanged to solver tolerance.

Suite: `make check-fast` 8/8 in the private release tree before and after the code-jury fix
and after the rebase onto 00b7972 (the final tip f6618b9 differs from it only in deleted example files and three documents); zero compiler warnings.

## Open Questions / residuals (for Christian)

- `gauge { nodes }` fixes its dofs after the free/fixed split (pre-existing; both auditors,
  high confidence; no shipped deck uses it). Either pre-fix it in the factory like a bearing,
  or refuse it. Not a row — Christian's rule; his call whether one is opened.
- Two disjoint symmetry planes in one component are both set to 0, which already fixes the MMF
  between them (Grok) — pre-existing, worth one sentence in the usage guide.
- `find_autopins` seeds `flag_all_nodes` without walking `original()` (Grok, selection side,
  unverified, low confidence that any standard plane hits it). Selection is unchanged here.
- No periodic fixture: a contact reached through partner adjacency is mirrored from
  `find_autopins` line by line, not tested.
- The factory-level tests hardcode `library : mumps` (existing pattern).

## Files Updated

`src/fem/maxwell/cl_MaxwellFactory.cpp`, `src/fem/maxwell/cl_MaxwellFactory.hpp`,
`tests/fem/test_AutoPinGauge.cpp`, `tests/fem/CMakeLists.txt`, `doc/input_file_reference.md`,
`doc/input_schema.yaml`, `examples/README.md`, `CHANGELOG.md`,
`todo/closed/autopin_symmetry_gauge.md`, `todo/README.md`, `devlog/README.md`, this file.
Exchange: `tmp/ai_exchange/review_autopin_symmetry_gauge.md` (mathematics),
`tmp/ai_exchange/review_autopin_gauge_fix.md` (plan, code, light pass; diffs v1–v3).
