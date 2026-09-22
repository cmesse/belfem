# Test-gap Tier 1: the unattended overnight run

**Date:** 2026-09-22
**Purpose:** Session record of working off Tier 1 of `tmp/devlog_sweep/triage.md` v2 with the
`code-with-jury` loop, unattended, under the standing approvals recorded in
`todo/test_gap_tier1_execution.md`
**Module:** `tests/` (mesh, fem, linalg, io, physics); `src/mesh`, `src/fem/kernel`,
`src/physics/materials` for jury-confirmed local fixes
**AIs involved:** Claude Fable 5.1 (loop, implementation agents, plan-drafting agents), Codex
(`gpt-6-astra`/high on critical rows; `gpt-5.6-sol`/medium plan and `gpt-5.6-terra`/high code
juries otherwise), Grok `grok-4.7`/high

## Summary

Tier 1 of the test-gap triage (22 rows) ran unattended overnight through the standing loop:
plan → plan jury (Codex + Grok) → Opus implementation → gates in both trees → code jury → one
commit per row. **21 rows landed on `devel`**, T-08 is carried with its plan v3 rulings. Nine of
the rows fixed a live source defect on the way (T-01 facet node lists, T-02 thin-shell
selection, T-26 constant-n dispatch, T-32 periodicity refusal, T-38 curve refusal, T-44 cut
coefficient refusal, T-51 factory flag, DR-178 the prism side-face surface increment — half the
area in every build — and N-19 the user-plugin dependency check); the rest pinned contracts
that had no test, several with `DISABLED_` pins on defects found by the new tests (DR-173 lead
cp, DR-174 YBCO raw λ). Twenty-two register rows were filed (DR-161 to DR-182), one closed
(DR-178). Every source fix was consumer-swept after Christian's mid-run rule (8/8 LOW, five
by-catch findings). Christian answered the 26 open rulings in question form; twelve follow-up
rows go to the next session. Final gate: `make check-fast` green in both trees. Vendor
incidents: two Codex usage windows, two Grok turn exhaustions, one two-makes collision, one
jury reading a scratch mutant; a Copilot wrapper was added as a reserve seat.

## Rows

Order follows `tmp/devlog_sweep/triage.md`; each row went plan → plan jury → implement → gate →
code jury → commit, with the exchange thread at `tmp/ai_exchange/review_tgap_<id>.md`.

### T-01 — `get_bottom_nodes` QUAD9 midside and HEX20 face size (`c41c1c0`)
Live defect: the HEX20 case sized its list to 9 and read node 21 of a 20-node element; QUAD9
and QUAD9TS read node 5 (the right edge) instead of the bottom midside, node 4. Fix in
`src/mesh/meshtools.cpp`; `tests/mesh/test_MeshTools.cpp` builds every supported element at its
Lagrange parameter coordinates and asserts the positional bottom/top pairing. Red on the unfixed
source for exactly the three predicted cases, green after; release 61/61 → 74/74 after T-13.
Both juries approve; five P3 nits applied. Residue: `get_top_nodes` accepts three types the
bottom helper refuses (DR-161).

### T-12 — `B * inv( A )` formed as one expression (`90c7bfe`)
Test-only pin of the lazy-inverse shape defect fixed on 2026-08-03. The named mutation could
not be applied to the tree (the wrapper header is included library-wide), so it ran as a
standalone probe linked against the debug library: lazy 2×3, wrapper 3×2. Both juries flagged
the comment for stating the bug as a requirement; reworded.

### T-11 — HDF5 read contract across the file type (`7d3dfff`)
Same-class widths convert (scalar between canaries, `Vector`, `Cell`, `Matrix`; the widening
target prefilled with ones so a short write is visible), other classes are refused in every
build. The plan's bool canary was dropped as a null gate (the loader reads into a local
`hbool_t`). Residue: the loaders leak two handles per refused read (DR-165); the philosophy's
"no automatic numerical conversions" has an undocumented I/O exception (O6).

### T-02 — thin-shell facet selection in `DofMgr::BlockData` (`b5a18cf`)
Live defect: the collector counted by selected shells and filled from every shell. One public
static predicate now serves count and fill; the serial linker skips unselected layer blocks.
`tests/fem/test_ThinShellSelection.cpp` (two ghost stacks, every id space continued). Named
mutation: two of three tests red, control green. Grok's first run exhausted its 30-turn budget
twice on this row; 60 turns delivered. Residue: the serial side-connector loop stays unguarded
(DR-162, a ruling on partial selection).

### T-10 — RRR identity and `Tmatrix` threshold (`81ffe4e`)
Test-only. The plan jury caught that public `rho_i()` is the 4 K spline on seven of nine metals
while `set_RRR` zeros `rho_i_custom`; the test reads the latter. The incident's literal
spelling (`abs` for `std::abs`) resolves to the double overload under these flags and stays
green; the integer-truncation mutants (`( int ) std::abs`) go red. Residue: ferromagnets carry
`rho_mag` outside the RRR identity (O7).

### T-13 — counter widths (round 2 in progress)
Sources uint16, dofs uint8, node/edge/face/element capacity uint8, facet capacity uint16 (the
plan had generalized 256 to all five), graph vertex counter uint32. Every refusal asserts the
`BELFEM_ERROR` substring and the unchanged count; the graph test asserts 70 000 fatally before
allocating. Residue: the graph counter's guards are assert-only in release (DR-166); an
allocated-but-unfilled source container leaks in `~Basis` (DR-169).

### T-26 — material derivative contracts (`99ac9d0`)
The plan jury found three live defects: the piecewise `dJ` tangent carries the parallel factor
below `j1` while the residual is the raw power law (DR-163), the factory Kim law has no field
derivative so `drho_powerlaw_dB` is zero (DR-164), and the constant-n three-argument overloads
dereferenced a null n-function on exactly their precondition (fixed at six sites; the named
mutation segfaults). The finite-difference oracle needed three test-side changes: a Richardson
estimate with a 3e-5 agreement gate, 1e-4 for the field overloads whose analytic tangent is
itself a fixed-step difference (iron at 7.3 K, 5 T: 1.7e-5 off, O5), and sample ratios ≥ 0.9
`J_cross` so the T-dependence is not roundoff. Residue: DR-170 (`n = 4` knot: `dJ` aborts,
`dT` falls back).

### T-24 — ghost interface coefficients (`73f8aa4`)
The code is symmetric (C12 holds); symmetry and the constant-field residual were already
pinned. A one-sided constant field makes the consistency and penalty rows hand-computable, so
the harmonic form of the mean, `1/h` and `α = η · k_pen` are now pinned on both flips. The plan
jury rebuilt the mutation map (the one-sided field annihilates the primal terms) and withdrew
the coercivity question. The formulation question (jump of the outward slopes vs SWIP's
weighted average) and the three-way disagreement on the zero-resistivity floor (code `η·k_reg`,
doc `2·η·k_reg`, kernel comment `2·k_reg`) are DR-168, marked `[F]`.

### T-32 — periodicity node match (`4871e04`)
Live R32 defect: `match_nodes` ended in `BELFEM_ASSERT` then `return true`; release ran on with
null partners and died later elsewhere. Promoted to a `BELFEM_ERROR` naming the unmatched
count (the count is hits, not distinct partners, so the message says "matched k of n"). Seven tests on a hand-built HEX8 bar: bijection with positional pairing, a twisted cap,
facet-stage and node-stage refusals, sub-tolerance acceptance, the kd-tree's named candidate on
a rejection, conductor edges through the master fallback. Mutation red by message text.
Residue: leaks on every refusal path (DR-167); off-axis rotation is not a supported transform
(O2).

### T-13 — counter widths (`e2f66ab`)
Three named mutants red (facet-guard deletion, source-guard deletion, graph counter narrowed
to 16 bits); the second and third chains collided in one tree and were repaired by a clean
serial rebuild. Round 2 of the code jury asked for fatal prerequisites, a scope owner for the
heap entities, `memory()` deltas as the release oracle for the capacity stores, and the face
parent in an owner; all applied.

### T-33 — pinning a hanging dof fixes its source (`8124a70`)
The row's "throws in debug" was wrong: both guards in `pin_dirichlet_dof` are `BELFEM_ERROR`.
The real gap was DR-135's untested residue. A hand-built two-square mesh with single-source
ties covers the reroute, the near-unit weight (where the division is observable — Codex
refuted "dead arithmetic" with a binary64 probe), the non-unit refusal, the post-freeze
refusal with its control, and the two-source fallback (which locks the helper's documented
behaviour, not enforcement — DR-172).

### T-51 — materials contract through the deck path (`2076ec2`)
The deck constructor hardcoded table building, so a defaulted `aBuildTables` seam was added
(production unchanged). The five documented forms are pinned; two silent mis-selections the
plan jury found (an in-name RRR token beating the deck's `RRR`, a valued `builtin` beside a
`curve` yielding an `Iron` under the key `copper`) and two accessor gaps are DR-171, marked
`[F]` because a fix changes what a shipped deck form means.

### T-34 — face creation on two TET10s (`abd598f`)
Seven faces, one shared with the lower id as master whichever element the factory meets
first, connectivity counts, and the sideset-facet id. Both auditors refuted the backlog's
"midsides sort below the corners" collision story; the discriminating layout is two corners
plus their edge midside. A third owner is silently dropped by the factory (O9).

### T-23 — 2-D generator sense through the suggest path (`5834e8b`)
The plan agent, Grok and Codex all derived +1 for the CCW inner-ring pairing on today's tree
against the row's −1 and the closed core's own comment; the test pins the sign the binary
produces as a named calibration constant, and the first run gave +1. The disagreement with
the comment is Gregory's (O10); the core was not touched.

### T-36 — relinked master node: thin-shell facet kept, plain facet re-derived (`892f862`)
The "restore" is a domain-type skip in `update_facet_nodes`; the release heap write the row
names is unobservable from `tests/mesh` (the incident's write is now a `BELFEM_ERROR` in
`ThinShellFactory`). The plan's relink-before-first-finalize case would have thrown on a
correct mesh (orientation computed after the relink) and was dropped.

### T-38 — `sort_end_nodes` on two chains (`f58524c`)
Live defect: the degree count passes two disjoint chains and the permutation indexes the
adjacency with `gNoIndex`. A complete-traversal `BELFEM_ERROR` naming the curve id was added
after the forced writes, the three sorters made public with their preconditions stated, and
six polyline tests cover refusals, the monotone sort, reversal invariance and a single loop.

### T-25 — builtin material values on a 0.5 K grid (first gate)
Two findings on the first run: lead's public specific heat is negative below 3 K (DR-173) and
YBCO's raw Callaway conductivity is NaN below about 1 K, hidden by the public spline (DR-174).
The envelope threshold pair is a named policy constant (O8). Both defects are pinned as
disabled cases through a known-defect table.

### T-44 — `collect_coefficients` ±2, k-hole unimodular pairing (implemented, gate pending)
The row's first contract was not in the tree: a ±2 silently became weight 0 in the open file
`cl_CutData.cpp`. An else-branch `BELFEM_ERROR` naming the cut and the edge was added there;
ten tests cover the refusals, the bitset pairing, the cleaned support and six k-hole cases.

### T-60 — facet normals, surface increment, measures (gated; one live defect, DR-178)
No defect on the normals; the plan round rebuilt the mutation map, and Codex derived two new
suspects: 2-D QUAD4TS ghost master normals point inward on the production clockwise winding,
and the PENTA side-face `dS` appears short by a factor of two. The probe settled it: in both
builds the side face integrates to exactly half its area (ratio 0.5000000000000001), the top
face and the ghost normals to machine precision. Mechanism in `normal_penta` cases 0-2: the
barycentric in-plane rows of the master Jacobian carry the full edge, the axial row half the
thickness, so the cross product is half the area and the quadrilateral factor 0.25 is half of
what Σ w = 4 needs. Filed as DR-178 with the one-token fix (0.25 → 0.5 in three cases) in a
plan under its own jury; whether any production sideset reaches a side face is Christian's
ruling. The probe is a disabled pin until then; the twelve other new cases are green in both
trees. Christian's note that PENTA6TS and HEX8TS facet numbering differ from the plain types
went into the brief and into memory; the agent found the two facet tables identical today and
`create_unity_nodes` without a HEX8TS arm (DR-176).

### N-17 — integer fidelity, bitset assignment, node coordinates (`c2490a4`)
All three claims were already contracts of the tree (the integer fill was typed since
`8465284`); the pins are regression guards, and two named mutations show they bite. Grok's
code round read a scratch mutant mid-chain and reported the tree broken — refuted, incident P5:
no mutation chain while a jury reads the tree.

### N-18 — dense algebra siblings (`a05c080`)
GEMM mixed transposes with odd stored row counts, `eigen_sym` on the upper triangle with the
lower one scribbled and the input compared afterwards, the Bezier third derivative on the
production-shaped non-affine curve, and the tensor Database's per-axis derivative steps. Both
juries rejected the plan's mutation map and the implementation agent rewrote it; four of the
corrected mutations ran and each reddened exactly its predicted tests. The agent also checked
every HEX27, QUAD9, HEX8 and QUAD4 table entry off-line against the Kronecker property.

### N-19 — constitutive dispatch (`f2d34c2`)
One live defect: the two-argument user `jc`/`n` registration ignored its declared dependencies,
so a plugin declaring `( T, normB )` was silently called as `( normB, angleNxB )`. Christian
ruled the refusal in and the shipped template corrected; the constructor now closes the
`dlopen` handle when a plugin's initializer throws. A MODULE plugin fixture pins the acceptance
by sentinel scales and the four refusals; the YBCO conductivity origin slope, the Callaway
optical channel against a closed form, and the convenience spline's end slopes are pinned too.

### T-08 — deck parser (carried)
The row grew twice under review. Round 1 rejected v1's hard errors on citations and format
strings; Christian then ruled the format deliberately failure tolerant (warn and proceed,
example decks flawless, a Python checker errors), which became v2; round 2 at xhigh found two
load-bearing holes in v2 — the `;`-statement-before-`{` shape is the factory's live third
builtin form, and an EOF close re-parents every later section — and Christian ruled both and
chose to carry the row. Plan v3 is a list, not a plan, at the end of the plan file.

## The quiz

Christian answered every open ruling in question form ("the quiz"): O1–O21 and the register rows
DR-163, 168, 170, 171, 172, plus the two T-08 rulings. The decisions are recorded inline in the
tracker's §4 and in the register rows; the twelve follow-up rows they produced are in
`todo/test_gap_followups.md`. Three answers changed the night's course: the input format is
failure tolerant by design (T-08 became a warning design), a blanket refusal of multi-source
Dirichlet dofs would break the cohomology cut duplicates (DR-172 narrowed), and `n = 4` was
never illegal (DR-170 is only about the abort tier).

## Tooling during the run

Codex's usage window closed twice (01:45 and 02:42). Both times the missing voice was reported
in the thread rather than substituted; both times Christian topped up and the Codex halves were
re-run from a script that rebuilds the jury prompt per row. Grok's turn budget went from 30 to
60 and, for one fixture-heavy diff, to 120 with a brief naming the three functions to read.

Christian's side quest from the same message: `.claude/scripts/ask_copilot.sh`, a third seat
mirroring the Codex wrapper (model and effort knobs, exchange slug, a `# COPILOT` entry). The
GitHub Copilot CLI runs read-only through its own permission rules: every writing tool denied,
a shell allow list of read-only commands, the GitHub MCP server off, the checkout never
trusted. The smoke test read a file and had an append refused. On this account only `auto`
(which routes to gpt-5.6-luna and ignores the effort knob) is policy-enabled; claude-fable-5.1,
claude-opus-5, gpt-6-astra, gpt-5.6-terra and the rest are in the catalog but policy-disabled,
refused before billing, so Copilot is a luna-class voice until the github.com toggles are
flipped.

## Consumer-impact sweep (Christian's rule, applied retroactively)

At 03:15 Christian asked that a fix which can change a value or reject an input be preceded by
an inventory of its consumers and of the shipped decks that reach them, and that the rule be run
retroactively over every source fix of the night. Eight inventories (`tmp/ai_exchange/
consumer_impact_<id>.md`, one read-only agent each, T-51 by hand) came back LOW without exception,
each for a structural reason rather than luck: T-01's two branches sit behind the Maxwell
factory's refusal of any non-TRI3/LINE2 shell; T-02 is a strict narrowing on the one funnel every
field passes, and both equations select by domain type so no deck can partially select a shell;
T-26's changed overloads have no production caller and a deck reaching them would have crashed,
not miscomputed; T-32's old release build segfaulted on its null targets and all eight periodic
decks are node-for-node images inside the 1e-9 gate; T-38's cut factory splits into components
before sorting; T-44's refusal stands behind two release-active unit gates inside `clean()`, and
the old 2-D release path ran a wrong-sign cut silently; DR-178's side faces are reachable only
through a prism volume block, which the calculator refuses at initialization. The sweep also
found five things nobody was looking for: DR-180 (a null sideset dereference on second-order
2-D curves), DR-181 (`background_phi` reads the surface increment at point 0 for every point),
DR-182 (a periodic stamp silently removes a kernel-bearing sideset), a dead bitset pair in the
cut processor, and a `Material::n()` without the constants fallback its sibling `jc()` has. The
reach questions that reading cannot settle — Christian's own decks with a hand-built prism tape
on a symmetry plane, the second-order thin-shell intent, the T-less power-law family — are in
the tracker's ruling list.

## Defects found by the run

(D1–Dn from the tracker, with the register rows they became)

## Process notes

- The two-hour clock skew in the first exchange stamps (~01:12 for a 23:12 entry) came from
  writing times by hand; every later stamp is `date`.
- gtest parameter names: `to_string( ElementType )` folds thin-shell variants onto the base
  name, so a `TEST_P` over element types needs the enum value in its name generator.
- The G4 named mutation of a header included library-wide (`fn_BZ_inv.hpp`) was run as a
  standalone probe linked against the debug library rather than as a tree edit.
