# Ferro Undulator: D1 mu Short-Circuit Fix and InterfaceFerroAir Backstop

**Date:** 2026-08-12
**Purpose:** Apply the audited fixes R1 and R6 from
`todo/ferro_undulator_mu_and_interface_bugs.md`.
**Module:** `src/fem/kernel`, `src/fem/maxwell`

## What was done

Two source edits, both pre-audited in the 2026-08-12 review round (Codex confirmed, Grok leg
failed; see the todo file §7 for the audit trail):

1. **R1 — `cl_FEM_Calculator.cpp` (MaxwellData constructor).** The unconditional
   `constant_property( mu )` call tripped its `BELFEM_ASSERT` for any material with a B-H curve,
   killing every ferro run in a debug build. Replaced with a short-circuit:

   ```cpp
   bool tIsConstantMu0 = tIsConstantMu
           && mMaterial->constant_property( MaterialProperty::mu ) == constant::mu0 ;
   ```

   `constant_property` is now called only when its assert would pass, the flag is always
   initialized, and the value matches the previous release-build behaviour on every path
   (release was accidentally correct: `NaN == mu0` is false). Gregory's local workaround —
   declaring the flag and assigning it only inside `if (tIsConstantMu)` — leaves it
   uninitialized on exactly the ferro path and is read in the bulk branch, i.e. undefined
   behaviour that can silently make the iron yoke behave as vacuum. That variant must not
   reach `BELFEM_Template` (todo step R2).

2. **R6 — `cl_IWG_Maxwell.cpp`, `link_to_group`, `case DomainType::InterfaceFerroAir`.** The
   kernel assignment was guarded by `if ( mFormulation == maxwell::Formulation::HPhi )` with no
   `else`, so any other formulation would silently inherit the previously linked group's
   `mFunMKF` — the same stale-kernel mechanism as removing the coil-interface errors (D2' in
   the todo). Added an `else` with a `BELFEM_ERROR`, matching the admission-tier error policy
   and the HPhi-only solving coverage recorded in `maxwell_usage_guide.md`.

## Verification

Both files pass `g++ -fsyntax-only` with the exact `CXX_FLAGS`/`CXX_DEFINES`/`CXX_INCLUDES`
from the build tree's `flags.make`. **Reviewed, not verified:** no full build or benchmark run —
the undulator debug-build gate (todo step R8) is still owed, and the build is run by Christian.

## Round 2 (same day): jury audit of the applied fix

Both auditors returned (the earlier Grok failures were brief-specific; a narrowed subject file
scoped to the two hunks ran clean) and **both confirm both edits — no P0, no P1**. The round's
sharpest contribution was Grok refuting the earlier framing: the ferro assembly kernels
(`phi_ferro_newton/picard`) call `Material::dmudH` directly and never consumed the two flags —
R1 unblocks the *constructor*, and independently sets the right MaxwellData helper pointer for
postproc/thermal consumers. Also corrected: round 1 misnamed the whitelist function
(`set_block_types_in_magnetic_equation`, not `select_sidesets`). Pre-registration, verification
pass, and reconciliation table: `tmp/ai_exchange/review_ferro_undulator_bugs.md`.

## R4 traced in-house (Christian's 2D hint) — and a new defect

The D2 admission machinery was then traced end to end: admission is by ID filtered on
type-at-snapshot (`set_block_types_in_magnetic_equation`), but the type `link_to_group` consumes
is re-read live from the mesh (`cl_FEM_SideSet.cpp:49`) and re-stamped unconditionally
(`cl_MaxwellFactory.cpp:733-738`). No writer stamps the coil/cond-air trio after admission, deck
strings cannot produce interface types, and git history shows the trio was never whitelisted —
so **the three "must be disabled!" errors are unreachable in this tree**, and Gregory's abort
implies version skew or a different message than assumed.

**D3, found while tracing:** `InterfaceCondFerro` *is* whitelisted, dof-carrying, and fully
active, but has **no case** in `IWG_Maxwell::link_to_group` — a conductor block touching iron
aborts `"Not implemented"` at link time. A 2D undulator whose coils are declared `conductor`
and sit on the iron poles hits exactly this. New steps R9 (fix D3, Christian's call on kernel
case vs deactivation) and a sharpened R3 ask (exact message text, commit hash, coil declaration,
coil-iron adjacency) are in the todo.

## Evening: Gregory answered R3 — D2 was misidentified, the real defect is D4

Gregory's error is **`fn_check_facet_orientation.hpp:78`** ("Could not determine orientation of
facets."), and **his model has no coils** — 2D, HTS thin-shell lines in the poles touching only
air, blocks 2/3 ferro. The whole coil-interface reconstruction of D2 is refuted by the reporter;
the round-1/round-2 code findings (D2' staleness, D3, R6) stand as code analysis but were never
his bug.

**D4, diagnosed and source-traced:** `check_facet_orientation` implements only the 3D
shared-edge test — it needs a shared node *pair*, which two distinct 2D line facets can never
have (adjacent segments share one node). Its only caller is the BFS in
`MaxwellFactory::fix_facet_masters` (`cl_MaxwellFactory.cpp:1298`, cuts path), whose 2D
facet-to-facet adjacency is node-based (`cl_Mesh_ConnectivityCalculator.cpp:375-406`). Seeds =
cross-block-type facets (his ferro-air interfaces); flagged = same-type facets (his air|air tape
lines); the first propagation step aborts. Structural: **any 2D cohomology model with a
same-type sideset facet node-adjacent to a cross-type facet dies here.** The node-based
adjacency crossing chains at a point contact is exactly Gregory's "blocks that were not even
neighbors" observation.

**His comment-out is silently wrong:** the fall-through `return false` makes the BFS flip every
flagged facet it reaches, which *preserves* the element-ID-based relative orientation —
consistency along each tape is inherited from gmsh element numbering, not established. Same
it-works-by-accident shape as his D1 patch.

**Proposed fix (R10, not applied — Christian's call):** dimension-aware
`check_facet_orientation` — keep the shared-edge loop for ≥3-node facets, add the 2-node chain
rule on `original()->id()`: consistent iff `a1==b0 || a0==b1`, inverted iff `a0==b0 || a1==b1`,
error otherwise. The twin cases reproduce the 3D semantics exactly. Filed alongside: R12 (the 3D
loop compares raw node ids, not originals — potential duplicate-blindness, unverified) and O4
(should 2D BFS propagate across chains at point contacts — semantics, deferrable).

## R10 applied and jury-audited (round 3)

On Christian's go-ahead, `check_facet_orientation` became dimension-aware: the ≥3-node
shared-edge loop is unchanged; two-node facets take the chain rule on `original()->id()` —
consistent iff one traversal continues the other (`tA1==tB0 || tA0==tB1`), inverted iff both
start or both end at the shared node, error preserved when no node is shared. Coincident twins
reproduce the old rule exactly. The load-bearing assumption — facet `node(0)→node(1)` is the
master element's boundary traversal, reversed by `flip()` — was verified in-house
(`cl_Facet.cpp:36-55`, `cl_Element_TRI3.hpp:50-68`) and then independently by both auditors.

**Round 3 (blind jury, both legs returned): both auditors accept the hunk, no P0.** The round's
catch — raised independently by Codex and Grok — was that my O4 side-claim ("chain consistency
holds regardless of which seed arrives first") is wrong for a seed touching an **interior**
node of a same-type chain: the BFS orients both arms from itself and they never compare against
each other — the chain tears. Endpoint contacts (Gregory's geometry) are safe. Reworded as a
latent P1 in the BFS + node-based adjacency (pre-existing, exposed rather than introduced by
R10), routed to Christian. Post-round polish, labelled: `tA0/tA1/tB0/tB1` rename (both auditors
flagged the `a`-prefix collision) and a comment wording fix; syntax gate re-run green.

**Reviewed, not verified** — static trace plus `-fsyntax-only`; no run. R8 is the gate.

## D5: the R8 run got further and broke differently — and the O4 warning went live

The debug undulator run cleared `fix_facet_masters` (R10 worked) and died in the cut pipeline:
"No surface found on tape/boundary" (`cl_CutFactory.cpp:2210`). Both round-3 auditors had named
this failure class hours earlier.

**Mechanism:** the 2D thin-shell cut pipeline silently requires a **uniform master side per tape
sideset** — never established anywhere, only inherited from gmsh element numbering (lower
element ID becomes master, and gmsh groups elements per physical surface).
`relink_slave_elements_with_duplicate_nodes` collects slave-side *blocks* over a tape's facets:
uniform masters give one block, so one side gets duplicates. **A single flipped facet puts both
blocks in the set, both sides get duplicates, and the original nodes are stranded** — so
`orient_terminal_curves_2D`, which searches the live sideset for a facet carrying the terminal
segment's original nodes, finds none. R10 made 2D propagation *reachable*, and it flips
selectively, so mixed sides are guaranteed once a tape touches seeds at more than one point —
the normal case in an undulator.

**Correction to the earlier entry:** Gregory's comment-out was *not* harmful in 2D. With the
error removed the check returns false for every distinct 2D pair, so the BFS flips *every*
reached facet — a uniform reversal per component, which **preserves** the invariant. That is why
his run passed the cut pipeline while a pairwise-*correct* rule broke it: flip-all and flip-none
both preserve uniformity; only selective flips destroy it. His workaround's real cost is
3D-scoped plus the lost diagnostic. R11's message to him carries the corrected rationale.

**Fix (R13, Christian chose option B):** `fix_facet_masters` gains
`tPropagate = ( number_of_dimensions() == 3 )` guarding four sites — connectivity build,
same-type flagging, BFS loop, connectivity teardown. Cross-type seed normalization still runs in
2D, which is all any 2D consumer ever used. **3D is bit-identical**; 2D returns to its
historically-working behaviour and skips a connectivity build it never needed. The R10 chain
rule stays correct, just unreached in 2D. No jury round: protocol §11 points at the executable
gate, not another review, for a four-site guard that restores prior behaviour in a
already-traced mechanism the auditors had already flagged.

**Left deliberately open (O5):** the uniform-master invariant is still implicit and still rests
on gmsh numbering. Enforcing or asserting it in CutFactory would have turned today's cryptic
downstream assert into a one-line diagnostic at the true cause — recommended as a follow-up with
its own round, not folded into a hotfix.

## Still open (in the todo file)

- **R2** — tell Gregory why his uninitialized D1 variant is unsafe (matters more than R1).
- **R11** — Gregory reverts his comment-out (with the §9.2-corrected rationale).
- **O5** — make the per-tape uniform-master invariant explicit in CutFactory.
- **R9** — D3 fix direction (kernel case vs post-condensation deactivation): Christian's call.
- **O4** — interior-junction chain tear in the BFS (latent P1, no known deck hits it): fix
  direction is Christian's call.
- **R7 / O1 / O2 / R12** — design decisions and the 3D raw-id follow-up, not silently resolved.
- **R8** — undulator benchmark in a debug build; gates D1 and D4 together.
