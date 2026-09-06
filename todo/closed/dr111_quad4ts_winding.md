# DR-111: The 2-D Thin-Shell Layer Elements Are Wound Clockwise

**Date:** 2026-08-30
**Purpose:** `ThinShellFactory` builds every 2-D thin-shell layer element (`QUAD4TS`, and the
dormant `QUAD9TS`) with a negative Jacobian determinant, because the 2-D branch derives its
extrusion normal by rotating the facet tangent *clockwise* while the element winding assumes a
counter-clockwise one. Reverse the in-plane node pair so the element is wound CCW, leaving every
quantity the magnetic path consumes bit-identical.
**Module:** `src/fem/kernel` (+ `src/fem/interpolation/nedelec`, `src/mesh`)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice)
**Status:** ✅ **CLOSED 2026-08-30 — DR-111 STRUCK on Christian's ruling ( DR-42/DR-49 pattern ).
Struck is not verified: one gate is owed and is named below.**
What landed: the fix moved from the factory to the CONSUMER on Christian's condition that the
factory not be touched and all changes stay within QUAD4TS — `Calculator::dV_quad4ts` and
`Pipette::measure_quad4ts`, both `|·|` wrappers wired at a single dispatch site each, both shared
bodies left byte-identical. Two plan-audit rounds ( Codex terra/high→xhigh, Grok high→xhigh ) plus
a code round returned unanimous functional CONFIRMED; Grok's round-2 catch that abs-on-the-
assignment would ship SILENTLY INERT ( shared `mDetJ` cache, k=0 hit ) is what fixed the shape.
The row's own root cause and severity were both refuted and rewritten in place, not merely struck.
Verified by execution: `make check` green; fix present in the binary by SYMBOL; 750+ coupled
timesteps on a 2-D thin-shell quench deck.
**Residual ( R4′ ):** that run used a RELEASE binary with the guard compiled out, so it cannot
discriminate fixed from broken-but-silent — and by the severity correction the pre-fix binary would
have given the same temperatures. An assertions-ON run of the same deck, on a binary first confirmed
to contain the string `Negative Jacobian determinant`, is the only thing that would make this
verified. By-catch filed as DR-151 / DR-152 / DR-153 — all three since struck ( DR-151/152 on 2026-08-30, DR-153 on 2026-08-31 with the orphan deleted and its captive pyramid fix ported ).

> **Scope guards:**
> - IN scope: the `LINE2 → QUAD4TS` winding, the `LINE3 → QUAD9TS` sibling (dormant, see §5 O1),
>   and the `Calculator::dV` tier demotion Christian ruled on 2026-08-30.
> - OUT of scope: the 3-D branches (`TRI3 → PENTA6TS`, `TRI6 → PENTA18TS`) — verified correct
>   by construction, §1. The `MeshChecker` swap-table entry for `QUAD4TS`/`QUAD9TS` is
>   diagnosed here but **not** changed (see D2).
> - No change to which side of the sideset the tape is extruded to — the geometry is untouched,
>   only the order the four corners are listed in.

---

## 1. Current Behaviour and How It Fails

`ThinShellFactory::process_nodes_line2` builds each LINE2 facet normal as
(`src/fem/kernel/cl_ThinShellFactory.cpp:923-924`):

```cpp
tN( 0 ) = tB( 1 ) - tA( 1 );      //  n = (  ty, -tx )
tN( 1 ) = tA( 0 ) - tB( 0 );      //  = tangent rotated CLOCKWISE
```

`compute_distances` (`:1315-1347`, centring at `:1347`) returns offsets ascending from `-t/2` to `+t/2`, so the "top"
layer always lies at `+h·n` from the "bottom" layer, `h > 0`. `create_elements_on_blocks_line2`
then winds the quad as `( bottom0, bottom1, top1, top0 )` (`:1543-1546`), whose edges are `t`
then `h·n`. Its signed area is therefore

    A = t × (h·n) = h · (t × n) = -h·|t|  <  0        ( n is normalized at `:926` )

for every facet at every orientation, for any physically meaningful thickness — `n` is *derived from* `t`, so reversing the facet
reverses both and the handedness survives. `QUAD4TS` interpolates with the plain CCW-referenced
`QUAD4` Lagrange functions (`cl_IF_InterpolationFunctionFactory.cpp:136-140`,
`cl_IF_QUAD4.hpp:49-64`), so `det J < 0` at every integration point.

The 3-D branch does the same thing correctly: `process_nodes_tri3` takes
`n = (B-A) × (C-A)` — the right-hand normal (`:1101`) — and `create_elements_on_blocks_tri3`
stacks `(bottom 0,1,2)` then `(top 0,1,2)` (`:1635-1641`), so the extrusion always agrees with the
reference ζ direction.

| Failure | Mechanism | Evidence |
|---|---|---|
| 2-D thin-shell thermal assembly aborts at the first element | `Calculator::dV` → `dV_ts` → `det( J )` on a CW quad | `Assertion adV >= 0.0 failed`, `adV = -5.000000000022252e-12`, gdb transcript at `cmake-build-claude/gate_dr72/gdb_probe.log` |
| Silent negative mass/volume once the guard is demoted | same det, consumed unguarded in release | `cl_FEM_Calculator.hpp:2050-2059` (the check itself at `:2055`) |
| `Pipette::measure` returns a negative area for these elements | `measure_linear` integrates the signed det | `cl_Pipette.cpp:144-148`, `measure_linear` body |
| Exported meshes carry CW quads | `QUAD4TS` is written out as `QUAD4` | `cl_Mesh.cpp:2639-2642` |

**Numeric reproduction (independent of the build).** Replicating the four functions above in
isolation — facet normal, `compute_distances`, the winding, and `cl_IF_QUAD4.hpp`'s `dNdXi` — on a
20 µm facet with a 1 µm layer gives `detJ = -5.000000e-12` for a tape along `+x`, along `+y`,
at 37°, **and for the reversed facet**; the same replication of the 3-D branch gives `+2.0e-16`
for CW, CCW and tilted triangles alike. That reproduces the measured gdb value to every digit it
prints. Confidence: **high** (arithmetic, re-derived twice, matches the instrument).

**Bottom line:** the 2-D extrusion normal has the wrong handedness for the winding the factory
then applies, unconditionally and since the element was introduced — and nothing noticed because
the magnetic path never computes a signed determinant.

### 1.1 Why this survived every audit and every run

`EF_QUAD4TS` — the only consumer on the magnetic path — is *structurally blind* to it:

- its volume weight is `mDetJ = 0.25 * mThickness * mLength`, a product of two magnitudes
  (`cl_EF_QUAD4TS.cpp:128`), never a determinant;
- its in-plane tangent comes from the **LINE2 facet**, not from the layer element's own nodes
  (`cl_EF_QUAD4TS.cpp:38-42`), and `BlockData::link_thin_shell_facets_serial` attaches that facet
  to each layer element (`cl_FEM_DofMgr_BlockData.cpp:484`);
- its through-thickness direction is read from the actual node positions with a comment saying
  outright that "the stacking side is not guaranteed to be `+mN`" (`cl_EF_QUAD4TS.cpp:79-85`);
- its edge signs come from `compute_edge_directions_thinshell`, which reads
  `mElement->physical_tag()` and never looks at node order at all
  (`cl_FEM_Element.cpp:1420-1436`).

The thermal path is the first and only consumer of the raw scalar-Lagrange determinant on these
elements, which is why the first-ever 2-D thin-shell thermal run failed at step one.

---

## 2. Architecture: Reverse the In-Plane Pair, Not the Through-Thickness Pair

Three candidate repairs. The deciding constraint is a contract that is easy to miss:
`link_elements_with_edges` assigns local edge 0 to the **bottom** layer's mesh `Edge` and local
edge 1 to the **top** layer's, positionally and independently of the node list
(`cl_ThinShellFactory.cpp:1934-1957`), while `Element_QUAD4TS::get_nodes_of_edge` declares edge 0
to be nodes `{0,1}` and edge 1 to be nodes `{3,2}` (`cl_Element_QUAD4TS.hpp:104-134`). **Nodes
`{0,1}` must remain the bottom curve and `{2,3}` the top curve**, or the two disagree.

| Candidate | Winding | det J | Verdict |
|---|---|---|---|
| **A — reverse the in-plane pair** | `( b1, b0, t0, t1 )` | `+5e-12` | **CHOSEN.** `{0,1}` still bottom, `{2,3}` still top; both edges keep the parallel-tangent property their header comment requires |
| B — swap bottom/top | `( t0, t1, b1, b0 )` | `+5e-12` | **REJECTED.** Node pair `{0,1}` becomes the *top* curve while local edge 0 is still linked to the *bottom* layer's mesh `Edge` — the two contracts disagree, and `mNablaEta` flips |
| C — flip the normal handedness in `process_nodes_line2` | unchanged | `+5e-12` | **REJECTED.** This moves the layer stack to the mirror side of the sideset. A geometry change, not a renumbering; out of the scope guard |

Under **A**, every quantity `EF_QUAD4TS` derives is provably unchanged:

| Quantity | Source | Under A |
|---|---|---|
| `mJ`, `mNablaXi`, `mN`, `mLength` | the LINE2 **facet** (`cl_EF_QUAD4TS.cpp:38-42`) | untouched — the facet is not renumbered |
| `mNablaEta`, `mCoeffs`, `mGrad` | `0.5*(n2+n3-n0-n1)` (`:79-85`) | `0.5*(t0+t1-b1-b0)` = `0.5*(t0+t1-b0-b1)` — **identical**, confirmed numerically (`(0,-1e-6)` before and after) |
| `mS` | `physical_tag()` (`cl_FEM_Element.cpp:1420-1436`) | untouched — node order is not an input |
| local edge ↔ mesh `Edge` | positional (`cl_ThinShellFactory.cpp:1934-1957`) | untouched, and still consistent with `get_nodes_of_edge` |

So the magnetic solve should be **bit-identical**, and the thermal determinant flips sign. That
prediction is the falsifiable claim this plan stands on (§4, R4). Confidence: **medium-high** —
it is a source trace, not a run.

### 2.1 D3 — the `mS` flip (Claude, self-caught; Grok CONFIRMED independently)

> **Downgraded 2026-08-30 from "fix A is refuted" to "fix A is unproven"** after the round: the
> flip is real, but it is compensated (§2.3). The original text is kept because the mechanism it
> documents is load-bearing.

**Severity: CRITICAL to the plan.** §2 claimed the reversal leaves everything `EF_QUAD4TS`
consumes unchanged. That is wrong on `mS`.

The layer block element takes the **generic** `Element::compute_edge_directions`
(`cl_FEM_Element.cpp:1337-1400`, called from the block-element constructor at `:47`), not the
`physical_tag`-based thin-shell variant — that one runs on the *sideset* element
(`cl_FEM_Element.cpp:936-939`). The generic routine derives each sign by comparing the mesh
`Edge`'s node IDs against `get_nodes_of_edge( e )`, i.e. **against the element's node order**. The
layer `Edge` objects are built with node copies in the *original facet edge order*
(`cl_ThinShellFactory.cpp:1858-1867`). So today, element edge 0 = nodes `{0,1}` = `(b0,b1)`
matches the facet direction and `mS[0] = +1`.

Reversing the in-plane pair flips both `mS` entries, while `mNablaXi` keeps following the **facet**
(`cl_EF_QUAD4TS.cpp:38-42`). The product `mS[k] * mNablaXi` is the Nédélec basis direction, so the
basis inverts for every layer element — and the interface coupling assembled on the sideset element
derives its signs from `physical_tag` instead, so it would **not** follow. The two sides would
disagree.

**The unstated invariant this exposes:** *the layer element's local edge node order must follow the
facet edge direction*, because `mS` is read from node order while `∇ξ` is read from the facet.
Fix A violates it; keeping it forces the winding `( b0, b1, t1, t0 )`, which is clockwise. **No
renumbering of the four corners can be both CCW and invariant.** Confidence: **high** — source
trace through four files, each re-opened.

This is also the second sighting of the same coupling: the 2026-08-04 layer-alternation work
recorded "`process_nodes_line2` stacks along (dy,−dx), the EF assumed (−dy,dx)" and resolved it by
making the *edge function* follow the factory
(`devlog/dl20260804_2d_thinshell_layer_alternation.md:96-101`), classifying the residue as
"factory-vs-EF normal convention (global sign only)".

### 2.3 D4 — why the flip is survivable (Grok, adjudicated by Claude 2026-08-30)

`hang_thinshell_edges_on_nodes_bottom` (`cl_MaxwellFactory.cpp:1776-1834`) labels the shell's
bottom nodes **by position** (`tNodesOnThinShell( k )->set_index( k )`) and then reads
`tNodesOnVolume( tEdge->node( k )->index() )`. The label is a permutation applied to both sides, so
for the 2-node LINE2 edge the source **set** is unchanged and only its **order** reverses. With
`q_edge = φ_src0 − φ_src1` (`cl_FEM_DofMgr_DofData.cpp:3639-3647`) the dof is exactly negated.
`mS` is negated too (D3). `h_t = mS · q` is therefore unchanged, and `EᵀE` / `CᵀC` were already
sign-invariant.

Codex read the same swap as a geometric cross-wire and rejected the plan on it. That reading would
be right if the tie were per-node; it is per-edge. **Adjudicated in Grok's favour on a first-hand
trace — but at ~70 %, which is why R4 exists.**

**Not assumed for LINE3/QUAD9TS:** that edge has a third node, and a reversal there is a
permutation whose sign story needs its own trace. R2 must not ride on this argument (see O1).

### 2.4 Round-2 defects — both changed the implementation, neither changed the design

**D5 — the `mDetJ` cache would have made the fix silently inert (Grok, round 2; Claude verified).**
`invJ2D3D` (`cl_FEM_Calculator.hpp:1713-1727`) writes the same `mDetJ` / `mDetJIndex` cache that
`dV_ts`'s fallback reads, and the thermal kernels call `B( k )` **before** `dV( k )`
(`mt_thermal_h.cpp:25`, then `:45-49`). QUAD4TS is linear, so `mIsCurved == false` and `dV_ts`
uses `tIndex = 0`. At k = 0 that is a **cache hit**: an implementation that took `abs` on the
inner `mDetJ = det( ... )` assignment would return the signed value and never execute the `abs`.
Gauss point 0 would still be negative and would still abort in debug — the fix would have looked
conservative and been inert exactly where it is measured. **Taking `abs` on the returned value
(the wrapper form) is immune.** This is why `dV_quad4ts` wraps `dV_ts` rather than copying it.

**D6 — `EF_QUAD4TS::det_J()` is not unconditionally positive (Codex, round 2).** It is
`0.25*mThickness*mLength` (`cl_EF_QUAD4TS.cpp:128`) and deck thickness carries **no positivity
check** (`cl_MaxwellFactory.cpp:3083-3094`; `cl_FEM_Kernel.cpp:987` asserts only `abs(t) > 0`).
The wrapper therefore also guards the edge-function branch. For valid (positive) thickness
`std::abs` on a positive double is bitwise identity, so the magnetic path is unchanged; for the
invalid case it is strictly better than the status quo. Filed as by-catch (R6′).

**Round-2 qualifications accepted, neither load-bearing.** Codex: G2's wording "the *only*
orientation-sensitive quantity" is too broad — thin-shell resistivity can consume a sideset normal
via `bn_angle`, but that normal comes from the *magnetic* calculator and is untouched here.
Grok: `dV_ts` evaluates one determinant at `tIndex = 0` for linear elements, so F′ corrects the
sign of the existing constant-J approximation rather than introducing per-point `|det J(ξ_k)|`
integration — exact for a parallelogram, unchanged in character for the averaged-normal trapezoid.

### 2.5 F′ — the chosen fix (Christian's ruling, 2026-08-30)

Christian rejected A: "I am reluctant to fix a Factory that is working; flipping the mS entries
seems a HUGE risk." On re-read that ruling has a principled footing, not just a risk one:

1. **The winding is forced, not accidental.** D3 established that the edge-linking invariant
   (local edge node order must follow the facet edge direction, because `mS` is read from node
   order at `cl_FEM_Element.cpp:1361-1372` while `∇ξ` is read from the facet at
   `cl_EF_QUAD4TS.cpp:38-42`) admits exactly one corner order — `( b0, b1, t1, t0 )` — and that
   order is CW given the CW extrusion normal. The factory therefore implements a *consistent
   left-handed convention*, and every magnetic consumer (edge signs, hang T-matrices,
   `get_bottom_nodes`/`get_top_nodes`, ghost facets) was built and validated around it.
2. **The genuine defect is in the consumer.** For an orientation-reversing element map, `N` is
   J-independent and `B = J⁻¹·dN/dξ` (`Bscalar`, `cl_FEM_Calculator.hpp:1833-1845`) gives exact
   physical gradients regardless of the map's handedness; the *only* orientation-sensitive
   quantity in `T_h_picard`/`T_h_newton` is the volume weight. The mathematically correct weight
   for any diffeomorphic map is `|det J|` — standard change-of-variables. So taking the absolute
   value in `dV_ts`'s no-edge-function branch (`cl_FEM_Calculator.hpp:2099-2108`) yields exactly
   the right M, K and f, not an approximation.
3. **Zero magnetic risk by construction.** The magnetic solve enters the *other* branch of
   `dV_ts` (`mEdgeFunction->det_J()`, `:2094-2097`) and never executes the changed line. Only
   `QUAD4TS`/`PENTA6TS`/`HEX8TS` without an edge function route here (`cl_FEM_Calculator.cpp:
   1437-1439, 1478-1481, 1490-1493`), i.e. the thermal / Nédélec-free path — precisely the path
   that is broken today. PENTA6TS's det is already positive, so `abs` is the identity there.

F′ supersedes A, C and F-as-originally-stated (F proposed `abs` *plus* declaring the mirror a
wart; under 1. the mirror is a documented convention, which is what makes `abs` honest).
`Calculator::dV`'s general path is NOT given an `abs` — for ordinary volume elements a negative
det remains a real defect that the assert must keep catching.

### 2.2 Options if F′ is refuted in round 2 — HISTORICAL, superseded by §2.5

| # | Option | Cost | Risk |
|---|---|---|---|
| **G** | Source `∇ξ` in `EF_QUAD4TS::link` from the element's **own** edge 0 instead of `aElement->facet()`, removing the invariant; then fix A becomes self-consistent | 2 lines in the EF + R1 | **Not** a no-op on curved tapes: the layer curve is not parallel to its facet when adjacent node normals differ, so `mLength` and the tangent change slightly. Changes the magnetic path for real |
| **C** | Flip the normal handedness in `process_nodes_line2` | 2 lines | Reverses which physical side the stack is extruded to, so the **material layer order mirrors** relative to the sideset's master/slave. A physics change, not a renumbering. Previously rejected in §2 and still the most dangerous |
| **F** | Leave the geometry; give the thermal path `std::abs( det J )` for thin-shell types only | small, contained | Papers over a genuine mirror rather than fixing it; leaves `Pipette` and mesh export still negative. Honest only if the mirror is declared *intended* |

Recommendation withheld pending the audit round: G is the only one that removes the coupling
rather than relocating it, but its "not a no-op on curved tapes" property means it needs its own
magnetic-parity gate on a curved deck, not just `2D_Tapestack`.

---

## 3. Gap Table

| # | State / behaviour | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | `QUAD4TS` corner order | thermal `det J` | **no** — CW | (c) | `cl_ThinShellFactory.cpp:1543-1546` |
| 2 | `QUAD9TS` corner order | same, when higher-order shells land | **no** — same defect, dormant | (c) | `:1589-1592`; unreachable, `cl_FEM_Calculator.cpp:1447-1453` errors out first |
| 3 | node pair `{0,1}` = bottom curve | edge↔dof consistency | yes | (a) | preserved by A |
| 4 | `EF_QUAD4TS` η vector | magnetic curl + gradient | yes | (a) | invariant under A, §2 |
| 5 | `mS` edge signs | Nédélec continuity | yes | (a) | from `physical_tag`, node-order independent |
| 6 | `Pipette` measure of a layer element | `MeshChecker` only | signed, currently negative | (a) | becomes positive; no other caller — the kernel measures **facets**, `cl_FEM_Kernel.cpp:946-953` |
| 7 | `Calculator::dV` guard tier | release cost | `BELFEM_ERROR` inside a debug-only `#if` | (c) | `cl_FEM_Calculator.hpp:2050-2059` (the check itself at `:2055`); Christian's ruling 2026-08-30 |
| 8 | mesh export node order | ParaView / Exodus well-formedness | CW today | (a) | improves as a side effect, `cl_Mesh.cpp:2639-2642` |

### 3.1 Cross-cutting findings

- **D1 — the winding is unconditional, not orientation-dependent.** The DR-111 register row says
  the winding "is CCW only when the tape normal points LEFT of the LINE2 facet direction — a
  right-hand pairing that nothing enforces." That is wrong: `n` is constructed *from* `t`, so no
  sideset orientation can produce a CCW quad. Every 2-D thin-shell layer element BELFEM has ever
  built has a negative determinant. The row must be corrected when it is closed.
- **D2 — `MeshChecker`'s swap table would corrupt a `QUAD4TS`, not repair it.** `link_to_block`
  maps `QUAD4TS → swap_quad4` (`cl_MeshChecker.cpp:93-97`), which swaps nodes 1↔3
  (`cl_MeshChecker.hpp:207-211`). Applied to a thin-shell element that makes node pair `{0,1}`
  the *tangential* pair, so `EF_QUAD4TS` would read the tape **length** as its thickness:
  `2e-5` instead of `1e-6`, a factor of 20. It is harmless today only because the checker runs at
  mesh load and these elements are born later. **Therefore: do not "fix" DR-111 by running the
  checker over the layer blocks** — an approach this plan explicitly rejects. Whether the swap
  table should refuse thin-shell types outright is O2.
- The `#if !defined( NDEBUG ) || defined( DEBUG )` wrapper in `Calculator::dV` is
  character-for-character the definition of `BELFEM_ASSERTIONS_ACTIVE` (`assert.hpp:26-33`), and
  `BELFEM_ERROR`/`BELFEM_ASSERT` expand to identical bodies (`assert.hpp:243-276`). The current
  spelling is therefore already an assert; R3 makes that honest rather than changing behaviour.

---

## 4. Ordered Steps

> **Steps rewritten 2026-08-30 for F′, then CONSTRAINED the same day on Christian's condition:
> "the changes are constrained within the QUAD4TS element" (+ explicit permission to touch
> `Pipette`). The original R1/R2 winding steps are withdrawn — the factory is not touched.**
>
> The constraint changed the shape of the fix. `dV_ts` is shared by `QUAD4TS`
> (`cl_FEM_Calculator.cpp:1439`), `PENTA6TS` (`:1481`) and `HEX8TS` (`:1493`), and
> `Pipette::measure_linear` is shared by `QUAD4TS` (`cl_Pipette.cpp:146-147`) and the 3-D
> `default:` branch (`:184-187`). Editing either body in place would touch the 3-D thin shells.
> Both fixes are therefore **new QUAD4TS-only functions plus a one-line dispatch change**, leaving
> the shared bodies byte-identical.

- [x] **R1′ — QUAD4TS-only volume weight.** *(landed 2026-08-30)* Add `Calculator::dV_ts_quad4ts` (declaration beside
      `dV_ts` at `cl_FEM_Calculator.hpp:1440`, body beside `:2080`): identical to `dV_ts` except
      `std::abs` on the no-edge-function branch. Change **one** line,
      `cl_FEM_Calculator.cpp:1439`, to assign it. `dV_ts`'s body is untouched, so `PENTA6TS`
      (`:1481`) and `HEX8TS` (`:1493`) are unaffected by construction.
- [x] **R2′ — QUAD4TS-only Pipette measure.** *(landed 2026-08-30)* Add `Pipette::measure_quad4ts` (an `abs` wrapper on
      the `measure_linear` integral) and assign it at `cl_Pipette.cpp:146-147` only.
      `measure_linear`'s body is untouched, so the 3-D `default:` branch (`:184-187`) is
      unaffected. Removes the negative layer-element area.
- [ ] ~~**R3′ — demote the `dV` guard.**~~ **HELD 2026-08-30** — `Calculator::dV` is the general
      path for every element type, so it breaches the QUAD4TS-only condition. Christian's earlier
      ruling stands and it is a semantic no-op in every build configuration (Codex Q6); it lands
      as its own separate change on his word, not inside this one.
- [x] **R3′′ — document the convention.** *(landed: `nedelec_thinshell.md` §2a)* Short note in the thin-shell doc: QUAD4TS layer elements
      are left-handed **by construction** (the forced-winding invariant of §2.5), `N`/`B` are
      orientation-independent, and the weight is therefore `|det J|` — so the next person who meets
      the negative det finds the ruling instead of re-deriving fix A.
- [ ] **R4′ — the gate. THE ONE RESIDUAL.** Debug **`belfem`** ( NOT `hphiTrun`, retired 2026-08-30 ) on the preserved deck at `cmake-build-claude/gate_dr72/`:
      completes step 1 instead of aborting (the register row's own fixed-gate). Magnetic parity
      on `2D_Tapestack` is expected **bitwise** — the changed lines are not on that path at all —
      so any difference refutes the branch analysis outright rather than needing interpretation.
- [x] **R5′ — regression.** `make check` green on the rebuilt binary (Christian, 2026-08-30).
- [x] **R6′ — close-out.** *(done 2026-08-30: row struck + root cause and severity corrected, DR-151/152/153 filed, devlog, check_doc_claims 37/37)*  Correct the DR-111 root cause per D1 and close the row; file the
      **three** by-catch defects as register rows or strike with reasons — (i) no positivity check
      on layer thickness (`cl_MaxwellFactory.cpp:3083-3094`), (ii) the `MeshChecker` swap-table
      trap (D2), (iii) `src/fem/postproc/cl_MeshChecker.{hpp,cpp}` is an **orphaned duplicate**,
      in no `CMakeLists.txt` and with a dispatch that does not list `QUAD4TS` — only
      `src/fem/kernel/cl_MeshChecker.cpp` is compiled; devlog; `todo/README.md`;
      `scripts/check_doc_claims.py`.

---

## 5. Open Design Questions

- **O1 — does R2 land in this change?** `QUAD9TS` carries the identical defect but is
  unreachable: `Calculator`'s 2-D dispatch raises "Higher order thin shells are not implemented!"
  for anything but `QUAD4TS`/`QUAD4` (`cl_FEM_Calculator.cpp:1447-1453`). Options: (a) fix now,
  so a known-bad winding is not left beside a fixed one; (b) leave it, keeping the diff minimal.
  Recommendation: (a) — it is unreachable, so it carries no run risk, and the trap it sets for
  whoever implements higher-order shells is real. **Christian's call.**
- **O2 — should `MeshChecker` refuse thin-shell element types?** Per D2 its swap would corrupt
  them. It cannot fire today, so this is a latent trap, not a defect. Options: remove the
  `QUAD4TS`/`QUAD9TS` cases from `link_to_block`, hard-error on them, or leave with a comment.
  Not decided here; out of scope for the fix itself.

---

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row mapped to a step or an open question.
- [ ] Each claimed gap backed by a citation re-read at its source this session.
- [ ] Plan audited by Codex **and** Grok before any source edit.
- [ ] R4 magnetic-parity prediction stated **before** the run, with the refutation branch named.
- [ ] `make check` / `check-fast` green.
- [ ] DR-111 row corrected per D1, not merely struck.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/dr111_quad4ts_winding.md`.
- Numeric replication scripts (scratchpad, not committed): 2-D winding, 3-D winding,
  candidate comparison.
