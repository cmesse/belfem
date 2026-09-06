# Thin-Cut Generation for Non-Unit Thick-Cut Coefficients {#homology_thin_cut_nonunit_rectification}

**Date:** 2026-06-22 (expanded 2026-07-01 with §1.1–1.3 on the geometry of the non-unit pushed object)
**Purpose:** Why coarse 3D meshes (e.g. CCT magnets, where a conductor loops under itself)
produce cohomology generators with edge coefficients `|c(e)| ≥ 2`, why the thin-cut
pipeline cannot consume them, and the graph algorithm that rectifies the generator to unit
coefficients — or certifies, constructively, that no unit representative exists on the given
mesh.
**Module:** src/homology

---

## Summary

BELFEM realizes the scalar-potential (`φ`) jump of the mixed h-φ formulation by converting
a **thick cut** (the integer 1-cochain `c` produced by the cohomology engine) into a **thin
cut** on element faces. The conversion machinery — cut cases, duplication, relinking, static
condensation — assumes every edge coefficient is in `{−1, 0, +1}`. The companion reference
[thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md) works this
out in full on a reference tetrahedron and a 2D hexagon.

When a conductor loops under itself and the mesh is coarse, the engine instead returns a
generator with `|c(e)| ≥ 2` on some edges.

> **Status, 2026-08-31: the rectifier described in §7 has shipped.** `Cohomology::clean()` now
> dispatches to `clean_spfa()` (`cl_Cohomology.cpp:191,358`), which solves the difference
> constraints through `graph::spfa_difference_constraints` (`:455`) and emits an obstruction
> certificate when no unit representative exists; pockets are removed afterwards by
> `remove_cut_pockets()`. The paragraph below describes the **superseded** greedy behavior and is
> kept because the failure it documents is what motivated the replacement.

The old code could not handle this: the greedy
`Cohomology::clean()` tried to remove non-unit coefficients by local node firing with no
global view and no termination guarantee, and the downstream `CutData` machinery is built only
for unit coefficients — a non-unit edge is flagged as support but recorded in neither `±1`
bitset (its weight becomes 0), so `determine_cut_case_3d()` ultimately rejects it as an
invalid cut pattern or non-uniform coefficient, rather than diagnosing the real cause. The
unit-only restriction is structural, not a bookkeeping accident: §1.1 explains why the push
construction has no place to store a second sheet.

The resolution has three parts:

1. `|c(e)| ≥ 2` has a precise geometric meaning — **sheet multiplicity** of the dual cut
   surface — and is *not* always an algorithmic artifact.
2. Whether a unit-coefficient representative exists on the **fixed mesh** is decided exactly
   by a **system of difference constraints**, solvable by **Bellman–Ford**. Feasible inputs
   yield an integer node correction `θ` directly; infeasible inputs yield a **negative-cycle
   certificate** that localizes the under-resolved region.
3. After rectification to unit coefficients the existing thin-cut pipeline applies unchanged
   (the admissible per-tetrahedron patterns are exactly the seven cut cases already coded).

This splits non-unit generators into **two regimes** (§4): *rectifiable* (a unit representative
exists — the thin-cut setup is used unchanged, only the rectifier upstream is fixed) and
*genuinely non-unit* (no unit representative exists on the mesh — the setup is insufficient and
the throat must be refined or duplicated multi-level).

---

## Notation

| Symbol | Meaning |
|--------|---------|
| `K` | the air-region cell complex (mesh of the conductor complement) |
| `e = (i, j)` | an oriented mesh edge, from node `i` to node `j` |
| `c` | the **thick cut**: an integer 1-cochain on edges (a generator of `H¹(K)`); `c(e)` is its coefficient on edge `e` |
| `θ` | the **integer node potential** — an integer 0-cochain on nodes (a discrete gauge / branch function, valued in ℤ). *Distinct from the physical potential `φ`* |
| `d` | the **coboundary** operator (discrete exterior derivative / discrete gradient); `(dθ)(e) = θ(j) − θ(i)` on edge `e = (i, j)`. `c` is *closed* (`dc = 0`); `dθ` is *exact* |
| `c − dθ` | the **rectified** cochain — a different representative of the *same* cohomology class (changing `θ` is a gauge move) |
| `z` | a closed edge-loop (1-cycle) in the mesh; `length(z)` = number of edges in it |
| `⟨c, z⟩` | the pairing of `c` with `z` = number of conductor passes linked by `z` = the transport current it threads |
| `φ` | the **physical scalar potential** of the h-φ formulation (real-valued); jumps by `I` across the cut |
| `I` | the transport current = the imposed `φ`-jump across the cut surface |
| `p`, `q` | an original node DOF and its duplicate across the thin cut; `q = p + I` (multi-level: `q_k = p + k·I`) |
| `w_e` | a per-edge weight (e.g. dual-cell area) in the optimal-representative objective |
| `λ` | an integer branch label per element (covering-space view of multi-level duplication) |
| `pattern` | the per-tetrahedron 6-bit edge bitmask, `Σ_k 2ᵏ·[edge k flagged]` (see §6) |

---

## 1. What `|c(e)| ≥ 2` means geometrically

Under Poincaré–Lefschetz duality, `H¹(K) ≅ H₂(K, ∂K)`: the thick cut **is already** a cut
surface, expressed in the **dual** complex. Each primal edge `e` corresponds to a dual
2-cell, and `c(e)` is the algebraic number of times the cut surface crosses that cell.

Hence `|c(e)| = 2` means **two parallel sheets of the cut surface are squeezed through the
same element layer**. In a CCT magnet this is the throat between adjacent turns: the spanning
surface of the circuit passes through the inter-turn gap once per turn, and if that gap is
one element across, two sheets share the same dual cells.

The two sheets are ordinary, non-intersecting surfaces: in normal-surface terms, two
parallel normal disks stacked inside the same tetrahedron (for example, two vertex caps at
the same vertex crossing the same three edges at different heights). The issue is
resolution: within that element layer, the mesh cannot keep them apart.

### 1.1 Why the push cannot emit a second facet copy

The thin-cut construction sends the interior disk of a cut element to a **face** of that
element (the conjugated face; the diagonal quadrilateral cases emit no face of their own,
§6 — the vertex-cap picture suffices for the argument here). This map has no multiplicity
slot, for two independent reasons:

- a doubled vertex cap has only **one** opposite face — both sheets are pushed onto the same
  conjugated face; and
- the mesh stores exactly **one** facet between an element pair; a facet is a combinatorial
  object, not something the data model can instantiate twice.

The pushed object is therefore not "two thin-cut facets"; it is **one facet with
multiplicity 2**, carrying a potential jump of `2I`. The jump across the thin-cut surface
becomes **non-uniform**: `2I` through the throat, `I` where only one sheet passes. The
duplication scheme — one duplicate per node, `q = p + I`, and a single abstract current
node — hard-codes a *uniform* jump. That is the structural reason the pipeline is
unit-only. What multiplies
in the general case is **node instances**, not facet instances: a node touched by potential
levels `0, I, 2I` needs two duplicates, `q_k = p + k·I` (§7, option C).

### 1.2 Watertightness: a chain, not an embedded surface

The multi-sheet cut is not "torn". Cochain closedness (`dc = 0`, zero oriented sum on
every face) **is** the watertightness condition, and it holds with multiplicities: the
boundary of the multiplicity-weighted facet chain `Σ_f m_f · f` vanishes on interior edges.
The pushed object is watertight **as an integer 2-chain**, not as an embedded 2-manifold.
Where the sheet count drops from 2 to 1 at the rim of the throat, the multiplicity-2 region
ends along a **branch line inside the surface** — the locus where the two sheets separate
onto different facets; the multiplicity balance around each such edge is exactly what
`dc = 0` enforces. Combinatorially consistent, but the two sheets have merged there and
cannot be separated on this mesh. The covering-space picture of
§7 makes this precise: a per-element integer branch label `λ` exists globally (up to the
winding around the generator) *because* `c` is closed, and the thin cut is exactly the facet
set where `λ` jumps, with jump size equal to the sheet count.

With multiplicities, the per-tetrahedron inventory grows beyond the seven unit cases of
§6: non-unit closed patterns decompose into *stacked* normal disks, possibly of different
types within one tet and not always uniquely. Haken's admissibility condition (at most one
quadrilateral type per tetrahedron) is the per-tet embeddability test; globally, the
matching conditions across faces must hold as well.

### 1.3 Can the class be resolved into a thin cut at all?

Split by the regimes of §4:

- **Regime 1 (feasible):** yes — into an ordinary **unit** thin cut. Re-gauging
  `c ← c − dθ` slides the second sheet onto a different element layer: every dual cell is
  crossed at most once, and no multi-sheet object reaches the thin-cut stage.
- **Regime 2 (infeasible):** **no embedded unit thin cut exists on this mesh**; the
  violating loop of §3 is the proof. The class is still resolvable in the generalized sense
  of §1.2: a facet chain with multiplicities, imposed algebraically by multi-level
  duplication `q_k = p + k·I` (§7, option C). That resolution is **topologically exact** —
  every Ampère loop threads the correct total current — but **locally blind where the mesh
  is under-resolved**: the intermediate potential level, i.e. the field between the
  two conductor passes, is squeezed into zero volume, so no element layer remains to carry
  the inter-turn gradient. This is the accuracy argument for preferring
  certificate-guided refinement (option B): it pulls the sheets apart geometrically,
  restoring Regime 1 and with it both topology and local accuracy.

## 2. Existence is guaranteed for the manifold, not the fixed mesh

Poincaré–Lefschetz guarantees the homology *class*, and (constructively, via Kotiuga's
cut algorithm — Gross & Kotiuga, Ch. 6) an **embedded** representative surface in the
*manifold*. It does **not** guarantee a unit-coefficient representative in a *fixed* mesh.

The obstruction is elementary. For any closed edge-loop `z` in the air mesh, the pairing
`⟨c, z⟩` is a topological invariant of the class — the number of conductor passes linked by
`z`. A unit cochain contributes at most `±1` per traversed edge, so

```text
|⟨c, z⟩| ≤ length(z)
```

is **necessary** for a unit representative to exist. If some loop violates it, no
cohomologous cocycle has all coefficients in `{−1, 0, +1}`.

**Example.** A loop around the bore of an N-turn CCT pairs to `N`. If the bore circle is
meshed with fewer than `N` edges, `|c(e)| ≥ 2` is forced for *every* representative. At that
resolution the obstruction is real, not an algorithmic failure — the mesh genuinely cannot
carry a unit cut.

## 3. Feasibility criterion: difference constraints and Bellman–Ford

Changing representative within the class means `c′ = c − dθ` for an integer node function
`θ` — exactly what `clean()` attempts, one node at a time. Requiring `|c′(e)| ≤ 1` on every
oriented edge `e = (i, j)` gives, per edge, the pair of constraints

```text
θ(j) − θ(i) ≤ 1 + c(e)        and        θ(i) − θ(j) ≤ 1 − c(e)
```

This is a **system of difference constraints**, decided by **Bellman–Ford** on the
constraint graph (one shortest-path computation):

- **Feasible ⇔ no negative cycle.** The data are integral and the node–edge incidence matrix
  is totally unimodular, so the shortest-path distances give an **integer** `θ` directly —
  no rounding. Then `c ← c − dθ` is a unit cochain in the same class.
- **A negative cycle is exactly a violating loop** with `|⟨c, z⟩| > length(z)`. It is a
  **constructive certificate** that no unit thick cut exists on this mesh: an **obstructive
  mesh-graph loop** whose edges are the under-resolved region (the too-coarse throat).

This is the same machinery as 2D phase unwrapping (the multivalued angle `φ` ↔ the cut, the
branch cut ↔ the thin cut, the branch labeling ↔ `θ`), where it is solved as a minimum-cost
network flow; see Costantini (1998).

### Optimal representative (optional refinement)

Feasibility (above) is the hard L∞ constraint `‖c − dθ‖∞ ≤ 1`, decided by Bellman–Ford.
*Choosing among* the feasible (unit) representatives is a separate, optional optimization:
pick the geometrically nicest one by minimizing the weighted L1 mass

```text
min over integer θ   Σ_e w_e · |c(e) − (dθ)(e)|     subject to  |c − dθ| ≤ 1
```

which on unit representatives equals the weighted support (each surviving edge contributes
`w_e`). This is an exact **min-cost circulation** problem (LP dual: maximize `⟨c, y⟩` over
circulations `Ay = 0` with `|y_e| ≤ w_e`); weights `w_e` can encode dual-cell area to bias
toward a short, smooth cut. Note this L1 *minimization* is **not** the feasibility test — run
unconstrained it can return non-unit coefficients on an infeasible mesh; it only ranks
already-feasible representatives. Integrality of `θ` comes from the totally-unimodular
incidence matrix and integer `c` regardless of the weights; integer weights additionally make
the dual circulation integral, real area weights need not.

(The object being optimized is, strictly, the cohomology **cocycle** `c − dθ`; under the
duality of §1 it corresponds to a relative-homology cycle, which is the form in which Dey,
Hirani & Krishnamoorthy 2011 state the TU/LP result.)

This places the problem in the tractable corner of homology localization, and it is worth
being precise about why, because the surrounding problem is hard:

- **over ℤ, with a totally-unimodular boundary matrix** (equivalently torsion-free relative
  homology — which holds for a codimension-1 cut in an orientable 3-manifold), the optimal
  representative is polynomial (Dey, Hirani & Krishnamoorthy 2011) — *this is our case*;
- **over ℤ in general** (with torsion) the optimal-cycle / bounding-chain problem is NP-hard
  (Dunfield & Hirani 2011);
- **over ℤ₂** it is NP-hard even to approximate within any constant factor (Chen &
  Freedman 2010/2011).

The codimension-1, integer, totally-unimodular setting we land in is exactly the lucky
intersection where the LP is both correct and cheap.

## 4. Two regimes: rectifiable vs. genuinely non-unit

The feasibility test of §3 partitions every non-unit generator into two regimes with very
different consequences. The distinction is the practical heart of this note, because it decides
whether the thin-cut machinery has to change at all.

**Regime 1 — rectifiable (feasible; the typical case).** A `|c(e)| ≥ 2` here is merely an
unlucky *choice of representative* by the cohomology engine: a cohomologous unit representative
exists on the same mesh, and `c − dθ` reaches it. The rectification happens entirely
**upstream of the thin cut**, on the cochain — it is exactly the job of `Cohomology::clean()`.
The `CutProcessor` thin-cut setup (cut cases, single node duplication `q = p + I`, relinking,
static condensation) then runs **completely unchanged** and never sees a non-unit coefficient.
So non-unit *generators* in this regime are already implementable with the current thin-cut
setup; the response — call it **(A) rectify upstream** — is just to replace the greedy
`clean()` with the global solve of §3. Once the representative is unwound the jump across the
cut is `I` everywhere and no dual cell is crossed twice.

**Regime 2 — genuinely non-unit (infeasible).** Here the negative-cycle certificate proves that
**no** unit representative exists on this mesh (e.g. an N-turn CCT bore meshed with fewer than
N edges). The jump across the throat genuinely is `2I` or more, and a single duplicate per node
(`q = p + I`) cannot represent two distinct potential levels. The current thin-cut setup is
therefore **insufficient as-is**; the supported responses are (B) certificate-guided refinement
— which lengthens the violating loop until the generator becomes rectifiable, returning to
Regime 1 — or (C) multi-level duplication, which extends the duplication+condensation machinery
to `q_k = p + k·I`.

| | Regime 1 — rectifiable | Regime 2 — genuinely non-unit |
|---|---|---|
| Feasibility (§3) | feasible (no negative cycle) | infeasible (negative-cycle certificate) |
| Meaning of `\|c\| ≥ 2` | bad representative choice | true sheet multiplicity forced by the mesh |
| Unit representative on this mesh | exists | does not exist |
| Physical jump | `I` (single sheet) | `2I` or more (multiple sheets) |
| Thin-cut setup | **unchanged** — rectify upstream, then run as-is | **insufficient** — needs B (refine → Regime 1) or C (multi-level) |
| Action | replace greedy `clean()` with the global solve | refine, or extend to `q_k = p + k·I` |

The takeaway: in the common case the `|c(e)| ≥ 2` "problem" is **not a thin-cut problem at
all** — it is a representative-selection problem solved before the thin cut ever runs. Only
Regime 2, the coarse-throat case that originally motivated this analysis, requires changing the
thin-cut machinery itself.

## 5. Relation to the source-field (h-formulation) foundation

The rectification and its objective are the discrete image of classical results in the
h-formulation source-field literature, which is worth making explicit because it grounds the
whole construction in established theory.

- **`c − dθ` is the discrete curl-kernel freedom.** In the continuous h-formulation the source
  field is non-unique: any `h_s + grad χ` has the same curl (Dular et al. 1997, §III.A — the
  kernel of the curl operator, with `h = −grad φ` where the field is curl-free). Choosing the
  integer node function `θ` is exactly the discrete analog of choosing `χ`; both move the
  representative within the cohomology class without changing the physics.
- **The minimal-support objective is theirs too.** Confining `h_s` to a thin neighborhood of
  the cut to save work (Dular et al. 1997, §II.B) is the continuous form of our minimal-weighted-
  support representative; the min-cost circulation of §3 is its algorithmic generalization.
- **The jump is the transport current.** The discontinuity the thin cut imposes is the global
  current `I` threading the cut surface (Dular, Geuzaine & Legros 1999): the cut field is
  `c_i = −grad q_i` (Eq. 9), where `q_i` is a discontinuous scalar equal to **1 on one side of
  the cut and 0 on the other**. That binary `q_i` *is* a unit-coefficient cut, and `⟨c, z⟩` is
  the current `I_i` linked by the loop `z` (Ampère; Gross & Kotiuga §3C).

Seen this way, `|c(e)| ≥ 2` is precisely the regime where Dular's binary `q_i` no longer
suffices — a single sheet cannot carry the jump, and `q_i` would have to step by more than one.
The unit case after rectification reproduces the classical construction exactly; the
multi-level duplication of §7 (option C, `q_k = p + k·I`) is the direct generalization of `q_i`
to a `q_i` that jumps by the sheet count.

## 6. Connection to normal surface theory

The case table of `CutData::determine_cut_case_3d()` — four vertex-cap cases `±1…±4` plus
three diagonal cases `±5…±7` — is exactly the **normal-disk inventory** of a tetrahedron in
Haken's normal surface theory: 4 normal triangle types + 3 normal quadrilateral types. A
unit-coefficient thick cut is precisely a normal surface with at most one disk per
tetrahedron; `|c(e)| ≥ 2` is a normal surface carried with **multiplicities** — in the simplest
case parallel copies of one disk, but in general several normal-disk levels in the same tet —
standard in that theory, merely outside what BELFEM's facet-pushing can represent.

This identification also pins down what rectification buys us. The code labels each tet's cut
by an integer **6-bit edge bitmask** — `pattern = Σ_{k=0}^{5} 2ᵏ · [edge k flagged]` over the
six local edges in the companion document's order (`cl_CutData.cpp:447-475`). Enumerating all
per-tetrahedron unit assignments subject to the cochain closedness condition (zero oriented sum
on every tet face) yields **exactly** the supports below, each nonzero pattern with exactly the
two `±` sign assignments the code admits:

| pattern | edges in support | case | normal disk |
|--------:|------------------|:----:|-------------|
| 0  | ∅ (no cut)        | —  | — |
| 38 | e1, e2, e5        | ±1 | vertex-cap triangle |
| 13 | e0, e2, e3        | ±2 | vertex-cap triangle |
| 19 | e0, e1, e4        | ±3 | vertex-cap triangle |
| 56 | e3, e4, e5        | ±4 | vertex-cap triangle |
| 43 | e0, e1, e3, e5    | ±5 | diagonal quadrilateral |
| 30 | e1, e2, e3, e4    | ±6 | diagonal quadrilateral |
| 53 | e0, e2, e4, e5    | ±7 | diagonal quadrilateral |

The enumeration is short: face closure forces every nonzero unit support to have size 3 or 4;
the only closed 3-edge supports are the four vertex stars (the triangle cases) and the only
closed 4-edge supports are the three edge-pair diagonals (the quadrilateral cases). This
matches the cut-case table of the companion reference exactly. Therefore, **once the generator
is rectified to unit coefficients, the existing thin-cut pipeline (cut cases, duplication,
relinking, condensation) applies unchanged.**

*Scope:* this enumeration assumes the TET4 six-corner-edge convention; it should be rechecked
if higher-order elements ever contribute non-corner edges to the cochain.

## 7. The rectification design ( the approach shipped; this listing is the design )

> The **approach** described here landed: `clean()` dispatches to `clean_spfa()`, which solves the
> difference constraints and emits an obstruction certificate (`cl_Cohomology.cpp:191,358,455`),
> with `remove_cut_pockets()` afterwards (`:629`). The listing below is the design as written,
> and still reads as a proposal — the shipped path is SPFA plus a greedy pass and a dense
> fallback, not every option enumerated here. Read it for the reasoning, not as a description of
> the code.

Replace the greedy `Cohomology::clean()` body with a global rectify-or-certify step (keeping
its signature and call sites):

```text
rectify_to_unit_or_certify( c ):
  1. Build the difference-constraint graph on the mesh node set.
       (On periodic meshes this must be the quotient graph — identify θ across each
        periodic node pair — or the certificates are wrong.)
  2. Run Bellman–Ford (or min-cost circulation with weights w_e for an optimal cut).
  3. Feasible   → apply c := c − dθ; assert closedness and |c| ≤ 1 globally;
                  run the existing thin-cut pipeline unchanged.
  4. Infeasible → extract the negative cycle (parent pointers) and fail loudly with the
                  certificate loop (edge IDs), recommending targeted refinement of the
                  throat it identifies.
```

The feasible branch — **(A) rectify upstream** (§4) — is the common case (Regime 1) and
needs no change to the thin-cut machinery. The infeasible branch (Regime 2) turns a crash (or a
silent hang in the current greedy loop) into an **actionable, physics-aligned adaptivity
instruction** — the obstruction loop is where field gradients are large anyway. Two follow-on
options handle it automatically when refinement is acceptable:

- **(B) Certificate-guided local refinement.** Split edges/elements along the certificate
  loop. `⟨c, z⟩` is invariant under refinement while `length(z)` grows, so guided splits drive
  that loop back below the obstruction threshold; a sufficiently fine — ultimately global —
  refinement is the guaranteed-terminating fallback, since the continuous embedded
  representative always exists. This is the pragmatic general solution.
- **(C) Multi-level duplication (covering-space view).** Assign an integer branch label `λ`
  per element (a lift to the ℤ-cover of the air region); the thin cut is the set of facets
  where `λ` jumps, with jump size equal to the sheet count, and each cut node gets one
  duplicate per level it touches (`q_k = p + k·I`). The hanging-DOF static condensation
  supports this algebraically, but the relinking and case bookkeeping become substantially
  harder. This is the heavy fallback, kept only for when refinement is unacceptable.
  (An in-element enrichment — XFEM / phantom-node, Hansbo & Hansbo 2004 — could represent the
  multi-sheet jump directly, but BELFEM's static-condensation assembly carries no enrichment
  infrastructure, so it is declined in favor of B/C.)

## References

The nine below are in the BELFEM `literature/` library (proprietary, internal use only — cite
the published work and DOI below, never the local `.txt`). Navigation entries:
`literature/papers/topology/index.md` (topology cluster, incl. Haken), `literature/books/`
(Gross & Kotiuga, CLRS), and `literature/papers/fem/index.md` (the two Dular papers).

| Reference | DOI / ISBN | Relevance |
|-----------|------------|-----------|
| Gross & Kotiuga, *Electromagnetic Theory and Computation: A Topological Approach*, Cambridge UP / MSRI Publ. 48, 2004 | ISBN 0-521-80160-5 | The definitive theory of branch cuts and (co)homology for EM; the cut algorithm (Ch. 6) and linking-number/cut foundations (§3C). Grounds the whole thick→thin-cut enterprise. |
| Cormen, Leiserson, Rivest & Stein, *Introduction to Algorithms* (CLRS) | book — cite §"Difference constraints and shortest paths" (§24.4, 3rd ed.) | The canonical source for the difference-constraint ⇒ Bellman–Ford reduction — i.e. the rectification step itself. |
| Dey, Hirani & Krishnamoorthy 2011, *SIAM J. Comput.* **40**(4):1026–1044 | 10.1137/100800245 | *Optimal Homologous Cycles, Total Unimodularity, and Linear Programming* — optimal codim-1 homologous cycle over ℤ is polynomial via LP when the boundary matrix is totally unimodular (⇔ torsion-free). The mathematical backbone of the optimal-representative step. |
| Dunfield & Hirani 2011, *Proc. SoCG '11*, pp. 135–144 | 10.1145/1998196.1998218 | Optimal Bounding Chain Problem; least-area spanning surfaces polynomial when H₂ = 0, NP-complete (over ℤ) in general — the general-ℤ hardness boundary. |
| Chen & Freedman 2011, *Discrete & Comput. Geom.* **45**(3):425–448 (conf. SODA 2010) | 10.1007/s00454-010-9322-8 | Homology localization over **ℤ₂** is NP-hard to approximate within any constant factor — the field that makes the problem hard, and why we stay over ℤ. |
| Haken 1961, *Theorie der Normalflächen*, *Acta Math.* **105**:245–375 | 10.1007/BF02559591 | Normal surfaces and normal coordinates with multiplicities — the cut-case ↔ normal-disk identification. Light identification only; not load-bearing. |
| Costantini 1998, *IEEE Trans. Geosci. Remote Sens.* **36**(3):813–821 | 10.1109/36.673674 | Phase unwrapping as minimum-cost network flow — the exact 2D twin (branch cut ↔ thin cut, branch labeling ↔ θ; a TU-structured LP). |
| Dular et al. 1997, *IEEE Trans. Magn.* **33**(2):1398–1401 | 10.1109/20.582518 | Generalized source field for h-formulations; the **curl-kernel freedom** (`h_s + grad χ`, §III.A) and **minimal-support** construction (§II.B) — the continuous origin of `c − dθ` and the optimal-representative objective. |
| Dular, Geuzaine & Legros 1999, *IEEE Trans. Magn.* **35**(3):1626–1629 | 10.1109/20.767308 | Global current/voltage and circuit coupling in h-formulations; the cut field `c_i = −grad q_i` (Eq. 9) with binary `q_i` — identifies the thin-cut jump as the transport current and fixes the unit-coefficient case our work extends. |

**Not in the library:**

| Reference | DOI | Relevance |
|-----------|-----|-----------|
| Hansbo & Hansbo 2004, *Comput. Methods Appl. Mech. Engrg.* **193**(33–35):3523–3540 | 10.1016/j.cma.2003.12.041 | Element-doubling / XFEM-equivalent in-element discontinuity (solid mechanics) — the in-element enrichment declined under §7 option C. |

## See also

- [thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md) —
  thick→thin conversion, cut-case convention, worked tetrahedron and hexagon examples.
- [cohomology_theory_and_implementation.md](cohomology_theory_and_implementation.md) —
  the Poincaré–Lefschetz reduction in prose.
- [cohomology_algorithms.md](cohomology_algorithms.md) — how the thick cut is computed.
- `literature/papers/topology/index.md` — computational-topology / optimal-cut references
  (Gross & Kotiuga, Dey et al., Dunfield & Hirani, Chen & Freedman, Costantini).
