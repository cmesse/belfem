# Buffer Layer Cut Topology for Transport Current {#fem_maxwell_buffer_cut_topology}

**Date:** 2026-04-09
**Purpose:** Document how the cohomology cut graph must be modified when a φ-formulation buffer layer splits the conducting region of a thin-shell stack into two electrically isolated halves.
**Module:** fem/maxwell
**Related:** `todo/thinshell_hphi_formulation.md` (design note, §3.4 three-way coupling), `literature/papers/fem/schnaubelt2023.txt` (NI coil cuts, lines 337–351)

## 1. The problem

Adding a φ-formulation buffer layer between the substrate (hastelloy) and the HTS film (ybco) electrically disconnects the upper and lower conducting halves of the tape. The buffer enforces ∇×H = 0 (no current), so no galvanic current can cross it.

The existing cohomology cut was designed for a single connected conductor. It imposes one transport current I through the entire tape. With two disconnected conductors, the single cut cannot correctly distribute current between the two halves, and the solver produces non-physical results.

## 2. Literature context

Schnaubelt et al. 2023 (NI coils, `schnaubelt2023.txt` lines 337–351) face the same topology: multiple conducting turns separated by an insulating T2TCL. Their solution uses **two automatically generated thick cuts**:

- **IC1 = Isrc** (imposed): the source current entering the current leads.
- **IC2 = Iθ,total** (free): the total azimuthal current, determined by the solver based on the T2TCL contact resistance.

For an insulated coil (no T2TCL shunt), IC2 is deterministic (all current azimuthal, none radial). For an NI coil (finite T2TCL), IC2 is a solver output.

Our case is the "insulated coil" limit: the buffer is a perfect insulator, and there is no rint contact resistance to shunt current between the halves. So the current partition is determined by the terminal geometry and impedances.

## 3. The single-cut rerouting approach

Instead of generating a second cut (Schnaubelt's approach), we **reroute the existing cut** through the buffer node. This exploits the fact that the buffer is now part of the φ region: the cut surface, which previously passed only through external air, now also passes through the buffer layer inside the tape stack.

### 3.1 Original cut graph (no buffer)

At the point where the cut crosses the tape sideset, the φ mini-graph is:

```
    3 (−)
   / \
  1   2
   \ /
    4 (+)
```

- **Node 1**: air node on one side of the cut (e.g. top).
- **Node 2**: air node on the other side of the cut (e.g. bottom).
- **Nodes 3, 4**: duplicate pair at a conductor material interface inside the thin shell. Node 3 is on the minus side of the cut, node 4 on the plus side.

The cut constraint is:

```
φ₄ = φ₃ + I
```

Node 4 is hanging with this constraint. The circulation of H around the entire tape equals I.

### 3.2 Modified cut graph (with buffer)

The buffer layer introduces a new φ-DOF node (5) at the buffer's position in the stack, between hastelloy and ybco. Node 5 is connected to both air endpoints because it belongs to the φ region:

```
  1 --- 3 --- 2
  |           |
  1 --- 5 --- 2     ← buffer node (φ DOF)
  |           |
  1 --- 4 --- 2
```

Modifications:

1. **Release node 4** — it becomes a regular free unknown, no longer carrying the cut constraint.
2. **Make node 5 hanging** with the cut constraint, choosing the correct partner based on orientation (see §3.3).

The circulation of H is now imposed around **only the HTS-side conductor**, not the entire tape. The substrate side has no cut constraint; its current is freely determined by the solver.

### 3.3 Orientation rule

The cut jump must enclose the HTS conductor. The correct pairing depends on which side of the buffer the HTS sits:

| Stack order at cut | HTS is between | Constraint |
|---|---|---|
| ... Has → **buffer(5)** → YBCO(3) ... | nodes 5 and 3 | φ₅ = φ₃ + I |
| ... Has(4) → **buffer(5)** → YBCO ... | nodes 4 and 5 | φ₅ = φ₄ + I |

For the standard REBCO tape (top-to-bottom: Cu/Ag/Has/buffer/YBCO/Ag/Cu), the HTS sits below the buffer. The actual constraint depends on the normal direction of the sideset and the layer indexing convention in `ThinShellFactory`. **This must be verified by tracing the node indices at the cut location in the debugger before hardcoding.**

### 3.4 Physical interpretation

- **HTS loop** (between buffer node 5 and the YBCO-side interface): carries the imposed transport current I. This is the superconducting path.
- **Substrate loop** (between hastelloy-side interface and buffer node 5): carries whatever current the physics dictates. For hastelloy at 77 K (ρ ≈ 1.2e-6 Ω·m) vs. YBCO in the superconducting regime (effective ρ ≈ 1e-12 Ω·m), the substrate current is negligible (~10⁻⁶ × I).
- **Total terminal current**: I_total = I_HTS + I_substrate ≈ I. Conservation is automatically satisfied because the terminal BC constrains the surface integral of J.

### 3.5 Relation to Schnaubelt's two-cut approach

This single-cut rerouting is equivalent to Schnaubelt's two-cut approach in the limit where the second cut's coefficient is known (the insulated-coil limit):

| Schnaubelt (two cuts) | Our approach (one rerouted cut) |
|---|---|
| IC1 = Isrc (imposed) | Cut jump at buffer node = I (imposed) |
| IC2 = Iθ,total (free / solver-determined) | Substrate current = free (solver-determined) |
| Two cut surfaces in the φ region | One cut surface, extended through the buffer |

The advantage of the single-cut approach is that it requires no changes to the terminal BC machinery (still one cut, one current). The CutFactory must learn to route the cut through the buffer node, but the rest of the pipeline stays the same.

## 4. Why the CutFactory needs modification

The existing `CutFactory` generates cuts in the volume φ region (air blocks). It does not traverse thin-shell layer nodes, because thin-shell layers are surface-embedded, not volume blocks. When the buffer is classified as a φ block (after the Phase 2 plumbing), the `CutFactory` still does not see it because:

1. `CutFactory` builds its topology graph from volume elements and their facets.
2. Thin-shell PENTA6TS elements are not volume elements.
3. The buffer node (5) does not appear in the volume topology graph.

**Required modification:** The `CutFactory` (or the post-cut wiring in `MaxwellFactory`) must explicitly insert the buffer node into the cut graph at every point where the cut crosses a thin-shell sideset that contains a buffer layer. The steps are:

1. **Identify** which cut edges cross a thin-shell sideset that has a buffer layer.
2. **Find** the buffer layer's node at the crossing position (node 5).
3. **Rewire** the cut constraint: release the old hanging node (4) and make the buffer node (5) hanging with the correct orientation (§3.3).

This is a targeted modification to the cut-sideset intersection handling, not a redesign of the CutFactory's cohomology algorithm.

## 5. Relation to Milestone C (contact-impedance / rint)

When the rint contact-impedance formulation (§9.A.3 of the design note) is eventually implemented, the current partition between the two halves becomes non-trivial. In that case:

- The rint provides a finite-resistance bridge between the two halves.
- The single-cut approach still works: the cut imposes I through the HTS loop, and the rint allows some current to leak to the substrate loop.
- The substrate current is no longer ≈ 0; it depends on the rint impedance.
- Schnaubelt's formula applies: Ir = Isrc − Iθ (radial = source − azimuthal).

No additional cuts are needed for the rint case. The single rerouted cut correctly imposes the source current, and the contact-impedance kernel handles the current transfer.

## 6. Implementation checklist

- [ ] Trace the actual node ordering at a cut-sideset crossing with a buffer layer in the debugger to confirm the orientation rule in §3.3.
- [ ] Modify `CutFactory` (or `MaxwellFactory` post-cut wiring) to detect buffer layers at cut crossings.
- [ ] Insert the buffer node into the cut constraint graph with the correct orientation.
- [ ] Release the old hanging node at the buffer-adjacent material interface.
- [ ] Verify with a test case: single tape, Cu/Ag/Has/buffer/YBCO/Ag/Cu, ramped transport current. Expected: YBCO carries ~100% of the current, substrate carries ~0%, B-field matches the no-buffer reference.
- [ ] Verify that the approach degrades gracefully when rint is later added (Milestone C).

## 7. References

- Schnaubelt, E. et al. (2023). "Electromagnetic simulation of no-insulation coils using H–φ thin shell approximation." IEEE Trans. Appl. Supercond. **33**(5), 4900906. Section IV, lines 337–351 of local text.
- Pellikka, M. et al. (2013). "Homology and cohomology computation in finite element modeling." SIAM J. Sci. Comput. **35**(5), B1195–B1214. The thick-cut algorithm used by BELFEM's CutFactory.
- Messe, C. et al. (2023). "BELFEM: ..." Supercond. Sci. Technol. **36**, 114001. Section 2 (H-φ formulation and cohomology cuts).
