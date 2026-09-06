# Mathematical Theory Instruction File for Advanced Finite Element Code in Electromagnetism {#homology_cohomology_theory_and_implementation}

> **Status — read before relying on this page.** The reduction algorithms it describes are in the
> build and cataloged in [cohomology_algorithms.md](cohomology_algorithms.md). The
> **manifold-filtering post-processing** described under "Post-Processing Cleanup" is **not**: its
> functions do not exist in the tree, and the Tarjan pocket detection survives only as an archived
> prototype outside the build. The design note in that section has the details; its status is
> Gregory Giard's to confirm.

## Overview
This document provides the mathematical and physical foundation behind an advanced finite element code for magnetostatics using the h-φ formulation in non-conducting domains. It covers key concepts from topology, homology, and cohomology, their computational algorithms, discretization techniques for imposing current boundary conditions, and post-processing for manifold surfaces. This aids in understanding and debugging the implementation, particularly for handling multiply connected domains without modeling conductors explicitly. The theory draws from recent advancements in computational (co)homology tailored to electromagnetic applications.

## Physical and Mathematical Background
In the context of magnetostatics within a non-conducting domain, the h-φ formulation employs the magnetic scalar potential φ as the primary degree of freedom, defining the magnetic field intensity as H = -∇φ, which satisfies the curl-free condition ∇ × H = 0. To ensure uniqueness, φ is gauged by setting its value to zero at an arbitrary reference point. When imposing an electric current I through a conductor that pierces this domain, Ampere's circuital law dictates that the line integral of H around any closed loop encircling the conductor equals I. Topologically, such loops represent non-trivial homology classes in multiply connected domains, necessitating the introduction of a cohomology cut—a surface across the domain—to render it simply connected. Along this cut, a discontinuity [φ] = φ⁺ - φ⁻ = I is enforced, leveraging the Ampere unit of φ to directly link the jump to the imposed current, thereby satisfying the global circulation condition without modeling the conducting domain explicitly.

## Computational Algorithms

The reduction algorithms (CCR, Pellikka, Generalized Pellikka, Smith Normal Form) and their
literature are cataloged with exact code locations in
[cohomology_algorithms.md](cohomology_algorithms.md). The narrative below covers only the
conceptual *why*.

### The Generalized Pellikka Algorithm

To implement topological handling in a discretized finite element framework, the generalized Pellikka algorithm, a novel enhancement proposed by Giard et al., combines the chain complex reduction (CCR) method of Kaczynski et al. with Pellikka's reduction techniques to efficiently compute homology and cohomology groups on tetrahedral meshes. It performs iterative local operations—elementary collapses (removing a k-simplex and its free (k-1)-face), interior face reductions (merging two k-simplices sharing a face), and general internal reductions (updating boundaries after removing a shared face among multiple neighbors)—wrapped in global loops: first a downward pass of pReduce(k) for k from dimension d to 1, followed by a pass interleaving pGeneralizedCombine(k) and pReduce(k-1). This preserves topological invariants while fully reducing the complex to minimize the subsequent Smith normal form computation for Betti numbers and generators, achieving complexity that the source paper does not state. **No bound is quoted here**: the Giard preprint
is in preparation and not in the local literature tree, and the implementation carries no
complexity annotation, so there is nothing to check a figure against. (Pellikka's own reduction
techniques are given as a *range*, O(n log n) to O(n²), applied cheapest-first —
`pellikka2013.txt:196-198`.) Ask Gregory Giard before quoting a bound for the generalized variant. For cohomology, dimensions are reversed (k from 0 to d-1), with an initial removal of a 0-simplex to initiate reductions.

The key insight is that the computed first cohomology generators H¹(K), essential for imposing currents in multiply connected domains, manifest as "thick cuts"—sparse sets of directed edges in the mesh with ±1 coefficients, defining oriented surfaces across which the scalar potential φ jumps by the current I (in Amperes). These edge-based representations, dual to homology cycles encircling conductors, enable straightforward enforcement of discontinuities in the FEM assembly while ensuring φ remains single-valued elsewhere, thus satisfying Ampere's law globally without conductor modeling. The algorithm's emphasis on collapses often yields visually interpretable thick cuts with fewer edges, aiding debugging, though post-processing can ensure binary coefficients if needed.

## FEM Discretization
> **Parity vs integer coefficients — flagged, not resolved.** The "odd number of times" phrasing
> below is the **ℤ₂** picture. The pipeline works over the integers: `CutData` carries
> `mCohomologyPlus`/`mCohomologyMinus` as signed membership, and `clean_spfa()` solves
> difference constraints to rectify representatives. The tree's own
> [`thin_cut_nonunit_rectification.md`](thin_cut_nonunit_rectification.md) exists precisely
> because **non-unit** coefficients occur — a notion with no meaning over ℤ₂ — and the debt
> register records that behavior as a mathematical property rather than a defect. So the
> statement below is at least narrower than what the code does. **Left as written pending Gregory
> Giard's ruling** (§7.1): the parity condition may be the intended reading for the
> single-imposition case even where the general pairing is integer-valued.

To discretize the non-conducting domain in the finite element method while enforcing the cohomology-derived discontinuity in φ, a key property is that every non-trivial homology cycle must intersect the cut an odd number of times, ensuring an odd parity in the additions and subtractions that manifest the jump condition. Although an extended finite element method (XFEM) could intuitively model this discontinuity through a "transition region" of enriched elements spanning the thick cut (the set of directed edges from the cohomology generator), a more efficient approach exploits the Poincaré-Lefschetz duality theorem stating that any thick cut can be equivalently reduced to a thin cut—a minimal surface where homology cycles cross oddly. This thin cut shifts the discontinuity precisely to element faces, implemented by duplicating nodes on one side of the cohomology surface (e.g., the positive side) and relinking the adjacent elements to these new nodes, thereby creating distinct degrees of freedom for φ⁺ and φ⁻. The jump [φ] = I is then imposed via static condensation: a transformation matrix T is constructed such that φ⁺ = T [φ⁻; I], and the local stiffness matrix undergoes a change of basis through Tᵀ K_local T to yield the condensed global system K_global, incorporating the current constraint without additional enrichment overhead.

For a concrete worked example of this push — thick cut, thin cut, conjugate edges, and DOF duplication on a single reference tetrahedron — see [thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md).

## Post-Processing Cleanup

> **Design note — partly prototyped, not wired in.** Of the manifold-filtering scheme described
> below, `manifold_filter_3d()` and `check_surface()` do **not** exist anywhere in the tree
> (searched `*.{cpp,hpp,f90,cmake}` across the repository with the ignore-file behavior disabled).
> **Tarjan articulation detection does exist, as an archived prototype**:
> `archive/graph/fn_Graph_tarjan.{hpp,cpp}` implements pocket detection by articulation point,
> with the size / cycle-density / compactness scoring the paragraph describes. `archive/` is
> outside the build — nothing in `CMakeLists.txt` references it — so none of it runs today.
> What the code does today is `Cohomology::remove_cut_pockets()` after SPFA rectification.
>
> The description is kept because it may be the design of record; it is marked so that no reader
> goes looking for the functions. **Its status is Gregory Giard's to confirm** (`doc/ai_collaboration_protocol.md` §7.1).

To ensure the computed cohomology cut forms a valid 2-manifold surface suitable for FEM assembly, a post-processing cleanup step is essential to resolve non-manifold artifacts, such as "double pockets"—small clusters of tetrahedral faces (typically 4-5) attached to the main cut via a single articulation face, leading to edges shared by more than two faces and causing deadlocks in iterative removal algorithms. Drawing from graph-theoretic approaches in the dual mesh graph (where faces are vertices and shared edges are graph edges), a hybrid detection strategy identifies and removes these pockets: Phase 1 employs Tarjan's DFS algorithm (O(V + E) complexity) to detect articulation points (neck faces) connecting unicyclic or near-unicyclic components (|E| ≈ |V|); Phase 2 uses region growing from suspicious areas to monitor genus changes and detect embedded pockets; Phase 3 validates via BFS layering, assessing metrics like cycle density (>0.3), compactness (<2.0), and size (<20 faces) before collective removal. This manifold filtering, integrated into functions like `manifold_filter_3d()` and enhanced with non-manifold edge pre-filtering from `check_surface()`, eliminates linear dependencies in cohomology constraints (as highlighted in the paper's numerical experiments), ensuring topological consistency and numerical stability in electromagnetic simulations without altering the global homology class.

## Implementation in BELFEM

The reduction algorithms and their exact code locations (Smith Normal Form, CCR, Pellikka,
Generalized Pellikka, BeltedTree) are cataloged in
[cohomology_algorithms.md](cohomology_algorithms.md); the main classes are summarized in the
[module index](@ref homology_index) Quick Reference.

The thick-to-thin cut conversion described above (Poincaré–Lefschetz duality, node duplication,
element relinking) is implemented in `cl_CutProcessor.cpp`, with static condensation performed
later in the FEM DOF manager. See
[thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md) for a worked
example with the full code map.

## Related Documentation

- [cohomology_algorithms.md](cohomology_algorithms.md) - Algorithmic overview with code locations
- [thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md) - Thick/thin cut worked example
- `../../fem/maxwell/` - Electromagnetic field solvers using these cuts
- `literature/papers/fem/index.md` - Research papers applying cohomology to HTS modeling
