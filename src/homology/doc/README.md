# Homology Module Documentation {#homology_index}

**Module:** src/homology
**Purpose:** Index of documentation for BELFEM's computational topology module

---

> **AI edit ban — cohomology core.** `cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`,
> `cl_Chain`, `cl_Cochain` and `fn_Smith` are closed to AI edits (Christian and Gregory Giard,
> 2026-08-31). The reason is that this is a **modified** Pellikka reduction — unpublished work
> (Giard et al., in preparation) — so an AI's nearest prior is textbook Pellikka and every
> deliberate modification reads to it as a defect. The attested case: `coreduceOmit`
> (`cl_SimplicialComplex.cpp:1236-1257`) removes one 0-cochain per stall and re-runs the whole
> `pCoreduce` cascade each time; batching it was proposed and struck as unsafe. The loop looks
> trivially hoistable and is not. Reading, tracing and reporting stay open; source changes need
> Gregory's named authorization, and comment-only edits are unruled — treat them as banned.
> Full policy: `doc/ai_collaboration_protocol.md` §7.1. The rest of the module, including this
> documentation, follows ordinary edit-safety rules.

---

## Overview

The `homology` module provides computational topology algorithms for computing homology and cohomology groups, primarily designed for electromagnetic field simulations in multiply-connected domains using finite element meshes.

**Key capabilities:**
- Cohomology and homology computation for mesh topology analysis
- Automatic generation of topological cuts for multiply-connected domains
- Smith Normal Form decomposition for integer matrices
- Multiple reduction algorithms (Pellikka, CCR, generalized Pellikka)
- Thick-to-thin cut conversion for FEM assembly
- Support for thin shell structures and terminal handling

---

## Documentation Files

### User Guides

- **[homology_usage_guide.md](homology_usage_guide.md)** - Comprehensive usage guide for the homology module
  - Main classes and workflow
  - CutFactory pipeline
  - Cohomology and homology computation
  - Cut algorithm selection
  - Practical examples
  - Common patterns and pitfalls

### Theory and Algorithms

- **[cohomology_theory_and_implementation.md](cohomology_theory_and_implementation.md)** - Mathematical foundation
  - Physical and mathematical background
  - h-φ formulation for magnetostatics
  - Generalized Pellikka algorithm
  - FEM discretization techniques
  - Manifold cleanup post-processing

- **[cohomology_algorithms.md](cohomology_algorithms.md)** - Algorithmic reference
  - Smith Normal Form algorithm
  - Kernel-image decomposition
  - Pellikka reduction techniques
  - Generator extraction and cleaning
  - Literature references

- **[thick_thin_cuts_and_conjugate_edges.md](thick_thin_cuts_and_conjugate_edges.md)** - Worked examples
  - Thick cut, thin cut, and conjugate edges on a single reference tetrahedron
  - EXODUS local node ordering note
  - Poincaré–Lefschetz push to the opposite face and DOF duplication (p0,p1,p2 → q0,q1,q2)
  - **Cut-case convention** as implemented in `CutData::determine_cut_case_2d/3d()` (2D + 3D tables, signs, diagonal cases)
  - **2D hexagon worked example**: thick cut → cases → conjugated facets → duplication/relinking (band vs conjugated-facet membership) → abstract node and static condensation
  - Code map into `cl_CutProcessor.cpp` and FEM static condensation
  - Relevance to periodic boundary conditions (asymmetric seam duplication; conjugated facets in the seam)

- **[thin_cut_nonunit_rectification.md](thin_cut_nonunit_rectification.md)** - Non-unit thick-cut coefficients
  - What `|c(e)| ≥ 2` means geometrically (sheet multiplicity of the dual cut surface; CCT inter-turn throat)
  - Fixed-mesh existence obstruction `|⟨c, z⟩| ≤ length(z)` (manifold vs. mesh)
  - Rectification `c′ = c − dθ` as a **difference-constraint system** solved by **Bellman–Ford** (feasible → integer θ; infeasible → negative-cycle certificate localizing the throat)
  - Optimal representative via min-cost circulation (total unimodularity); normal-surface-theory identification of the cut cases
  - Source-field foundation: `c − dθ` ↔ curl-kernel freedom and minimal support (Dular et al. 1997); jump = transport current, binary `q_i` = unit cut (Dular et al. 1999)
  - Proposed `rectify_to_unit_or_certify()` algorithm and refinement / multi-level-duplication fallbacks

---

## Quick Reference

### Entry Point Class

| Class | File | Purpose |
|-------|------|---------|
| **`CutFactory`** | cl_CutFactory.{hpp,cpp} | **Main orchestrator** for cut generation |

### Core Computational Classes

| Class | File | Purpose |
|-------|------|---------|
| `Cohomology` | cl_Cohomology.{hpp,cpp} | Computes cohomology groups H^k |
| `Homology` | cl_Homology.{hpp,cpp} | Computes homology groups H_k |
| `SimplicialComplex` | cl_SimplicialComplex.{hpp,cpp} | Simplicial complex reduction engine |
| `CutProcessor` | cl_CutProcessor.{hpp,cpp} | Converts thick cuts to thin cuts for FEM |

### Data Structures

| Class | File | Purpose |
|-------|------|---------|
| `Chain` | cl_Chain.{hpp,cpp} | Formal sums of k-simplices (homology) |
| `Cochain` | cl_Cochain.{hpp,cpp} | Cochains (cohomology) |
| `CutData` | cl_CutData.{hpp,cpp} | Per-cut topology and metadata |
| `Topology` | cl_Topology.{hpp,cpp} | Mesh topology analysis |
| `Protoshell` | ../mesh/cl_Protoshell.hpp (header-only, namespace `belfem`) | Thin shell configuration |

### Algorithms

| Function/File | Purpose |
|---------------|---------|
| `smithForm()` (`fn_Smith.{hpp,cpp}`) | Smith Normal Form for integer matrices |
| `reduce_complexPellikka()` | Pellikka reduction algorithm |
| `reduce_complexCCR()` | Chain Complex Reduction (Kaczynski) |
| `reduce_complexPellikkaGeneralized()` | Enhanced Giard algorithm — the production default (`cl_MaxwellFactory.hpp:48`) |

---

## Typical Usage Pattern

```cpp
// 1. Prepare mesh and topology
mesh::Topology tTopology(aMesh);
tTopology.run();

// 2. Create CutFactory
mesh::CutFactory tFactory(
    aMesh,                          // Mesh to process
    &tTopology,                     // Topology analysis
    mProtoshells,                   // Thin shell configurations
    mesh::CutAlgorithm::PellikkaGeneralized,  // Algorithm choice
    false                           // Use enrichment
);

// 3. Set up terminals (if needed)
tFactory.set_terminals(mTerminals, mThinShellIndices);

// 4. Run pipeline
tFactory.run();

// 5. Retrieve results
Cell<mesh::Node*>& abstractNodes = tFactory.abstract_nodes();
Cell<mesh::SideSet*>& cuts = tFactory.cuts();

// 6. Use in FEM assembly
for (mesh::SideSet* cut : cuts) {
    cut->set_domain_type(mesh::DomainType::Cut);
    // Apply discontinuity constraints in FEM system
}
```

---

## Cut Algorithms

Four reduction algorithms available via `CutAlgorithm` enum:

| Algorithm | Value | Performance | Use Case |
|-----------|-------|-------------|----------|
| **PellikkaGeneralized** | 3 | no measured comparison recorded | **The production default** (`cl_MaxwellFactory.hpp:48`) |
| Pellikka | 0 | Good | Original algorithm, widely tested |
| CCR | 1 | Moderate | Classical Kaczynski method |
| BeltedTree | 2 | Alternative | Spanning tree approach |

---

## Physical Context

The homology module is designed for **electromagnetic field simulations** in multiply-connected domains, specifically:

- **Magnetostatics** using h-φ formulation
- **Transport current constraints** in superconducting magnets
- **Thin shell modeling** for layered conductors
- **Multiply-connected geometries** (toroids, solenoids, cables)

### Key Concept: Cohomology Cuts

In multiply-connected domains (e.g., a torus), the magnetic scalar potential φ must satisfy Ampere's circuital law:

```
∮_C H · dl = I
```

For any closed loop C encircling a conductor carrying current I. This requires introducing a **cohomology cut**—a surface across the domain where φ experiences a discontinuity [φ] = φ⁺ - φ⁻ = I.

**Workflow:**
1. Compute cohomology generators (thick cuts = sets of directed edges)
2. Convert to thin cuts (minimal surfaces for FEM assembly)
3. Duplicate nodes on cut surface to create distinct DOFs for φ⁺ and φ⁻
4. Enforce jump condition via static condensation in FEM system

---

## Source Code

**Module location:** `../../`

**Key source files:**
- Entry point: `cl_CutFactory.{hpp,cpp}`
- Computation: `cl_Cohomology.{hpp,cpp}`, `cl_SimplicialComplex.{hpp,cpp}`
- Algorithms: `fn_Smith.{hpp,cpp}`, reduction methods in SimplicialComplex
- Data structures: `cl_Chain.{hpp,cpp}`, `cl_Cochain.{hpp,cpp}`
- Post-processing: `cl_CutProcessor.{hpp,cpp}`, manifold filtering

---

## External References

> The canonical citation list for the whole project is
> [`doc/literature_references.md`](../../../doc/literature_references.md) — check there first,
> and add new citations there rather than here. This module keeps its own list because the
> algorithms below are the module's subject matter and the mapping from algorithm to paper is
> part of understanding the code; where the two disagree, the project list is authoritative.

### Primary Literature

1. **T. Kaczynski, K. Mischaikow, and M. Mrozek**, "Computational Homology," Springer, 2004.
   - Chain Complex Reduction (CCR) method

2. **M. Pellikka et al.**, "Homology and cohomology computation in finite element modeling," *SIAM J. Sci. Comput.*, 2013. [DOI: 10.1137/130906556](https://doi.org/10.1137/130906556)
   - Pellikka reduction algorithm

3. **G. Giard et al.**, "Generalized Pellikka algorithm," in preparation (drafted 2025, resumed 2026).
   - Enhanced algorithm combining CCR with Pellikka (implemented in BELFEM)

### Application Papers

4. **V. Lahtinen et al.**, "Cohomology basis functions for H-oriented formulation," *J. Supercond. Novel Magnetism*, 2015.
5. **B. Alves et al.**, "3D Thin-Shell Model for HTS Tapes," *IEEE Trans. Appl. Supercond.*, 2022.
6. **A. Riva et al.**, "H-φ with Domain Decomposition for HTS," *IEEE Trans. Appl. Supercond.*, 2023.
7. **N. Schnaubelt et al.**, "No-Insulation Coils Using H–φ Thin Shell," *IEEE Trans. Appl. Supercond.*, 2023.

See `cohomology_algorithms.md` and `cohomology_theory_and_implementation.md` for detailed references.

---

## Related BELFEM Modules

- **Mesh** (`src/mesh/`): Provides Mesh, Node, Element, SideSet classes used by homology
- **FEM Maxwell** (`src/fem/maxwell/`): Electromagnetic solvers that consume cohomology cuts
- **Linear Algebra** (`src/linalg/`): Matrix operations for boundary/coboundary matrices

---

## Development Notes

### Adding New Cut Algorithms

1. Add enum value to `en_CutAlgorithm.hpp`
2. Implement reduction method in `cl_SimplicialComplex.cpp`
3. Add case to the `switch( mAlgorithm )` in `CutFactory::compute_cohomologies()` (`cl_CutFactory.cpp`, currently `:280`)
4. Benchmark against existing algorithms and document performance characteristics

### Debugging Cut Generation

**Tip:** Debug builds write one `cut_<index>.vtk` per cut automatically
(`CutFactory::compute_thin_cuts_and_duplicate_interface_nodes()`, `cl_CutFactory.cpp:477-479`,
compiled out under `NDEBUG`); release builds do not. The call behind it is public:

```cpp
CutProcessor* processor = ...;
processor->save_debug_meshes();  // one mesh per cut
```

Output files:
- `cut_<index>.vtk` - one debug mesh per cut (`cl_CutProcessor.cpp:132-140`)
- `curve_<id>.vtk` - written by the public `CutFactory::save_curve_debug_meshes()`
  (`cl_CutFactory.hpp:184`)

Earlier revisions of this file listed `thick_cut_*.vtk`, `thin_cut_*.vtk` and `manifold_*.vtk`.
No code emits those names; `manifold` does not appear anywhere in the module's sources.

### Common Pitfalls

- **Non-manifold cuts**: post-processing cleanup may fail for complex geometries. Note that `check_surface()` does **not** exist in the tree — see the design note in `cohomology_theory_and_implementation.md`. The implemented cleanup is `Cohomology::remove_cut_pockets()` after SPFA rectification.
- **Algorithm selection**: `PellikkaGeneralized` is the production default (`cl_MaxwellFactory.hpp:48`). No timing comparison against the other three is recorded here, so treat "fastest" claims as unmeasured.
- **Terminal orientation**: Ensure terminals are oriented consistently to avoid current sign errors.

---

## See Also

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`
- **General Documentation**: `../../../doc/README.md`
