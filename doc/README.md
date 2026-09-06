# BELFEM Documentation {#doc_index}

Index of general project documentation.

## Guidelines and Philosophy

- [documentation_guidelines.md](documentation_guidelines.md) - How to organize and name documentation files
- [coding_philosophy.md](coding_philosophy.md) - **BELFEM coding philosophy** (nomenclature, memory management, container selection, performance patterns)
- [literature_references.md](literature_references.md) - Reference guide to literature used in BELFEM development
- [ai_collaboration_protocol.md](ai_collaboration_protocol.md) - Two-AI review protocol, exchange format and confidence calibration
- [ai_workflow_best_practices.md](ai_workflow_best_practices.md) - **Experience report on the multi-AI method**: which practices were adopted, what changed after each, what they cost, and a five-item minimum adoption set — written to be readable outside this project
- [lessons_learned.md](lessons_learned.md) (repository only, not rendered on the site) - **The tripwire layer**: operating rules distilled from the incident catalog, indexed by the activity you are starting and by the symptom you are staring at. Layer 1 is what a session loads; Layer 2 carries one evidence card per rule. The file states its own rule and incident counts, and the reconciliation between clustered incidents and catalog rows, in its header
- [lessons_learned_evidence.md](lessons_learned_evidence.md) (repository only, not rendered on the site) - The cataloged incidents behind those rules, one row each, every `INC-NNN` citation resolving to the dated `devlog/` entry it was mined from — a reference appendix, not a document to read through

## User Reference

- [getting_started.md](getting_started.md) - **Build BELFEM and run your first simulation**: configure/build, the example workflow (gmsh mesh generation, then `belfem`), restart and material-database gotchas, deck validation with `belfem-conf`
- [input_file_reference.md](input_file_reference.md) - **Living reference of the `input.conf` contract**: syntax, every parsed section/key with types, units, defaults and parse-site citations, aliases, dead keys, and pitfalls. Extend it in the same session as any input-feature change.
- [mpi_support.md](mpi_support.md) - **Open MPI is the only supported MPI**: why MPICH and Intel MPI are refused at configure time, what the `-DALLOW_UNTESTED_MPI=ON` override costs, and why the Open MPI link flags in the MUMPS and MKL configs must not be "made portable"
- [parallel_execution.md](parallel_execution.md) — choosing MPI ranks, OpenMP threads and allocator settings for a production run: assembly scales with ranks (measured 3.6–3.8×), the factorization does not, and memory is set by a rank-independent factor

## Module-Specific Documentation

### Infrastructure

- [Core](../src/core/doc/README.md) - Fundamental utilities (logging, timing, types, constants, string tools)
- [Containers](../src/containers/doc/README.md) - BELFEM container classes usage guide
- [Linear Algebra](../src/linalg/doc/README.md) - Backend-agnostic linear algebra API and LAPACK wrappers
- [Sparse](../src/sparse/doc/README.md) - Sparse matrix storage and solver interfaces (UMFPACK, MUMPS, STRUMPACK, PARDISO, PETSc)
- [I/O](../src/io/doc/README.md) - File input/output (HDF5, ASCII, CSV, configuration files, XML)
- [Communication](../src/comm/doc/README.md) - MPI communication abstraction layer
- [Mesh](../src/mesh/doc/README.md) - Mesh data structures, I/O (Gmsh, HDF5, Exodus, VTK), and parallel partitioning
- [Visualizer](../src/visualizer/doc/README.md) - Optional VTK rendering of meshes and curves (`USE_VTK`, off by default)

### Mathematics

- [Graph](../src/math/graph/doc/README.md) - Graph algorithms and partitioning (BFS, DFS, RCM, METIS, SCOTCH)
- [Tensor](../src/math/tensor/doc/README.md) - Fourth-order tensor helper
- [Quaternion](../src/math/quaternion/doc/README.md) - Quaternion value type for 3D rotations
- [Spline](../src/numerics/spline/doc/README.md) - Spline interpolation
- [Optimizer](../src/numerics/opt/doc/README.md) - NLOPT-backed bound-constrained minimization (`USE_NLOPT`, default ON; not used by the FEM path)
- [Homology](../src/homology/doc/README.md) - Cohomology algorithms, cut generation and theory

### Finite Elements

- [FEM](../src/fem/doc/README.md) - Cross-cutting FEM notes, including the time-stepping strategy
- [FEM Kernel](../src/fem/kernel/doc/README.md) - DOF management, hanging nodes, assembly, time stepping and the nonlinear controller
- [Interpolation](../src/fem/interpolation/doc/README.md) - Shape functions, integration points, and Nédélec elements
- [IWG](../src/fem/iwg/doc/README.md) - Integral Weak Form physics module (Poisson, heat conduction, elasticity, time-stepping)
- [Maxwell](../src/fem/maxwell/doc/README.md) - Electromagnetic physics (h-φ formulation, thin shells, HTS, cohomology cuts)
- [Thermal](../src/fem/thermal/doc/README.md) - Transient heat conduction and the magneto-thermal coupling `belfem` selects for a deck with a thermal solver section
- [Postprocessing](../src/fem/postproc/doc/README.md) - Mesh-level postprocessing: volumes, surfaces, normals, gradients, mesh checks

### Physics and Circuits

- [Materials](../src/physics/materials/doc/README.md) - Material properties framework (pure metals, HTS, alloys, user-defined materials)
- [Database](../src/physics/database/doc/README.md) - Precomputed lookup tables on tensor grids (projection, interpolation, HDF5 persistence)
- [Gas Models](../src/physics/gasmodels/doc/README.md) - Equation-of-state and transport models
- [Gas Tables](../src/physics/gastables/doc/README.md) - Tabulated thermophysical property data
- [Circuit](../src/circuit/doc/README.md) - Lumped-element circuit simulator coupled to the FEM problem
- [Executables](../src/executables/doc/README.md) - Solver application (`belfem`) and the supporting tools (`material`, `gas`, `db2exo`, `msh2exo`), with their command-line interface

---

**See Also:**

- The generated API reference (`make doc`) covers classes, files and namespaces.
- `README.md` in the repository root - project overview and license.
