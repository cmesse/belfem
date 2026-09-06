# Cohomology Computation Algorithms in BELFEM {#homology_cohomology_algorithms}

## Overview

The homology module in BELFEM implements computational topology algorithms for computing cohomology groups, primarily designed for electromagnetic field simulations using finite element meshes.

## Core Algorithms

### 1. Smith Normal Form Algorithm
- **Location**: `fn_Smith.cpp:142-218`
- **Purpose**: Computes the Smith normal form of integer matrices to analyze chain complex structure
- **Method**: Transforms matrices into diagonal form revealing topological invariants and torsion coefficients
- **Output**: Matrices Q, Q⁻¹, R, R⁻¹ and ranks s, t such that QAR = D (diagonal)

### 2. Kernel-Image Decomposition
- **Location**: `fn_Smith.hpp:334-368`
- **Purpose**: Computes kernel (cycles) and image (boundaries) of boundary/coboundary operators
- **Method**: 
  - Row echelon form reduction
  - Matrix transposition techniques
- **Formula**: H^k = Ker(δ^k) / Im(δ^{k-1})

### 3. Quotient Group Computation
- **Location**: `cl_Cohomology.cpp:108-141`
- **Purpose**: Computes cohomology groups as quotient spaces
- **Method**:
  - Solves integer linear systems `W\V` to find generators
  - Extracts free and torsion parts via Smith decomposition
  - Identifies both infinite order (free) and finite order (torsion) generators

### 4. Pellikka Reduction Algorithm
- **Location**: `cl_SimplicialComplex.cpp:693-962`
- **Reference**: Pellikka et al., SIAM J. Sci. Comput., 2013
- **Purpose**: Reduces simplicial complex via topological operations
- **Key Operations**:
  - `pReduce()` (lines 693-755): Elementary collapses - removes k-simplex and its free (k-1)-face
  - `pCombine()` (lines 759-852): Interior face reductions - merges neighboring simplices
  - `reduceOmit()` (lines 937-962): Initialization pass before combine operations
- **Performance**: Efficiently reduces complex while preserving homology

### 5. Generalized Pellikka Algorithm (Giard Enhancement)
- **Location**: `cl_SimplicialComplex.cpp:994 (reduce) and :1291 (coreduce)`
- **Reference**: Giard et al., in preparation (drafted 2025, resumed 2026)
- **Purpose**: cohomology computation combining CCR with Pellikka techniques (the characterization is the authors'; no cost claim is made here)
- **Method**:
  - Downward pass: `pReduce(k)` for k from dimension d to 1
  - Second pass (again from d down to 1): interleaves `pGeneralizedCombine(k)` (lines 856-933) with `pReduce(k-1)`
  - For cohomology: Reversed dimensions (k from 0 to d-1) with initial 0-simplex removal
- **Complexity**: Pellikka et al. give a **range**, not a bound: the chain equivalences run from
  O(n log n) to O(n²) in the size of the chain-group basis being modified, and are applied in
  succession so that the cheapest do most of the work (`pellikka2013.txt:196-198`). The
  implementation carries no complexity annotation, so treat any single figure as indicative
- **Implementation**: Algorithm selection via `enum CutAlgorithm` in `en_CutAlgorithm.hpp`

### 6. Simplicial Complex Construction
- **Location**: `cl_SimplicialComplex.cpp`
- **Purpose**: Builds chain/cochain complexes from mesh topology
- **Features**:
  - Creates boundary/coboundary matrices from flagged mesh entities
  - Supports dimensions 0-3 (nodes, edges, faces, elements)
  - Handles both chains and cochains simultaneously
  - Data structures: `mChainsMap` and `mCochainsMap` storing k-chains/cochains

### 7. Generator Extraction and Cleaning
- **Location**: `cl_Cohomology.cpp` (`generatorsOfCohomology`, `clean_spfa`)
- **Purpose**: Identifies and cleans cohomology generators
- **Method**:
  - Extracts generators from last columns of U matrix after Smith decomposition
  - Cleans coefficients to ensure they are in {-1, 0, 1}
  - Removes illegal coefficients (|coeff| > 1) by subtracting node coboundaries.
    `clean()` dispatches to `clean_spfa()`: an SPFA difference-constraint solve
    proves a unit representative exists (or aborts with an edge-ID certificate
    naming the too-coarse loop), greedy sweeps rectify the feasible generator,
    and `remove_cut_pockets()` strips Tier-A pockets afterwards. See
    `thin_cut_nonunit_rectification.md` for theory and the certificate.

#### Cleaning runs twice on the factory path (by design)

On the production path (`CutFactory::compute_cohomologies()`, all cut
algorithms), `clean_spfa()` runs **twice** for each problem. Both runs are
required:

1. **Constructor run.** Every `Cohomology` constructor ends with `this->clean()`,
   so every constructed object holds unit-coefficient generators. Direct
   constructor users (the homology test fixtures, and any future non-factory
   caller) rely on this invariant. This run also raises the "no thin cut exists
   on this mesh" certificate for each raw Smith-form generator, before any
   recombination.
2. **Post-recombination run** (`cl_CutFactory.cpp`, after
   `updatekGeneratorsFromHomology()`). The update replaces the H¹ generators
   with integer linear combinations chosen so each cut pairs with one suggested
   homology generator. Adding unit cochains generally creates non-unit
   coefficients on shared edges, so the combined generators must be rectified
   again. This run is required for correctness: `CutData` stores coefficients as
   two ±1 bitsets and cannot represent |c| ≥ 2 (3D decks abort in
   `determine_cut_case_3d`; a 2D release build would silently mis-orient), and
   this is the only `remove_cut_pockets()` pass over the final generators.

Do not remove either run as an optimization. The first run would be redundant
for the final cohomology classes only if every suggested homology generator
were an exact cycle, a property the code never checks (suggested chains are
built with boundary tracking suppressed), and even then the concrete cut
representatives and the placement of the coarse-mesh certificate would change.
The cost of keeping both runs is one additional certify-and-rectify pass,
which is acceptable in practice.

### 8. Algorithm Selection
- **Location**: `en_CutAlgorithm.hpp`; the `switch( mAlgorithm )` in `CutFactory::compute_cohomologies()`
- **Available Algorithms**:
  - `Pellikka` (value 0): Original Pellikka algorithm
  - `CCR` (value 1): Chain complex reduction (Kaczynski)
  - `BeltedTree` (value 2): Spanning tree approach (alternative method)
  - `PellikkaGeneralized` (value 3): enhanced Giard version — **the production default**
    (`cl_MaxwellFactory.hpp:48`)

## Key Features

### Dual Computation
- Supports both homology (chains) and cohomology (cochains)
- Maintains duality between boundary and coboundary operators

### Physical Integration
- **Electromagnetic Applications**: Designed for Maxwell equation simulations
- **Terminal Handling**: Special support for:
  - Input/output terminals
  - Voltage conditions
  - Thin shell terminals
  - Conductor blocks and interfaces

### Mesh Integration
- Works directly with finite element meshes
- Supports multiple mesh formats (HDF5, Exodus, VTK)
- Handles both 2D and 3D geometries

### Advanced Features
- **Belted Tree Algorithm**: Alternative cohomology computation method (`cl_Cohomology.cpp:41-58`)
- **Orientation Detection**: Computes generator orientations relative to specified directions
- **Field Visualization**: Creates mesh fields to visualize cohomology generators

## Implementation Details

### Matrix Operations
- Integer matrix arithmetic for exact computations
- Specialized solvers for integer linear systems
- Efficient sparse matrix representations

### Parallel Support
- ~~MPI-aware implementations for distributed computing~~ — **not implemented.** Cohomology is computed serially on an undistributed mesh; `homology_usage_guide.md` states the same ("parallel cohomology is planned but not yet implemented")
- Scalable to large meshes and complex geometries

## References

The algorithms implemented in BELFEM follow standard computational topology approaches from the literature:

### Primary References

1. **T. Kaczynski, K. Mischaikow, and M. Mrozek**, "Computational Homology," Applied Mathematical Sciences, Springer, 2004.
   - **Referenced in code**: `fn_Smith.cpp:12-15`, `cl_SimplicialComplex.cpp:478,533,856`
   - **Algorithm**: Chain Complex Reduction (CCR) and Smith Normal Form

2. **M. Pellikka, S. Suuriniemi, L. Kettunen, and C. Geuzaine**, "Homology and cohomology computation in finite element modeling," *SIAM J. Sci. Comput.*, vol. 35, no. 5, pp. B1195–B1214, 2013. DOI: [10.1137/130906556](https://doi.org/10.1137/130906556)
   - **Referenced in code**: `cl_SimplicialComplex.cpp:693-962` (pReduce, pCombine, reduceOmit)
   - **Algorithm**: Pellikka reduction via elementary collapses and interior face reductions
   - **Widely cited**: Used in Sparselizard, GetDP, and various HTS electromagnetic modeling codes

3. **G. Giard et al.**, "Generalized Pellikka algorithm for cohomology computation," in preparation (drafted 2025, shelved, resumed 2026).
   - **Referenced in code**: `SimplicialComplex::reduce_complexPellikkaGeneralized()`
     (`cl_SimplicialComplex.cpp:994`) and `coreduce_complexPellikkaGeneralized()` (`:1291`)
   - **Algorithm**: a novel enhancement combining CCR with Pellikka's techniques. **No complexity
     figure is claimed here.** The preprint is in preparation and is not in the local literature
     tree, so the earlier "optimal O(n log n)" could not be checked against a source; the
     implementation carries no complexity annotation either. Ask Gregory Giard before quoting a
     bound for this variant
   - **Implementation**: algorithm selection switch in `CutFactory::compute_cohomologies()` (`cl_CutFactory.cpp:280`)

### Application Papers

The cohomology algorithms are applied in BELFEM for electromagnetic field simulations, as described in:

4. **V. Lahtinen, A. Stenvall, F. Sirois, and M. Pellikka**, "A finite element simulation tool for predicting hysteresis losses in superconductors using an H-oriented formulation with cohomology basis functions," *J. Supercond. Novel Magnetism*, vol. 28, no. 8, pp. 2345–2354, Aug. 2015.
   - **Concept**: Cohomology basis functions for transport current constraints

5. **B. d. S. Alves, V. Lahtinen, M. Laforest, and F. Sirois**, "3-D Finite-Element Thin-Shell Model for High-Temperature Superconducting Tapes," *IEEE Trans. Appl. Supercond.*, vol. 32, no. 6, 2022, Art. no. 6900310. DOI: [10.1109/TASC.2022.3165543](https://doi.org/10.1109/TASC.2022.3165543)
   - **Concept**: Thick cuts from relative homology H₂(Ω,Ωₛ); boundary of relative homology basis yields cohomology basis dual to resulting homology basis

6. **A. Riva, J. Dular, C. Geuzaine, and B. Vanderheyden**, "H-φ Formulation with Domain Decomposition for Magnetodynamic Simulations of HTS," *IEEE Trans. Appl. Supercond.*, vol. 33, no. 5, 2023, Art. no. 4900105. DOI: [10.1109/TASC.2023.3260666](https://doi.org/10.1109/TASC.2023.3260666)
   - **Concept**: Cohomology basis functions for multiply-connected domains in parallel computing framework

7. **N. Schnaubelt, M. Zhang, and M. Spenko**, "Electromagnetic Simulation of No-Insulation Coils Using H–φ Thin Shell Approximation," *IEEE Trans. Appl. Supercond.*, vol. 33, no. 5, 2023, Art. no. 4900605. DOI: [10.1109/TASC.2023.3260295](https://doi.org/10.1109/TASC.2023.3260295)
   - **Concept**: Automatic thick cut generation for multiply-connected no-insulation coil geometries

## Usage Context

These algorithms are primarily used in BELFEM for:
- Computing cuts in multiply-connected domains
- Ensuring gauge uniqueness in electromagnetic simulations
- Topological analysis of conductor geometries
- Handling periodic boundary conditions in Maxwell solvers