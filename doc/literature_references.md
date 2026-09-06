# BELFEM Literature References {#doc_literature_references}

**Purpose:** Reference guide to the literature used in BELFEM development.

**Note:** The `./literature/` directory contains proprietary content and is maintained as an independent repository (in `.gitignore`). This document provides references so users can locate the cited works independently.

**Last Updated:** 2026-08-11

---

## Paper Alias Index (retired — decoder for older records)

**The `paperN` aliases were retired on 2026-08-11.** New citations use author and year
("Messe et al. 2023, §2.7"); see `CLAUDE.md`. They were converted out of the source comments and
the module documentation in the same pass, so the aliases now survive only in dated devlog entries,
which are kept as written.

This table remains the authoritative alias → source file → citation map, so those older records stay
readable and auditable. Source `.txt` files live in the proprietary `literature/` repository; the
citations below let anyone locate the works independently.

| Alias | Source file | Citation |
|---|---|---|
| paper0 | `dular2021.txt` | Dular et al. 2021, "Stability of Mixed FE Formulations for HTS", IEEE TASC, DOI:10.1109/TASC.2021.3098724 |
| paper1 | `messe2023.txt` | Messe et al. 2023, "BELFEM: a special purpose FE code for magnetodynamic modeling of HTS tapes", SUST, DOI:10.1088/1361-6668/acf7f9 |
| paper2 | `messe2022.txt` | Messe 2022, "A Special Purpose Finite-Element Framework for HTS Applications", HTS 2022 Conference, HAL:hal-03791404 |
| paper3 | `arsenault2023.txt` | Arsenault et al. 2023, "Magnetodynamic H-φ Formulation", IEEE TASC, DOI:10.1109/TASC.2023.3293449 |
| paper4 | `riva2023.txt` | Riva et al. 2023, "H-φ in Sparselizard with DDM", IEEE TASC, DOI:10.1109/TASC.2023.3240389 |
| paper5 | `alves2022a.txt` | Alves et al. 2022a, "3D FE Thin-Shell Model for HTS Tapes", IEEE TASC, DOI:10.1109/TASC.2022.3143076 |
| paper6 | `alves2022b.txt` | Alves et al. 2022b, "Thin-shell approach for modeling superconducting tapes in H-φ formulation", SUST, DOI:10.1088/1361-6668/ac3f9e |
| paper7 | `arsenault2021.txt` | Arsenault et al. 2021, "Implementation of H-φ Formulation in COMSOL Multiphysics", IEEE TASC, DOI:10.1109/TASC.2020.3033998 |
| paper8 | `alves2024.txt` | Alves et al. 2024, "2-D Thin-Shell Model Based on H-φ-Formulation in COMSOL", IEEE TASC, DOI:10.1109/TASC.2024.3473850 |
| paper9 | `schnaubelt2023.txt` | Schnaubelt et al. 2023, "Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA", IEEE TASC, DOI:10.1109/TASC.2023.3258905 |
| paperA | `schnaubelt2023.txt` | Same work as paper9 — July 2026 devlogs use `paperA` for it (alias drift); both resolve here |

Newer papers cited by author-year (no `paperN` alias assigned):

| Source file | Citation |
|---|---|
| `arsenault2026.txt` | Arsenault et al. 2026, Erratum to "Magnetodynamic H-φ Formulation" (corrects Arsenault et al. 2023's air-domain coupling), IEEE TASC, DOI:10.1109/TASC.2026.3686487 |
| `denis2026.txt` | Denis et al. 2026, "Simultaneous Multi-Scale Homogeneous H-Phi Thin-Shell Model for Stacked HTS Coils", IEEE TASC, DOI:10.1109/TASC.2026.3652981 |
| `dular1997.txt` | Dular et al. 1997, "A Generalized Source Magnetic Field Calculation Method for Inductors of Any Shape", IEEE Trans. Magn., DOI:10.1109/20.582518 |
| `dular1999.txt` | Dular et al. 1999, "A Natural Method for Coupling Magnetodynamic H-Formulations and Circuit Equations", IEEE Trans. Magn., DOI:10.1109/20.767308 |
| `luccini2025.txt` | Lucchini 2025, "Evaluating Magnetization Losses in 3-D CORC Tapes With Integral and Finite-Element Methods", IEEE TASC, DOI:10.1109/TASC.2025.3544512 (filename `luccini` — cite as Lucchini) |
| `wozniak2025.txt` | Wozniak et al. 2025, "Influence of Critical Current Defect on Operation, Quench Detection and Protection of a Conduction-Cooled Pancake REBCO Coil", IEEE TASC 35(5):4604006, DOI:10.1109/TASC.2025.3532246 |
| `badel2021.txt` | Badel et al. 2019, "Modeling of 'quench' or the occurrence and propagation of dissipative zones in REBCO high temperature superconducting coils", SUST 32(9):094001, DOI:10.1088/1361-6668/ab181f (filename `badel2021` is the HAL deposit year, HAL:hal-02509494 — cite as Badel et al. 2019) |

FVM papers (aliases F0–F4):

| Alias | Source file | Citation |
|---|---|---|
| F0 | `aavatsmark2002.txt` | Aavatsmark 2002, "An Introduction to Multipoint Flux Approximations for Quadrilateral Grids", Comput. Geosci. 6:405-432, DOI:10.1023/A:1021291114475 |
| F1 | `agelas2008.txt` | Agélas, Di Pietro, Masson 2008, "A Symmetric and Coercive Finite Volume Scheme for Multiphase Porous Media Flow", FVCA V, DOI:10.1515/JNUM.2008.006 |
| F2 | `klausen2006.txt` | Klausen & Winther 2006, "Robust Convergence of Multi Point Flux Approximation on Rough Grids", Numer. Math. 104:317-337, DOI:10.1007/s00211-006-0023-4 |
| F3 | `ingram2010.txt` | Ingram, Wheeler, Yotov 2010, "A Multipoint Flux Mixed Finite Element Method on Hexahedra", SIAM J. Numer. Anal. 48(4):1281-1312, DOI:10.1137/090766176 |
| F4 | `wheeler2011.txt` | Wheeler, Xue, Yotov 2011, "A Family of Multipoint Flux Mixed Finite Element Methods for Elliptic Problems on General Grids", Procedia CS 4:918-927, DOI:10.1016/j.procs.2011.04.097 |

---

## Overview

BELFEM development is informed by four complementary literature categories:

1. **FEM Fundamentals** — Authoritative textbooks and lecture notes on finite element theory, implementation, and numerical methods
2. **BELFEM FEM Specialization** — Research papers specific to HTS modeling, h-φ formulations, and thin-shell approximations
3. **BELFEM FVM Methods** — Research papers on finite volume methods, specifically MPFA for diffusion problems
4. **Computational Topology and Optimal Cuts** — Algebraic-topology foundations and optimization theory behind cut generation in multiply-connected domains (the `homology/` module)

---

## Quick Reference: Which Literature Do I Need?

### For BELFEM-Specific Questions

Start with **research papers** if your question relates to:
- h-φ formulation and interface coupling
- Thin-shell models for HTS tapes
- Transport current constraints and cuts (thin/thick)
- Cohomology and multiply-connected domains
- HTS material models (E-J power law, Kim model)
- Validation benchmarks for superconductor applications

### For BELFEM FVM Questions

Start with **FVM research papers** if your question relates to:
- Multipoint Flux Approximation (MPFA) methods
- Anisotropic diffusion on general grids
- Transmissibility matrix computation
- Physical vs reference space evaluation
- Enhanced BDDF₁ spaces for 3D hexahedra
- Coercivity and convergence conditions

### For Cuts, Cohomology, and Topology Questions

Start with the **topology references** if your question relates to:
- Why cuts exist and what they are topologically
- How a cohomology basis is computed on an FE mesh (BELFEM's thick cuts)
- Whether a *minimal* cut is computable, and at what cost
- Reduction/coreduction preprocessing before Smith normal form
- The graph and LP algorithms underneath any of the above

### For General FEM Theory

Start with **textbooks** if your question relates to:
- Weak formulations and Galerkin methods
- Element technology (isoparametric, mixed, reduced integration)
- Locking phenomena (volumetric, shear, membrane)
- Inf-sup stability and mixed formulations
- Solver strategies (direct, iterative, eigenvalue)
- Nonlinear solvers (Newton-Raphson, arc-length)
- Time integration and dynamic analysis

---

## FEM Textbooks

### Core References

**Bathe, K.-J.** *Finite Element Procedures*, 2nd ed., 4th printing. Prentice Hall, 2016.
- **Use for:** Robustness, reliability, production code implementation
- **Key chapters:** Ch. 4 (formulation), Ch. 5 (elements), Ch. 8 (solvers), Ch. 9 (time integration)
- **Personality:** Disciplined engineer — emphasizes proven, robust methods

**Hughes, T.J.R.** *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*. Dover, 2000 (reprint of 1987 original).
- **Use for:** Clear linear FEM implementations, algorithm design
- **Key chapters:** Ch. 1 (weak forms), Ch. 4 (locking), Ch. 5 (isoparametric elements)
- **Personality:** Algorithm designer — clean, systematic approach

**Zienkiewicz, O.C., Taylor, R.L., and Zhu, J.Z.** *The Finite Element Method: Its Basis and Fundamentals*, 7th ed. Butterworth-Heinemann, 2013.
- **Use for:** Comprehensive reference, adaptivity
- **Personality:** Encyclopedia — covers everything

**Zienkiewicz, O.C., and Taylor, R.L.** *The Finite Element Method for Solid and Structural Mechanics*, 7th ed. Butterworth-Heinemann, 2014.
- **Use for:** Nonlinear implementation, exact shell elements
- **Key chapters:** Ch. 5 (nonlinear theory), Ch. 8 (contact)

### Advanced References

**Belytschko, T., Liu, W.K., Moran, B., and Elkhodary, K.** *Nonlinear Finite Elements for Continua and Structures*, 2nd ed. Wiley, 2014.
- **Use for:** Nonlinear mechanics, XFEM, stability theory
- **Not for:** First-time learning (dense presentation)

**Brenner, S.C., and Scott, L.R.** *The Mathematical Theory of Finite Element Methods*, 3rd ed. Springer, 2008.
- **Use for:** Rigorous proofs, convergence theory
- **Key chapters:** Ch. 4-5 (convergence), Ch. 8 (mixed methods)

**Bronshtein, I.N., Semendyayev, K.A., Musiol, G., and Mühlig, H.** *Handbook of Mathematics*, 6th ed. Springer, 2015.
- **Use for:** Quick formula verification and lookup

### Lecture Notes

**Evans, J.A.** *Mathematical Foundations of the Finite Element Method*, Lecture Notes, University of Colorado Boulder, 2017.
- **Use for:** The functional-analysis layer under everything else — Sobolev spaces, distributions, Lax-Milgram, Babuška-Brezzi, interpolation error estimates
- **Cite as:** "Evans 2017, Ch. N" or "Evans 2017, §N.M"
- **Pairs with:** Brenner & Scott (proofs in book form), Boffi et al. (mixed extensions)

**Felippa, C.A.** *Introduction to Finite Element Methods (IFEM)*, Course Notes, University of Colorado Boulder.
- **Use for:** Step-by-step learning, clear derivations; the penalty-method and Lagrange-multiplier tutorials, and static condensation
- **Cite as:** "Felippa IFEM, Ch. N, §X.Y"
- **Available:** Online as course materials

### Specialized Advanced References

**Arnold, D.N., Falk, R.S., and Winther, R.** *Finite Element Exterior Calculus*. SIAM, 2018.
- **Use for:** FEEC theory, topological foundations of FEM spaces
- **Key topics:** de Rham complex, structure preservation, why FEM spaces must exist
- **Personality:** Topologist — explains the mathematical structure underlying FEM

**Boffi, D., Brezzi, F., and Fortin, M.** *Mixed Finite Element Methods and Applications*. Springer, 2013.
- **Use for:** Mixed FEM theory, saddle-point problems, inf-sup stability
- **Key chapters:** Ch. 3-5 (inf-sup conditions), Ch. 11 (Maxwell eigenvalues)
- **Personality:** Mixed FEM theorist — comprehensive treatment of stability

**Monk, P.** *Finite Element Methods for Maxwell's Equations*. Oxford University Press, 2003.
- **Use for:** Maxwell equations, edge elements (Nédélec), electromagnetic scattering
- **Key chapters:** Ch. 5-6 (edge elements), Ch. 7 (convergence), Ch. 9-12 (scattering)
- **Personality:** Maxwell specialist — definitive source for electromagnetics FEM

### Element Distortion on Quadrilateral/Hexahedral Meshes

Why QUAD/HEX elements in the magnetic solve must be perfectly rectangular: see
`src/fem/interpolation/doc/nedelec.md` §6.6 and the maxwell usage guide §1.7. On
non-affine (bi-/trilinearly mapped) quads and hexes, the mapped H(curl)/H(div) spaces lose
completeness, and the lowest-order elements can lose convergence entirely. Monk 2003
develops the hexahedral edge-element theory only for parallelepipeds (§6.1) and cites this
line of work as the reason (§8.2–8.3). Boffi et al. 2013 summarize it in §2.2.4/§2.5.5.

**Arnold, D.N., Boffi, D., Falk, R.S., and Gastaldi, L.** "Finite element approximation on
quadrilateral meshes". *Communications in Numerical Methods in Engineering*, 17(11):805-812,
2001. DOI:10.1002/cnm.450
- Necessary and sufficient conditions for approximation order under bilinear maps

**Arnold, D.N., Boffi, D., and Falk, R.S.** "Approximation by quadrilateral finite
elements". *Mathematics of Computation*, 71(239):909-922, 2002.
DOI:10.1090/S0025-5718-02-01439-4
- Scalar theory showing that distorted quads already break serendipity elements

**Arnold, D.N., Boffi, D., and Falk, R.S.** "Quadrilateral H(div) finite elements".
*SIAM Journal on Numerical Analysis*, 42(6):2429-2451, 2005. DOI:10.1137/S0036142903431924
- The vector-valued case: RT[0] divergence does not converge on general quad meshes

**Falk, R.S., Gatto, P., and Monk, P.** "Hexahedral H(div) and H(curl) finite elements".
*ESAIM: Mathematical Modelling and Numerical Analysis*, 45(1):115-143, 2011.
DOI:10.1051/m2an/2010034
- Extends the negative results to 3D hexahedral H(curl), the case that bites h-φ directly

---

## BELFEM Research Papers

### Tier 1: Core BELFEM

**Messe et al. 2023** — "BELFEM: a special purpose FE code for magnetodynamic modeling of HTS tapes" (primary reference)
- BELFEM architecture, h-φ formulation, thin-shell implementation
- Static condensation preferred over Lagrange multipliers
- Solver strategies: STRUMPACK vs MUMPS
- Nonlinear iteration: hybrid Picard/Quasi-Newton/Newton
- Checkerboarding prevention: tight tolerances (ε < 10⁻¹¹) — the paper's recommendation, tightened per deck; BELFEM's default `tolerance` is `1e-6`

**Messe 2022** — "A Special Purpose Finite-Element Framework for HTS Applications"
- Project motivation and scope

### Tier 2: h-φ Thin-Shell Focus

**Alves et al. 2022b** — "Thin-shell approach for modeling superconducting tapes in H-φ formulation"
- Full derivation with interface conditions
- Cohomology basis theory (thick cuts)
- N=1 vs N>1 analysis
- **Key finding:** N>1 needed for stacked tapes to capture top/bottom losses

**Alves et al. 2024** — "2-D Thin-Shell Model Based on H-φ-Formulation in COMSOL"
- Thin cut implementation
- Alternative approach to Alves et al. 2022b

**Schnaubelt et al. 2023** — "Electromagnetic Simulation of No-Insulation Coils Using H-φ TSA"
- Practical interpretation for coils
- Transport current handling

### Tier 3: Related Methods

**Arsenault et al. 2023** — "Magnetodynamic H-φ Formulation"
- **Recommended variant:** Magnetodynamic (H-φ/D) over quasi-static
- Interface coupling methodology
- Includes ∂(μh)/∂t term for correct transients

**Alves et al. 2022a** — "3D FE Thin-Shell Model for HTS Tapes"
- Adaptable to h-φ
- Relative homology background

**Arsenault et al. 2021** — "Implementation of H-φ Formulation in COMSOL Multiphysics"
- Performance comparison with pure H-formulation
- Element order recommendations (quadratic for 2D)

### Tier 4: Supporting Theory

**Dular et al. 2021** — "Stability of Mixed FE Formulations for HTS"
- Inf-sup conditions for h-φ
- Hierarchical basis functions
- HTS-specific stability analysis

**Riva et al. 2023** — "H-φ in Sparselizard with DDM"
- Domain decomposition methods
- **Note:** Different software (NOT BELFEM)

### Tier 5: HTS Material Models (E-J laws)

**Rhyner 1993** — "Magnetic properties and AC-losses of superconductors with power-law
current-voltage characteristics"
- Physica C 212, pp. 292–300
- The E-J power law underlying all three BELFEM resistivity laws
- Cited in `src/physics/materials/powerlaws.hpp`

**Plummer & Evetts 1987** — "Dependence of the shape of the resistive transition on composite
inhomogeneity in multifilamentary wires"
- IEEE Trans. Magn. 23 (2), pp. 1179–1182
- Phenomenology of the n-value characterizing the transition sharpness

**Duron et al. 2004** — "Modelling the E–J relation of high-Tc superconductors in an arbitrary
current range"
- J. Duron, F. Grilli, B. Dutoit, S. Stavrev — Physica C 401 (1-4), pp. 231–235
- The parallel combination ρ = ρ_PL·ρ_n/(ρ_PL+ρ_n) used by the `powerlaw` and `riva` laws
- **DOI:** 10.1016/j.physc.2003.09.044

**Riva 2021** — "Quench behavior of high-temperature superconductor tapes for power applications:
a strategy toward resilience"
- N. Riva, EPFL doctoral thesis no. 8754 (2021)
- Measured overcritical-regime resistivity of REBCO tapes 77–90 K (Ch. 3–5); the ρ_ηβ
  overcritical model and its continuity analysis (§5.1); impact of the law choice on simulated
  quench speed (Ch. 6)
- Namesake of the `resistivity type : riva` law (which implements the Duron parallel form the
  thesis builds on, not the fitted ρ_ηβ model itself)
- Local copy (ephemeral): `./tmp/EPFL_TH8754.pdf`

See `src/physics/materials/doc/resistivity_laws.md` for how these map onto the three laws.

---

## BELFEM FVM Papers

### Primary References

**Aavatsmark 2002** — "An Introduction to Multipoint Flux Approximations for Quadrilateral Grids"
- **DOI:** 10.1023/A:1021291114475
- **Use for:** MPFA O-method fundamentals, transmissibility matrices (start here)
- **Key topics:** Interaction volumes, T = CA⁻¹B - D, K-orthogonality, monotonicity

**Ingram et al. 2010** — "A Multipoint Flux Mixed Finite Element Method on Hexahedra"
- **DOI:** 10.1137/090766176
- **Use for:** Enhanced BDDF₁ for 3D hexahedra (critical for implementation)
- **Key topics:** 4 DOF/face, curl enrichment, trapezoidal quadrature, superconvergence

### Convergence and Stability

**Agélas et al. 2008** — "A Symmetric and Coercive Finite Volume Scheme"
- **DOI:** 10.1515/JNUM.2008.006
- **Use for:** Coercivity conditions, when your mesh will converge
- **Key topics:** Local coercivity coer(D,Λ), general meshes, L∞ coefficients

**Klausen & Winther 2006** — "Robust Convergence of Multi Point Flux Approximation on Rough Grids"
- **DOI:** 10.1007/s00211-006-0023-4
- **Use for:** Physical vs reference space, robust convergence on rough grids
- **Key topics:** Jacobian at cell center, stability condition (27), broken RT

### Theory and Connection

**Wheeler et al. 2011** — "A Family of Multipoint Flux Mixed Finite Element Methods"
- **DOI:** 10.1016/j.procs.2011.04.097
- **Use for:** MFMFE overview, variational framework
- **Key topics:** BDM₁/BDDF₁ with quadrature, FEM-FVM connection, h²-perturbations

---

## Computational Topology and Optimal Cuts

The references behind the `homology/` module: why cuts exist, how BELFEM computes them, and what is
and is not tractable if you want them minimal.

### Theory and Method

**Gross, P.W., and Kotiuga, P.R.** *Electromagnetic Theory and Computation: A Topological Approach*. MSRI Publications, vol. 48, Cambridge University Press, 2004. ISBN 0-521-80160-5.
- **Use for:** Why cuts must exist at all — cohomology for electromagnetics, branch cuts for the magnetic scalar potential, duality theorems
- **Cite as:** "Gross & Kotiuga, §X" or "Gross & Kotiuga, Ch. N"
- **Answers:** existence and construction in principle; not minimality

**Pellikka, M., Suuriniemi, S., Kettunen, L., and Geuzaine, C.** "Homology and Cohomology Computation in Finite Element Modeling", *SIAM Journal on Scientific Computing*, 2013.
- **DOI:** 10.1137/130906556
- **Use for:** **BELFEM's actual method** — chain complex → reduction → Smith normal form → thick cuts (the Gmsh solver)
- **Key topics:** cohomology basis on an FE mesh, thick cuts vs thin cuts

**Giard, G., et al.** "Generalized Pellikka algorithm for cohomology computation", **in preparation** (drafted 2025, shelved, resumed 2026).
- **DOI:** none — unpublished draft, not in `./literature/`
- **Use for:** the algorithm BELFEM actually runs by default — `reduce_complexPellikkaGeneralized` / `coreduce_complexPellikkaGeneralized`, selected at `cl_MaxwellFactory.hpp:48`
- **Key topics:** generalized Pellikka reduction and coreduction; the cohomology-side 0-cochain omit loop
- **Note:** this is the reason the cohomology core is closed to AI edits (`doc/ai_collaboration_protocol.md` §7.1). A modified algorithm with no published implementation gives a model no correct prior — its nearest neighbor is textbook Pellikka, so deliberate departures read as defects. Ask Gregory for the draft; do not reconstruct the method from the code.
- **Author's name is Giard.** The module documentation carried "Giarda" in nine places until 2026-08-31; if it reappears, it is a typo.

**Mrozek, M., and Batko, B.** "Coreduction Homology Algorithm", *Discrete & Computational Geometry*, 2009.
- **DOI:** 10.1007/s00454-008-9073-y
- **Use for:** Why the reduction step is fast — linear-time coreduction collapses the complex before Smith normal form

### Optimality and Complexity

**Dey, T.K., Hirani, A.N., and Krishnamoorthy, B.** "Optimal Homologous Cycles, Total Unimodularity, and Linear Programming", *SIAM Journal on Computing*, 2011.
- **DOI:** 10.1137/100800245
- **Use for:** When a *minimal* cut is computable in polynomial time — over **Z**, with a totally unimodular boundary matrix, via LP

**Chen, C., and Freedman, D.** "Hardness Results for Homology Localization", *Discrete & Computational Geometry*, 2011 (also Proc. ACM-SIAM SODA 2010).
- **DOI:** 10.1007/s00454-010-9322-8
- **Use for:** The negative result — over **Z₂** the problem is NP-hard to approximate within any constant factor. Read before proposing a "just minimize it" cut algorithm

**Dunfield, N.M., and Hirani, A.N.** "The Least Spanning Area of a Knot and the Optimal Bounding Chain Problem", *Proc. 27th Annual Symposium on Computational Geometry (SoCG '11)*, 2011.
- **DOI:** 10.1145/1998196.1998218
- **Use for:** Minimal spanning surfaces — NP-complete in general, polynomial when H₂ = 0

**Costantini, M.** "A Novel Phase Unwrapping Method Based on Network Programming", *IEEE Transactions on Geoscience and Remote Sensing*, 1998.
- **DOI:** 10.1109/36.673674
- **Use for:** The concrete grid analogue — branch cuts as min-cost network flow, i.e. a TU-structured LP solved at scale

**Haken, W.** "Theorie der Normalflächen: Ein Isotopiekriterium für den Kreisknoten", *Acta Mathematica*, 1961 (German).
- **DOI:** 10.1007/BF02559591
- **Use for:** Foundational — normal surface theory, decidability of unknot recognition. Historical context for why surface-finding is hard

### Algorithms Reference

**Cormen, T.H., Leiserson, C.E., Rivest, R.L., and Stein, C.** *Introduction to Algorithms*, 3rd ed. MIT Press, 2009. ISBN 978-0-262-03384-8.
- **Use for:** The graph and optimization algorithms underneath the above — BFS/DFS, spanning trees, shortest paths, min-cost flow, LP duality, total unimodularity, NP-completeness
- **Cite as:** "CLRS, §X.Y" or "Cormen et al., Ch. N"

---

## Key Implementation Guidelines from Literature

### FVM Implementations (from FVM papers)

**MPFA O-Method (Aavatsmark 2002 - Aavatsmark 2002)**
- Use interaction volumes at vertices with continuity points at subedge midpoints
- Compute transmissibility T = CA⁻¹B - D for local flux-pressure relationship
- Check K-orthogonality for TPFA feasibility

**Physical Space Evaluation (Klausen & Winther 2006 - Klausen 2006)**
- Evaluate Jacobian Dₖ = D(x̂ₖ) at cell CENTER, not at integration points
- Accept non-symmetric system matrix for robustness on rough grids
- Check stability condition det(Λ + Λᵀ) ≥ γ₀ > 0

**Enhanced BDDF₁ for 3D (Ingram et al. 2010 - Ingram 2010)**
- Standard BDDF₁ has only 3 DOF/face → insufficient for MPFA elimination
- Must use enhanced BDDF₁ with 6 curl terms → 4 DOF/face
- Requires h²-parallelepipeds for full O(h) convergence
- Superconvergence O(h²) at cell centers on regular grids

**Convergence Check (Agélas et al. 2008 - Agélas 2008)**
- Compute coercivity condition coer(D,Λ) ≥ θ > 0 before running
- Condition is computable for any mesh + tensor combination
- Predicts convergence on general polygonal/polyhedral meshes

### FEM Implementations (from FEM papers)

**Static Condensation (Messe et al. 2023)**
- **Preferred method** for hanging DOFs
- Avoids zero diagonal blocks in system matrix
- Better conditioning for direct solvers (MUMPS/STRUMPACK)

### Magnetodynamic Coupling (Arsenault et al. 2023)
- Use H-φ/D formulation (not quasi-static)
- Must include ∂(μh)/∂t term
- Better convergence than pure quasi-static

### Thin-Shell Element Order (Alves et al. 2022b, Arsenault et al. 2021)
- **N=1:** Linear through-thickness distribution (faster, similar to T-A)
- **N>1:** Needed when tangential field penetration matters
- Use N>1 for closely-packed tapes

### Convergence Tolerance (Messe et al. 2023)
- Drive ε < 10⁻¹¹ to prevent checkerboarding
- Tighter than typical FEM tolerances
- Critical for nonlinear HTS problems

### Solver Selection (Messe et al. 2023)
- **STRUMPACK:** Typically 2× faster
- **MUMPS:** More robust for difficult problems
- Try STRUMPACK first, fall back to MUMPS if needed

---

## Question Type → Literature Mapping

| Question | Primary Source | Secondary Sources |
|----------|----------------|-------------------|
| How does BELFEM implement X? | Messe et al. 2023, Messe 2022 | Specific papers |
| Theory behind h-φ formulation? | Arsenault et al. 2023, Arsenault et al. 2021 | Bathe Ch. 4, Hughes Ch. 3 |
| Why isn't this converging? | Messe et al. 2023 Section 2.7, Bathe Ch. 8 | Dular et al. 2021 (stability) |
| How do thin-shell models work? | Messe et al. 2023, Alves et al. 2022b | Alves et al. 2022a, Alves et al. 2024, Schnaubelt et al. 2023 |
| What are cuts for? | Alves et al. 2022b, Alves et al. 2024, Schnaubelt et al. 2023 | Alves et al. 2022a (homology) |
| How does MPFA work? | Aavatsmark 2002 | Wheeler et al. 2011 |
| Will my FVM mesh converge? | Agélas et al. 2008, Klausen & Winther 2006 | Aavatsmark 2002 (fundamentals) |
| How to implement 3D MPFA? | Ingram et al. 2010 | Aavatsmark 2002, Klausen & Winther 2006 |
| Physical vs reference space? | Klausen & Winther 2006 | Ingram et al. 2010 (3D details) |
| How to implement element X? | Hughes Ch. 5, Bathe Ch. 5 | Zienkiewicz Vol. 1 |
| What's the weak form for Y? | BELFEM papers | Bathe Ch. 3, Hughes Ch. 1 |
| Why is my element locking? | Hughes Ch. 4, Bathe §4.3-4.5 | Zienkiewicz Ch. 10 |
| How to prove convergence? | Brenner Ch. 4-5 | Bathe Ch. 4 |
| Which solver should I use? | Messe et al. 2023, Bathe Ch. 8 | Hughes Ch. 11 |
| How much Ic defect can a coil tolerate? | Wozniak et al. 2025 | Badel et al. 2019 |
| Will quench detection see a defect in time? | Wozniak et al. 2025 §IV-C | Badel et al. 2019 §5 |
| How do I control the time step in a nonlinear transient? | Badel et al. 2019 §4.3.2 (integral controller) | Bathe Ch. 9 |
| How do I warm-start Newton-Raphson between steps? | Badel et al. 2019 §4.3.1 | Messe et al. 2023 §3 |
| How is current shared between the tape and its stabilizer? | Badel et al. 2019 §3.1 | Wozniak et al. 2025 §III |

---

## Common Pitfalls and Solutions

| Pitfall | Symptom | Solution | Reference |
|---------|---------|----------|-----------|
| Loose convergence tolerance | Checkerboarding | Drive ε < 10⁻¹¹ | Messe et al. 2023 |
| N=1 for stacked tapes | Missing losses | Use N > 1 | Alves et al. 2022b |
| Identical FE orders at interface | Oscillations | Hierarchical enrichment | Dular et al. 2021 |
| Missing cuts | Wrong current flow | Add cuts (thin/thick) | Alves et al. 2022b, Alves et al. 2024 |
| Lagrange multipliers | Zero diagonal | Use static condensation | Messe et al. 2023 |
| QS formulation + linear elements | Poor convergence | Use magnetodynamic | Arsenault et al. 2023 |
| Neglecting ∂(μh)/∂t | Incorrect transient | Include time derivative | Arsenault et al. 2023, Alves et al. 2022b |
| TPFA on non-K-orthogonal grids | O(1) error | Use MPFA | Aavatsmark 2002 |
| Reference space on rough grids | Divergence | Use physical space Jacobian | Klausen & Winther 2006 |
| Standard BDDF₁ for hexahedra | Only 3 DOF/face | Use enhanced BDDF₁ | Ingram et al. 2010 |
| Violating coercivity | No convergence | Check coer(D,Λ) ≥ θ | Agélas et al. 2008 |

---

## Citation Guidelines

### When Referencing Papers in Code/Documentation

**DO:**
- Cite by **author and year**, with the exact section or equation: "Messe et al. 2023, §2.7", "Bathe, §8.6"
- Use full academic citations with author names, journal, volume, pages, year in reference lists
- Include DOIs as clickable links: `DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)`
- Reference the published work directly
- Example: "Following the approach of Arsenault et al. 2023 (IEEE Trans. Appl. Supercond.)..."

**DO NOT:**
- Write a new `paperN` or `FN` alias — they were retired on 2026-08-11 and are decoded above for older records only
- Reference local file paths (e.g., "alves2022a.txt")
- Include extracted text content
- Assume literature files will be available in public repository

**Rationale:** Proper citations ensure academic attribution and provide public access via DOIs.

---

## How to Access These References

### Research Papers
Most BELFEM-related papers are published in:
- IEEE Transactions on Applied Superconductivity
- Superconductor Science and Technology
- COMPEL - The International Journal for Computation and Mathematics in Electrical and Electronic Engineering
- SIAM Journal on Scientific Computing

Search for papers using DOIs (when provided) or author names and keywords like:
- "h-phi formulation superconductor"
- "thin-shell HTS tape"
- "cohomology cuts transport current"

### Textbooks
Standard references available from:
- University libraries
- Publishers: Springer, Wiley, Dover, Prentice Hall
- Some older editions available as Dover reprints (Hughes, Zienkiewicz)
- Felippa's IFEM notes available online as course materials

---

## Workflow Examples

### Implementing a New Feature
1. Understand physics → BELFEM papers for HTS-specific, books for general FEM
2. Find weak form → Research papers
3. Check element technology → Hughes Ch. 5, Bathe Ch. 5
4. Verify stability → Dular et al. 2021 (h-φ specific) or Bathe/Brenner (general)
5. Implement → Cross-reference papers (equations) with books (algorithms)
6. Validate → Papers contain benchmark cases

### Debugging Convergence Issues
1. Identify symptom → Oscillations (Dular et al. 2021, Hughes Ch. 4) or Divergence (Bathe Ch. 8, Messe et al. 2023)
2. Check formulation → Interface coupling (Messe et al. 2023, Arsenault et al. 2023, Alves et al. 2022b)
3. Verify implementation → Compare with cited equations
4. Adjust solver strategy → Messe et al. 2023 (nonlinear), Bathe Ch. 8 (linear)

---

## Special Topics

### h-φ Formulation
Consult in this order:
1. Messe et al. 2023 — BELFEM's implementation choices
2. Arsenault et al. 2023 — Magnetodynamic coupling (recommended)
3. Arsenault et al. 2021 — Performance comparison
4. Alves et al. 2022b — Interface conditions and cuts
5. Bathe Ch. 4-5 — Mixed formulation theory
6. Dular et al. 2021 — Stability analysis

### Thin-Shell Implementation
1. Messe et al. 2023 — BELFEM's macro-element approach
2. Alves et al. 2022b — Full derivation
3. Alves et al. 2022a — 3D extension
4. Alves et al. 2024 — Alternative COMSOL approach
5. Hughes Ch. 5 — Isoparametric element foundations

### Cuts and Cohomology
1. Alves et al. 2022b — Cohomology basis theory (thick cuts)
2. Alves et al. 2024 — Thin cut implementation
3. Schnaubelt et al. 2023 — Practical interpretation for coils
4. Alves et al. 2022a — Relative homology background

---

## Tags for Searching

### BELFEM FEM
`#h-phi-formulation` `#thin-shell` `#cuts` `#cohomology` `#transport-current` `#hts-material` `#belfem-architecture` `#magnetodynamic` `#checkerboarding`

### BELFEM FVM
`#mpfa` `#transmissibility` `#coercivity` `#physical-space` `#enhanced-bddf1` `#k-orthogonal` `#mfmfe` `#anisotropic-diffusion`

### Computational Topology
`#cohomology-computation` `#thick-cuts` `#coreduction` `#smith-normal-form` `#total-unimodularity` `#np-hardness` `#min-cost-flow` `#normal-surfaces`

### General FEM
`#weak-form` `#isoparametric` `#locking` `#mixed-methods` `#inf-sup` `#newton-raphson` `#time-integration` `#eigenvalue` `#convergence`

---

## See Also

- [Documentation Guidelines](documentation_guidelines.md) - How to cite literature in BELFEM docs
- Project README: `../README.md`
- CLAUDE.md: `../CLAUDE.md` (includes literature directory information)

**Note:** If you have access to the full `./literature/` repository, see `./literature/README.md` for detailed navigation guides, routing tables, and search patterns for the complete text extractions.
