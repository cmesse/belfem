# Test Gap Backlog

**Date:** 2026-09-21
**Purpose:** The ranked list of tests the suite is missing, distilled from the incident record.
Every row names the mechanism a test would pin, the suite it belongs in, and what would have gone
red before the fix. It is a list, not code: nothing here has been written yet.
**Module:** `tests/` (all suites); enabling fixtures under `tests/fem/support`, `tests/mesh`, `tests/comm`
**AIs involved:** Claude Fable 5.1 (five blind slice sweeps with tree access, reconciliation,
verification), Claude Fable via web API (two pre-sweeps over the raw devlogs, one with tree access
at `ef6d47f2`), Christian (corpus preparation, rulings). Codex and Grok were not used on the sweep;
they take the jury on each test as it is written.
**Status:** OPEN — list compiled and citation-checked 2026-09-21; no test written; rulings O1–O9 pending

> **Scope guards:**
> - This file ranks tests. It does not fix the defects the sweep found alive; those are listed
>   in §2 so they can be filed or fixed, and their tests are the red-then-green gates.
> - The cohomology core (`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`,
>   `cl_Cochain`, `fn_Smith`) is closed to AI edits. Tests may exercise it through its public
>   interface; a test that needs a change inside it is listed under §4 for Gregory.
> - Infrastructure (sanitizer, coverage, CI legs) is in §6, separate from tests.
> - Priority follows the criteria the sweep was briefed with: a test pins a mechanism, not a
>   symptom; it would have gone red on the pre-fix tree; cheapest tier that reaches the
>   mechanism; a defect found once and left without a test can be found again
>   (Martin 2025 Ch. 15, Ch. 17, Ch. 32; Thomas & Hunt 2020 Topic 51).

---

## 1. Method and what was verified

**Inputs.** The incident catalog `doc/lessons_learned_evidence.md` (601 rows, 20 threads, locked
2026-09-01), the 49 devlogs dated September 2026 that postdate the lock, the seven open rows of
`todo/debt_register.md`, and two web-API pre-sweeps over the merged devlogs. The catalog was split
into seven thread groups; each group, the September devlogs and the register rows went to one
blind reviewer with read-only tree access and the brief in `tmp/ai_exchange/test_gap_sweep_prompt.md`.

**Verification performed (Claude, 2026-09-21).**

- Every incident id, register id and devlog file name cited in the seven answers resolves
  (script over all files: 0 missing of 352 incident cites, 73 register cites, 50 devlog cites).
- The "red today" claims that set the top of the list were re-read in source by the reconciler,
  not taken from the reviewers: §2 marks which.
- Cross-reviewer duplicates were merged by mechanism; where two reviewers disagreed on cost or
  reachability the more conservative reading is kept and the disagreement noted.

**Not verified.** Cost estimates are the reviewers' judgment. "Red pre-fix" is an argument from
the incident record, not a run against the pre-fix tree. Whether each proposed fixture can be
built at all is the open question of §3 and is the single largest uncertainty in the ranking.

---

## 2. Defects the sweep found still in the tree

These are not test gaps; they are open defects. Each needs a ruling (fix, file, or accept) and
its test is the red-then-green gate. Confidence is on the defect being present, not on severity.

| # | Defect | Where | Verified by | Test (see §3 id) |
|---|---|---|---|---|
| L1 | `get_bottom_nodes`/`get_top_nodes` HEX20 case sizes 9 and reads `node( 21 )` on a 20-node element (copied from the HEX27 case) | `src/mesh/meshtools.cpp` `get_bottom_nodes` HEX20 branch and the matching top branch | Claude, read | T-01 |
| L2 | `BlockData::collect_thin_shell_facet_ids` sizes the output from shells that own a selected block, then fills it from every shell: an out-of-bounds write whenever a thin shell has no selected block | `src/fem/kernel/cl_FEM_DofMgr_BlockData.cpp` `collect_thin_shell_facet_ids` | Claude, read | T-02 |
| L3 | `mesh::interpolation_order()` has no QUARTIC branch; LINE5/TRI15/TET35 report UNDEFINED although their shape-function classes declare QUARTIC | `src/mesh/meshtools.cpp` `interpolation_order` switch | Claude, read | T-03 |
| L4 | `hdf5::group_exists` and `hdf5::dataset_exists` are both bare `H5Lexists`; neither can tell a group from a dataset | `src/io/hdf5_tools.hpp` both helpers | Claude, read | T-04 |
| L5 | `Map::operator()` release path aborts with "Key not found in map." and no key; the debug path names it | `src/containers/cl_Map.hpp` release branch | Claude, read | T-05 |
| L6 | `.bfm` loader reads the writer's git hash into metadata and prints it, never compares it with the running binary (DR-160) | `src/mesh/cl_Mesh_BfmFile.cpp` load path | Claude, read | T-06 |
| L7 | No plugin ABI stamp exists; a user-material `.so` built against an older vtable registers nothing and the layer becomes Buffer (INC-316) | `src/physics/materials` (grep for any ABI guard: none) | Claude, grep | T-49 (needs the guard first) |
| L8 | `Ascii` resolves a relative path against `getenv( "PWD" )`, not the process cwd; under ctest the two differ | `src/io/cl_Ascii.cpp` constructor | Claude, read | T-07 (O1) |
| L9 | A deck line `key : value // comment;` is dropped silently: comments are stripped at `//` and `Section` registers only lines that still carry `;` (INC-361) | `src/io/cl_Input_Section.cpp` registration | reviewer, read | T-08 |
| L10 | `SolverParameters::set_relative_tolerance` accepts 0, negative and non-finite values; the absolute sibling validates (INC-228 residual) | `src/sparse/cl_SolverParameters.cpp` | reviewer, read | T-09 |
| L11 | `CutData::collect_coefficients` tests `== ±1` only; a ±2 coefficient is dropped to weight 0 at this layer, caught today only by the upstream post-rectification error (INC-479) | `src/homology/cl_CutData.cpp` coefficient loop | Claude, read | T-44 |
| L12 | Nickel thermal conductivity coefficient set gives about half the handbook value at 300 K; flagged for refit in a devlog, tracked nowhere (INC-129) | `src/physics/materials/cl_Material_Nickel.cpp` | reviewer, from the record; not evaluated | T-30 (O3) |
| L13 | `CurveFactory::sort_end_nodes` accepts an intersection that is two disjoint chains, then indexes adjacency with the `gNoIndex` sentinel | `src/mesh/cl_CurveFactory.cpp` | reviewer, read | T-38 |
| L14 | Two mesh globals created with the same explicit id collapse the id map; only the call site moved to auto-ids | `src/mesh/cl_Mesh.cpp` `create_global_variable` | reviewer, read | T-39 |
| L15 | HEX8TB wall reload: `reconstruct_edge_connectivity` takes the FIRST twin for wall slots ≥ 2, `normalize_edge_hangs` re-points 1:1 hangs to the LAST twin (INC-066, "pending Christian") | `src/mesh/cl_ProtoMesh.cpp` | reviewer, read | T-52 (O4) |
| L16 | `DISABLED_Gas_Entropy_RefState` is still disabled; the enabled `Cubic_Entropy` gate is an r² test blind to a constant offset, which is exactly the defect class (INC-419, INC-426) | `tests/physics/gasmodels/cl_GM_Gas_Consistency.cpp` | reviewer, read | T-27 |
| L17 | `Vector<T>::operator=( expression )` resizes the target silently in both build types (`b = bt + bn` grew 2 → 3 with no diagnostic, INC-364) | `src/linalg/armadillo/cl_AR_Vector.hpp` (Blaze twin unchecked) | reviewer, read | T-53 (O5) |
| L18 | Grouped terminal curves: `updatekGeneratorsFromHomology` discards `s` from `smithForm`, so a k-crossing group carries k× the imposed current with no refusal | `src/homology/cl_Cohomology.cpp` (closed core) | reviewer, read; reported to Gregory 2026-09-19 | §4 |
| L19 | `compute_edge_directions_thinshell()` tests `physical_tag() == 1` as a direction flag while the tag is set as "first material" (recorded "looks accidental", never ruled) | `src/fem/kernel/cl_FEM_Element.cpp` | reviewer, read | O6 |

L1 and L2 were recorded in the catalog as "confirmed real but latent" (INC-084, INC-072) and are
still in the tree: a latent-defect note that nobody re-greps is a documentation failure, not a
test failure. Both deserve a register row or a fix this week.

---

## 3. The backlog

Tiers: **1** = serial gtest suite, **2** = MPI suite at 2 and 4 ranks, **F** = needs a fixture
from §3.5 first. Cost: S = an afternoon on an existing fixture, M = a day, L = new fixture or
suite. Sources: `INC-n` catalog rows, `dl…` September devlogs, `pre` = web pre-sweep proposal
numbers, `s<k>#n` = slice `k` proposal `n` in `tmp/ai_exchange/test_gap_sweep_claude_*.md`.

### 3.1 Write first: red today, cheap, mechanism verified

| id | suite | test asserts | sources | cost | |
|---|---|---|---|---|---|
| T-01 | fem, extend `TopBottomNodeAlignment` | for every case label of `get_bottom_nodes`/`get_top_nodes`: size equals the face node count, every index `< number_of_nodes()`, bottom at min ζ (2-D: min η), top at max, positionally paired | L1; INC-084, INC-367, INC-063; s4#2 | S | [ ] |
| T-02 | fem (DofSeeding fixture + two thin shells) | with two thin shells and a block selection covering one, the returned facet-id list has exactly that shell's facet count and no id of the other, and no bounds assert fires | L2; INC-072; s4#1, dl20260904_thinshell_volume_shell_align.md (s8#17) | M | [ ] |
| T-03 | mesh | for every element type the interpolation factory builds, `interpolation_order()` equals the shape class's declared order and is never UNDEFINED; `interpolation_order_numeric()` returns | L3; INC-404; s7#1 | S | [ ] |
| T-04 | io (HDF5Test fixture) | for a file holding dataset `x` and group `g`, `group_exists( "x" )` and `dataset_exists( "g" )` are both false | L4; INC-028 residual; s1#8 | S | [ ] |
| T-05 | containers | the release-path error of `Map::operator()` names the key (assert under `NDEBUG`, or compare both branches' formats) | L5; INC-308, INC-528; s6#10 | S | [ ] |
| T-06 | mesh (BfmFile fixture, §3.5 F3) | writing a `.bfm`, overwriting its `git` dataset with a foreign hash and reloading is refused (or raises the flag the ruling chooses) unless the override is set; a matching hash loads | L6; DR-160; s9#1 | M | [ ] |
| T-07 | io | `Ascii( "rel.txt" )` opens the file in the process cwd when `PWD` names another directory — polarity depends on O1 | L8; dl20260901_corc_solder_hn_recovery.md; s8#14 | S | [ ] |
| T-08 | io (Deck helper from `test_GhostSwitch`) | a section line `name : value` with no `;` outside a comment makes the parser throw, or at minimum `key_exists` is false AND a diagnostic is emitted | L9; INC-361; s5#14 | M | [ ] |
| T-09 | sparse | `set_relative_tolerance( 0 / -1e-8 / NaN / Inf )` throws and leaves the stored value untouched | L10; INC-228; s3#6 | S | [ ] |

### 3.2 Pin the fix: P0, existing fixture, an afternoon each

| id | suite | test asserts | sources | cost | |
|---|---|---|---|---|---|
| T-10 | physics + fem | (a) for each built-in pure metal at RRR 50, `( rho_i( T_ref ) + rho_0 ) / rho_0 == RRR` to 1e-9; (b) a `Tmatrix` built from a matrix whose nonzeros are all 0.25 keeps every nonzero and reproduces `T·x` | dl20260913_ybco_test_int_abs.md; pre#5, s8#2, s8#3 | S | [ ] |
| T-11 | io (HDF5Test) | loading an int64 dataset into a 32-bit `Vector<int>` throws in BOTH build types; a 4-byte dataset read through `load_bool_from_file` into a canary-guarded bool leaves the canary untouched | INC-330, dl20260831 DR-116; pre#3, s6#19 | M | [ ] |
| T-12 | linalg | `Matrix<real> B(3,2), A(2,2); C = B * inv(A)` is 3×2 with hand-computed entries; `2×1 · inv(1×1)` does not abort; same on both backends | INC-271, INC-272; pre#4, s4#6 | S | [ ] |
| T-13 | mesh + math | narrow-counter family: a `Basis` given 300 sources reports 300 and returns the last intact, above the member width it throws rather than wraps; `graph::Vertex` with 70 000 insertions reports 70 000; `Vertex::allocate_element_container( 256 )` throws before allocating; reset twice then destroy is clean | DR-108, DR-13, INC-036, dl20260904_dlist_code_fixes.md; pre#1, s1#3, s8#21 | S | [ ] |
| T-14 | core | a message containing `%s`, `%d` and `50%` written through `Logger` arrives byte-identical | dl20260913_ubuntu_format_security.md; s8#1 | S | [ ] |
| T-15 | comm or fem | after the test main's communicator init, `gRhoMax > gRhoMin` and `gRhoMax ≥ 1e6`; the clamp is the identity on `[1e-12, 1e2]` | INC-261; s4#4 | S | [ ] |
| T-16 | physics (test database) | an `Alloy` reports `depends( rho, normB )`, `depends( rho, angleBxJ )` and the same for `lambda`; a T-only metal reports neither | INC-270; s4#5 | S | [ ] |
| T-17 | sparse (`BELFEM_STRUMPACK`) | the option setter on a default `SolverParameters` yields rel_tol == class default, abs_tol == 1e-14, maxit == 50, Krylov == PREC_GMRES, compression NONE at 100 000 rows, matching on; with BLR the cutoff equals `compression_cutoff()`; `set_krylov_method( CG )` throws on the first solve | INC-207, INC-208, INC-209, INC-226, INC-228, INC-229; s3#8, s3#7, pre#20 | M | [ ] |
| T-18 | sparse (`BELFEM_MUMPS`) | `warning_message( 9 )` decodes to exactly the +1 and +8 lines with only +1 fatal; `error_message` table over `(INFO(1), INFO(2))` pairs decodes the right size and ICNTL number with no literal `%i` | INC-206, INC-214, INC-215, INC-499; s3#1, s3#3 | S | [ ] |
| T-19 | linalg | a degree-8 polynomial sampled at 60 points on `[0,1]` is recovered by `polyfit` with max relative coefficient error < 1e-6 on the built backend | INC-454; s3#14 | S | [ ] |
| T-20 | sparse (`BELFEM_PETSC`) | a solve capped at one iteration on a system that needs more returns a non-success status or throws, never status 0; every `KrylovMethod` and `ReorderingMethod` enum value including AUTO/AUTOMATIC solves the tridiagonal without a PETSc error | INC-319, INC-355; s6#15, s5#13 | S | [ ] |
| T-21 | fem (`GhostSwitch` cache-tag fixture) | decks differing only in `…,11` vs `…,13`, only in a circuit terminal pair's `input curves`, or with the ids of two unlabelled pairs swapped, produce different `mesh_config_tag` values | INC-383, INC-386, INC-518; s6#21 | S | [ ] |
| T-22 | fem + physics | `bn_angle( b, n ) + bn_angle( b, −n ) == π` for random `b`; `JcFunctionDatabase` at θ = π returns the θ = π row, not θ = 0, in value and all three derivatives on a deliberately asymmetric table | INC-322, INC-323, INC-324; s6#16 | S | [ ] |
| T-23 | homology (annulus fixture, public interface) | a generator built through the suggest path and passed through `reorient_generators` pairs with the CCW inner ring as exactly −1 (the existing test takes `abs()` and cannot see the 2-D mirror) | dl20260828 2-D current sign; pre#10 | M | [ ] |
| T-24 | fem (`InterfaceOrientation` PENTA6TS ghost fixture) | on one ghost facet with random inputs `Kms == trans( Ksm )` and the assembled `K` is symmetric to 1e-14; for a field constant across the interface the ghost residual is zero; for a unit jump the consistency term equals `rho_harm·[h]/h` with the harmonic mean, finite at `ρm == ρs` | INC-090, INC-093, INC-094, INC-095, INC-100; s6#1, s6#2 | M | [ ] |
| T-25 | physics | for every builtin label, cp, lambda and rho sampled on a 0.5 K grid over the declared range are finite, positive, and no adjacent relative jump exceeds a few percent, including at every switch temperature | INC-103, INC-123, INC-125, INC-126, INC-127, INC-130, INC-132; s2#1 | M | [ ] |
| T-26 | physics | for every builtin material, `dcpdT`, `drhodT`, `dlambdadT` agree with a central difference of the value accessor to 1e-6 relative at B = 0 and 5 T; every `rho_*` / `drho_*_dJ` overload pair agrees across the three regimes and on both sides of the `mRhoMin` floor; `drho_piecewise_dT` matches an FD of `rho_piecewise` for n ∈ {2, 4, 7, 20} at nff = 3 with no NaN | INC-114, INC-115, INC-120, INC-121, INC-122, INC-125, INC-127, INC-292, INC-296; s2#2, s2#3, s5#15 | M | [ ] |
| T-27 | gasmodels | methane: `p( Tc, ρc ) = pc` to 1e-5, `∂p/∂ρ|_T > 0` for every T in `(Tc, Tc+10 K]`, ancillary gives 101 325 Pa at the NBP; cubic air: absolute `s( T_ref, p_ref )` equals the mixture standard entropy within 1 % (re-enable `Gas_Entropy_RefState`) | INC-416, INC-423, INC-424, INC-426, INC-427; s2#13, s2#14 | S | [ ] |
| T-28 | math | for the Silver-high control net and random nets, `xi_by_x( x_by_xi( ξ ) ) == ξ` and `xi_by_y( y_by_xi( ξ ) ) == ξ` on a ξ sweep of `[−1, 1]`; `x` outside the range returns ±1 | INC-123, INC-124; s2#9 | S | [ ] |
| T-29 | physics (needs `share/material/sp-ap.hdf5` reachable, check) | on the shipped REBCO table, `jc( B )` is monotone non-increasing in B over the whole (T, θ) grid, `jc( 0 ) > jc( B_min )` everywhere, `jc( 0, θ1 ) == jc( 0, θ2 )`, and `d jc/dB` has no jump at the former 10 mT clamp; `n` non-increasing in log₁₀B | DR-128, INC-549, INC-553, INC-554, dl20260904_rebco_tables_n_monotone_b16.md; pre#28, s2#12, s8#25 | M | [ ] |
| T-30 | physics | `Nickel::lambda( 300 K )` within 15 % of 90.7 W/(m K) and `Iron::lambda( 300 K )` within 15 % of 83.2; Nickel reduced magnetization at 600 K in `[0.42, 0.52]` and 0 at Tcurie — land `DISABLED_` until O3 | INC-128, INC-129; s2#10, s2#11 | S | [ ] |
| T-31 | fem (DofSeeding fixture + `IWG_Timestep`) | with two dof fields and two history steps, `qold` on a block calculator returns the value written for (step 1, field 0), not (step 0, field last); changing the timestep method afterwards still returns the right value | INC-264, INC-373, INC-514; s3#13 | M | [ ] |
| T-32 | mesh (tensor mesh, two opposite sidesets) | after `create_periodicity()`: every source node has exactly one partner at source + period (bijection); a target perturbed by 1e-6 m is rejected; on a mesh scaled to 10 m every pair equals the brute-force nearest neighbour; with edges finalized on a conductor crossing the plane the rebuild ties > 0 conductor edges | INC-006, INC-007, INC-009, INC-022, INC-031; s1#1, s1#2, s1#4, s8#11 | M | [ ] |
| T-33 | fem | pinning a hanging dof with one unit-weight source fixes the SOURCE to the value; a non-unit single-source weight throws in debug | DR-135; pre#11 | M | [ ] |
| T-34 | mesh | a face-creation test: on two TET10s whose midside indices sort below the corners every created face pairs exactly two coincident facets and no key is shared by three owners; after connectivity every face reports its type's edge count and each face neighbour appears once | DR-102, INC-233; pre#6, s7#2 | M | [ ] |
| T-35 | fem | `ElementMapper` on TRI3/QUAD4/TET4/HEX8/PYRA5: xi→x→xi round-trips on affine and distorted elements, a non-convergent Newton returns false with finite output, `inside_pyra` rejects points outside the cross-section | INC-238, INC-239; s7#4 | M | [ ] |
| T-36 | mesh | after relinking a master element's node and re-finalizing, a `ThinShell`-typed sideset facet keeps its restored node list while a plain sideset facet is re-derived | dl20260825 facet clobbering; pre#32 | M | [ ] |

### 3.3 Pin the fix: P1 and P2, cheap

| id | suite | test asserts | sources | cost | |
|---|---|---|---|---|---|
| T-37 | fem | a subclass constructed straight from `IWG` with default arguments reports `SymmetryMode::Unsymmetric`; so does `IWG_MaxwellPostproc` | INC-544, INC-552; s5#12, s7#11 | S | [ ] |
| T-38 | mesh | an open LINE2 polyline whose lowest element id is not the start sorts to a monotone node list identical to the reverse sort; two disjoint chains handed to `sort_end_nodes` abort with an error naming the curve (L13) | dl20260830 CurveFactory, dl20260903_pancake_deck_curve_crash.md; pre#31, s8#12 | M | [ ] |
| T-39 | mesh | two auto-id globals get distinct ids and are retrievable by label and id with the value passed; two explicit globals with the same id are refused or both retrievable (L14) | INC-313, INC-314, dl20260901_heatloss_and_iv_globals.md; s6#12, s8#7 | S | [ ] |
| T-40 | comm (Tier 1) | at comm size 1, `allreduce` and `allreduce_min` copy send→recv for counts 0, 1, n into a distinct buffer; a count above `INT_MAX` is refused | INC-565, INC-566; s7#10 | S | [ ] |
| T-41 | sparse | `solver_type( to_string( k ) ) == k` for every `SolverType` below UNDEFINED, names unique, `matrix_type` and `preferred_matrix_format` do not throw, `solver_type( "no-such" )` throws; `check_naupd( 3 )` mentions NCV, `check_naupd( 1 )` warns, `check_neupd( 1 )` is the Schur message, `check_saupd( 3 )` warns and returns | INC-204, INC-212, dl20260904_dlist_code_fixes.md; s3#4, s3#11, s8#20 | S | [ ] |
| T-42 | sparse (`BELFEM_SUPERLU`, `BELFEM_ARPACK`) | SuperLU solves the tridiagonal in CSR and CSC to the dense reference, a re-solve after in-place value mutation reflects the new values, an unlinked build throws like the others; ARPACK `nev = 5` on a tridiagonal returns the five extreme eigenvalues to 1e-10 with `info(3) == 5` and `info(7) ≥ 11` | INC-203, INC-204, INC-211, INC-213, INC-409, INC-536, INC-546; s3#5, s3#10, s7#7, s7#8 | M | [ ] |
| T-43 | io + math | a vlen string dataset written `H5T_CSET_UTF8` loads verbatim; `Spline::save` then `load` reproduces knots, coefficients and evaluations | INC-240, INC-550; s7#9, s7#16 | S | [ ] |
| T-44 | homology (public interface) | a cochain carrying a ±2 coefficient handed to `CutData` is rejected naming the edge (L11); on a disk with 2 and 3 holes every algorithm returns exactly k generators with a unimodular pairing matrix; `Homology` and `Cohomology` generator counts agree per dimension on every fixture; after `clean()` on the perturbed annulus the support is no larger than before perturbation | INC-406, INC-433, INC-437, INC-479; s1#13, s1#14, s1#15, s7#18 | M | [ ] |
| T-45 | mesh | gmsh code → `ElementType` round-trips for every type, 36 → QUAD16, unknown throws naming the number; the same geometry read as 2.2 and 4.1 yields identical physical tags; a truncated `$Entities` line throws; N point elements yield N vertices with ids equal to tags | INC-308, dl20260904_dlist_code_fixes.md; s6#11, s8#9, s8#10 | M | [ ] |
| T-46 | mesh + fem | `load_fields` sizes an empty target to the dump, aborts on a sized target of different length, loads equal lengths; the number of BDF history levels equals the order whether `set_timestepping_method` is called before or after `create_old_dof_fields()` | INC-333, INC-334, DR-112, dl20260902_memdump_history_shell_length_guard.md; pre#24, s5#11, s8#15 | M | [ ] |
| T-47 | mesh (Exodus in build) | saving a mesh with one empty block yields a block count equal to the non-empty count and `num_elem` equal to the element sum; reading it back succeeds | INC-304; s6#9 | S | [ ] |
| T-48 | mesh | element topology tables: `get_nodes_of_facet( k )` has the facet's node count, each `get_edges_of_facet` edge is a cyclic-consecutive node pair of the facet, and the set of facets that throw is exactly the documented lateral set; `get_nodes_of_edge` identical across HEX8/20/27/64; a mesh with k node-reversed elements makes `MeshChecker` flip exactly k, and running it after `create_edges` throws; PYRA13/14 reorientation keeps midsides at edge midpoints | INC-080, INC-081, INC-084, INC-237, INC-446, INC-516, DR-153; s4#3, s7#3, s2#18, pre#8 | M | [ ] |
| T-49 | physics (plugin fixture from `tests/circuit/plugin`) | a plugin whose exported ABI stamp differs from the library's is refused with a named error; a matching plugin's registrations set `have( rho )` — needs the guard first (L7) | INC-316; s6#13 | M | [ ] |
| T-50 | linalg + fem | a 2×2/3×3 Gram scaled by `h = 1e-4` inverts with relative error < 1e-12 and a sliver throws in debug; `static_assert` that `Controller::time()` and `timestep()` return lvalue references (bound to `MaxwellData::mTime`); `Kernel`, `Solver`, `ElectricalCircuit` are not copy-constructible | INC-259, INC-273, DR-129, DR-117, DR-141; s4#7, s4#8, pre#27 | S | [ ] |
| T-51 | physics | `Alloy( "pb40sn60", NaN )` throws the intended message, `Alloy( label, 5 )` and `Alloy( "pb40sn60rrr5" )` construct identically; a deck with the three `builtin` forms and an alloy with `RRR : 5` yields four materials with RRR applied; after `set_constant( jc )` and `set_constant( n )`, `jc()`, `n()`, `rho_powerlaw` are finite and equal the closed form | INC-021, INC-104, INC-108, INC-117, INC-118; s2#6, s2#7, s2#8, pre#19 | M | [ ] |
| T-52 | mesh (BfmFile fixture) | save (tag T) → load → save → load yields tag T at every stage; a base mesh saved, loaded and saved again reproduces connectivity, facet, face and control-point datasets byte-for-byte; two edge records with identical node pairs both survive a load; on the hand-wired `PeriodicBar` fixture every periodic back-pointer survives save→load and the plane matrix is 3×3; HEX8TB wall slot 2/3 equals its hang target (L15, O4) | INC-028, INC-066, INC-382, INC-385, INC-387, INC-388, INC-395; s1#7, s6#22, s6#23, s6#25, s4#12 | M | [ ] |
| T-53 | linalg | `Vector<real> b( 2 ); b = bt + bn` with a 3-long right side throws in debug once the assert exists (L17, O5); `lapack::leading_dimension` on a shrunk matrix gives the value a `gesv` round trip needs on both backends; after grow-then-shrink `broadcast`, `distribute` and `collect` of a `Matrix` deliver the current shape on every rank (commmpi) | INC-364, dl20260809 Matrix-MPI; pre#15, pre#29 | S | [ ] |
| T-54 | core + fem | with `set_throw_on_error( true )` a failed `BELFEM_ERROR` throws `std::runtime_error` carrying the message in both build types and the default equals `BELFEM_ASSERTIONS_ACTIVE`; `residual_string( 5.8e4 )` renders in scientific notation, `residual_string( 9.5 )` is not "9.000000"; a Controller (or its reader) with no `anderson` key reports Anderson off | INC-254, INC-255, INC-275, INC-287; s7#15, s5#21, s5#19 | S | [ ] |
| T-55 | mesh + fem | `unflag_all_edges` on a mesh with edges on one block completes and `edge( 0 )` on an edge-less element throws in debug; after `allocate_element_container( 3 )` then `( 0 )`, `number_of_elements()` is 0 and the container null; two 2-D line facets sharing one node report consistent/inconsistent correctly for both orderings and a 3-D pair's verdict is unchanged; `select_sidesets( {2, 4} )` creates exactly two sideset dof tables; the ghost `number_of_dofs_per_element( SideSet* )` count sums master and slave | INC-014, INC-032, INC-089, INC-098, INC-448, INC-451, INC-452; s1#10, s2#19, s2#20, s6#4, s6#5, pre#9 | S | [ ] |
| T-56 | io + circuit + homology | `42:48` → 7 id groups of 1 and `[42:48]` → 1 group of 7; a `terminal pair` carrying two bracket groups is refused naming the section; a cohomology call with more conditions than generators throws a named error in debug | INC-317, INC-318, INC-443; s6#14, s5#17 | S | [ ] |
| T-57 | sparse | `rcond( diag( 1..n ) ) ≈ 1/n` and `[[1, 1e3], [0, 1]]` returns O(1e-6) not O(1); a graph with a duplicate column index makes the debug `SpMatrix` constructor throw; `multiply( x, y )` with `y` pre-filled with NaN equals the 5-argument overload on a zero-filled `y` (MKL and Fortran branches); `Integrator::time()` equals the sum of accepted steps; `faces_are_finalized()` tracks face finalization | INC-213, INC-240, INC-328, dl20260901_dr155_plan_jury.md, dl20260904_dlist_code_fixes.md; s3#12, s6#18, s7#17, s8#4, s8#8 | S | [ ] |
| T-58 | physics | at `|j| = 10³·jc`, n = 40, every HTS law family returns a finite `ρ ≤ ρₙ` and `rho_riva` returns exactly `ρₙ` at its cap; every species in `share/fluid` constructs in SPLINE and POLY mode with cp, μ, λ continuous across interval junctions; `trans.inp` keeps the NIST low-T intervals and `gasdata.inp` has unique labels | INC-418, INC-421, INC-422, INC-515, dl20260901_dotq_rename_and_docs.md; s2#16, s2#17, s8#26 | M | [ ] |
| T-59 | fem | for TRI3/TRI6/QUAD4/QUAD9 volumes every edge facet's mapped integration points coincide with direct 1-D integration on master and every slave orientation; for every (facet, orientation) of a TS prism the slave integration index is in range and distinct; the rule `IntegrationData::populate` selects for order 0 is exact to `auto_integration_order()` for every geometry | INC-097, INC-234, dl20260903_intpoints_exactness.md; s7#6, s6#29, s8#24 | S | [ ] |
| T-60 | fem + physics | facet normals equal the normalized tangent cross product, unit length, outward from the master; `Pipette` volumes of a unit tet/hex match; `MaxwellData` over a `Metal` with a B-H curve neither throws nor selects the constant-μ branch; a default `ThinShellFactory` reports its fuse/connector flags false; a non-`Metal` throws on `set_bh_curve` | INC-057, INC-058, INC-253, INC-445, INC-450; s7#5, s2#21, s2#23, s4#15 | M | [ ] |
| T-76 | sparse (`BELFEM_PARDISO`, MKL build only) | after symbolic factorization the matrix's entry base is the caller's base and `operator()` reads the same entries as before (`PARDISO::initialize` left the matrix 1-based) | DR-101 strike; pre#22 | S | [ ] |
| T-77 | homology (public interface) | a generator with a non-unit coefficient that no thin cut can absorb makes `clean_spfa` throw a named error in debug rather than return a silently unrectified cochain (the annulus fixture never exercises an infeasible system) | DR-23 stage 1, dl20260905_thick_cut_abort_wording.md; pre#23 | M | [ ] |

### 3.4 Blocked on a fixture (P0 and P1 behind §3.5)

| id | fixture | test asserts | sources | cost | |
|---|---|---|---|---|---|
| T-61 | F1 | at np 2 and 4, with rank 1 owning shell 2 but no element of shell 1, every layer element's volume equals its own shell's facet area × thickness; with a shell block whose owned elements are interleaved across ranks every `_Volumes` weight equals its own facet | dl20260904_thinshell_volume_shell_align.md, INC-327; pre#14, s8#16, s6#17 | L | [ ] |
| T-62 | F1 | `DofManager::solve` and `SolverData::collect_matrices` at 2 ranks pass their `comm_drain_check` chokepoints with zero aborts (DR-113 gate a); the global free-dof count is partition-independent; a stray born inside `collect_matrices` is named at the next boundary | DR-113, DR-105 (2), INC-218, INC-533; s9#2, s3 6.3 | L | [ ] |
| T-63 | F1 | after `synch_source_field` on an order-2 field every aura edge holds both dofs equal to the owner's (DR-31); when the thermal peer's linked element id differs, the temperature used is `gTbulk` (DR-30); the sideset integral returned on rank 1 equals rank 0's; after recovery no interface-duplicate row is 0.0 on any rank and every duplicate's original is present with resolving links | DR-30, DR-31, INC-217, INC-299, INC-302, INC-303, dl20260916 F1; pre#2, pre#16, s6#7, s6#8, s8#23 | L | [ ] |
| T-64 | F1 | after a warm restart from a history-carrying dump at 2 and 4 ranks every rank's history fields at the verified depth are non-empty and equal rank 0's on shared entities; serial face: the label set the restore distributes equals the numbered history labels the IWG registers | INC-341, INC-343, INC-534, INC-535, DR-78; s5#9, s5#10 | L | [ ] |
| T-65 | F1 | after distribution of a small enriched mesh every owned hanging entity's sources exist locally on every rank and the halo is one layer; `set_currents( aI )` with no abstract-node dofs on a worker returns without throwing | INC-023, INC-025, INC-026, INC-517; s1#11, s1#12 | L | [ ] |
| T-66 | F2 | after writing NaN into a free dof and `reset_timestep()`, every free dof equals its field value and the next residual is finite; after the loop returns a fresh `compute_residual()` of the committed state is below tolerance; after any number of resets `mDeltaTime ≥ mDeltaTimeMin` and `mIteration2` counts thermal iterates | INC-540, INC-346, INC-353, INC-293, INC-277, DR-127, DR-91; pre#12, pre#25, s5#1, s5#2 | L | [ ] |
| T-67 | F2 | fed the logged sequences: a ±0.05 dB plateau over 5 samples trips the stall detector, `+8.91 → +5.05 → −65.9` does not cut, a −67.6 dB attempt with ω recovering is not cut by the watchdog; with `absolute tolerance : 1e-6` and a residual floored at 1e-7 the attempt is accepted; on a magnetic-only Controller the moved-baseline branch never accepts; a deck with `tolerance switch` below `tolerance` throws naming both keys; a stall exit is taken only within the ruled floor (O2) | INC-279, INC-286, INC-291, INC-293, INC-319, INC-356, INC-541, DR-159; s5#5, s5#6, s5#7, s5#20, s9#6 | L | [ ] |
| T-68 | F2 + HDF5 | `save_memdump` then `load_memdump` on a fresh Controller restores `mTime` and history, re-enters at `min( dumped, initial )` Δt, restores nothing with `restart : false`, and throws when the dump carries a circuit group while none is attached | INC-295, INC-339, INC-340, INC-351, INC-441; s5#8, s5#18 | L | [ ] |
| T-69 | DofSeeding fixture + solver | after a Newton solve on a linear system `residual()` equals `‖A·x_new − b_original‖/‖b_original‖` computed independently (≈ 0), not 1.0 and not the pre-solve value; with Anderson on, the bootstrap-iterate residual is O(1) against a current-assembly recomputation, not ≤ 1e-12; one thermal step of a uniform-T adiabatic problem leaves every nodal T equal to T0 to 1e-12 independent of the linear rtol (O7) | INC-281, INC-288, INC-357, INC-360, INC-531, INC-532; s5#3, s5#4, s5#16 | M–L | [ ] |
| T-70 | F5 (IWG element assembly) | after one element assembly `dotQ` from the Picard and Newton μ₀ kernels agree and equal `∫ρ|j|²dV` from the calculator's own points; a ferro element's `dMdX_times_h` block changes sign when the history field is reversed and vanishes as h → 0; after `Calculator::allocate()` every conductor-block calculator owns a `MaxwellData` with its dispatch set and a non-null thermal peer on a coupled pair; after a vector solve a hanging dof's value equals the weighted sum of its sources; the three FD derivative helpers scale by their own element steps | INC-064, INC-260, INC-262, INC-265, INC-266, INC-449, INC-512, dl20260901_heatloss_and_iv_globals.md; pre#17, s2#22, s4#9, s4#10, s4#11, s8#6 | L | [ ] |
| T-71 | F7 (Postprocessor) | for a field that jumps across an interface the recovered value at the original equals its own side's value and at the duplicate the other side's; on a one-layer TS fixture the written layer B has a non-zero tape-normal component `μ·h_n` and JJC's `b` equals the solve's; a PENTA6TS block reclassified to `Buffer` is still skipped by the aura pass | INC-298, INC-300, INC-331, INC-332; s6#6, s6#20, s6#28 | L | [ ] |
| T-72 | F4 (ThinShell/CutFactory serial) | after save→load of a thin-shell-enriched mesh `number_of_facets()` is unchanged and no facet appears in two sidesets unless it did before; after `create_thin_shell_cuts` every terminal-curve facet carries an orientation and no node has `index() == gNoIndex`; on an all-air 2-D tape stack tips are not duplicated and each interior tape node once; an empty thin cut is refused with the named error; a closed LINE3 loop yields 2n nodes | INC-315, INC-363, INC-366, INC-391, INC-392, INC-393, dl20260902_corc_strip_diagonal_and_d2_gate.md; s6#24, s6#26, s4#13, s4#14, s8#18 | L | [ ] |
| T-73 | F6 (python `belfem_conf`) | every key literal consumed via `key_exists`/`get_*`/`section_exists` in `src/` (case-normalized, multi-line aware) has an anchor in `doc/input_schema.yaml` and every anchor resolves (the executable form of the Input Contract, `drift.py` exists and nothing runs it); a deck with a comment between two statements yields exactly `key1`, `key2`; the validator rejects a deck missing `mesh { file }`, rejects `crc`/`crank-nicolson`/`galerkin`, rejects an unknown BC key and a dimensionless value where a unit is required | INC-493, INC-494, INC-495; s7#12, s7#13, s7#14 | L | [ ] |
| T-74 | 3-D homology fixture | a TET4 solid torus yields exactly one generator under all three algorithms with no abort (the constructed dead-face cascade of DR-29 would surface here); a two-conductor periodic mesh whose period map sends set A's front face to set B's back face ties every duplicate across sets | DR-29, INC-029; s9#4, s1#6 | L | [ ] |
| T-75 | new suite `physicsmpi` / 2 | building a metal ρ database at 2 ranks completes and agrees with the serial build to 1e-6 relative (not bit-identical) | INC-309, INC-311, INC-529; s6#27 | L | [ ] |

### 3.5 Enabling fixtures

Each unlocks a cluster above. Their feasibility is the ranking's largest unknown; the first
session on each should start by proving the fixture can be built at all.

| id | fixture | unlocks | open question |
|---|---|---|---|
| F1 | **Tier 2 FEM fixture**: a tensor mesh, enriched, partitioned and distributed, with a `Kernel` and a `DofManager`, registered through `TESTRANKS 2 4` as a new suite (`femmpi`) | T-61 – T-65, DR-113 gate (a), DR-105 (2) | whether the `Kernel` can partition an in-memory mesh (production reads from file), and whether a METIS partition of a 2×2×2 hex fixture is deterministic on the two-core runner; the deleted test-hardening plan's O5 recorded `mesh::Partitioner` as METIS-only with no injection point |
| F2 | **Serial Controller fixture**: a `Controller` on the `test_DofSeeding` Kernel with a linear IWG and a stub second kernel, or the accept/reject and guard predicates extracted into pure functions taking `(ε_ref, ε_trial, ω, backtracks)` | T-66 – T-68, T-54 in part, headroom gates | the constructor needs live Kernels and a solver today; the seam route is cheaper and Martin's "timely" point applies: code written without a test resists one |
| F3 | **BfmFile round-trip fixture**: save → load → save on a tensor mesh, dataset-by-dataset compare | T-06, T-52 | none known; `BfmFile` is constructible from a `Mesh` |
| F4 | **Serial ThinShellFactory/CutFactory fixture** on a two-tape mesh without a `MaxwellFactory` or an `input.conf` | T-72, T-02 in part | whether `ThinShellFactory` can run without the deck-driven material map (`TS_TestStack` bypasses it) |
| F5 | **IWG element-assembly fixture**: `EF_TestVolume` extended with a Maxwell IWG and a database material | T-70, T-24 in part | whether `IWG_Maxwell` can be built on the hand-made TRI3 mesh with a database material |
| F6 | **Python suite** for `python/belfem_conf`, registered with ctest as `python -m unittest` | T-73 | whether `belfem_conf` is release tooling or a one-off audit aid (O8) |
| F7 | **Serial Postprocessor fixture**: two-block strip with duplicated interface nodes | T-71 | the SPR base carries MPI collectives; a serial write-order seam may be needed |

---

## 3a. State after the Tier 1 run (2026-09-22)

The triage `tmp/devlog_sweep/triage.md` (ephemeral) crosswalked all 365 Astra rows: 216 matched
the other voices, 36 were excluded with a reason, 113 became the 32 mechanisms N-17 … N-48.
Tier 1 (22 rows, three of them Astra-only: N-17, N-18, N-19) is done except T-08 (carried,
`todo/test_gap_followups.md`). What the run changed for the Astra rows still queued:

- **N-12 (Tier 3, PENTA normal dispatch)** — largely covered: T-60 pins the ghost normals as unit,
  outward and opposite through the public facet seam and the side-face surface increment on all
  three quads (`041c1dd`); the side-face SIGN and a curved side face stay open, as does Q2 (2-D
  QUAD4TS ghost master normals inward on the clockwise winding).
- **N-42 (Tier 3, cut on a Dirichlet boundary keeps the imposed jump)** — answered by reading in
  the DR-172 consumer sweep: the cut duplicate's alias fix is deliberately inert
  (`phi_dup = phi_org + Σ I_k`), so the jump survives; a `gantry`/`costheta` run counting
  multi-source hanging nodes per Dirichlet sideset is the gate. Follow-up R10 writes the
  narrow refusal and the comment.
- **N-21 (Tier 2, default contracts)** — partly ruled: DR-171 settles the material-deck
  contradictions (refuse a contradiction, allow agreeing redundancy); the `eta`/gauge/penalty
  defaults are untouched.
- **N-20 (Tier 2, admission that silently selects another problem)** — the T-08 warning design
  now governs the parser half; the sideset/key half is unchanged.
- **N-07 (Tier 2, integer validation)** — unchanged; the parser warnings (T-08 v3) do not cover it.

Next: Tier 2 (48 rows, 13 Astra-only) through the same loop, consumer inventory before any fix.
No Astra finding is outside a tier; the 36 excluded rows keep their reasons in the triage §5.

## 4. For Gregory (closed core)

Read-and-report only; no BELFEM-side session starts these.

- L18: `updatekGeneratorsFromHomology` discards `s` from `smithForm`; the refusal form (`s < t`)
  lives inside `cl_Cohomology`. The pairing form is testable outside (T-44 could carry it once
  the refusal exists). Reported 2026-09-19.
- DR-29 step 0: the debug invariant `mD·mD = 0` and live↔live row/column agreement need access
  inside `cl_Cohomology`/`fn_Smith`. T-74 is the public-interface proxy.
- DR-71 (2): the unbounded `++it2` in `pGeneralizedCombine`; a test would hang, not fail.
- `SimplicialComplex` has no guard against a complex built from a mesh without connectivity
  (INC-437); T-44's Homology/Cohomology count agreement is the outside proxy.
- `pCoreduce` pairs on `size() == 1` without the ±1 coefficient check (INC-477).
- The `Matrix<real>` → `Matrix<int>` change in the Smith transform turned floor into truncation
  for negative quotients (dl20260902_greedy_cleaner_removal_trace.md); which quotient is intended
  is his ruling.

---

## 5. Rulings needed (logged, not decided)

- **O1** `Ascii` relative paths: is `$PWD` the intended contract or the process cwd? T-07's
  polarity follows. The `GhostSwitch` test was patched to absolute paths; the library was not.
- **O2** DR-159 stall exit: is an accepted step at −20 to −30 dB a silent failure or intended?
  T-67's stall assertion follows.
- **O3** Nickel λ refit (L12): fix the coefficient set, or land T-30 `DISABLED_` with the
  handbook value as the target.
- **O4** HEX8TB wall reload twin rule (L15): FIRST twin at reconstruct or LAST twin at
  normalize. T-52's wall assertion follows.
- **O5** `Vector::operator=( expression )` silent resize (L17): add a length assert in the
  expression assignment, or document the resize as the contract. T-53 follows.
- **O6** `compute_edge_directions_thinshell` reads `physical_tag() == 1` as a direction flag
  (L19): the comment or the code is wrong; no test until one is chosen.
- **O7** Thermal system in absolute form drifts by `rtol·T` with zero source (INC-360, parked
  P3): keep parked, or move to the increment form and let T-69's adiabatic assertion gate it.
- **O8** Is `python/belfem_conf` release tooling? If yes, F6/T-73 become the top three of the
  list as one suite, because the code→schema drift check is the only thing that has ever
  caught an Input Contract violation.
- **O9** DR-160 refusal polarity: refuse a foreign-hash `.bfm`, or warn and load. T-06 follows.

---

## 6. Infrastructure, not tests

One line each; these are CI or build-system work and are listed so they are not lost.

- A sanitizer leg in the nightly matrix (no `-fsanitize` in any CMake file; INC-378, DR-116,
  DR-129, L2 are the memory defects found only by hand).
- An MKL leg (INC-464, INC-467, INC-525, T-57's MKL branch are compiled nowhere automatically).
- A configure smoke over the optional TPLs: `USE_NLOPT`, `USE_SUITESPARSE`, `BELFEM_INT64`
  (INC-402, INC-408, INC-240, INC-527).
- Label `sparse` and `physics` `fast` (0.58 s measured, INC-375; `physics` has no label so
  `check-fast` never runs it).
- A configure-time check that every ctest-registered target is a `check` dependency
  (`sparsempi` was registered and absent from the list; dl20260904_parmetis_d1_owner_guard.md).
- A count-parity gate: the release suite registers the same number of tests as the debug suite
  (46 bare `#ifndef NDEBUG` blocks remain in `tests/`, each vanishing silently in release).
- Run `python -m belfem_conf drift` under `make check` (the executable form of T-73).
- `-Wformat-security` on every platform, and a tripwire for unqualified `abs`/`sqrt`/`pow` on
  floating point (`-Wall -Werror` does not catch it).
- Valgrind in CI, or at least a scripted deck-gate hygiene step: delete `*.bfm`,
  `memdump.hdf5`, `iv_results.csv` and assert "Creating cuts" in the log (DR-160 interim rule).

---

## 7. Definition-of-Done for this backlog

- [x] Every proposal cites a catalog row, a devlog or a register row that exists (checked by script).
- [x] Every "red today" claim at the top of the list was re-read in source by the reconciler,
      or is marked as the reviewer's reading.
- [x] Cross-reviewer duplicates merged by mechanism.
- [ ] O1–O9 ruled.
- [ ] L1–L6 filed or fixed.
- [ ] Each fixture F1–F7 has a feasibility session recorded before its cluster is scheduled.

---

## 8. Audit trail

- Brief: `tmp/ai_exchange/test_gap_sweep_prompt.md`.
- Slice answers: `tmp/ai_exchange/test_gap_sweep_claude_slice_{1_2,3_4,5_6,7_9,8}.md`.
- Web pre-sweeps: `tmp/devlog_sweep_claude/fable_pre{1,2}.md` (Claude Fable via web API; the
  second pass had tree access at `ef6d47f2` and withdrew two rows of the first).
- Reviewers' own doubts about the record (catalog rows that contradict the tree or themselves)
  are in each answer's §6.5; the ones that matter for this list are folded into §2 and §5.
- What each reviewer caught that the others did not: the slice-4 reviewer found L1 and L2 alive;
  slice 7 found L3 and the unrun `drift.py`; slice 1 found L4 and L11; slice 6 found L5 and L7;
  slice 8 found L8, L13, L14; slice 3 found L10; slice 5 found L9; the web pre-sweep found the
  `reorient_generators` sign gap (T-23) and the Blaze `inv` gap (T-12) before the slices did.
