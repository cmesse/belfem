# 0.9.3 pre-check fixes: rho cache revision, hanging-dof double free, relative μ0 band

**Date:** 2026-10-04
**Purpose:** Close the three code findings that Christian ruled must land before the 0.9.3 tag
(`devlog/dl20261004_release_precheck_093.md`): a cached resistivity table that survives a change of
the law it samples, a double free on the `create_field` refusal path, and an absolute μ = μ0 band
that refuses a rounded μ0.
**Module:** `src/physics/materials`, `src/fem/kernel` (DofData), `src/fem/maxwell` (MaxwellFactory); tests in `tests/physics`, `tests/fem`
**AIs involved:** Claude Opus 5.5 (plan, implementation); Codex `gpt-6-astra`/high + Grok `grok-4.7`/high (plan and code juries; the critical-code row, because DofData and the Maxwell factory are on it)
**Status:** COMPLETE — R1–R7 landed 2026-10-04 and committed the same day; code jury closed (Codex astra/high on the pre-P3 diff, Grok 4.7/high re-run after a credit top-up on the post-P3 diff; no P0/P1, three P3s fixed); gates: `make check` 17/17 debug/reference-BLAS and 19/19 release/MKL/gasmodels in private trees, positive control red on HEAD (SIGSEGV), valgrind 0 errors on every new test. CHANGELOG and devlog written; CITATION with the tag; R4b filed as test_gap_followups R13

> **Scope guards:**
> - Out of scope: the alloy's 4 K cp grid, the exact-equality μ0 test in `cl_FEM_Calculator.cpp:100-101`,
>   the absolute band on alloy components (`cl_Material_Alloy.cpp:261`), the multi-source Dirichlet
>   refusal (ruled: keep), the jc-function database caches.
> - No `input.conf` key changes; the deck contract is untouched.
> - The CHANGELOG / CITATION edits are a separate record step after the code jury.

---

## 1. Current behavior and how it fails

| Failure | Mechanism | Evidence |
|---|---|---|
| Silent stale result — a warm run keeps the pre-fix Kohler law | the cache is named `<label>_RRR<int>.hdf5`, and "current" means only that the file has an `RRR` dataset, a marker already written by 0.9.2; the loaded table replaces `rho_kohler` as the field path | `fn_rho_database_is_current.hpp:47`, `cl_Material_Metal.cpp:893-940`, `cl_Material_Alloy.cpp:516` |
| Double free in throw mode | `collect_hanging_dofs` copies pointers of `mDOFs` into `mHangingDOFs` and can refuse before the two containers are separated; the field owner then runs `DofData::reset`, which deletes both | `cl_FEM_DofMgr_DofData.cpp:1957` (alias), `:1959` → `:3333` etc. (refusals), `:2014` / `:2336` (split), `:78-88` (reset); `cl_FEM_Kernel.cpp:596` (owner) |
| Hard refusal of a rounded μ0 | `abs( mu - mu0 ) < BELFEM_EPSILON` is an absolute band of about 2.2e-15 H/m ( ten machine epsilons, `typedefs.hpp:89` ); `1.256637e-6` misses by 6.2e-14 | `cl_MaxwellFactory.cpp:2924-2925` |

**Bottom line:** the cache needs a content revision, `DofData::reset` must delete each `Dof` exactly
once whatever state the setup left, and the μ0 test must be relative.

## 2. Mechanism

**A — rho cache revision.** A `constexpr uint gRhoDatabaseRevision = 1` in
`fn_rho_database_is_current.hpp`, documented as "bump whenever a law sampled into the table changes"
(kohler, rho_0T, the grid). Both savers (`Metal::save_rho_database` `cl_Material_Metal.cpp:1035`,
`Alloy::save_rho_database` `cl_Material_Alloy.cpp:770`) write a `revision` dataset with that value.
`rho_database_is_current` returns true only when `RRR` exists **and** `revision` exists as a
one-element dataset, reads as an unsigned integer through the raw HDF5 API ( `H5Dopen2`, extent check,
`H5Dread` with `H5T_NATIVE_UINT` ), and equals the constant; every probe error, a missing dataset, a
wrong extent or a different value means rebuild (the recoverable outcome, as today). The probe never
throws ( a rank-0 throw before `broadcast( tUsable )` would hang the other ranks ) and closes every
identifier on every return. The saver uses `save_data( const uint & )` ( `cl_HDF5.hpp:209-210` ), the
portable scalar overload. The header comment ( `fn_rho_database_is_current.hpp:20-32` ), the two
call-site comments ( `cl_Material_Metal.cpp:897-901`, `cl_Material_Alloy.cpp:508-511` ) and
`src/physics/materials/doc/alloy_homogenization.md:138` change in the same edit. The rebuild warning text changes from "( no RRR entry )"
to name the revision. Every 0.9.2 and pre-fix cache rebuilds once, for every material that builds a table ( a metal
with angle-dependent Kohler and tables on, `cl_Material_Metal.cpp:286-288`; an alloy with tables on,
`cl_Material_Alloy.cpp:240-242` ), copper included. The master-probes-and-
broadcasts structure is unchanged, so the decision stays rank-uniform. Alternative rejected: hashing
the sampled law (heavy, and a hash of floating-point output is not stable across backends).

**B — delete each Dof once.** `DofData::reset` (`cl_FEM_DofMgr_DofData.cpp:71-95`) first unflags every
non-null hanging dof, flags every non-null dof in `mDOFs`, deletes the non-null hanging dofs that are
**not** flagged, then deletes `mDOFs` ( `delete nullptr` is defined ). The order is load-bearing, and
null slots are skipped by every dereference: `create_dofs` pre-fills `mDOFs` with nullptr
( `:399` ) and an allocation failure in its fill loop unwinds through `~DofData` ( Codex, Grok ). On the normal path the containers are disjoint and nothing changes; on a refusal
inside the alias window the shared pointers are deleted once, through `mDOFs`. This uses the
`graph::Vertex` flag that `Dof` inherits (`cl_FEM_Dof.hpp:27`), the house idiom for marking
(flags over maps). Alternative rejected: a member "aliased" bool set and cleared around the window,
which would have to be maintained at three exits (serial return, parallel master, workers) and is
silently wrong if a future edit adds a fourth.

**C — relative μ0 band.** `std::abs( mu / constant::mu0 - 1.0 ) < tMu0RelativeTolerance` with a
named `constexpr` of `1e-6` local to the factory check ( no `g` prefix: it is not a framework
global ). It is a **rounding allowance**, not a physics budget: it admits μ0 written to seven
significant digits ( `1.256637e-6`, relative 4.9e-8; `4π·1e-7`, relative 5.4e-10 ) and refuses
six-digit `1.25664e-6` ( 2.3e-6 ), five-digit `1.2566e-6` ( 2.95e-5 ) and the susceptibilities of
real conductors ( |χ| ~ 1e-5 ). A value inside the band does not take the kernel's exact-μ0
shortcut ( `cl_FEM_Calculator.cpp:100-101` ): the bulk uses the deck's constant through
`compute_mu_const`, while the shell recovery keeps assuming μ0 — an error of the order of the
relative deviation, which the band bounds at 1e-6. Open: O1.

## 3. Steps

- [x] **R1** A: constant + saver writes + reader check + warning text (Metal and Alloy).
- [x] **R2** A test (`tests/physics/test_MaterialFactoryDeck.cpp`, beside `ARhoDatabaseForAnotherRRRIsRefused`):
  build `Pb40Sn60` RRR 5 with tables in a `WorkDir`; `rho_database_is_current` is true; overwrite
  the file's `revision` with `gRhoDatabaseRevision + 1` (raw HDF5) → false; delete the `revision`
  dataset → false; a revision dataset of the wrong extent → false; a second build rebuilds and leaves a
  current file. The Metal saver is covered separately ( the alloy builds its components with tables
  off, `cl_Material_Alloy.cpp:293` ): `lead` RRR 100 with tables writes a current file ( Codex ).
- [x] **R3** B: the flag pass in `DofData::reset`.
- [x] **R4** B unit test — landed in `tests/fem/test_DofSelfContainment.cpp` (`DofDataReset.*`, on its `HangingStrip` field) rather than in `test_DirichletPin.cpp`, whose subject is Dirichlet pins. Planned as: (`tests/fem`, on the `GluedStrip` field of `test_DirichletPin.cpp`): a
  `dofmgr::DofData` constructed on that field, populated with a **counting `Dof` subclass**
  ( `~Vertex` is virtual ) in four shapes — one dof in both containers, disjoint dofs, null slots in
  both containers, and dofs entering with their flag already set — then destroyed; the oracle is
  "each counting dof destroyed exactly once", which also catches a leak.
- [x] ~~**R4b**~~ filed as `todo/test_gap_followups.md` R13 (O2). Full-path refusal regression: **not reachable with the existing node fixtures**. Probed
  2026-10-04 on unfixed HEAD in the debug tree: a node chained on a hanging node ( single and
  two-source ) never reaches the `create_dofwise_t_matrices_master` refusals, because
  `Mesh::finalize` flattens node source chains ( `cl_Mesh.cpp:2012` ). The refusals inside the window
  are on hanging **edges** ( weights other than ±1, `:3214`, `:3221`, `:3333` ), which needs an h-curl
  fixture with a hanging edge. Decision for Christian ( O2 ). ~~Build now~~ — **RESOLVED 2026-10-04 → filed as a test-gap row**, not built for 0.9.3 (Christian).
- [x] **R5** C: relative band, named constant, message unchanged.
- [x] **R6** C test — no fixture reaches the factory check, so the predicate became `MaxwellFactory::mu_is_mu0` and `Mu0Band.*` in `tests/fem/test_AdmissionRefusals.cpp` tests it. Planned as: if a cheap fixture exists for the factory check, μ = 1.256637e-6 passes and
  1.2566e-6 refuses; otherwise a unit test of the predicate after extracting it to a static helper
  (decided at implementation, reported in the code round).
- [x] **R7** Gate ( required, not optional ): `make check` in the private release/Armadillo tree and in
  the non-MKL debug tree (`/opt/scls/debug`, built 2026-10-04); valgrind on the R4 test ( and R4b if
  approved — it was not, O2 ) in the debug tree, with a positive control: the same test on unfixed HEAD must show the
  invalid free.
- [x] **R8** Code jury on an own-files diff (done); CHANGELOG (done, swept); devlog (done); commit (done, 2026-10-04); `CITATION.cff` waits for the tag date.

## 4. What could break

- A: a rank that is not rank 0 never reads the file (unchanged). A cache directory shared between a
  0.9.2 and a 0.9.3 binary does **not** ping-pong: 0.9.3 rebuilds the old file once, and 0.9.2 then
  loads the new file silently, because it only checks `RRR` ( Codex, Grok ). Old binaries are out of
  reach; the CHANGELOG says so. The existing
  `ARhoDatabaseForAnotherRRRIsRefused` still loads a current cache and must still see "Found RRR".
- B: `reset` is also called at the start of `set_equation` via `DofManager::reset` and from the
  destructor; flags on dofs carry no meaning after reset. A dof that is in `mHangingDOFs` twice would
  still be deleted twice — not a shape any path produces (confidence medium).
- C: a deck with a deliberately non-μ0 constant within 1e-6 now passes the refusal; the bulk uses that
  constant, the shell recovery assumes μ0 ( error ≤ 1e-6 relative ). Alloy components keep the
  absolute band ( `cl_Material_Alloy.cpp:261`, out of scope ), so a rounded μ0 on a user-defined
  alloy component can still be refused there.

## 5. Open questions

- **O1** Tolerance value for C: 1e-6 (proposed, seven-digit rounding) vs. 1e-5 ( admits six-digit
  rounding, and also |χ| up to 1e-5 ). Design value, Christian's. **RESOLVED 2026-10-04 → 1e-6** (Christian).
- **O2** R4b: build an h-curl hanging-edge fixture for the full-path refusal regression now ( more test
  code, a day-scale fixture risk ), or land R4 + valgrind now and file R4b as a test-gap row.

## 6. Exchange

`tmp/ai_exchange/review_release093_fixes.md` (plan and code rounds).
