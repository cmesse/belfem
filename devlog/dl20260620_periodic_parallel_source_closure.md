# Devlog 2026-06-20 — Parallel CORC Resolved: Source-Closure Fixpoint + Buffer/MUMPS Red-Herrings

**Date:** 2026-06-20
**Topic:** The z-periodic CORC reproducer ( cohomology thin cuts + periodic BCs + thin-shell layers incl. a buffer/φ layer ) now runs the full pipeline on 4 MPI procs and tracks the serial result. Closes the parallel-distribution problem; isolates the remaining convergence stagnation as a serial-reproducible buffer/ghost-stabilization formulation issue.
**Module:** src/mesh ( Distributor ), src/fem/kernel ( Kernel ), src/sparse ( SolverMUMPS )
**AIs involved:** Claude ( diagnosis, implementation, probes ), Codex + Grok ( independent audits ).
**Claude Confidence:** high on the distribution fix ( empirical 0-miss + parallel-tracks-serial + tri-AI audit ).
**Codex/Grok Audit Confidence:** high ( both: changes sound for CORC ).
**Literature References:** N/A ( MPI distribution + solver-wrapper mechanics ). MUMPS 5.9 user guide consulted for the INFO sign convention.

Continues `dl20260618_periodic_parallel_distribution.md` ( the facet-87025 and first element-0 crashes ); supersedes that devlog's element-0 root cause ( re-diagnosed below ) and its volume-exchange sizing fix ( reverted ).

---

## The arc ( each stop was a different bug )

The serial reproducer solved earlier; the parallel run hit a chain of distinct failures, each masking the next.

1. **`facet 87025` crash** ( `create_t_matrices` ) — a periodic slave facet hangs on its master facet ( a FACET-typed source ), but the Distributor's owned-facet loop never called `select_sources( tFacet )`, so the master facet was never ghosted. **Fix:** add that call ( `cl_Mesh_Distributor.cpp`, owned-facet loop ). Codex+Grok audited.

2. **`element id 0` crash** ( `compute_element_volumes` ) — **re-diagnosed.** Not a distribution defect ( rank 0's full mesh shows it too ): the **buffer thin-shell layer ( Block 23 ) is tagged `DomainType::Buffer`, not `ThinShell`**, so loop A ( "skip ThinShell" ) does not skip it while loop B ( `tShell->blocks()` ) counts it — a double-count that over-sizes the volume-exchange arrays → trailing zero ids → `element(0)`. Diagnostic confirmed `owned-in-UNTAGGED` exactly matched the mismatch ( 1946/1874/1906/1856 ) and that the buffer block is the *only* Buffer block ( no standalone buffer φ region ). **Fix ( Christian ):** loop A skips `ThinShell || Buffer`. ( The earlier "size arrays from the fill" fix from 0618 was reverted as a band-aid in favour of this root-cause fix. )

3. **`MUMPS -9` then `+8`** — red herrings, see below.

4. **The actual remaining issue:** the Picard iteration stagnates at timestep 4 ( ~−55 dB ) with the buffer + ghost-stabilization config. A **serial ( np=1 ) baseline stagnates identically** — so this is a **formulation/conditioning problem of the buffer φ thin-shell layer, not a distribution bug**. Out of scope for this devlog; serial-reproducible, so debuggable without MPI.

## The real distribution fix ( 5b — source-closure fixpoint )

**Root cause:** in `Distributor::select_entities()` the pre-aura edge/face `select_sources` loops gated on `tElement->is_flagged()`, but no element is flagged until *after* the aura — so they were **dead code**. Only NODE sources were closed; hanging EDGE/FACE/FACET sources were never ghosted. A debug check measured **29 hanging edges on procs 1/3 missing their NODE sources** ( the conductor↔buffer interface edges ).

**Fix:** replace the node-only closure with a **source-closure fixpoint** — call `select_sources` on every set node ( +duplicates ), edge, face, and facet until the summed bitset population stops growing. Delete the dead loops.

**Key correction to the ( tri-AI-audited ) plan:** the plan called for the geometric aura *inside* the fixpoint. That is **wrong** — it grows the halo one element-layer per iteration and walks the whole mesh ( observed: distribution hung ). A source entity only needs to **exist as a basis** ( plus its sub-nodes, which `select_sources` adds ), **not** its element neighborhood. So the geometric aura stays a **single pass before** the fixpoint. Codex + Grok confirmed this is correct for CORC ( a general off-element edge/face source discovered post-aura is a documented limitation, not exercised here ).

**Result:** the debug check reports **0 missing sources on all procs**; parallel CORC runs many timesteps ( 58+ ) and tracks serial.

## MUMPS red-herrings ( and the hardening that came out of them )

With the matrix complete, the parallel MUMPS solve still failed — but for solver reasons, not distribution:

- **`INFO(1) = -9`** ( internal work array too small, on a worker ) — `MemoryRelaxation` ( ICNTL(14) ) was `0`. The Fortran bridge only forwards it when `> 0`, so `0` left MUMPS's own default ( ~20–35% ), too small for the parallel factorization's dynamic pivoting. **Hardening:** force 30%.
- **`INFO(1) = +8`** ( iterative refinement exceeded ICNTL(10) steps — a *warning* ) — the wrapper aborted on `INFO != 0`, treating the warning as fatal. MUMPS convention: `< 0` error, `> 0` warning. **Hardening:** abort only on `< 0`; log `> 0` via `cl_Logger` and continue ( all three INFO-check sites ).

These let MUMPS reach and complete the solve ( raising ICNTL(14) to 50 had it run 3 timesteps converging to −110 dB ). But the **serial baseline** then revealed the stagnation is config/physics ( the buffer ), not the solver or the distribution — so the MUMPS issues were never the cause, just obstacles. ( Christian switched to STRUMPACK; the MUMPS changes are latent hardening for its robustness-fallback role. )

## Audit + cleanup

Tri-AI audit ( Codex + Grok, independent, separate scratch files ) found both changes **sound for CORC**, verified against the MUMPS 5.9 manual and the Fortran bridge. Applied their corrections: accurate ICNTL(14) comment, and the convergence metric now includes the control-point bitset. The temporary `#DIAG-5A` source-closure check was **converted to a debug-only `BELFEM_ERROR` invariant** ( `#if !defined(NDEBUG)` ) — a zero-release-cost permanent guard against a source-closure regression. All `#DIAG-B*` thin-shell probes stripped.

## Changes ( files )

- `src/mesh/cl_Mesh_Distributor.cpp` — `select_sources(tFacet)` for owned periodic facets; source-closure fixpoint ( replaces dead loops + node-only closure ); `#DIAG-5A` → debug invariant.
- `src/fem/kernel/cl_FEM_Kernel.cpp` — loop-A `ThinShell || Buffer` skip ( Christian ); `#DIAG-B*` probes removed.
- `src/sparse/cl_SolverMUMPS.cpp` — ICNTL(14) 0→30; warnings logged ( not fatal ) at all three INFO checks.

## Open / deferred

- **Buffer + ghost-stabilization convergence** ( serial-reproducible ) — the live problem now; a formulation/IWG issue, not distribution. Likely the higher-leverage next target.
- **Latent ( no CORC trigger ):** closing the sources *of* elements / control points ( fixpoint does not iterate them; `select_sources` lacks a `CELL` case — `Element::entity_type()` returns `CELL` ).
- **MUMPS follow-ups:** decide whether `+1` ( out-of-range entries silently dropped ) should escalate to a hard error; the wrapper reads rank-0 local `INFO`, so worker-only *warnings* can be missed ( errors are still caught via the `-1` relay ).
- A formal parallel-vs-serial comparison of a global quantity ( total current / AC loss ) once the buffer convergence is fixed — the clean final check on 5b.

## Acknowledgements

A genuinely joint effort. **Christian ( Prof. Messe )** drove the physics and the key calls — notably insisting that "it solves" ≠ "it's right" ( the unusually slow parallel solve ), which redirected from the MUMPS rabbit hole back to the matrix; the element-0 `Buffer`-skip fix; and the STRUMPACK switch that isolated solver-vs-matrix. **Claude** ran the diagnosis, the `#DIAG` probe instrumentation, and the implementation/audit loop ( and owned the one real implementation bug — geometric aura in the fixpoint — caught by its own timeout and the audit ). **Codex** and **Grok** gave independent, file:line-grounded audits that confirmed soundness, corrected the ICNTL(14) rationale ( verified in the Fortran ), flagged the control-point convergence gap, and checked the MUMPS sign-split against the 5.9 manual. The discipline that kept it honest: measure before fixing, cross-check non-trivial claims, and reverse conclusions when a probe contradicted them ( the "invalid distribution → MUMPS failure" theory was walked back by the serial baseline ).

## Files Updated

- src/mesh/cl_Mesh_Distributor.cpp
- src/fem/kernel/cl_FEM_Kernel.cpp
- src/sparse/cl_SolverMUMPS.cpp
- todo/periodic_parallel_source_closure_fix.md ( plan, checkboxes updated )
- devlog/dl20260620_periodic_parallel_source_closure.md ( this file ) + devlog/README.md
