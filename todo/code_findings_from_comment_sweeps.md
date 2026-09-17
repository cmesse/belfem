# Code findings surfaced by the comment sweeps

**Date:** 2026-09-16
**Purpose:** Defects and dead declarations that the comment sweeps' auditors found while reading for contracts. None is a comment; each needs the plan + audit → code + audit jury loop before it is touched, so they are parked here for a session that focuses on them. Christian's ruling 2026-09-16: behavior-altering changes are not made in a comment-cleaning session.
**Module:** `src/fem/postproc`, `src/fem/iwg`, `src/sparse`, `src/mesh`, `tests`
**AIs involved:** found by Codex, Grok and a Claude subagent during sweep 2's second-tier round; each verified by Claude against the cited lines
**Status:** OPEN — parked; nothing applied.

## 1. Findings

- [ ] **F1 — `fn_Mesh_integrate_scalar_over_sidesets.cpp:167-171` broadcasts the master rank, not the value.** The integral is computed inside `if ( rank == aMasterRank )`, then `proc_t tMasterRank = aMasterRank; broadcast( tMasterRank );` and `return aValue;`. Every non-master rank returns its untouched `aValue`. Callers that read the value off rank 0 get zero silently. Behavior-altering fix (`broadcast( aValue, aMasterRank )`), needs a test that checks the value on a non-master rank. Confidence high (lines read). Found by Codex.
- [ ] **F2 — `cl_IWG.hpp:820-821` declares `const Matrix< real > & N( const uint & aIntegrationPoint );` and nothing defines it** (grep over `src/fem/iwg/*.cpp`). A public declaration that cannot link; the working equivalent is `Calculator::N`. Fix: delete the declaration. Trivially safe, but it is code; do it with F3 in one audited commit. Found by the Claude subagent.
- [ ] **F3 — `cl_SolverWrapper.hpp:74-75` `mX`, `mY` are never assigned** anywhere in `src/sparse`; `x()`/`y()` (`:438,446`) dereference null. Either dead (delete members and accessors) or a lost contract (a subclass should set them). Grep for callers of `x()`/`y()` first. Found by the Claude subagent and Grok.
- [ ] **F4 — `IWG::collect_node_coords` (`cl_IWG.cpp:1624`) loops `i <= mNumberOfSpatialDimensions`** and writes `nDim + 1` columns into `aX`. No in-tree caller (the sideset integrator calls the free `collect_node_coords` of `geometrytools.hpp`), but `src/fem/iwg/doc/iwg_usage_guide.md` advertises it, so a user with an `nDim`-column matrix writes out of bounds. Fix: `<` and a bounds assert, or retire the function with the guide entry. Found by the Claude subagent.
- [ ] **F5 — `cl_Element_PENTA6TS.hpp` carries two include guards** (`CL_ELEMENT_PENTA6TS_HPP` and `BELFEM_CL_ELEMENT_PENTA6TS_HPP`), harmless, carried over from sweep 1's O6 pass. Obvious fix; can ride with F2.
- [ ] **F6 — `tests/` was outside both sweeps:** `tests/physics/backendfree/test_MaterialBackendFree.cpp:18` cites a date, `tests/fem/test_EdgeFunctions.cpp:32` cites `tmp/tet10/tet4_circulation_probe.py`. Comment-only; a small batch of its own after sweep 2.
- [ ] **F8 — abstract nodes filled past the owner.** `Mesh::set_abstract_nodes()` adopts nodes into `mNodes` (`cl_Mesh.cpp:2196-2205`) and nothing in `src/` calls it; `CutFactory` moves freshly allocated nodes (`cl_CutProcessor.cpp:817`) straight into the mesh's abstract-node list (`cl_CutFactory.cpp:439`), and `~Mesh` deletes `mNodes` only. Either those nodes leak, or another owner deletes them; find out, then route them through the adopting setter or document the owner. Found by Codex in the sweep-2 C3 check. Confidence high on the call graph, unknown on the leak.
- [ ] **F9 — chained kernel allocates a placeholder mesh nobody deletes.** On a worker chained onto a parent kernel, `distribute_mesh()` does `mSubMesh = mMesh; mMesh = new Mesh( dim );` (`cl_FEM_Kernel.cpp:697-699`); `~Kernel` deletes `mMesh` only under `mOwnMesh`, which nothing sets, and `mSubMesh` only under `mOwnSubmesh`, set only on the fresh-distribution path (`:674-680`). The empty placeholder leaks. Found by Codex in the sweep-2 C4 check. Confidence high.
- [ ] **F7 — closed cohomology core narration for Gregory** (report only, no AI edit): `cl_Homology.cpp:26` narration, `:45` a commented-out `delete tHomology` fused into a banner, `:82,103,142,147` and `cl_Cohomology.cpp:92,113,151,159` "loop over all dimensions", `cl_SimplicialComplex.cpp:44-100` "delete pointers"/"delete map" ownership facts in narration form; plus sweep 1's §6 TODO list. Found by Grok.

## 2. Steps

- [ ] R1 — F2 + F5 as one obvious, non-behavioral commit with a Codex check (declaration deletion, guard deletion; `make check`).
- [ ] R2 — F1 with a reproducer on 2 ranks, plan + jury, then code + jury.
- [ ] R3 — F3 and F4 after a caller grep decides delete vs. contract; jury.
- [ ] R4 — F6 comment batch; F7 handed to Gregory.
