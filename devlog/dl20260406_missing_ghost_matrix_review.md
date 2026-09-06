# Devlog 2026-04-06 — Missing Ghost Matrix Review

**Date:** 2026-04-06
**Topic:** Read-only investigation of missing `K++` in ghost calculator matrix map
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** N/A

## Summary

Investigated the `Key K++ not found in map` abort during `h_ghost()` assembly. The failure is caused by initialization order: the FEM sideset calculator allocates its custom matrix map before the sideset is marked as `DomainType::Ghost`, so ghost-only matrices are never created.

## Key Findings

- `SideSet` initializes its calculator during construction (`src/fem/kernel/cl_FEM_SideSet.cpp:70`).
- `Calculator::set_integration_order()` triggers `allocate_memory()` (`src/fem/kernel/cl_FEM_Calculator.cpp:149`), which calls `IWG::create_custom_vectors_and_matrices()` (`src/fem/kernel/cl_FEM_Calculator.cpp:330`).
- `IWG_Maxwell::create_custom_vectors_and_matrices()` creates `K++`, `K+-`, `K-+`, `K--`, `D+`, and `D-` only for `DomainType::Ghost` (`src/fem/maxwell/cl_IWG_Maxwell.cpp:599`).
- FEM groups default to `DomainType::Default` (`src/fem/kernel/cl_FEM_Group.hpp:71`), and `SideSetData::create_sidesets()` does not propagate the mesh sideset domain type during sideset construction (`src/fem/kernel/cl_FEM_DofMgr_SideSetData.cpp:143`).
- MaxwellFactory copies the mesh sideset domain type onto FEM sidesets only later (`src/fem/maxwell/cl_MaxwellFactory.cpp:523`), after the calculator map has already been allocated.
- `h_ghost()` later accesses `K++` and aborts because the key was never inserted (`src/fem/maxwell/matrices/mt_maxwell_h.cpp:2091`).

## Changes Made / Proposed

- No source changes made.
- Appended the investigation to `todo/ai_exchange.md`.
- Proposed likely remedies:
  - propagate sideset domain type before calculator setup, or
  - rebuild calculator custom matrices after domain type propagation.

## Open Questions

- Whether there are any other custom matrix/vector branches in other IWGs that also depend on late domain-type propagation.

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260406_missing_ghost_matrix_review.md
