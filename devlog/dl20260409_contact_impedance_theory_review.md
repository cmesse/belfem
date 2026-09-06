# Devlog 2026-04-09 — Contact Impedance Theory Review

**Date:** 2026-04-09
**Topic:** Read-only review of `src/fem/maxwell/doc/contact_impedance_theory.md`
**AIs involved:** Codex
**Claude Confidence:** N/A
**Codex Audit Confidence:** high
**Literature References:** Schnaubelt et al. 2023 §III-IV, Alves et al. 2024, Juntunen & Stenberg 2009

## Summary

Reviewed the contact-impedance theory note against the current Maxwell/thin-shell implementation, local literature, and online paper metadata. The thin-resistive-interface idea remains plausible, but the note overstates how close BELFEM already is to implementing it and still carries one stale comparison against the current `h_ghost()` kernel.

## Key Findings

- `src/fem/maxwell/doc/contact_impedance_theory.md:122-127` says the current `h_ghost()` path needs `alpha > max(rho/h)` and therefore blows up at the rint/Ag contrast, but `src/fem/maxwell/matrices/mt_maxwell_h.cpp:1811-1849` now uses a regularized harmonic-mean penalty that deliberately bounds `alpha`.
- `src/fem/maxwell/doc/contact_impedance_theory.md:208-216` understates the amount of framework plumbing needed. `DomainType::InterfaceTsCond` exists in `src/mesh/en_DomainType.hpp:56-63`, but there is still no parser/domain/topology/IWG/field-list support for such an interface in `src/mesh/en_DomainType.cpp:97-166`, `src/fem/kernel/cl_FEM_Domain.cpp:28-69`, `src/homology/cl_Topology.cpp:270-375`, `src/fem/maxwell/cl_IWG_Maxwell.cpp:558-567`, and `src/fem/maxwell/cl_Maxwell_FieldList.cpp:411-428`.
- The proposed "store rint material and thickness on the interface sideset" path is not present today. Maxwell material assignment is block-based in `src/fem/maxwell/cl_MaxwellFactory.cpp:2226-2440`, while thin-shell material/thickness metadata is currently indexed one-to-one with existing layer blocks in `src/mesh/cl_ThinShell.hpp:45-50`, `src/mesh/cl_ThinShell.hpp:80-89`, and `src/mesh/cl_ThinShellFactory.cpp:109-113`.
- The claim that no consistency / adjoint terms are needed is plausible but still not fully proven for BELFEM's H(curl) interface context. The note states it as settled fact at `src/fem/maxwell/doc/contact_impedance_theory.md:99-100` and `src/fem/maxwell/doc/contact_impedance_theory.md:155-161`, while the design note still lists that proof as Milestone C's precondition in `todo/thinshell_hphi_formulation.md:691-702`.

## Changes Made / Proposed

- Added a CODEX audit entry to `todo/ai_exchange.md`.
- Added this devlog entry.
- No source-code changes; read-only review only.

## Open Questions

- Should BELFEM repurpose the dormant `DomainType::InterfaceTsCond` enum or introduce a fresh interface type with explicit material/thickness attachment semantics?
- Does the intended collapsed-rint implementation keep the existing duplicated conductor traces explicitly enough in the note, or should the document state that requirement up front?
- Is there a short derivation note that can settle the "no consistency terms needed" claim specifically for the H(curl) conductor-conductor interface, rather than by analogy to general Robin/Nitsche papers?

## Files Updated

- todo/ai_exchange.md
- devlog/dl20260409_contact_impedance_theory_review.md
