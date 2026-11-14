# tapestack3d: severe J/Jc oscillation with `edge coating : on` — candidate mechanisms

**Date:** 2026-09-10
**Purpose:** Read-only diagnosis of Christian's report that switching the copper side-coating
walls on in `examples/tapestack3d` produces severe oscillations in the YBCO J/Jc field, smooth
with the walls off. Claude pass first, then a blind Codex pass (gpt-6-astra, high) on the same
brief. No file under `src/` or `examples/` modified, nothing compiled or run.
**Module:** fem/kernel (ThinShellFactory, MaxwellData), fem/maxwell (wall kernels, postprocessor)
**Exchange:** `tmp/ai_exchange/sidecoating_jjc_oscillation.md` (ephemeral; this devlog is the
distillation)
**Status:** reviewed, not verified. The pattern (along z, across the width, in time) was not
available at the time of writing; the discriminators below are written so that one look at the
output settles which candidate applies.

## One confirmed defect (Codex F1, source-verified by Claude; confidence high on the defect,
medium on causality)

The buffer (MgO, `DomainType::Buffer`) layer gets a copper wall like every other layer
(`cl_ThinShellFactory.cpp:562`, `cl_FEM_DofMgr_BlockData.cpp:512`). The wall's field recovery
`MaxwellData::compute_h_side_connector` (`cl_FEM_Calculator.cpp:662-690`) evaluates the binormal
component as `Em * mReferenceCalc->q()` with `Em` the PENTA6TS **edge** basis (3×6,
`cl_EF_PENTA6TS.cpp:25`). For a buffer master, `q()` (`cl_FEM_Calculator.cpp:2995`) returns the
element's six **nodal φ** values, because Buffer blocks use the Air dof table
(`cl_Maxwell_FieldList.cpp:281-286`). The sizes match by accident, so nothing asserts. The result
is gauge-dependent (a constant φ offset changes it) and enters the wall's copper
ρ(T, |B|, β) through `compute_rho_metal` → `compute_b` (`cl_FEM_Calculator.hpp:3129`); copper
carries Kohler magnetoresistance with its `rho` field-dependency flags intact
(`cl_Material_Copper.hpp:34-48`, `dl20260903_copper_lambda_field_flag.md`). The normal component
`compute_hn` dispatches per side (φ region vs. h conductor) and is not affected. The postprocessed
wall H/B cannot show this: `copy_seam_fields` overwrites them with tape values
(`cl_MaxwellPostprocessor.cpp:656`).

Consequence: the buffer-level wall's resistivity varies with the local φ level along the tape.
That wall segment is the only through-thickness path between the YBCO half and the Hastelloy half
of each tape, so a station-to-station texture in its ρ is a station-to-station texture in the
inter-half coupling current. Caveat against over-reading: the defect is not September-specific
(buffer interfaces never carried duplicates, so it was live in the 2026-08-15 run that showed only
0.06–0.7 % YBCO ripple with the walls on — at tolerance 1e-11), and the buffer-level wall's ζ
resistance is small (0.2 µm tall) unless Kohler drives ρ far up.

Not fixed (read-only). The fix is local: for a Buffer master, build the binormal component from
the nodal gradient of φ (or skip it, the buffer carries no in-plane h of its own) instead of the
edge basis; a `BELFEM_ASSERT` on the master's dof type would have caught the size coincidence.

## Ranked candidates (Claude + Codex, reconciled)

| # | mechanism | signature | cheapest check | conf. |
|---|---|---|---|---|
| 1 | **Loose tolerance unmasks Messe 2023 §4 checkerboarding.** Deck runs `tolerance : 1e-7`; Messe et al. 2023 need 1e-11 to stop "current density per element oscillating between ~0 and ±2·jc in adjacent elements" (`messe2023.txt:556-563`). The walls add stiff copper terms at every rim station, so the relative residual can pass 1e-7 while the YBCO front has not converged. Codex adds: the wall's Newton tangent freezes the master-field dependence of its ρ (`mt_maxwell_h.cpp:74,105`), a convergence handicap that appears only with the walls on. | element-alternating J/Jc at the flux front; more iterations, relaxation cuts, Δt cuts with coating on | rerun the onset window with **only** `tolerance : 1e-11`; log prints 10·log10(ε), −70 dB = 1e-7, −110 dB = 1e-11 | medium (both) |
| 2 | **Buffer-master field recovery defect** (above) | texture in inter-half coupling current and YBCO rim J/Jc, pinned to stations, gauge-sensitive | plot per-element `rho` on the coating blocks across the MgO level; in-memory φ-offset probe must leave the physics unchanged | high defect / medium cause |
| 3 | **Rim redistribution through free outer traces.** The wall's outer edges are free unknowns with no source (`cl_ThinShellFactory.hpp:285-310`, both fuse flags false by the 2026-08-13 ruling), so the wall's outer face is an open boundary: all planes' rim ξ-edges are tied together through per-station copper resistance and a costless surface sheet. The YBCO rim edges (planes 2/3) are thereby coupled to the metal planes. Codex adds the real slab ratios: (w/d)² = 10⁴ at the MgO level, ≈156 at YBCO, not the w = d of the Sept 5 note. HEX8TB carries no longitudinal component (`cl_EF_HEX8TB.cpp:240`), so it smooths nothing along z. | standing z-mesh texture near the width edges, growing with current, surviving tighter tolerance, changing under refinement | J/Jc and wall transverse current vs. z at fixed width positions across saves; halve `tapeResolutionTip` for a short onset run | medium (modeling consequence, not a defect) |
| 4 | **Ghost OFF + coating ON never exercised as a combination.** Walls were built under ghost ON (default until 2026-09-01). Code handles `hasDuplicates == false` (`cl_ThinShellFactory.hpp:355`, `.cpp:572-580`); no September devlog shows a run. Codex found HEX8TB primitive tests only (`test_InterfaceOrientation.cpp:533`), none of the factory-built seven-layer φ-buffer path. | none specific | add `nitsche ghost penalty { eta : 4 ; k_reg : 1e-3 Ohm ; }` with coating on (the August configuration) | risk medium, mechanism low |
| 5 | **Orientation slip on one connector sign** — the wall reference's pending handedness gates (`side_coating_wall_element.md:435-440`). | station-alternating sign in the wall's **native** J = C·q, one side only | plot J on the coating blocks (not seam-copied H); Codex's proposal: constant-field patch on factory-built walls, signed exact edge integrals must give Cq = 0 on both signs | low |
| 6 | **Temporal ringing, BDF5 + adaptive Δt + Anderson.** Free-outer-trace subspace τ ≈ μ₀w²/ρ ≈ 0.3 µs; BDF is stiffly damped there. Δt adaptation targets iteration count, not temporal error (`cl_FEM_Controller.cpp:2904`). | alternating J/Jc at fixed elements across accepted steps, correlated with Δt changes | save every accepted step around onset; repeat with `method : bdf1` at fixed tolerances | low (both) |

Ruled out statically (Claude, Codex agrees): wall nodes are decoupled duplicates, so nodal SPR
does not mix wall J into YBCO J/Jc; wall targets exclude `JJC*` and seam copy excludes J
(`cl_MaxwellPostprocessor.cpp:117,612`); the buffer-level wall's inner edges hang on the buffer
φ nodes with the same sign convention as their source edge (`cl_FEM_DofMgr_DofData.cpp:3615-3637`);
the wall's M and K are a positive Gram and a curl stiffness with only the common-mode null
vector, so no intrinsic checkerboard mode (Sept 5 devlog, re-read against source). Codex's
caveat stands: none of that clears the coupled system, only the wall in isolation.

## Reading the picture

Element-alternating at the flux front → 1. Station-pinned texture at the rims, surviving 1e-11
→ 2 or 3 (φ-offset probe separates them). One-sided station-alternating wall J → 5. Step-to-step
at fixed elements → 6. First two runs to make: tolerance 1e-11 alone, then ghost ON alone.

## Where Claude and Codex differed

Claude ranked the tolerance first and had noted the buffer-master `q()` question in passing
without verifying it; Codex ran the check and rated the defect high, the causality medium. Codex
also corrected the Sept 5 note's w = d scaling with the true slab ratios. No disagreement on
direction; the reconciled ranking above puts the defect second because it predates the symptom.
