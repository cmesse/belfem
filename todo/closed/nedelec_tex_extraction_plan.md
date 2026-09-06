# Nédélec LaTeX Notes: Extraction into Module Docs + Drift Record

**Date:** 2026-08-14
**Purpose:** The pre-BELFEM theory notes in `tmp/nedelec/` (LaTeX, Lagrange-multiplier era) contain
a weak-form derivation and a Nédélec-element derivation that are still correct and worth keeping.
Extract those two parts into the module documentation (`src/fem/maxwell/doc/`,
`src/fem/interpolation/doc/`), and record in `tmp/nedelec/drift.md` where the notes have drifted
from the current code — most notably the absence of the cohomology cut machinery and of the
hanging-node/hanging-edge concept.
**Module:** `src/fem/maxwell` (+ `src/fem/interpolation`), source material in `tmp/nedelec/` (separate git repo)
**AIs involved:** Claude (exploration + plan + extraction), Codex (audit), Grok (third voice)
**Status:** ✅ COMPLETE (2026-08-14). Both extractions landed and passed the file-stage
Codex + Grok audit round with all findings applied: `src/fem/interpolation/doc/nedelec_derivation.md`
(now including the full TET10 set with the E4 fix), `src/fem/maxwell/doc/maxwell_weak_forms.md`,
`tmp/nedelec/drift.md`, plus both module README indexes and the stale maxwell README rows.
Everything is **reviewed, not verified** (static + hand algebra, no executable gate).
Residual follow-ups: **D1 (confirmed `EF_TET4::E` edge defect) and D2 (suspected `EF_TET10`
counterpart) → `todo/nedelec_edge_function_defects.md`, awaiting Christian's ruling**; the
exchange thread sweeps after the devlog. Session record: `devlog/dl20260814_nedelec_tex_extraction.md`.

> **Scope guards:**
> - **Documentation only.** No C++ source, no CMake, no `input.conf` keys are touched — the
>   Input Contract rule is not triggered.
> - The LaTeX sources in `tmp/nedelec/` are **not modified** (separate repository); the only file
>   added there is `drift.md`.
> - The `discrete/` chapter (h-a and h-φ discretizations, Lagrange-multiplier interfaces, manual
>   cuts, thin-shell, lumped-mass thermal) is **not extracted** — it is where the drift lives, and
>   its current-state truth is already carried by the module docs and by Gregory's forthcoming
>   paper and thesis (cohomology). It is summarized in `drift.md` only.
> - `methods/` (θ-stepping, Newton/Picard, adaptive relaxation) and `materials/powerlaws.tex` are
>   likewise drift-summarized, not extracted.

---

## 1. Source Inventory and Current Behaviour

The document (`tmp/nedelec/main.tex`, "Finite-Element Discretization of the Quasi-Magnetostatic
Maxwell Equations") predates the moves to static condensation and cohomology cuts. Inventory, with
disposition:

| tex source | Content | Disposition |
|---|---|---|
| `introduction/notation.tex` | N/E/B/C operator concept; 2D and axisymmetric specializations | extract → interpolation doc |
| `triangle/triangle_lagrange.tex` | barycentric coordinates, TRI3/TRI6 Lagrange N and B, geometry Jacobian + the Jᵀ pitfall | extract → interpolation doc |
| `triangle/triangle_nedelec.tex` | Nédélec TRI3/TRI6 (edge + face functions), curl operator `C = 2/detJ [s1 s2 s3]`, edge-sign convention s_k, TET4/TET10 edge/face functions, face ownership, edge/face generation algorithm | extract → interpolation doc |
| `weakform/intro.tex`, `fundamental_lemma.tex` | fundamental lemma of variational calculus; Gauss and Stokes divergence-theorem corollaries | extract → maxwell doc |
| `weakform/leastsquares.tex` | L2 projection of edge fields onto nodes | extract → maxwell doc (still live: `matrices/mt_maxwell_l2_h.cpp`, `mt_maxwell_l2_b.cpp`, `mt_maxwell_l2_phi.cpp`) |
| `weakform/heat.tex` | thermal-conduction warm-up example, ρ₀ mass-conservation remark | extract → maxwell doc (O1) |
| `weakform/maxwell.tex` | Maxwell equations, MQS simplification, b-conform vs h-conform taxonomy | extract → maxwell doc |
| `weakform/bconform.tex` | a-formulation weak form | extract → maxwell doc, marked "derived, not implemented" |
| `weakform/bconform_voltage.tex` | a-v formulation (orphan — not `\input` by `main.tex`) | extract → maxwell doc (O2) |
| `weakform/hconform.tex` | h-conform weak form — the implemented core | extract → maxwell doc |
| `discrete/intro.tex`, `ha.tex`, `hphi.tex` | h-a and h-φ discretizations, λ-multiplier interfaces and symmetry/antisymmetry elements, manual cuts (node doubling + 0-d λ element), Alves thin-shell, lumped-mass thermal | NOT extracted → summarized in `drift.md` |
| `methods/theta.tex`, `newton.tex` | θ-timestepping stability table, Newton/Picard, arctan adaptive relaxation, adaptive Δt | NOT extracted → `drift.md` |
| `materials/powerlaws.tex` | E-J, E-J-T, E-J-B, E-J-B-T power laws | NOT extracted → `drift.md` (code has a superset) |
| `methods/intro.tex` | two-line chapter header | NOT extracted (absorbed; added 2026-08-14, was missing — Codex/Grok A3) |
| `title/{title,nomenclature,intro}.tex` | front matter; `title/intro.tex` is the prose overview of the four chapters | NOT extracted; overview content informs the new docs' intros (added 2026-08-14) |
| `triangle/tri{3,6}.tex`, `*.eps/.ai`, `matlab_curl/*.m`, `trash/*` | TikZ figures, graphics, the notes' own MATLAB curl check, discarded drafts | NOT extracted; `matlab_curl/fragment.m` is a useful independent check for R1 (added 2026-08-14) |

**Reviewed drift anchors** (static review, no executable gate; per-row confidence; reconciled
2026-08-14 against both audits):

| Document says | Code does today | Evidence |
|---|---|---|
| h-φ interface coupled via Lagrange multiplier block `L∥` (`discrete/hphi.tex:73-87`, Eq. `khphiparallel`) | Conductor-air / conductor-ferro interfaces are coupled by hanging-edge condensation; the interface sidesets are hard-errored off as weak-form groups [high] | `src/fem/maxwell/cl_IWG_Maxwell.cpp:311-319, 341-353`; hanging-edge creation in `cl_MaxwellFactory.cpp:1359-1473`; `cl_Maxwell_TMatrix.{hpp,cpp}` |
| Symmetry as saddle-point λ elements (`discrete/ha.tex:163-191`, `discrete/hphi.tex:293-304`) | Penalty-style quadratic form `K += (n×B)ᵀ(n×B)` resp. `(n×E)ᵀ(n×E)`, no λ dofs, for **both** φ- and h-side symmetry [high] | `matrices/mt_maxwell_symmetry.cpp:23-49, 52-79` (`symmetry_phi_2d/3d`), `:83-135` (`h_symmetry_2d/3d`) |
| Antisymmetry as saddle-point λ elements (`discrete/ha.tex:195-223`, `discrete/hphi.tex:306-317`) | **Not** replaced by a penalty form: all anti-symmetry groups are hard-errored off ("must be disabled"). Drift is removal, not substitution [high, Grok A1.2b] | `src/fem/maxwell/cl_IWG_Maxwell.cpp:364-368, 388-397` |
| Manual thin cuts: hand-placed, duplicated φ nodes, 0-d λ element with imposed current I (`discrete/hphi.tex:230-283`, Eq. `cutgov`) | Cohomology generates the cut **automatically** (thick cut), then `CutProcessor` pushes it to a thin cut (node duplication + hanging condensation). Two distinct drifts: (1) automatic generators vs hand-drawn cuts, (2) condensation vs the 0-d λ element. Thick-**then**-thin, not thick-instead-of-thin [high on mechanism; medium that `CutProcessorManual` is off the live path — it is still compiled, `new CutProcessorManual` not found in `src/`] | `src/homology/doc/thick_thin_cuts_and_conjugate_edges.md:11-15, 65-79`; `en_CutAlgorithm.hpp:21-28`; `cl_CutFactory.cpp:113-156, 209-358`; `cl_CutProcessor.cpp:21-104`; `cl_CutProcessorManual.hpp:26-54` |
| No mention of hanging nodes/edges anywhere | Hanging-edge condensation is the central interface mechanism [high] | kernel docs (`src/fem/kernel/doc/`), `cl_Maxwell_TMatrix` |
| h-a formulation fully derived (`discrete/ha.tex`) | Only HPhi is a solving formulation; `L2*` enum values are postprocessing projections [high] | `en_Maxwell_Formulations.hpp:22-29`, maxwell doc README |
| φ-domain via Gauss's law derived alongside Faraday's law — and the notes themselves already call Gauss ill-suited and Faraday "better on all points" (`discrete/hphi.tex:46-58, 93-120`) | Code follows the magnetodynamic Faraday coupling (Arsenault et al. 2023); the 2026 erratum correcting the air-domain coupling postdates the notes. Do **not** write "notes recommend Gauss" — they don't [high, Grok] | `doc/literature_references.md` routing |

**Bottom line:** the weak-form and element-derivation chapters are formulation-agnostic mathematics
that the module docs currently lack (the existing `nedelec.md` documents the *implementation
framework*, not the derivation); the discretization chapter is a historical snapshot of the
Lagrange-multiplier era and must be recorded as drift, not imported as truth.

## 2. Approach

Three new markdown files, house style (ASCII math in code fences, as in
`contact_impedance_theory.md`; `lowercase_with_underscores.md`; header block with Date/Purpose/
Module; facts stated inline, no devlog citations):

1. **`src/fem/interpolation/doc/nedelec_derivation.md`** — the mathematical derivation:
   interpolation operators N/E/B/C (3D, 2D, axisymmetric), barycentric coordinates, Lagrange
   TRI3/TRI6, geometry Jacobian (including the Jᵀ implementation pitfall), Whitney/Nédélec edge
   functions for TRI3/TRI6/TET4/TET10 with the s_k orientation rule, curl operators, face-function
   redundancy (`F1+F2+F3 = 0`) and face ownership, and the sort-unique edge/face generation
   algorithm. Complements (does not duplicate) `nedelec.md`: that file keeps the class/factory/
   lifecycle story; the new file carries the math the code implements. Cross-links both ways.
2. **`src/fem/maxwell/doc/maxwell_weak_forms.md`** — the weak-form primer: fundamental lemma,
   the two divergence-theorem corollaries, L2 projection example, heat-conduction example,
   Maxwell equations + MQS assumption, b-conform (a and a-v) and h-conform weak forms. b-conform
   sections carry an explicit "derived here for completeness; BELFEM solves h-φ only" marker.
3. **`tmp/nedelec/drift.md`** — the drift record (audience: future readers of the LaTeX repo):
   what was extracted where; what the notes lack relative to the code, led by (a) the cohomology
   cut machinery — deferred to Gregory's paper and thesis for the full treatment — and (b) the
   hanging-node/hanging-edge concept that replaced Lagrange-multiplier interface coupling; plus
   the smaller drifts (symmetry penalty form, h-a not implemented, Faraday-vs-Gauss φ domain with
   the Arsenault 2026 erratum, power-law superset, controller evolution, thin-shell evolution);
   and the transcription errata found during extraction (§3.1).

**Prose rules (Christian, 2026-08-14):** the extraction keeps Christian's voice. Codex runs a
sentence-level pass over the extracted prose to fix the German-English patterns ("differ" for
"distinguish", "over" for "via", article slips), under two hard constraints: (1) **no em-dashes**
anywhere in the produced documents; (2) the dry humor and the practitioner asides (the
transposed-Jacobian warning, "Reader's discretion is advised", the edge-ownership discussion) are
preserved, not flattened. The pass is read-only; Claude applies the returned edits.

Rejected alternative: porting the discretization chapter into `maxwell/doc` with HISTORICAL
markers. Rejected because the current-state equivalents are already documented per topic
(interfaces, thin shell, cuts) and a parallel outdated derivation would be a standing citation
hazard.

## 3. Content-Mapping Gap Table

Class (a) = mechanical transcription, (b) = open question, (c) = needs explicit handling
(fix/mark/verify during extraction).

| # | Item | Destination | Class | Notes |
|---|---|---|---|---|
| 1 | Operator concept N/E/B/C + 2D/axisymmetric | interpolation doc | (a) | |
| 2 | Lagrange TRI3/TRI6, Jacobian, Jᵀ pitfall | interpolation doc | (a) | keep the pitfall paragraph — it is the kind of thing only this document records |
| 3 | Nédélec TRI3 edge functions + C operator | interpolation doc | (a) | **reviewed, exact match** (Grok A4, hand algebra): `EF_TRI3.cpp:80-113` implements the notes' Whitney forms verbatim, `:64-68` the `C = 2/detJ [s]` operator. First-order Whitney is already unit-circulation; no scale fight |
| 4 | Nédélec TRI6 / TET4 / TET10 + face functions | interpolation doc | (c) | **three-way reconciliation, not "code wins"** (Grok): (i) TET4 — the notes' edge order (ξ→η, η→ζ, ζ→ξ, …) does not match the code convention (node0↔ξ, node1↔ζ, node2↔η; edges ξ→ζ, ζ→η, η→ξ, ξ→τ, ζ→τ, η→τ per `cl_Element_TET4.hpp:107-142`); the doc must present the code's table. (ii) TET4 `E()` edge 2 carries defect **D1** (see §4.0) — the doc presents the correct Whitney form (which `C` implements) and flags the defect. (iii) TRI6 — notes and `EF_TRI6.cpp:86-119, 311-330` agree with each other, but each parent-edge function has circulation `s/2`, not `δ_jk`, so both contradict the blanket unit-circulation policy in `nedelec.md:141-161`; TRI6/TET10 are tagged unvalidated PoC there (`:201-212`). Record the convention split, do NOT renormalize the formulas. (iv) reconcile `∫E·dl = δ_jk` (global) vs `s_k δ_jk` (local) phrasing — pick one sentence and keep it |
| 5 | Edge/face generation (sort-unique on `d = a + b·n`) | interpolation doc | (c) | mark as the *concept*; actual generation lives in `src/mesh` — cross-ref, don't claim it is the implemented algorithm |
| 6 | Fundamental lemma + divergence theorems | maxwell doc | (a) | |
| 7 | Least-squares projection | maxwell doc | (a) | anchor to the live `mt_maxwell_l2_*` kernels |
| 8 | Heat-conduction example | maxwell doc | (b) | O1 |
| 9 | Maxwell eqs + MQS + taxonomy | maxwell doc | (c) | fix Eq. (ampere) typo (stray `\label` splits `∂d/∂t`) |
| 10 | b-conform weak form (a-formulation) | maxwell doc | (c) | mark not implemented |
| 11 | a-v formulation (orphan file) | maxwell doc | (b) | O2 |
| 12 | h-conform weak form | maxwell doc | (a) | this is the implemented core; anchor to `mt_maxwell_h.cpp` |
| 13 | Drift: cohomology cuts | drift.md | (a) | defer full treatment to Gregory's paper + thesis |
| 14 | Drift: hanging nodes/edges | drift.md | (a) | |
| 15 | Drift: symmetry penalty, h-a unimplemented, Faraday-vs-Gauss + erratum, power-law superset, controller, thin shell | drift.md | (a) | one short section each, with code anchors |
| 16 | README index updates (both module docs) | both READMEs | (c) | while editing `maxwell/doc/README.md`: its Quick Reference still lists `mt_maxwell_interface.{hpp,cpp}` and `mt_maxwell_aphi.{hpp,cpp}`, which no longer exist in `matrices/` — fix the stale rows in the same pass (R4) |

### 3.1 Transcription errata found in the LaTeX source (fix silently in the extraction, list in drift.md)

E1-E5 (Claude, pre-audit) were confirmed by both Codex and Grok on independent hand algebra
[high]. E7-E18 were added by the audit round (Grok's numbering kept; overlapping Codex findings
merged in). "Extracted?" marks whether the erratum lands in a produced document (must be fixed)
or only in `drift.md` (recorded).

| ID | Where | What | Extracted? |
|---|---|---|---|
| E1 | `notation.tex:57-61` | 3D curl matrix (nodal a): sign errors in the **first** column block, rows 2-3 (`+N¹_{,x}` → `−N¹_{,x}`; `−N¹_{,x}` → `+N¹_{,x}`); last block correct | yes |
| E2 | `weakform/maxwell.tex:27` | stray `\label` splits `∂d/∂t` (typo only, math intact) | yes |
| E3 | `weakform/heat.tex:30-39` | `grad T ≈ Bᵀ T̂` spurious transpose; B is (dim × n), so `grad T ≈ B T̂`; the discretized `∫ Bᵀ k B` is the consistent form | yes |
| E4 | `triangle_nedelec.tex:194` | TET10 face function F¹¹: missing `+` before `16ηξ∇ζ` | yes |
| E5 | `discrete/ha.tex:130-134` | 3D `E×n` row 3, **last column only**: `n_x E_y^m − n_x E_y^m` is identically zero; must be `n_y E_x^m − n_x E_y^m`. First column is correct (refined by Grok) | no |
| E6 | throughout | prose typos ("strait edged", "quadliterals", "Lagrangiam multiplyer", …); extraction rewrites prose | yes |
| E7 | `hphi.tex:75-76` | interface functional `Π = ∫ λᵀ n×(h − h) dS` is identically zero as written; meant `h_SC − h_NSC` | no |
| E8 | `hphi.tex:261` | cut element variation carries `+ I` inside the δφ term; `I` is prescribed, `δI = 0`; the K/f actually stated are consistent with the correct variation | no |
| E9 | `ha.tex:18` | non-conducting domain written `Ω_SC` instead of `Ω_NSC` | no |
| E10 | `ha.tex:75` | interface dof vector lists `â_z³` in the node-1 triplet; must be `â_z¹` | no |
| E11 | `bconform.tex:24` | Stokes LHS written `∫_{∂Ω} δaᵀ curl h dV`: surface domain with volume measure; must be `∫_Ω` | **yes** |
| E12 | `notation.tex:29-32` | vector-N matrix row 3: misplaced `\zero`/`\hdots` in the column layout | yes |
| E13 | `notation.tex:102-106` | axisymmetric section says "magnetic **scalar** potential" then derives the vector potential `a_t` | yes |
| E14 | `notation.tex:113-116` | axisymmetric φ-domain current labeled `j_z`, conducting labeled `j_t`; stray `\kma` inside a matrix | yes |
| E15 | `notation.tex:20-22, 94-96` | E-matrix last column indexed `n` in row 1, `m` in the others | yes |
| E16 | `triangle_nedelec.tex:144-153` | notes' TET edge order/labels do not match the code convention (node0↔ξ, node1↔ζ, node2↔η; edges ξ→ζ, ζ→η, η→ξ, ξ→τ, ζ→τ, η→τ, `cl_Element_TET4.hpp:107-142`); the doc must present the code's table, not transcribe the notes' | **yes** |
| E17 | `hconform.tex:3,28`, `hphi.tex:36` | "applied to to"; leftover empty `{array}` environments | yes |
| E18 | `bconform_voltage.tex` | continuity-equation sign and a missing `\diff V` (Codex, medium-high; re-derive during R2 before publishing) | yes (if O2 stays include) |

Non-errata clarifications from the round: `ha.tex` 2D uses `n×E` while 3D uses `E×n` (operand
flip, sign trap, recorded in drift.md); TRI6 face functions in notes and code agree
(`ξ−η−1 = −(η−ξ+1)`, no erratum).

## 4. Ordered Steps

### 4.0 Audit round 1 — plan stage (2026-08-14, Codex + Grok, blind; reconciled by Claude)

Both auditors returned; every A1-A5 item answered. Summary: drift anchors upheld with two
corrections (antisymmetry = removed not replaced; cuts = thick-then-thin, automatic generation is
the drift), E1-E5 all confirmed, thirteen further errata added (E7-E18), inventory closed
(methods/intro, front matter, figures/matlab), TRI3 confirmed an exact Whitney match, and the
"code wins" rule replaced by a three-way reconciliation (gap row 4). Claude independently
re-derived the load-bearing findings before inclusion. Full threads in
`tmp/ai_exchange/nedelec_tex_extraction.md`.

**Defect tracker:**

- [ ] **D1 — `EF_TET4::E()` edge 2 uses `∇ζ` where `∇η` belongs.** Severity: **HIGH**
  (found by Grok as an algebraic mismatch, independently confirmed by Claude 2026-08-14).
  `src/fem/interpolation/nedelec/cl_EF_TET4.cpp:193-195` computes `η∇ξ − ξ∇ζ`; the Whitney form
  for edge 2 (node2→node0, i.e. η→ξ) is `η∇ξ − ξ∇η`, which is what the code's own comment
  (`:192`) and the separately-coded curl operator (`mC(:,2)` at `:122-124`, `= 2∇η×∇ξ`) both
  say. Consequences (hand-derived): the edge-2 basis has circulation 1/2 instead of 1 on its own
  edge and −1/2 instead of 0 on edge 0, so the interpolated `E` operator breaks tangential
  conformity on every TET4 wherever `E` is consumed (mass matrices, L2 projections), while the
  stiffness path through `C` is correct. No test covers edge-function circulation (the
  `todo/falsification_tooling.md` D1 battery, which would catch exactly this, is unimplemented).
  `EF_TET10` (237 Nabla lines) needs the same sweep — follow-up, not done. **Read-only session:
  not fixed. Needs Christian's ruling; the one-line fix is `mNablaZeta` → `mNablaEta` in the
  three edge-2 lines, followed by a rebuild and a 3D bulk-conductor regression run.**
  Confidence: high on the algebra (reviewed, not verified — no numeric probe was run).

### 4.0b Audit round 2 — file stage (2026-08-14, Codex + Grok on the produced files)

**Codex:** every equation block present in the two docs is clean against TeX and code; E1-E3 and
E11-E16 confirmed fixed; D1 warning correctly preserved. Four findings, all applied: (1) E4 was
recorded but not truly fixed because TET10 was only summarized → resolved by extracting the
full TET10 edge/face set with the E4 fix (per Christian's "copy my original text" instruction);
(2) drift.md "fixed in extracted markdown" was overbroad for the `discrete/` errata → reworded;
(3) drift.md listed contact impedance as implemented machinery although it is a theory note
with no kernel → corrected; (4) em-dashes in the weak-forms references violated the prose rule
→ stripped, plus three Germanism pairs applied ("interpolated with", "using the edge field",
"in terms of both potentials"). One suggestion rejected: dropping "volume elements" (that is
BELFEM's own terminology for the non-thin-shell family).

**Grok:** all load-bearing weak-form signs re-derive cleanly (Stokes corollary, `weaka`,
`hstokes`, `weakh`, axisymmetric and 2D curls); E1/E3/E11/E16 landed correctly with no new
slips in the big operator matrices; the D1 callout is accurate and the TRI6 circulation-1/2
integral was independently reproduced; every drift.md claim holds against the tree (cut
pipeline citation upgraded with `CutFactory.cpp:688-692`). Three findings, all applied:
(1) **TET10 was wrongly lumped into the circulation-1/2 convention** — the notes' TET10
polynomials are twice the TRI6 pair and integrate to unit circulation → narrowed to TRI6 in
both files; (2) the `h ≈ B phi_hat` (unsigned) vs `h = -grad phi` convention split between the
two docs → one clarifying note added to each, citing the kernel's deliberate sign drop
(`mt_maxwell_phi.cpp:103-115`); (3) the face-generation paragraph had silently rewritten the
notes' element-pair keying to node triples → original scheme restored. Grok also surfaced
**D2** (below) and the L2-interface quarantine clause.

**Defect tracker (continued):**

- [ ] **D2 — `EF_TET10` edge polynomials suspected of a TET4-style relabeling defect.**
  Severity: **MEDIUM** (suspected, unconfirmed; found by Grok, hand integral only). The
  implemented edge 0 pair (`cl_EF_TET10.cpp:247, 260, 602-604`) appears to be the notes'
  `(xi,eta)` polynomial with `grad eta` rewritten to `grad zeta` but the scalar factor left
  η-based, giving circulations 2/3 and −1 instead of 1 and 0 on the implemented edge 0.
  Confidence medium-high (~80 %) per Grok; nobody has executed a probe. TET10 is tagged
  proof-of-concept and is not on the production path, so this is lower urgency than D1.
  Tracked with D1 in `todo/nedelec_edge_function_defects.md`.

### 4.1 Steps

- [x] **R1** — Write `src/fem/interpolation/doc/nedelec_derivation.md` (items 1-5), verifying
  every element formula against `src/fem/interpolation/nedelec/EF_*.{hpp,cpp}` per the gap-row-4
  three-way rule (notes vs `EF_*` vs `nedelec.md` policy). TET4 section presents the code's edge
  convention (E16) and the correct Whitney forms, and flags D1 with a pointer to this plan.
  `matlab_curl/fragment.m` serves as an independent check on the TRI3 curl. State plainly that
  the doc covers TRI3/TRI6/TET4/TET10 only; LINE3, HEX8 and the TS family stay in `nedelec.md` /
  `nedelec_thinshell.md`.
- [x] **R2** — Write `src/fem/maxwell/doc/maxwell_weak_forms.md` (items 6-12) with
  implementation-status markers and code anchors.
- [x] **R3** — Write `tmp/nedelec/drift.md` (items 13-15, §3.1 errata, extraction map). (after: R1, R2 — the errata and convention-conflict lists feed it)
- [x] **R4** — Update the two module-doc indexes (`src/fem/interpolation/doc/README.md`,
  `src/fem/maxwell/doc/README.md`); fix the stale `mt_maxwell_interface`/`mt_maxwell_aphi` rows
  in the maxwell README (item 16; both files confirmed absent from `matrices/` by Grok).
- [x] **R5** — Codex + Grok audit of the three new files (equation transcription vs LaTeX and vs
  code; drift claims vs tree). Codex additionally runs the prose-polish pass under the §2 prose
  rules (no em-dashes, voice preserved). Apply confirmed findings. *(round 2 record in §4.0b;
  all findings applied 2026-08-14)*
- [x] **R6** — Devlog entry `devlog/dl20260814_nedelec_tex_extraction.md` + `devlog/README.md`
  index line; D1/D2 split into `todo/nedelec_edge_function_defects.md`; `todo/README.md`
  updated.

**Plan-stage audit (this session, before R1):** Codex and Grok audit *this plan* — the
inventory/disposition table, the drift-anchor table, and the errata list.

## 5. Open Design Questions

- **O1** — Keep the heat-conduction warm-up example in `maxwell_weak_forms.md`?
  **RESOLVED 2026-08-14 → keep**, with E3 fixed and a pointer to `src/fem/thermal`. Both
  auditors concurred: it is the only short Gauss-theorem walkthrough and the doc needs a scalar
  rehearsal before Maxwell.
- **O2** — Include the orphaned `bconform_voltage.tex` (a-v formulation)?
  **RESOLVED 2026-08-14 → include**, explicitly marked "not `\input` by `main.tex`, not
  implemented", with E18 re-derived before publishing. Both auditors concurred with the
  condition.
- **O3** — Math notation: house ASCII-in-code-fences vs LaTeX `$...$`.
  **RESOLVED 2026-08-14 → house style** (per `contact_impedance_theory.md`). Both auditors
  concurred; mitigation for the transcription risk Codex flagged: every equation block in the
  new docs keeps a comment naming its TeX source file and equation label, so slips remain
  traceable.
- **O4** — How to cite Gregory's paper/thesis in `drift.md`?
  **RESOLVED 2026-08-14 → descriptively** ("Gregory's forthcoming paper and PhD thesis"), no
  invented bibliographic metadata; for *current* cut mechanics, link `src/homology/doc/`
  (`thick_thin_cuts_and_conjugate_edges.md`). Both auditors concurred.

## 6. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step or an On. *(2026-08-14, post-reconciliation)*
- [x] Every extracted equation checked against the LaTeX source; element formulas additionally
  reviewed against `EF_*` implementations (or the divergence recorded in drift.md). Static
  review only — "reviewed", not "verified" (protocol §11). *(round 2, both auditors)*
- [x] Drift claims each carry a code citation. *(re-checked by Grok round 2, all hold)*
- [x] Both module-doc README indexes updated; stale rows in maxwell README fixed.
- [x] Codex + Grok audits ran on plan and on the produced files; findings reconciled and
  recorded in §4.0/§4.0b and the exchange thread.
- [x] D1 surfaced to Christian for a ruling (fix is out of scope for this docs task; tracked
  with D2 in `todo/nedelec_edge_function_defects.md`).
- [x] Devlog written, plan closed.

## 7. Audit Trail

- Exchange thread: `tmp/ai_exchange/nedelec_tex_extraction.md` (AI-only, ephemeral; distilled
  into this plan and the devlog before sweep).
- Plan-stage round (2026-08-14, both blind): **Codex** confirmed E1-E5, split the antisymmetry
  drift, caught the missing `methods/intro.tex` inventory row, the stale README rows, the
  `bconform.tex` `∂Ω` slip, the axisymmetric scalar/vector mislabel, the `â_z³` slip, and first
  raised the TET convention caveat. **Grok** confirmed E1-E5 with refinement of E5, refuted the
  antisymmetry-penalty and claim-6 wording, corrected the cut story to thick-then-thin, added
  E7-E17, established the TRI3 exact-match and the TRI6 circulation-1/2 findings, and surfaced
  the `EF_TET4::E` edge-2 mismatch (→ **D1**). Claude independently re-derived the TET4
  conventions, the edge-2 circulation values, and spot-confirmed the antisymmetry hard-errors
  and E10/E11/E13 against the sources before inclusion. Static review throughout; nothing here
  is "verified" in the protocol sense.
