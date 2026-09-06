# DR-103 Strike: The Second-Order Weights Todo Ruled Stale

**Date:** 2026-08-27
**Purpose:** Record the investigation that retired `// todo: fix weights for second order!`
(`cl_CutProcessor.cpp:97`), the ruling that closed it, and the strike + archive of DR-103.
**Module:** homology (cut pipeline), fem/kernel (T-matrix consumption)

## Context

DR-103 (the never-executed TET10 cut-pipeline crash) was fixed, gated, and committed on
2026-08-24 (`0db03e32`; the register's "Uncommitted" note had gone stale). During the
2026-08-27 register-cleaning pass, the one residue blocking a clean strike was the ancient
todo in the `CutProcessor` constructor — present since the very first CutProcessor commit
(`6ee74493`, 2025-02-05) — and tracked as C5 in `todo/deferred/2d_thinshell_todo.md`.
Christian asked whether the weights could be fixed, pointing at `src/numerics/integration`
as the source for any quadrature weights needed.

## Investigation (static, full path trace)

The question was: what weights, and are they wrong at second order?

**What the weights are.** The only weight-producing step in the entire cut pipeline is
`CutSet::create_duplicates` (`src/homology/cl_CutSet.cpp`): every duplicated node hangs on
`{abstract current nodes…, original}` with all weights 1.0, encoding φ_dup = φ_org + Σ I_c
as a nodal T-matrix tie. The DofManager consumes it as `TᵀKT`
(`cl_FEM_DofManager.cpp`, `consolidate_dofs`). Confidence: high (verified on disk).

**There is no weak form to weight.** `DomainType::Cut` sidesets never reach an IWG: the
sideset whitelist in `MaxwellFactory` (search `set_sidesets`) does not include `Cut`, the
`FieldList::Cut` doftable is never filled, and the `SideSetDofLinkMode::Cut` branch in
`cl_FEM_Element.cpp` is dead code. The historical `compute_jacobian_and_rhs_cut` (±1 lambda
rows, itself order-agnostic) was deleted long ago. No facet quadrature is ever formed
anywhere in the cut path. Confidence: high.

**Why 1.0 is exact at any order.** The physical jump across a thin cut is [φ] = I,
constant over the cut surface. φ is a nodal Lagrange basis at order 2 as well
(TRI6/TET10), so tying every node pair — corners and midsides alike — with weight 1.0
reproduces the constant jump exactly. There is no ½/⅓/edge-length structure to correct,
because the constraint is collocated, not integrated. Gauss tables from
`src/numerics/integration` would only enter if the constraint were reformulated weakly
(∫_Γ λ([φ]−I) dΓ) — a design change, not a repair. Confidence: high on the mathematics;
the todo's original intent is recorded as Christian's to confirm, and he ruled.

## Ruling and actions

**Christian's ruling: remove the todo so DR-103 can be struck.** Applied same session:

- `src/homology/cl_CutProcessor.cpp` — todo line removed (no code change).
- `src/homology/cl_CutSet.cpp` — three-line rationale comment added at the weight vector,
  recording why unit weights are exact at any Lagrange order.
- `todo/deferred/2d_thinshell_todo.md` — C5 struck with the ruling.
- `todo/debt_register.md` — DR-103 struck (ID + description; status cell carries the
  closure evidence per register convention) and moved to `debt_register_closed.md`;
  preamble counts refreshed ([P] 24→23, [W] 6→7). The stale "DR-31's gate remains
  blocked, now on DR-103" clause in DR-102's status cell struck in place.
- **DR-122 filed ([CODE][W], P3)** — the genuine order-2 residue split out at strike time,
  all unreachable behind the deliberate order-1 guard in
  `MaxwellFactory::create_hanging_edges_and_facets` (`cl_MaxwellFactory.cpp:1413`):
  1. diagonal thick-cut cases ±5/6/7 rely on chord-midside neighbor coverage that
     `thick_thin_cuts_and_conjugate_edges.md` documents as unenforced;
  2. `CutSet::create_duplicates` has no hanging-source cascade (unlike
     `InterfaceProcessor::duplicate_nodes` and the ThinShellFactory) — at order 2 that
     trips the `is_hanging` assert at `cl_FEM_DofMgr_DofData.cpp:3513` in debug and writes
     a wrong T-matrix row in release;
  3. no TRI6/TET10 homology fixture exists.

## Process note

No vendor audit round was run: the source edit is comment-only (zero executable change),
made on Christian's direct instruction; the plan+audit protocol targets behavioural fixes.
The path trace behind the ruling was a single-session static investigation — reviewed, not
verified, which is the appropriate rung for a comment removal. The claims about the live
tie mechanism were verified against the tree (file paths and cited asserts checked on
disk before writing).
