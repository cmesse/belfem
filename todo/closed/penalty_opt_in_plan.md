# Ghost and Gauge Penalties: Honour eta = 0, Make the Penalties Opt-In

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): opt-in ghost and gauge landed with parse tests 2026-09-01; R7/R8 gates owed. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-09-01
**Purpose:** Let a deck switch the Nitsche ghost coupling of thin-shell layers OFF — not just weaken
it — so that the duplicate interface dofs and ghost facets are never created, and demote the
Coulomb-gauge penalty (and, as a second step, the ghost itself) from default-on to opt-in. Trigger:
today's measurements on `tape_quench` — the conditioning of the magnetic system is set by Δt
(κ ≈ 2.1e15 × Δt[ms], constant to 5 % over 30× in Δt), `chi = 1e-4` does not move it at all,
`chi = 0.01` moves it by a factor 2–3, and a run at `eta = 4e-6` converges exactly like one at
`eta = 4`. Neither penalty is doing what its default promised on this deck.
**Module:** `src/fem/kernel` (ThinShellFactory, Controller), `src/fem/maxwell` (MaxwellFactory),
`doc/input_*`
**AIs involved:** Claude (plan), Codex (plan jury 2026-09-01), Grok (relay re-run 2026-09-01 after two max-turn failures)
**Status:** IN PROGRESS — **code landed 2026-09-01 (R1–R6, R9 unit fix), code jury done, fixes
applied; owed: R7 factory-level test, R8 gate (D6 two-run experiment first), R10 devlog numbers.**
Plan jury (Codex + Grok relay): D1–D6 folded in. Code jury (Codex + Grok): no P0; one P1 and four
P2 on this change fixed; four findings on Christian's co-resident hunks routed to him. Exchange:
`tmp/ai_exchange/review_penalty_opt_in.md`.

> **Scope guards:**
> - The ghost *formulation* (`mt_maxwell_h.cpp` ghost term, DG edge groups, `cl_ThinShellFactory`
>   duplicate machinery) is NOT touched. This plan adds a switch and moves defaults.
> - The `.bfm` reload path keeps whatever the saved mesh has; a mismatch between the deck's switch
>   and the saved mesh is refused, not repaired (gap row 5).
> - Side-connector walls (`edge coating`) are out of scope.
> - No physics claim about which coupling is *better*; the plan makes both reachable from the deck.

---

## 1. Current Behaviour

| what | where | today |
|---|---|---|
| duplicate edge/face dofs at a layer interface | `cl_ThinShellFactory.cpp:146-158` | created iff `mCreateGhostFacets && tA != tB && both have(rho)`; buffer interfaces never (`:2131`) |
| ghost facets + `DomainType::Ghost` sideset | `:245-270`, `:294`, `:308` | created iff `mCreateGhostFacets`; `nullptr` sideset is a supported path |
| `mCreateGhostFacets` | `cl_ThinShellFactory.hpp:335` | **hardcoded `true`, no setter, no ctor argument** — this is the flag Christian remembers as `mUseNitsche` |
| `eta`, `k_reg` | `cl_FEM_Controller.cpp:4459-4480` | `eta > 0` enforced ("eta must be positive"), `k_reg > 0` with unit; **absent block = defaults (4.0, 1e-3 Ohm)** — "load-bearing, must not switch off by omission" (`devlog/dl20260824_penalty_input_keys.md`) |
| `chi` | `:4437-4457` | absent block = **ON at 1e-4** (flipped from off on 2026-08-27, Christian's conditional ruling: "the ruling's basis will be whichever configuration traverses", `dl20260827_evening_quenchfront_campaign.md:71-73`); opt-out `chi : 0` |
| where the deck is first read | `cl_MaxwellFactory.cpp` (ctor, `create_thinshells()` rank 0) | the thin shells are built **before** `Controller::set_params()` parses the penalty block — the switch must be read in the factory |

**What "off" is physically.** With the ghost off, adjacent layers share their interface edges
again: tangential H is C⁰ across the interface, there is no interface resistance, and the
pre-2026-03-18 model is recovered (`devlog/dl20260318_ghost_thinshell.md` §1: the ghost was
introduced because sharing forced C⁰ continuity "across extreme resistivity contrast"). An
explicit resistive interface is then modelled the way `tape_quench` already does it — a thin `rint`
layer with its own ρ — not by the Nitsche term.

**Evidence from 2026-09-01 (this deck, builtin materials, MUMPS with the new cap):**
- κ(Δt): 1.2e13 at 0.005 ms … 4.5e14 at 0.23 ms; κ/Δt = 2.0–2.4e15 throughout. The mass term
  M/Δt sets the small end of the spectrum. Penalties are not where κ comes from.
- `chi = 0` vs `1e-4`: κ identical to two digits at five steps. `chi = 0.01`: κ 2–3.5× lower at
  steps 2–4 (2.0/1.5/1.3e13 vs 3.5/3.9/4.6e13). The estimator sees the gauge; the default is
  below its noise.
- `eta = 4e-6` (six decades below default): step 1 at 2.2e-16, Δt grown to 0.19 ms by 1 ms, one
  soft fail (a `-20`, unrelated). The ghost *magnitude* is not load-bearing for convergence here;
  whether the ghost *structure* (duplicates) is, nobody has measured — that is what the switch is for.
- The hastelloy|buffer J speckle (15 ms frame) sits on the one interface that has no ghost facet
  at all; the ghost is not its cause.

**Why the penalties existed, and why that reason is gone (Christian, 2026-09-01 evening).** Both
were stabilizers for a system that was *exactly singular*: an unpinned φ constant on every
connected φ component — the buffer layers enclosed between conductors float — gave κ ≈ 1e17
(DR-126). The Nitsche ghost (2026-03) and the Coulomb gauge (2026-08) were built to tame that
matrix. The actual fix landed 2026-08-28: `find_autopins` (`cl_MaxwellFactory.cpp:3154`,
commit `04c33664`) labels every φ component and fixes one safely pinnable node per component
to 0 at setup ("prescribed (fixed) : 3" in today's tape dof line is that). With the singularity
removed, κ is set by the timestep (2e15 × Δt[ms]) and neither penalty has a job left — which is
what the A/B measured. That is the causal reason for opt-in; §1's numbers are the confirmation.

## 2. Design

**S1 — the switch.** `nitsche ghost penalty { eta : 0 ; }` means OFF: `MaxwellFactory` reads the
block from `solver → nonlinear magnetic` before `create_thinshells()` and passes
`aCreateGhostFacets = ( eta > 0 )` to a new `ThinShellFactory` constructor argument (default
`true`). `Controller::set_params()` relaxes its guard to `eta >= 0` and, at 0, skips
`set_penalty( …, 0 )` — there are no ghost facets to penalise. `k_reg` is meaningless at
`eta = 0` and is ignored with a one-line notice if stated. Two readers of one block, one rule:
`eta == 0` ⇔ no ghost.

**S2 — gauge opt-in.** Absent `coulomb gauge penalty` block = `chi = 0` (the pre-2026-08-27
semantics, `dl20260824_penalty_input_keys.md:24`). Stated block keeps requiring `chi`. One line in
the controller (`:4451-4456`) plus the comment, schema and reference defaults, and the memory
note "never cite the campaign as a reason to default chi off again" is superseded by **this**
ruling and dated.

**S3 — ghost opt-in (the design call, O1).** Absent `nitsche ghost penalty` block = OFF. This is
what Christian leaned to ("demote all three"). Cost: **every thin-shell deck without the block
silently changes model** — examples, gate decks, and any test mesh that relies on duplicates
(`tests/fem/test_InterfaceOrientation.cpp` touches ghost facets). Recommended sequencing: land
S1 + S2 now (they change nothing for a deck that says nothing), sweep the decks, then flip S3 in
its own commit with the deck sweep recorded.

## 3. Gap Table

| # | State / behaviour | Class | Handling |
|---|---|---|---|
| 1 | `mCreateGhostFacets` has no setter | (c) | ctor argument, default `true` |
| 2 | penalty block parsed twice (factory + controller) | (c) | one helper `read_ghost_eta( const input::Section * aSolver )` returning `eta`, or `-1` for an absent block; used by both readers; **absent → OFF** (O1: S3 lands in this change). A present block **requires `eta`**; `k_reg` alone or an empty block is a setup error naming the rule (D4) |
| 3 | `eta = 0` rejected today (`:4468`) | (c) | `>= 0`; at 0 skip the penalty write |
| 4 | `k_reg` stated with `eta = 0` | (c) | ignored + notice (not an error: a deck flipping eta to test should not have to delete k_reg) |
| 5 | `.bfm` reload carries shells built with or without ghosts | (c) | the **create-ghost boolean goes into the mesh config tag** (`fn_mesh_config_tag.hpp`), so a flipped deck misses the cache and rebuilds; an explicitly named `.bfm` compares the boolean stored in the file against the deck and refuses on mismatch. NOT a sideset-null test (D3: `ThinShell` has `ghost_id()`/`ghost_facets()`, and creation ON legitimately yields zero facets) |
| 6 | no-ghost thin shell through DofManager / IWG / FieldList | (b) → O2 | the `nullptr` sideset path exists (`cl_IWG.cpp:650`, `cl_Maxwell_FieldList.cpp:413`); exercised today only by single-layer or all-buffer stacks. Needs the smoke run |
| 7 | multi-rank: factory runs on rank 0, the flag must not diverge | (a) | read from the deck, which every rank holds; the controller's `set_penalty` is already collective |
| 8 | memdump / warm start across a switch flip | (c) | **the existing checks do NOT refuse it (P0, D5):** `Mesh::compute_checksum` hashes node count + element node ids only (`cl_Mesh.cpp:2416-2456`), blind to duplicate edges, and `load_vector_from_file` resizes to the dump length (`hdf5_tools.hpp:672`). Fix: store the create-ghost boolean in the dump and refuse on mismatch; always-active length check in `load_fields()` naming the field |
| 9 | input contract | (c) | `eta` constraint `>= 0` with the 0 semantics; `chi` default; `nitsche ghost penalty` absent-block semantics (S3) — schema + reference same session |
| 10 | prior rulings | (c) | 2026-08-24 "must not switch off by omission" and 2026-08-27 "chi on by default" are reversed by Christian 2026-09-01 on the evidence in §1 — recorded in the devlog and the memory notes, not silently |

## 4. Ordered Steps

- [x] **R1** `ThinShellFactory( …, const bool aCreateGhostFacets = true )`; member set from it.
- [x] **R2** `read_ghost_eta()` helper (kernel-level, deck section in, `real` out; −1 = absent, absent = OFF; present block requires `eta`).
- [x] **R3** `MaxwellFactory::create_thinshells()` reads the switch and passes it (rank 0 only
  builds, but the read is rank-uniform).
- [x] **R4** Controller guard `eta >= 0`; `k_reg` notice at `eta = 0`; gauge absent → 0 (S2).
- [x] **R5** create-ghost boolean into the mesh config tag and the `.bfm`/memdump headers; mismatch refusals; `load_fields()` length check (gaps 5, 8 — D1, D3, D5).
- [x] **R6** Docs: schema, reference, `src/fem/maxwell/doc/ghost_penalty_stabilization.md` (§3
  says the ghost cannot be switched off), devlog, memory notes.
- [ ] **R7** Tests: a NEW factory-level case (the existing `test_InterfaceOrientation.cpp` drives a
  hand-built `TS_TestGhostStack`, not the factory) with the factory at `false`: no duplicates, no
  ghost sideset; a controller parse test for `eta : 0`, for a `k_reg`-only block (error), and for
  absent `chi`; a config-tag test that the boolean changes the tag.
- [ ] **R8** Gate: `tape_quench` with `eta : 0` vs `eta : 4` from cold, same Δt cap
  (`maximum timestep : 0.2 ms` so κ stays < 5e14): dof counts (170941 → fewer), κ, Picard
  counts, and the sheet-by-sheet J table at a matched frame.
- [◐] **R9 (S3, same landing per O1)** — default flipped and the `Ohm*m` unit fixed; per-deck sweep table in the devlog; no block added anywhere (O1) absent ghost block → OFF. Deck sweep with a per-deck
  *duplicates possible?* column: single-layer / same-label stacks (`undulator2d`, `gantry`,
  `tapestack2d_gregory`, `block3d`, `garber`, `corc`) are no-ops and get no block; the
  multi-metal decks (`tapestack2d_christian`, `circuit`, `tapestack3d_*`, `tape_quench_usermat`,
  `corc2` in the build tree) get an explicit block **only once D6 says what they run today**. The
  commented `tapestack3d_*` block carries `k_reg : 1e-3 Ohm*m` (illegal unit) — fix, do not
  uncomment. `make check` on the flipped default.

### 4.0 Defect tracker (plan jury round 1, 2026-09-01 — Codex jury + Grok relay re-run; all verified by Claude)

**Code round (2026-09-01, Codex + Grok jury):**
- [x] **D7 (P1, Codex + Grok).** A named `.bfm` from before the switch (no tag line) bypassed the
  refusal. Fixed: no line ⇒ built ghost-ON (the only layout the old factory produced).
- [x] **D8 (P2, Codex + Grok).** Reference prose and example still said on-by-default. Fixed.
- [x] **D9 (P2, Grok).** Gauge theory doc, maxwell README, IWG ctor comment stale. Fixed.
- [x] **D10 (P2, Grok).** Memdump/`.bfm` refusal overclaimed for pre-change files; `num_faces`
  read without its own existence check. Docs narrowed; check added.
- [x] **D11 (P2, Grok).** Tests littered the cwd; two cases missing. Fixed. Factory-level test
  (R7) still owed.
- [x] **D12 (P0, Claude, resolved D6 — 2026-09-01 evening).** **The `.bfm` reload collapses a DG
  stack.** Measured on the two caches on disk: the 2026-08-30 usermat cache (fresh build, flag ON)
  holds 526385 edges = 177673 + **14** × 24908 (9 interface sheets + 5 duplicated), a ghost sideset
  (`SideSet_25`) and `thinshells/ghost = 25`; today's flag-OFF cache holds 501477 = 177673 + **13**
  × 24908 for the 12-layer stack, ghost 0. Yet every run that *reloaded* the 08-30 cache reported
  170941 dofs: `BfmFile::load_edge_data()` ends in `ProtoMesh::reconstruct_edge_connectivity()`,
  which relinks elements to edges **by node pair** (`cl_ProtoMesh.cpp`, the twin-edge caveat in its
  own comment) — the duplicate interface edges share their node pairs with the originals, come
  back as orphans, get no dofs, and the ghost sideset then couples nothing. So: the flag was ON in
  every commit since 08-26 (checked), the factory did build the DG stack on 08-30, and **every
  cache-loading ghost run since has silently been a shared-edge run**; today's fresh runs had the
  flag hardcoded OFF. Guard landed: with the ghost ON the cache is never reused (rebuild with a
  message) and a named `.bfm` is refused; the real fix — a reload that preserves twin edges — is a
  separate item (mesh module, `reconstruct_edge_connectivity` needs the per-layer edge lists, not
  node pairs).
- **Routed to Christian (his hunks in the same files):** `I_*`/`U_*` uniqueness (Codex P1),
  thermal notices over-deleted (Grok), `T_max` rank-0 max (Grok), `temperature` label legal again (Grok).


- [ ] **D1 (P0, Codex).** The `.bfm` cache identity (`fn_mesh_config_tag.hpp:321-381`) hashes
  topology, homology, terminals and layers — never the solver section — so an unchanged `.msh`
  reuses a ghost-enabled `.bfm` after the deck flips `eta` to 0 (`cl_MaxwellFactory.cpp:401-418`).
  Fix: put the semantic boolean *create ghost facets* into the tag text (not raw `eta`); keep the
  refusal only for an explicitly named `.bfm`. Verified by Claude.
- [ ] **D2 (P0, Codex).** Skipping `set_penalty( 0, 0 )` at `eta = 0` leaves the IWG default
  `mPenalty = { 4.0, 1e-3, 1e-4 }` (`cl_IWG_Maxwell.cpp:51`) in place and the log would print
  `ghost eta 4` — contradicting O3. Fix: write 0 collectively on every no-ghost path (`set_penalty`
  broadcasts, `cl_IWG.cpp:1073-1077`). Verified.
- [ ] **D3 (P0, Codex).** Gap row 5 names `ghost_sideset()`, which does not exist (`ThinShell` has
  `ghost_id()` / `ghost_facets()`, `cl_ThinShell.hpp:117-130`), and its predicate is wrong: with
  creation ON a stack can legitimately have zero ghost facets (no duplicated interface, or all
  buffer-bounded), so "switch ⇔ non-null sideset" rejects valid meshes. Fix: persist the boolean
  in the `.bfm` (config tag, D1) and compare booleans. Verified.
- [ ] **D4 (P1, Codex).** Ghost-block grammar after S3 is underspecified: today `eta` and `k_reg`
  are independent (`input_file_reference.md:297-299`); the helper only sees eta present/absent.
  Fix: a present block **requires `eta`**; `k_reg` alone or an empty block is a setup error naming
  the rule. Verified.
- [ ] **D5 (P0 — Codex P1, promoted by Grok with the mechanism).** Warm restart across a flip:
  `Mesh::compute_checksum` hashes node count + element node ids only (`cl_Mesh.cpp:2416-2456`) —
  duplicate edges are invisible to it — so the memdump refusal at `:3053` cannot fire, and
  `load_vector_from_file` resizes to the dump length (`hdf5_tools.hpp:672`): a ghost-on dump
  silently seeds a ghost-off mesh. Fix: the boolean in the dump header + refusal; always-active
  length check in `load_fields()`. Verified.
- [ ] **D6 (P0, Claude, from the running experiment — NOT seen by the jury).** On this deck the
  switch is **inert**: Christian's hardcoded `mCreateGhostFacets = false`, fresh `.msh` build,
  gives 170941 free dofs — the same count as every flag-ON run today; only the 2026-08-30 usermat
  run 1 (fresh, ON) ever showed 295481. So no duplicate interface dofs were being created with the
  flag ON either. Grok refuted both of my candidates: plugin `rho` DOES set `have`
  (`cl_Material_UserDefined.cpp:55-65` → `set_custom`), and twin edges carry unique ids and are
  counted per id (`cl_ThinShellFactory.cpp:1853-1891`, `cl_FEM_DofMgr_DofData.cpp:1305-1337`).
  Grok's candidate: the `.bfm` short-circuit — `BfmFile` persists shells and `create_thinshells()`
  returns early on a shell-bearing mesh (`cl_MaxwellFactory.cpp:1256-1286`), so a run that loads a
  cache built without duplicates keeps that structure whatever the flag says. Today's arithmetic
  (295481 − 170941 = 5 × 24908; 270573 − 170941 = 4 × 24908, one edge sheet per interface) says
  the duplicates existed once. **Blocks R8.** Deciding experiment, no code: flag ON, `rm tape.bfm`,
  cold run → dof line; same deck again with the cache present → dof line.

## 5. Open Questions

- **O1 — flip the ghost default (S3) now or after the sweep?** **RESOLVED 2026-09-01, Christian
  ("I'm not afraid of no-ghost"): all three opt-in now, in one change.** The deck sweep (R9)
  becomes part of the same landing, not a later commit: every thin-shell deck in `examples/` and
  the gate decks gets an explicit block if it is meant to keep the ghost, and `make check` runs
  on the flipped default.
- **D6 — RESOLVED 2026-09-01 by D12** (the `.bfm` reload, not the factory, not `have(rho)`).
- **O2 — is the no-ghost path exercised anywhere with duplicated *materials*?** Today
  `mCreateGhostFacets = false` never happens with a multi-metal stack. R8 is the first such run;
  if the DofManager path needs the ghost sideset for something other than ghost dofs, R8 finds it.
- **O3 — should `eta = 0` also be the value the controller reports in the penalty-slots line?**
  Yes: `ghost eta 0 ( no ghost facets )`, so a log states the model in force.

## 6. Definition of Done

- [ ] `eta : 0` vs `eta : 4` on `tape_quench` from cold with `tape.bfm` deleted differ in the dof
  line by 5 × 24908 (the five duplicated interfaces); D6 resolved by that experiment first.
- [ ] A deck without a gauge block reports `gauge chi 0`.
- [ ] Both input-contract artifacts updated; `check_doc_claims.py` run.
- [ ] R8 numbers in the devlog; prior rulings marked superseded with the date.

## 7. Audit Trail

- To be dispatched after O1: plan jury (Codex `gpt-5.6-terra`/high, Grok `grok-4.6`/high).
- Records reversed: `devlog/dl20260824_penalty_input_keys.md` (absent ghost block = defaults),
  `devlog/dl20260827_evening_quenchfront_campaign.md` §3 (gauge on by default, conditional).
- Evidence: `cmake-build-debug/tape_quench/out.txt` launches 6–7 of 2026-09-01 (κ vs Δt, chi A/B),
  the 15 ms frame sheet table (this session).
