# Devlog 2026-08-31 — riva n-source setup gate; DR-105 status; the np=8 "wall" re-read

**Date:** 2026-08-31
**Topic:** A setup-time validity gate for the riva law's n source (three call sites, four bad
states closed); DR-105 re-answered against the rebuilt `tape_quench_usermat`; the np=8
floor-collapse re-measured from the surviving logs and reclassified.
**AIs involved:** Claude (Fable 5, then Opus 5), Codex `gpt-5.6-terra`/`xhigh`, Grok `grok-4.6`/`xhigh`
**Claude Confidence:** high on the gate's behaviour (probe-executed); high on the np=8 log
statistics; medium (~75%) on the np=8 *mechanism* ranking, which is not yet discriminated
**Verification:** focused regression (ladder level 2) — **`make check` green**, run by Christian
on branch `claude` @ `aecea892` + this working tree. Supported by: 12/12 standalone probe
(level 4) against a freshly compiled `cl_Material.cpp` linked to `build/lib/libbelfem.a`;
`-fsyntax-only` green on `cl_Material.cpp` and `cl_MaterialFactory.cpp` with the tree's own
flags; `check_doc_claims.py` 37/37; `input_schema.yaml` parses. Note what `make check` does and
does not establish here: it proves the gate breaks nothing that the suite covers, and the suite
contains no riva-gate case — the four closed states rest on the probe, not on the suite.

## Summary

Three threads, one session. (1) DR-105's status was asked and the answer changed mid-session:
the rebuilt `examples/tape_quench_usermat` moved from `piecewise` to `riva`, which removes the
abort *class* the row died on — independently of the 08-29 `critical temperature` remedy the row
already records. (2) Elaborating the row's "unexplained parallel wall" residue from the surviving
logs showed it is **not a wall and nothing happened at step 264**. (3) Christian asked for an
`n <= 1` checker in riva; riva's n<=1 handling is deliberate, so the check went in at setup
instead, and a three-vendor round found four holes in the first draft, all now closed.

## 1. The riva n-source gate

### Why it is not in riva

`riva_rho_pl` (`powerlaws.hpp:2481-2512`) takes `n <= 1` as the ohmic limit `ec/jc` and falls back
to fully-normal on any non-finite result; `n_eval` floors at 1 (`powerlaws.hpp:83-90`) and
`dn_eval_d{B,T}` return 0 while the floor binds, so value and tangent differentiate the same
clamped law. That totality is DR-105's remedy — a `BELFEM_ASSERT( n > 1 )` inside riva would
re-create the `rho_piecewise` abort (`powerlaws.hpp:780`) that riva exists to avoid. The check
therefore belongs at setup, where the failure is attributable and `BELFEM_ERROR` (always-active,
once-per-run) is the correct tier.

There was already partial coverage: `set_n_function` warns on rank 0 when a table's nodal minimum
drops below 1.01 (`cl_Material.cpp:673-689`). It is a `message`, not an error, and it is silent for
any source whose `min_value()` is NaN.

### What landed

`Material::check_riva_n_source()` (`cl_Material.cpp:790`), declared `cl_Material.hpp:395`, called
from **three** sites so the outcome is installation-order independent:

- `set_resistivity_law` (`:783`) — law set after the source;
- `set_n_function` (`:716`) — source installed after the law (plugin `_init` order);
- `set_constant` (`:254`) — a late constant n, the C5 hole both auditors flagged.

It refuses: a table whose `min_value()` is `<= 1.0`; a constant n that is not
`std::isfinite( tN ) && tN > 1.0`. NaN `min_value()` (analytic and plugin fits — no cheap bound)
fails the `<=` and passes unchecked, deliberately.

Plus a fourth guard in the factory (`cl_MaterialFactory.cpp:248-256`): a usermat that registers
`jc` but no n source, under `resistivity type : riva`, is now fatal at deck load. Previously a
release build ran it fully normal via NaN, silently.

### Round findings, and what they changed

Codex's first dispatch at `gpt-5.6-sol`/`xhigh` died vendor-side ("Selected model is at capacity",
twice, 156k tokens spent). Redispatched on `terra`. **The sol-vs-terra comparison §9.1 has been
waiting for is still unmeasured.**

Both auditors refuted C1. The first draft admitted `+inf` (passes `tN > 1.0`) while its own message
and both doc artifacts promised "finite"; it left the `set_constant` bypass open; it accepted a
source-less riva material; and a constant NaN died on `constant_property`'s unrelated
`"Property is not constant"` debug assert before ever reaching the riva message. All four fixed —
the last by reading the raw slot instead of the accessor.

**C3 was my error and both caught it.** I had claimed the shipped tables' B-spline control net keeps
the interpolant above 1.01. It does not: database eval is piecewise **Lagrange**
(`cl_Database.cpp:205,235`), not a Bernstein hull, so shape functions go negative and an interior
value can undershoot a nodal minimum. The control-net property is a *build-pipeline* invariant, not
an eval-time bound. Recorded as an accepted limit with the runtime ohmic branch as the safety net.

**Split verdict, ruled and open for ratification.** Both challenged the *global* table test —
Codex: values above `T_crit` are rejected though riva never reads n there; Grok: ~70% that this is
a misclassification, since the module's own prose calls softening through n=1 expected physics.
Kept global: the build pipeline floors *entire* tables at 1.02 (`/meta/n_floor` = 1.02 in
fujikura, fysc-sch04, sp-ap, sst-1, superox — read out of the HDF5 by the sol run before it died,
the one useful thing it produced), so any stored value `<= 1` anywhere marks a table that escaped
the pipeline. The gate tests the construction convention, not riva's read set. Documented as such
in all three artifacts.

Grok also found live doc drift neither the prompt nor Codex targeted: `resistivity_laws.md:157`
still listed riva's valid n range as "any table output (floored at 1)". Fixed.

**Grok breach check: clean.** Concurrent modifications to `fn_material_data_path.cpp`,
`doc/README.md` and `UserLibraryTemplate.cmake` appeared mid-round and are Christian's own
`search_data_file()` refactor in the shared checkout, outside audit scope — not touched.

## 2. DR-105 status

The row (`todo/debt_register.md:160`, `[MIXED] [P]`, P2) is open with three residues: a ~60 h
serial rerun, the np=8 collapse, and a physics ruling. The rebuilt shipped deck changes the row
without discharging the run debt:

- It is **not** the T5 deck — different mesh (`tape.geo`, single tape vs. the two-tape
  `tapestack3d.geo`), mumps not strumpack, no Anderson, timestep 0.5→1 ms vs 0.002→0.02 ms. The
  register's "two decks share this row's name" warning still holds, and now more strongly.
- The defect rewrite **is** the physics ruling, taken in code: depth 0.25 → 0.1, re-sized against
  an explicit ~0.95 mm MPZ argument, switched on at 20 ms with a C2 quintic inside the 5–155 ms
  plateau — both of the row's named physical drivers addressed.
- `piecewise` → `riva` removes the abort class independently of the `critical temperature` key,
  so the row now has two remedies where it records one.

Row text is stale in two ways (it describes the shipped deck as `piecewise`, and predates this
gate). Not amended this session — Christian's call.

## 3. The np=8 "wall" is a residual floor

Re-measured from `cmake-build-claude/tape_quench_dr105/`. **Both runs used the identical binary**
(Build Date `Aug 27 2026 03:25:34`, commit `dcca23cc`, DEBUG on, `hphiTrun`), verified from both
log headers — so the register's "debug binary" confound is **common-mode and eliminated**. Deck
identity rests on the 08-27 devlog's record ("identical deck, 1 rank"), not on the logs, which do
not echo solver settings. That is a logging gap worth closing.

| | serial | np=8 |
|---|---|---|
| accepted steps | 4191 (t = 9.14 ms) | 264 (t = 0.3749 ms) |
| rejections | 1010 | 265 |
| floor-retry warnings | 47 | 229 |
| floor warnings / accepted step | 1.1 % | 87 % |

The np=8 run was pinned at the 1 µs floor essentially throughout; step 264 is merely where it
failed 20 times consecutively and tripped `mMaxFloorRetries` (`cl_FEM_Controller.cpp:2267`).
Serial is *also* near the floor at its own step 264 (Δt = 1.1 µs). Nothing happened at step 264 —
**the register's "floor collapse at step 264" framing misleads.**

The discriminator is the achievable magnetic residual, at comparable Δt:

- serial, step 264: 1.5e-3 → 6e-6 → 4e-6 → 2e-6 → 1e-6 → … → 4e-7 (−28 → −64 dB), relax 1.0
  throughout. A clean 3.5-decade descent in 9 Picard iterations.
- np=8, step 264: flat 1.6e-4 – 2.2e-4 (−36.6 … −37.9 dB) across **every** iteration of **all
  twenty** retries. No descent. Relax cut to 0.5, 0.28, 0.25 — the residual does not respond.

Deck tolerance is 1e-6, so np=8 parks ~200× above it while serial reaches 4e-7: a **~500×
difference in achievable residual from rank count alone**. The abort message's own wording is
exactly right — the residual "does not respond to the timestep size", which by definition is not a
temporal-accuracy problem. Second symptom, same direction: np=8 thermal reads −156.54 dB (≈2e-16,
machine epsilon) bit-identically over and over, while serial thermal at the same phase is
−81…−114 dB and descending.

Mechanism ranking (**not discriminated — hypotheses, not results**): (1) rank-dependent
inconsistency in the residual assembly or its MPI reduction over owned/ghost dofs — 2e-4 is ~1e9×
machine epsilon, far too large for roundoff on a well-conditioned 2133-node problem, and the
thermal machine-zero is a second instance of the same shape; (2) STRUMPACK's distributed
factorization path; (3) partitioning of thin-shell layers / cohomology cuts at 267 nodes/rank;
(4) Anderson depth 3 — ranked low, the residual is already 2e-4 at Picard iteration 1, before
Anderson has history. The bearing was checked and cleared: `link_bearings_with_dofs`
(`cl_FEM_DofManager.cpp:1472`) broadcasts the vertex→dof table and each rank links wherever the
node exists.

Cheap discriminators, each ~15 min since the np=8 run walls fast: `compute conditioning : true` at
np=8 (the error's own hint — a modest κ rules conditioning out and points at assembly); np=8 with
`library : mumps`; np=2 and np=4 to see whether the floor scales with partition count or is a
cliff. **The new shipped deck cannot answer any of this** — it changes solver, mesh, geometry and
physics at once.

## Changes Made

- `src/physics/materials/cl_Material.cpp` — `check_riva_n_source()` (`:790`); calls at `:254`,
  `:716`, `:783`; `<cmath>` include.
- `src/physics/materials/cl_Material.hpp` — declaration + doc comment (`:389-395`).
- `src/physics/materials/cl_MaterialFactory.cpp` — usermat riva-without-n guard (`:248-256`);
  stale `cl_Material.cpp:799-802` line cite converted to a searchable anchor (`:147`).
- `src/physics/materials/powerlaws.hpp` — riva design comment: "no n > 1 precondition anywhere"
  qualified as runtime-only (`:2472`).
- `doc/input_file_reference.md`, `doc/input_schema.yaml` — the input contract, both artifacts,
  same session; schema anchors are the three function names, no line numbers.
- `src/physics/materials/doc/resistivity_laws.md` — riva's "valid n range" cell (`:157`).

## Open Questions

- **Ratify or overturn the global-table ruling** (§1). Windowing it to `T <= T_crit` is the
  auditors' alternative.
- **Residue recorded, not fixed:** NaN inside a table's stored values poisons `min()` → NaN → gate
  passes (belongs in database loading, not this gate); Lagrange interpolant undershoot with all
  nodes above 1; the pre-existing factory raw-pointer leak on a debug throw and the skipped
  `dlclose` on a ctor throw.
- **np=8 is diagnosed only to the level of "residual floor, not physics."** The mechanism needs the
  three discriminating runs above.
- ~~DR-105's row text needs a dated CURRENTNESS clause on both counts (§2).~~ **Done this
  session** — the row now carries a 2026-08-31 CURRENTNESS clause recording that the
  reconstructed deck is repaired (residue (1) unblocked), that the shipped deck is a different
  model on `riva` (two remedies where the row recorded one, and residue (3)'s ruling taken in
  code), and that the 08-31 riva gate is inert for a `piecewise` deck. **Residue (2)'s
  "floor-collapse at step 264" framing is retracted in the row itself**, not merely amended.
  DR-105 stays OPEN: residues (1) and (2) are undischarged.
- The solver-settings echo missing from `belfem.log` is why deck identity across the two DR-105
  runs cannot be proven from the artifacts alone.

## Files Updated

- todo/debt_register.md (DR-105 row: CURRENTNESS clause + residue (2) retraction/re-measurement)
- src/physics/materials/cl_Material.cpp
- src/physics/materials/cl_Material.hpp
- src/physics/materials/cl_MaterialFactory.cpp
- src/physics/materials/powerlaws.hpp
- src/physics/materials/doc/resistivity_laws.md
- doc/input_file_reference.md
- doc/input_schema.yaml
- tmp/ai_exchange/riva_n_guard.md (pre-registration + resolution; ephemeral)
