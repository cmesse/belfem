# DR-38 Sweep: the Thermal `1.000000` Entry Residual Is an Unseeded-Dof Artifact

**Date:** 2026-08-18
**Purpose:** Three-voice read-only sweep (Claude + Codex + Grok) into the question the
step-124 jury left open — why the thermal Picard-1 residual reads exactly
`1.000000 ( 0.00 dB )` on the first step after a warm restart, on a run whose banner
claims the thermal BDF history was restored.
**Thread:** `tmp/ai_exchange/dr38_thermal_entry_residual.md` (distilled here before sweep)
**Register:** new row DR-91; DR-38's ledger note now points at it

## What it is

`ThermalFactory` initializes the thermal DofManager at `cl_ThermalFactory.cpp:242`
and fills the temperature field with `gTbulk` only at `:247-248`. Initialization
runs `init_dof_values()`, whose free-dof branch copies field → `Dof::value()`
(`cl_FEM_DofMgr_DofData.cpp:2891-2893`) — at that moment the field is the zero
vector `create_fields()` just made. The two `fill()` calls afterwards, the second
fill in `hphiTrun.cpp:85`, and every later `load_fields()` touch the *mesh vector*
only. Nothing pushes them into the dof objects, because `DofManager::initialize()`
early-returns on `mInitializedFlag` (`cl_FEM_DofManager.cpp:330-334`, set at `:361`,
never cleared) and the public `init_dof_values()` wrapper (`:1245-1247`) has no callers.

The Picard residual is `norm( A x − b ) / norm( b )` evaluated against the
pre-update iterate (`SolverData.cpp:2287-2290`, `:2411-2425`). With `x = 0` that is
`‖−b‖/‖b‖` — exactly 1.

**The magnetic asymmetry is the same fact from the other side.** `MaxwellFactory`
never calls `initialize()`, so the call inside `Controller::load_memdump`
(`cl_FEM_Controller.cpp:4169`) is the first one for the magnetic manager and runs
*after* the fields are restored: magnetic dofs are seeded from the dump. The thermal
call at `:4178` hits the early return.

Cold-start control, found in two places independently: `sidecoatings/out.txt:116-117`
(Claude) and `tapestack3d_out_20260817.txt:108-109` (Grok) show **both** physics
printing `1.000000` at step 1 of a cold start, normal values at step 2. The artifact
fires once per process, warm or cold — it is not a restart defect at all.

## What it costs

Physics: nothing shown. `A` and `b` are assembled from mesh field data
(`Calculator::q()` `:2881-2889`, `qold()` `:2781-2805`), not from the zero dof
vector, and BCs are recomputed before the solve. Grok supplied the decisive
argument: at ω = 1 with an empty Anderson window the Picard update is
`x_new = x + β( G − x ) = G` **regardless of x**, so the unseeded zero writes the
same field a correctly seeded step would. The residual is a lie; the solution is not.

Controller: also nothing, on this deck. Claude's claim that the fake
`mEpsilon20 = 1.0` hands the thermal ω schedule its maximum arctan growth step was
refuted — the first coupled iterate skips the AIMD branch (`Controller:1571-1574`)
and thermal ω is already pinned at `max relaxation : 1.0`.

**The reason to fix it is that the safety is accidental.** A driver that starts the
first thermal iterate at ω < 1 would write `T ← ω·G` and genuinely pull the restored
temperature field toward zero. Nothing in the code records that ω = 1 is load-bearing
here. Open and not shown: thermal Anderson depth 1 stages the junk pair `( x = 0,
r = G )` on every process's first iterate.

## Consequence for the step-124 jury (unanimous)

That jury's V2 leg — "thermal entry states differ, 0.00 dB vs −82.99 dB" — is
**withdrawn**. It compared this artifact (`out.txt:433`, first step of a process)
against a real residual (`out.txt.old:13157`, live step); each log contains exactly
one `1.000000`, at its own first step. It never showed that the restored `T`
differed from the live `T`, and the ω = 1 argument says the field-space update was
the same map either way.

What survives of H3: fresh versus stale STRUMPACK matching, and the still-unexplained
11 dB magnetic Picard-2 split after a bit-identical Picard-1. "The commits made 124
converge" remains unproven, and so does "the restart state was worse".

## Fix (implemented 2026-08-18, plan+audit → code+audit with both vendors; reviewed, not verified)

Landed the same day once the campaign reached a memdump, through the full round —
plan audit (Codex high / Grok high, no blocker) and code audit (both high, no
blocker) on the diff:

- `DofData::init_dof_values( const bool aFreeDofsOnly = false )` — default in the
  header only; seeding mode skips ONLY the fixed-dof dof → field write, keeps
  `flag()` so the dedup invariant is untouched.
- `DofManager::seed_dof_values()` — free-dofs-only wrapper, guarded by
  `BELFEM_ASSERT( mInitializedFlag, ... )` (Grok's hardening).
- `ThermalFactory`: seed directly after the two `fill( gTbulk )` calls (cold-start
  half). `load_memdump`: seed after `distribute_fields()` for both kernels
  (restart half — the thermal one is the DR-91 write; the magnetic one is an
  idempotent repeat, kept for the single after-distribute-seed rule).

Audit corrections folded in: my magnetic-worker rationale was refuted by both
vendors (`init_dofs()` distributes fields *before* seeding, so an uninitialized
manager already seeds workers from restored rank-0 fields); the header comment's
"never runs again" softened per Grok (`reset()` via `set_equation()` clears the
flag — construction path, not this one). Syntax gate passed independently twice:
Claude with `g++ -fsyntax-only`, Codex with `mpicxx`, both on the targets' own
`flags.make` flags.

Deliberately NOT done, and still open: the fixed-dof back-write during a full
`initialize()` after a rank-0 restore (H-C, `cl_FEM_Controller.cpp:4170`) — out of
scope by decision, recorded on the DR-91 row.

## Gate: PASSED, verified by execution (2026-08-18)

Warm restart on the rebuilt binary, from the step-152 dump:

```
WARM RESTART from memdump.hdf5
resuming at t = 2150.0000 ms, timestep 152, delta t = 13.4137 ms
  Magnetic Picard 1, residual 0.024278 ( -16.15 dB )
  Thermal  Picard 1, residual 0.000000 ( -91.57 dB )   <-- was exactly 1.000000
```

−91.57 dB sits inside the campaign's −73…−99 dB continuation band, and the
magnetic entry is unchanged at −16.15 dB. So the seeded `x` is the restored `T`,
and the restart satisfies the thermal operator to the same depth a live step
does. **The entry state was never wrong — only the yardstick.** That retroactively
confirms the sweep's central claim and closes the step-124 jury's thermal leg for
good.

By-catch of the fix, as designed: the first-iterate Anderson column is now the
real `( x, G − x )` pair instead of `( 0, G )` — thermal depth 1 on this deck.

Process note worth keeping: the first gate reading was a **false fail** produced
by my own watcher, not by the code. `wc -l < out.txt || echo 0` returned 0 on one
poll, which reset the watcher's baseline and made it re-report line 433 — a line
from the *previous* run, still present because the new process appends to the same
log. The lesson is to anchor a log gate to a unique event marker (the new WARM
RESTART banner) rather than to a line count on a file another process is writing.

## Method note

Both vendors independently proposed the same discriminator Claude had already run
(grep a cold-start log), and both corrected Claude: Grok on the ω-schedule damage
claim and on "field → dof happens in exactly one place" (true for the campaign path,
too strong in general — `Dof::value()` has five other writers), Codex on the generic
case where a manager *is* first initialized after a field load, which is precisely
the live magnetic restart path.
