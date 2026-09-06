# Warm-Restart Loader Contract: Δt Consistency and the Missing-Triple Hole

> **CLOSED 2026-08-30 — OBSOLETE, moved out of the active set by the register/todo currentness sweep.**
> Both halves of this plan are gone. Half (a), the missing-BDF-triple hole, was fixed and **DR-112 was
> struck 2026-08-30** after its three-case restart gate ran by hand. Half (b), the Δt inconsistency, was
> **retracted as a false positive** on 2026-08-29 — the archived DR-112 row records the retraction and
> states in terms that this file "is superseded by that retraction". The Status line below ("OPEN — R1
> not started") predates both events. Kept for the record.

**Date:** 2026-08-29
**Purpose:** mitigate both halves of DR-112 — the permissive missing-BDF-triple acceptance, and the
Δt-inconsistency of the first post-restart timestep measured by the DR-124 R6 gate
**Module:** fem/kernel (warm restart), with a circuit-visible symptom

**Status:** OPEN — R1 (localization) not started. The defect is **measured and reproducible**, the
owning module is **not yet established**, and no fix is proposed until it is. Nothing in this plan
has been implemented. Scope covers DR-112 only; this plan exists because Christian ruled
(2026-08-29) that findings extend the row being worked rather than branch new IDs — a DR-149 was
filed for the Δt half and withdrawn the same session.

---

## 1. What is wrong

Two defects in one loader contract (`Controller::load_memdump`, `cl_FEM_Controller.cpp:4734+`).

**(a) The missing-triple hole (DR-112 as originally filed).** A dump carrying no
`bdf_h` / `bdf_step_count` / `bdf_last_dt` triple is accepted by cold-starting the BDF order ramp,
even when the restored run is configured for BDF2–5. A partial triple already hard-errors by name;
the wholly-absent case is the permissive hole. Christian ruled 2026-08-27 that breaking old files
is acceptable for a pre-release format, so this should become a loud loader error.

**(b) The first post-restart step is not Δt-consistent (measured 2026-08-29).**
`load_memdump` caps a restored `mDeltaTime` at the deck's `initial timestep` before the first solve
— the deliberate DR-92 cap, which exists because a full restored Δt poisoned STRUMPACK's lifetime
setup. That first step then produces the state the uninterrupted run reaches a **full dumped Δt**
later.

### Where the cap came from — read this before touching it

The cap is not an accident and must not simply be deleted. It is DR-92's fix (landed
2026-08-19, row struck 2026-08-28): a large-timestep warm restart poisoned STRUMPACK's
lifetime setup because the first factorization was built from the full saved Δt. Capping the
first re-entry step at the deck's `initial timestep` fixed that, with a production gate showing
the first preconditioned residual falling 10851 → 0.903286 and Krylov iterations 50 → 16.

**(b) is therefore a consequence of the DR-92 fix, not a recurrence of DR-92.** The solver-setup
symptom is genuinely fixed; what the cap left behind is that the clock and the integrated state
now advance by different amounts on that one step.

That also explains why it went unseen: DR-92's gate measured preconditioner quality and
convergence, both of which are blind to a step that advances the right amount of physics over
the wrong amount of time. Only a trajectory comparison against a never-dumped twin shows it.

DR-92's closure additionally records that a post-landing audit "confirmed the clamp ordering,
rank consistency, and uncontaminated BDF history". That may still hold — history VALUES can be
correct while being consumed at the wrong Δt — but R1 should re-read that audit against the
measurement below rather than treat it as settled.

### Evidence for (b)

`examples/RLC_Circuit`, `belfem` executable, np=4, 4 ms dump resumed to 8 ms, deck otherwise
byte-identical to an uninterrupted 8 ms baseline:

| quantity | value |
|---|---|
| clock advance on the first post-restart step | `1.0e-5 s` (= deck `initial timestep`) |
| dumped Δt | `2.0e-4 s` |
| first-row currents `IL1/IC1/IR1` | **bitwise identical** to the baseline row `2.0e-4 s` after the seam |
| that step's current rise | `144 A`, where neighboring steps rise `~4.6 A` |
| node voltage `V0` | `0.961389` = `sin(2*pi*50*0.00411258)` — correct for the NEW time |
| residual over the 4–8 ms overlap | `7.5e-3` relative on the currents, decaying |

So the **sources are evaluated at the new time while the reactive components integrate a full old
Δt**. Two controls make this admissible rather than suggestive:

- **Determinism.** The 4 ms leg and the 8 ms baseline are bitwise identical over all 27 common rows
  and all 7 channels, so nothing after the seam is run-to-run noise.
- **Not state loss.** Re-running the same dump with `initial timestep` raised to the dumped 0.2 ms
  makes the resumed run bitwise identical to the uninterrupted baseline over all 20 compared rows,
  every channel. That is what struck DR-124: the restored state itself is exact.

---

### Who else is exposed — DR-105 is the sharp case

The RLC measurement used a 20x ratio between the dumped Δt and the deck's `initial timestep`.
**DR-105's tape_quench deck has a 10x ratio and a far higher cost of being wrong:**
`initial timestep : 0.002 ms`, `maximum timestep : 0.02 ms`, `simulation time : 0.3 s`
( `cmake-build-claude/tape_quench_dr105/input.conf:51-54` ), projected at ~60 h serial. A run
that long is restarted, and every restart then advances the clock by 0.002 ms while integrating
0.02 ms of state.

That deck's whole open question is Δt behaviour near a contraction radius, which is the worst
possible place for a silently mis-advanced first step: it would present as deck physics, not as
a restart bug. Anyone tuning DR-105 across a restart should treat its post-restart steps as
suspect until R1 lands.

Two further notes for that deck specifically: it still uses the old `custom { }` subsection
spelling and will now hit the rename diagnostic, and its `T_crit` abort at the quench front is
what the `critical temperature` key added on 2026-08-29 exists to fix.

## 2. Gap table

| # | Gap | Status |
|---|---|---|
| G1 | Owning module of (b) is unknown — circuit vs FEM BDF history | **open, blocks any fix** |
| G2 | Is (b) reachable outside a circuit-coupled deck? Untested on a magnetic-only restart | open |
| G3 | Does (b) exist for BDF1, or only for multistep orders? The RLC deck's order is not pinned | open |
| G4 | The DR-92 cap's own justification has never been re-measured since STRUMPACK changed | open |
| G5 | No regression test exists for restart Δt consistency at any level | open |

---

## 3. Steps

- [ ] **R1 — Localize (b). Blocks everything else.** The controller *does* call
      `set_timestep -> shift -> compute_MNA_matrix` on that step
      (`cl_FEM_Controller.cpp:246-249`), so the circuit's own companions should already be rebuilt
      at the capped Δt. The stale Δt therefore more plausibly enters through the restored **magnetic
      BDF history** (`bdf_last_dt`) feeding the FEM-coupled terminal pair. Decide by instrumenting,
      not by reading: dump `mDeltaTime`, the circuit's `mDeltaTime`, and the IWG's BDF `h` at the
      first post-restart assembly and compare all three. Confidence that it is the FEM side rather
      than the circuit: **medium** — the circuit path was audited twice this week and the coupled
      terminal current is the only channel showing it.
- [ ] **R2 — Decide the contract (after: R1).** Two candidate semantics, and this is a design
      question, not a bug fix:
      (i) the first post-restart step must integrate exactly the Δt the controller advances — the
      cap stays and the history is re-derived for it; or
      (ii) the restored Δt is honoured for one step and the cap is applied from the second — which
      would preserve DR-92's protection only if STRUMPACK's setup is not touched on that step.
      Prefer (i) unless R1 shows the re-derivation is not well posed for the restored order.
- [ ] **R3 — Close the missing-triple hole (a), independent of R1/R2.** After `/meta` is read and
      fields restored, require the complete magnetic triple whenever a warm restart resumes a
      nonzero timestep under an equation whose `timestepping_order() > 1`; in coupled mode require
      the thermal triple too when a thermal equation exists at order > 1. BDF1 and non-multistep
      methods stay exempt. The error must be rank-symmetric, name the offending dump, and fire
      **before** `restore_history_state` so no half-restored BDF state exists.
- [ ] **R4 — Regression gate (after: R2).** Unit level cannot reach this — it needs the controller.
      Add a coupled restart A/B to the example suite: dump at a Δt above `initial timestep`, resume,
      and assert the first post-restart row matches the uninterrupted twin to truncation. The
      dt-matched variant already passes bitwise and makes a good companion assertion.
- [ ] **R5 — Answer G2/G3 (after: R1).** Repeat the measurement on a magnetic-only restart and on a
      deck with an explicit `scheme : bdf1`, to bound the blast radius.

---

## 4. Open questions

- [ ] **O1** — Is the DR-92 cap still needed at all? It was introduced for a STRUMPACK lifetime-setup
      failure. If that no longer reproduces, removing the cap dissolves (b) entirely and is a smaller
      change than re-deriving history. Needs a measurement, not an opinion. (G4)
- [ ] **O2** — Should a restart be allowed to change integrator order at all, or should the loader
      refuse a dump whose order differs from the deck's? Adjacent to R3.

---

## 5. Reproducer

```
# baseline: uninterrupted
examples/RLC_Circuit, simulation time : 8 ms, np=4        -> CircuitResults.txt
# restart pair: identical deck, truncated then resumed
simulation time : 4 ms   -> run, keep memdump.hdf5 AND CircuitResults.txt
simulation time : 8 ms   -> rerun in place, resumes from the dump
```

`init_output_file` truncates while `save_timestep` appends, so **copy leg 1's
`CircuitResults.txt` aside before resuming** or it is destroyed. Discriminator: the first
post-restart row's currents match the baseline row one full dumped-Δt later.

---

## 6. Audit trail

- `tmp/ai_exchange/dr124_circuit_restart_v2.md` — the v2 per-component restore this gate vindicated
- Grok predicted this interaction on the DR-124 row before it was run: "the controller caps the
  first post-restart Δt at `mDeltaTimeInitial` independent of circuit state"
- `devlog/dl20260829_r6_restart_and_phase_sign.md` — the gate run and its two controls

---

## 7. Refuted readings — do not re-derive

- **"The physics is permanently ahead of the clock."** REFUTED 2026-08-29 by measurement:
  re-comparing the capped run against a baseline shifted by the 1.9e-4 s lead makes the fit *worse*
  (452 A vs 137 A). It is one bad step that then decays, not an accumulating offset.
- **"The circuit state is lost on restart."** REFUTED the same day: the dt-matched rerun is bitwise
  identical to the uninterrupted baseline over all 20 rows.
