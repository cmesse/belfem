# Thermal Joule Source: Missing-Heat Investigation and Probe Verdict

**Date:** 2026-08-13
**Purpose:** Record why the tapestack3d quench run showed no temperature rise
after 5 s, what the in-kernel probe proved, and what remains open
**Module:** fem/thermal, fem/kernel (Calculator / MaxwellData), physics/materials

## The observation

`cmake-build-debug/tapestack3d`, 8-proc coupled h-φ + thermal, initial
condition 77 K. At t ≈ 5.2 s, with the postprocessor showing J/Jc up to 1.44
on the tapes, the entire temperature field sat in

    76.97982996530327 < T < 76.99522170776994

— no heating anywhere, and 20–50 mK **below** the initial state. Two
independent anomalies: **A1** no Joule source, **A2** cooling that an
adiabatic problem (this deck has no thermal boundary-condition section)
cannot produce.

## Three-AI round (`tmp/ai_exchange/thermal_joule_missing.md`)

Pre-registered five questions (rho provenance and field transfer, dV scaling,
clamps, the cooling mechanism, air dilution), then dispatched Codex and Grok
in parallel. Both returned independently, and converged: **no static defect
in the Joule path.**

Cleared on read, by all three tracers:

- `T_h_picard` / `T_h_newton` do assemble `f += w·Nᵀ·ρ·|j|²·dV`
  (`mt_thermal_h.cpp:38-49`, `:92-123`), and `IWG_Timestep` adds it to the
  residual as `+Δt·f` (`cl_IWG_Timestep.cpp:880`) — no sign inversion.
- The thermal-side `MaxwellData` binds the magnetic block of the same group
  id (`cl_FEM_Calculator.cpp:66-76`); `link_element_thermal_maxwell`
  (`:3019-3036`) relinks that peer to the same element before assembly, and
  `compute_j` is `C(k)·q()` off the magnetic calculator, reading mesh fields
  the magnetic solve has already written back.
- Layer thickness is carried on **both** sides: the magnetic edge functions
  fold it into `det_J` (`cl_EF_PENTA6TS.cpp:182-187`, `cl_EF_HEX8TS.cpp:365`),
  and the thermal path — which has no edge function — falls back to
  `det(J_Lagrange)` over the geometrically extruded layer nodes
  (`cl_FEM_Calculator.hpp:1998-2027`, `cl_ThinShellFactory.cpp:1358-1387`).
- Clamps cannot bite: `gRhoMin = 0`, `gRhoMax = 1e10`
  (`cl_Communicator.cpp:63-64`), material `mRhoMin = 0`.
- No default thermal BC is injected when the deck has no thermal section —
  `ThermalFactory::create_bc_factory` builds an empty factory
  (`cl_ThermalFactory.cpp:274-283`). The problem is genuinely adiabatic.
- Air is not in the thermal mesh at all (`cl_ThermalFactory.cpp:71-106`), so
  dilution into λ_air is not available as an explanation.

**Both auditors also caught the session reading a deck that had drifted.** The
`input.conf` on disk carried `resistivity type : power-law` and no
`density correction` key — neither of the settings the earlier sessions had
established. That correction is what re-scaled the whole expectation:
power-law dissipation near Jc is orders of magnitude gentler than piecewise
flux-flow.

## The probe settles it

Christian instrumented `mt_thermal_h.cpp` with a `#CHECK` print of
(element, k, T, ρ, |j|, λ, cp, density) at each integration point. The
numbers reconstruct the power law exactly:

| element | T [K] | ρ [Ω·m] | \|j\| [A/m²] | ρ·\|j\| [V/m] |
|---|---|---|---|---|
| P2399590 | 76.985 | 2.93e-15 | 3.310e10 | **0.97e-4 ≈ Ec** |
| P2444032 | 76.984 | 2.43e-15 | 3.389e10 | ≈ Ec |
| P2446058 | 76.998 | 2.63e-15 | 2.919e10 | ≈ Ec |
| P2434609 | 76.983 | 4.73e-37 | 5.677e9  | far below Ec |

Two independent confirmations in one printout:

1. The loaded elements sit at `ρ·J ≈ Ec = 1e-4 V/m` — the definition of
   J = Jc for a power law. They are exactly at critical current, not above it.
2. The cold element reconstructs the exponent: with jc ≈ 3.3e10 and n ≈ 29,
   `(5.68e9 / 3.3e10)^(n-1) · Ec/jc ≈ 5e-37` against the printed 4.73e-37.

**A1 is resolved as physics, not a defect.** The transfer chain is live and
correct. The real finding is that the **table** gives
jc(77 K, |B|, folded θ) ≈ 3.3e10 A/m² — over 3× the 1e10 the deck comment
assumed — so the stack is running *at* its true Ic, not 44 % above it. The
1.44 in the plot is the signed-angle lobe of the same table (DR-69); the
solve samples the folded lobe. Peak dissipation is then ρJ² ≈ 3.3e6 W/m³ in a
1.6 µm layer, and below the knee ρ ~ 1e-37 means the first ~4.9 s produced
**no heat at all**. The expectation of 77.1 K belonged to the piecewise law,
which the drifted deck was not running.

## What stays open: A2

With zero Joule input for the first 4.9 s, the −20 mK cannot be source
related. Codex read `memdump.hdf5` directly and found the band already inside
it (`/fields/T` ∈ [76.9805, 77.0001] at t = 5.23), and `hphiTrun` warm-restarts
from that dump by default (opt-out, `cl_FEM_Controller.cpp:3689-3703`; the deck
carries no `restart : false`). So the cooling is **inherited through the
restart chain**, not produced by the current binary — but its origin is not
identified. All five pre-registered candidates were refuted: no default BC,
savepoint/restore pairing is consistent, BDF adds `+Δt·f`, Anderson and
relaxation cannot invent a colder state from a flat field, and the
`density correction` was not even active.

Registered as **DR-70** — and the register already held its twin: **DR-06**,
open since 2026-07-09, is "thermal-coupled: −1.9 K dip not proven fixed",
with `dV(k) < 0` on wafer elements among its suspects. Same signature class:
a coupled thermal run going below its initial temperature with no sink
available. That suspect is directly testable here, because
`Calculator::dV` guards `adV >= 0.0` with a `BELFEM_ERROR` that is compiled
under debug only (`cl_FEM_Calculator.hpp:1969-1976`) — a release run with an
inverted layer element would integrate a negative mass silently, and a
negative-mass element cooling while its neighbors hold is exactly the
observed band. The two rows are now cross-linked and should be worked as one
investigation. The discriminating experiment is now running:
Christian restored `piecewise` and `density correction : 0.044`, deleted the
memdump, and restarted from scratch. Below the knee the source is ~1e-37, so
an adiabatic fresh run **must hold exactly 77.000** for the first several
seconds. Drift there means the leak is live in the current binary; a flat
77.000 closes A2 as damage inherited from the pre-instrumentation runs (the
era of the DR-08 defect, when diverged thermal solves were accepted silently).

**First reading from the fresh run, and a correction to the test design.**
At t = 50 ms the field reads 77.000 at three decimals, and the log confirms
the physics — `Thermal Picard 1, residual 0.000000 ( -70.64 dB )`, a flat
unforced field. But that reading is **not decisive**: a leak at the old rate
(20 mK over ~5 s ≈ 4 mK/s) would amount to 0.2 mK here and would display
identically. The gate needs full double-precision min/max, or ~1–2 s of
simulated time. Worse, a time-uniform test may miss the mechanism entirely if
the sink is **event-driven** — per rejection or per warm restart, of which the
old campaign had many — in which case the smooth sub-knee phase holds 77.000
and the leak only appears at the knee.

**The spread was under-read at first, and it is the sharper clue.** The old
field spanned 76.9798–76.9952: a 15 mK *range*, not a uniform offset. With
zero source and pure conduction from a uniform initial condition, the field
cannot develop structure at all — so the sink is **spatially localized**. That
fits DR-06's negative-`dV` suspect far better than any global integrator
drift: a few inverted layer elements would integrate negative mass, cool
locally, and conduct the deficit into their neighbors, producing exactly a
graded band.

A second, cheaper discriminator noted by Grok and not yet run: the exodus `T`
field is written over the whole mesh, but air nodes are not thermal dofs and
must sit at exactly 77 K. If the quoted min/max spans the whole domain and the
maximum is below 77, the write path touched non-thermal nodes — a different
bug entirely.

## Mechanism "verified" — then refuted the same session

**Read the retraction below before using anything in this section.** What
follows was correct about the residual floor and wrong about the cause of the
temperature dip.

### The measurement (still valid)

The rerun with `relative tolerance : 1e-10` in `linear thermal` (nonlinear
thermal 1e-7, three decades of headroom) confirms the diagnosis
quantitatively rather than by inference:

| | rtol 1e-8 | rtol 1e-10 | predicted |
|---|---|---|---|
| thermal residual floor | −75.9 dB | −95.6 dB | −20.0 dB shift |
| absolute | 2.57e-8 | 2.75e-10 | ×1/100 |
| floor / rtol | 2.6 | 2.8 | constant |

The floor moved 19.7 dB against the 20.0 dB a 100× tolerance change
predicts, and the floor-to-rtol ratio is unchanged at ~2.7. **The nonlinear
residual floor is the linear solver's tolerance, scaling 1:1 with it.**

The same log carries the control case: the magnetic side, solved directly by
STRUMPACK, reports **−156.5 dB = 2.2e-16** — machine precision. That is the
whole argument for direct solvers in one line, and it is why Messe et al.
2023 could specify a nonlinear ε_n = 10⁻¹¹ while no iterative solver at
default settings can approach it.

Residual drift under the accepted setting is 3.7e-8 K per step: ~0.9 µK at
t = 50 ms and a few tenths of a mK across a full multi-thousand-step run,
comfortably below the several-mK physical signal. The row stays open for the
architectural fix (increment form) only.

### Retraction — itself withdrawn (read to the end of this section)

The retraction below was written on a stale reading and is **wrong**. It is
kept because the sequence is the useful part of the record.

Rerunning at `rtol 1e-10` gave **min T = 76.99991423348413 at t = 50 ms —
bit-for-bit identical to the rtol 1e-8 run, all sixteen digits**, on a freshly
written output file. A hundredfold change in linear tolerance cannot leave a
solver-produced number bit-identical. The dip does not originate in the
linear solve.

What was actually established is narrower than claimed: the *nonlinear
residual floor* tracks rtol one-for-one. That floor was then assumed to set
the temperature error, and it does not. The warning sign was present from the
first run and was explained away rather than pursued — the drift was
3.7e-6 K/step against `rtol × T` = 7.7e-7, a factor of five too large, which
was attributed to condition-number amplification instead of being treated as
a discrepancy that falsified the model.

**The corrected signature is much sharper than the original one:**

    max = 77.0 exactly        min = 76.99991423348413

The bulk of the field is *untouched*. This is a **localized, deterministic,
bit-reproducible cold spot**, not a global drift. The old campaign's
whole-field depression (max 76.9952, min 76.9798) is what that spot becomes
after thousands of steps of conduction spreading it outward.

That signature fits **DR-06** closely: a few elements with inverted Jacobians
would integrate negative mass and act as a fixed sink at a fixed place —
indifferent to solver settings, identical from run to run. The `adV >= 0.0`
guard that would catch it is compiled only under debug
(`cl_FEM_Calculator.hpp:1969-1976`), so a release run integrates a negative
Jacobian in silence.

Next steps, cheap because the defect is bit-reproducible: locate the minimum
(a handful of nodes, or a region?); run this mesh under a debug build so the
`dV` guard can fire; and check whether the cold nodes sit on a thin-shell
layer interface or a particular geometric feature. DR-06 and DR-70 should now
be treated as one defect.

### Re-confirmed by direct file read — the original mechanism holds

Reading the Exodus frames directly (`scipy.io.netcdf_file`) rather than
trusting a value reported out of the viewer:

| t [s] | min T | deficit | max T |
|---|---|---|---|
| 0.05 | 76.99999974013564 | 2.60e-7 | 77.00000001 |
| 0.10 | 76.99999952930314 | 4.71e-7 | 77.00000002 |
| 0.20 | 76.99999918555523 | 8.14e-7 | 77.00000001 |
| 0.30 | 76.99999906458766 | 9.35e-7 | 77.00000000 |
| 0.40 | 76.99999876161793 | 1.24e-6 | 77.00000000 |

Against the old run's **8.577e-5 K at t = 0.05 s**, the new run's 2.60e-7 K is
a **330× improvement** — 100× from `rtol` and roughly 3× more from the
nonlinear tolerance moving 1e-6 → 1e-7. Growth is smooth and roughly linear
in step count; extrapolated to t = 5.2 s it is ~1.6e-5 K, three orders below
the physical signal.

**The bit-identical 50 ms value was the previous run's file still held open
in the viewer** (the filename `hphi_results.e-s.00001` is reused), which is
exactly why it matched to sixteen digits. The "localized cold spot" reading
built on it is withdrawn too: `#nodes below 77` is 374 310 at t = 0.05 s, so
the drift is global, as first described. The DR-06 link is downgraded
accordingly — a debug run checking for `det(J) <= 0` remains worthwhile on
its own merits, but it is not this row's blocking question.

**Methodological note, which is the durable lesson here.** A claim was called
verified on a partial check, retracted on a stale number, and re-confirmed
only when the file was finally read directly — three positions in one
session, none of them separated by new physics. The moment two readings
disagreed, the file should have been opened. Reported values from a viewer
holding a reused filename are not measurements.

The `relative tolerance : 1e-10` + nonlinear 1e-7 pairing is verified
sufficient in practice. The row stays open for the architectural fix
(increment form) only.

## A parser hazard found the same session

The 1e-10 setting was very nearly lost to a syntax trap worth recording,
because it would have falsified a correct hypothesis. Written as

    relative tolerance : 1e-10  // if PETSc;

the terminating semicolon sits *inside* the comment.
`InputFile::remove_comments` truncates each line at `//`
(`cl_InputFile.cpp:67`), and `Section` only registers a statement that
contains a `;` (`cl_Input_Section.cpp:92`) — so the key is **silently
ignored, with no diagnostic**, and the built-in default applies. The run
would have drifted at the old rate and read as "1e-10 did not help".

Since `doc/input_schema.yaml` already enumerates every legal key, a warning
for "line looks like `key : value` but was not registered" is cheap and would
close a whole class of silent deck drift — the same class that put a
`power-law` / no-`density correction` deck under this investigation in the
first place.

## The tolerance-headroom rule (found while checking the 1e-11 citation)

Christian asked whether the paper's ε < 10⁻¹¹ might actually have been the
*linear* solver tolerance. It was not — Messe et al. 2023 §4 is explicit:
`r = Aq − b` (Eq. 10) with `ε = ‖r‖/‖b‖` (Eq. 11, already relative, the same
definition the code uses today), Picard until `ε_p = 10⁻³`, then
Quasi-Newton-Raphson "until a convergence of `ε_n = 10⁻¹¹` is achieved",
required to suppress checkerboarding (current density oscillating between 0
and ±2·j_crit in adjacent elements).

But the question exposed a coupling that is nowhere written down:

**A nonlinear tolerance is only reachable if the linear solve is
substantially tighter than it.** The outer residual cannot go below the noise
floor of the inner solve. The paper's 10⁻¹¹ was free because its linear solves
were direct (~10⁻¹⁴, three decades of headroom). With an iterative solver at
the BELFEM default `rtol = 1e-8` (`cl_SolverParameters.hpp:100`), a 10⁻¹¹
nonlinear target is **unreachable**: the loop grinds to `max iterations` and
is rescued by a timestep cut, with nothing in the output attributing the
stall to the solver stack rather than to the physics.

This is the same family as DR-52 (the Anderson "convergence" that turned out
to be linear-solve roundoff). It is currently latent rather than active: the
magnetic side uses STRUMPACK (direct), so its 1e-7 is safe, and the thermal
side has two decades between its 1e-6 and PETSc's 1e-8. The rule to apply
when changing either: **keep at least two to three decades between the
nonlinear tolerance and the linear one, and tighten the linear tolerance
first.**

Note the two findings are distinct and should not be merged: the paper's
10⁻¹¹ addresses *checkerboarding*, a nonlinear-iteration artifact of the
highly nonlinear ρ(J) law; the drift documented above is a relative tolerance
applied to an absolute-temperature right-hand side. Only the second is a
defect.

## Follow-ups

- The `#CHECK` probe has answered its question and should come out at the next
  rebuild (probe-removal policy).
- The deck comment claiming `jc 1e10 A/m², constant` is now known to be wrong
  by 3.3× against the attached table; the Ic arithmetic in that comment block
  (512 A for the stack) understates the real value correspondingly.
