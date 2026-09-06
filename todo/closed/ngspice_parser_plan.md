# ngspice / SPICE Netlist Parser for the Circuit Module

**Date:** 2026-06-28
**Purpose:** Propose a parser that reads an ngspice (`.cir` / `.net`) netlist and builds a
`belfem::electronics::ElectricalCircuit` through the existing factory `create_*` calls;
minimise the `circuit{}` block in `input.conf`; and add a circuit-only runner that has no
FEM overhead.
**Module:** `src/circuit`
**AIs involved:** Claude (exploration + plan + reconciliation); Codex + Grok (plan audit round 1,
2026-08-27, blind parallel — exchange `tmp/ai_exchange/ngspice_parser_plan_audit.md`)
**Status:** ✅ **v1 COMPLETE 2026-08-27 — promoted from `deferred/`, designed, signed off,
implemented, audited and shipped in ONE DAY.** What landed: `fn_spice_number` +
`cl_NetlistParser` + `cl_NgspiceCircuitFactory` + the `circuit { file : ... }` hybrid path in
`ElectricalCircuitFactory` (both executables inherit it), 74 tests in the `fast` suite, the
Input Contract in both artifacts, and `examples/circuit` converted to `tapestack.cir` —
**running in production hphirun, verified by Christian's `./Allrun`**. Five audit rounds
(1 plan + 4 code, Codex + Grok each time); DR-115 fixed and struck along the way (its np=4
gate passed on the classic deck first); by-catch filed: DR-117, DR-118, DR-121 (switch-
before-V refused in v1 pending the MNA fix), DR-123, DR-124. Design decisions O5–O9 are
FROZEN — see §12. **Residual, post-v1:** Phase 4 (`circuitrun`; ~~+ shared `CircuitSolver`, discharging DR-118~~ — helper leg
superseded 2026-09-04, see §7), Phase 5 (6a one-file mode), Phase 6 (PULSE/PWL, ex-DR-39 item b; constraints carried across 2026-09-04), the
classic terminal-pair deck test, cross-simulation against a real ngspice binary (none on
this host), and the §12.5 open risks.

**Phase-0 decisions (2026-08-27, Christian):**
1. **FEM linkage = 6b** — the `.cir` holds only lumped elements; `terminal pair` blocks stay in
   `input.conf` beside `topology{curves{}}`. Note both variants keep the `.cir` loadable by stock
   ngspice (the extension directives are comments either way); 6b was chosen for separation of
   concerns once that was clear.
2. **BELFEM-only lumped elements = Option A** — `* belfem:` comment directives (superconductor,
   timed switch). Under 6b the `terminal_pair` directive form is not used.
3. **No `.subckt` in v1** — hard-error on `X` cards and `.subckt`; flattening is a later phase if
   a real deck needs it. This re-scopes the `.subckt` leg of DR-39 as the register suggested.
4. **v1 scope = Phases 1–3 + wiring `circuit { file : ... }` into `hphirun`**, converting
   `examples/circuit`'s lumped part to a `.cir`. `circuitrun` (Phase 4) and Pulse/PWL (Phase 6)
   follow later.

**Currentness (2026-08-27, corrected after audit round 1):** still no code — `src/circuit/` has
no `NetlistParser` / `NgspiceCircuitFactory` / `spice_number_to_si` and no `.cir` handling, and
there is no *file-driven* circuit runner. ~~Do not mistake `src/executables/electricalCircuit.cpp`
for one — it pulls in `cl_MaxwellFactory` / `cl_CutFactory`~~ **stale (struck 2026-08-27, both
auditors + tree check): those includes were removed and the executable's includes are circuit/core
only (`electricalCircuit.cpp:13-22`) — it IS a FEM-free driver, just hardcoded, not file-driven**
(the executables `CMakeLists.txt` still *links* the fat FEM library list, which is a link-bloat
note, not a header dependency). The reference PDF (`tmp/ngspice-manual.pdf`) lived in ephemeral,
gitignored `tmp/` — **re-fetch it when implementation starts rather than assume it survived**
(several §8/§10 semantics items were audited from vendor memory, not the manual).
Related debt: `debt_register.md` DR-39 (items (b) Pulse/PWL and (c) netlist parser; the `.subckt`
leg is re-scoped by decision 3 above) and **DR-115 — `examples/circuit` aborts at np=4 on the
first timestep rejection. A FIX gate before the showcase conversion, not waivable (§12.4.1
carries the verified root-cause diagnosis); its `Allrun` must survive its own demo.**

**DR-39 verification sync (2026-08-13 audit, folded in 2026-08-27; refined by audit round 1):**
the Newton loop exists in **two** copies, not the three §7 anticipated —
`Controller::solve_circuit()` (anchor: search `Controller::solve_circuit`; today
`cl_FEM_Controller.cpp:2104-2162`) and `electricalCircuit.cpp:126-156`, same algorithm and
relaxation update, differing in deck-driven vs hardcoded constants. The `tEpsilon0` assignment
sits at the loop TOP in the controller (`:2122`) and the BOTTOM in the executable (`:154`) — a
**structural** difference that Grok's unrefuted equivalence argument shows does not change the
ω update (both compare the current residual against the same previous value; `residual()` holds
no hidden state, `cl_ElectricalCircuit.cpp:900-903`). The Phase-4 `CircuitSolver` extraction
should still pick one site deliberately, but this is a code-unification concern, not a numerical
drift. All line citations in this file are advisory — §12.3 has the corrected anchor table.

---

## 1. Motivation and the three questions

The user asked three questions; short answers first, with details below.

1. **Can we pull the circuit definition from an ngspice file instead of `input.conf`?**
   **Yes, for the lumped part.** R, L, C, independent V/I sources and diodes map almost
   one-to-one onto existing BELFEM components and the existing `create_*` API. The only
   things that have *no* standard SPICE element are the FEM-coupling `terminal pair` and
   the lumped `superconductor`; those need a BELFEM-specific extension (Section 5).

2. **Can we minimise the circuit definition in `input.conf`?**
   **Yes — hybrid model.** Keep one line in `input.conf` pointing at the netlist
   (`circuit { file : magnet.cir ; }`), move all lumped topology into the `.cir`, and keep
   only the FEM-linkage that SPICE cannot express (terminal-pair ↔ curve/terminal mapping)
   either in `input.conf` or in a BELFEM extension block inside the `.cir`. See Section 6.

3. **Can we run a circuit model without the FEM overhead?**
   **Yes — and one already exists.** `src/executables/electricalCircuit.cpp` is a
   FEM-free driver that builds an `ElectricalCircuit` with hardcoded `create_*` calls and
   runs its own time loop (Newton + adaptive Δt). The remaining work is to make it
   **file-driven** (read the topology from a netlist or a trimmed `input.conf`) rather than
   recompiling for each circuit. Note: `Controller::solve_circuit()`
   (`cl_FEM_Controller.cpp:950-1004`) is only the *Newton loop*; the full per-step driver
   (`set_timestep → shift → compute_MNA_matrix → solve_circuit`) lives in
   `initialize_timestep()` (`:138-156`). See Section 7 for the corrected loop.

---

## 2. What we are mapping onto (recap of the existing API)

The factory builds components by calling these `ElectricalCircuit` methods with **raw SI
`real` values** (the unit parsing happens earlier, in the factory). A netlist importer can
call the same methods directly and skip the unit machinery entirely:

| BELFEM call (`cl_ElectricalCircuit.hpp`) | Args |
|------------------------------------------|------|
| `create_resistor(value, n+, n-, label)` | R [Ω] |
| `create_capacitor(value, order, n+, n-, label)` | C [F], BDF order |
| `create_inductor(value, order, n+, n-, label)` | L [H], BDF order |
| `create_voltage_source(SourceFunction*, n+, n-, label)` | excitation |
| `create_current_source(SourceFunction*, n+, n-, label)` | excitation |
| `create_switch(isClosed, switchTime, n+, n-, label)` | — |
| `create_diode(Is, Vt, n+, n-, label)` | — |
| `create_superconductor(Ic, n, Ec, length, n+, n-, label)` | — |
| `create_terminal_pair(n+, n-, label)` | + FEM BCs (factory only) |

Nodes are `index_t` indices `0 … number_of_nodes-1`, **last = ground**.

---

## 3. SPICE element → BELFEM component mapping

SPICE identifies an element by the first letter of its instance name (case-insensitive).
General forms (from `tmp/ngspice-manual.pdf`, Chapters 3–4):

| SPICE card | General form | BELFEM target | Notes |
|------------|--------------|---------------|-------|
| Resistor `R` | `Rxxx n+ n- value` | `create_resistor` | value may be `r=...`; ignore `ac=/m=/tc1=/temp=` initially |
| Capacitor `C` | `Cxxx n+ n- value [ic=...]` | `create_capacitor` (default order 1) | `ic` ↦ initial condition (optional, see §8) |
| Inductor `L` | `Lxxx n+ n- value [ic=...]` | `create_inductor` (default order 1) | `ic` ↦ initial condition |
| Voltage source `V` | `Vxxx n+ n- [[DC] val] [AC ...] [TRAN_FUN]` | `create_voltage_source` | transient function → `SourceFunction` (§4) |
| Current source `I` | `Ixxx n+ n- [[DC] val] [AC ...] [TRAN_FUN]` | `create_current_source` | same |
| Diode `D` | `Dxxx n+ n- mname` + `.model mname D(is=... )` | `create_diode(Is, Vt)` | `Vt = n·k·T/q`; from model `n`+`temp`, else default 0.026 |
| Switch `S`/`W` | voltage/current-controlled switch + `.model` | ~~`create_switch`~~ **hard-error (audit round 1)** | SPICE switches are *controlled*, BELFEM's is *timed* — a semantic mismatch, not a mapping. Per Phase-0 decision 2 the timed switch comes ONLY from a `* belfem: switch` directive; a real `S`/`W` card must hard-error |
| Subckt `X` | `Xxxx nodes subname` | ~~flatten or map~~ **hard-error (Phase-0 decision 3)** | no `.subckt` in v1 |

**No SPICE-standard equivalent** (need BELFEM extension, Section 5):
- `terminal pair` (FEM coupling — needs curve/terminal ids).
- `superconductor` (lumped E–J power law — has Ic/n/Ec/length).

### Control / analysis cards

| SPICE card | Maps to |
|------------|---------|
| `.tran tstep tstop [tstart [tmax]] [uic]` | ~~timestep config (initial = tstep, …)~~ **corrected (audit round 1, both auditors): `tstep` is the PRINT/plot increment in SPICE, not the integrator's initial step; `tmax` bounds the internal step.** v1 `hphirun` **ignores `.tran` with a log line** — `solver{timestep{}}` owns the timestep (O7). The Phase-4 standalone runner may interpret `.tran`, mapping `tmax`→max Δt, `tstop`→simulation time, and choosing its own initial Δt (never `tstep`) |
| `.ic v(n)=...` | initial node voltages (optional, §8) |
| `.model name TYPE(...)` | parameter table consumed by `D` (and future `R`/`C`/`L`/`SW` models) |
| `.end`, `.ends`, `.subckt` | structural |
| `.param`, `.include`, `.lib`, `.global` | **out of scope for v1** — log + ignore or error |
| First line | **title**, always ignored by SPICE (must be consumed, not parsed) |

---

## 4. SPICE transient functions → `SourceFunction`

`SourceFunction` (`numerics/sources/cl_SourceFunction.hpp`) already covers most SPICE
source shapes. Mapping:

| SPICE | Args (SPICE order) | `SourceFunction` |
|-------|--------------------|------------------|
| `SIN(Vo Va Freq Td Theta Phase)` | offset, amp, freq, delay, damping, phase | `set_periodic(Sine, Va, 1/Freq, Phase)` — **offset `Vo`, delay `Td`, damping `Theta` unsupported**, warn if nonzero. SPICE `Phase` is in **degrees** → convert to radians (`set_periodic` takes rad). Note BELFEM `sine` reads `phase()` rather than `time_offset()` (`cl_SourceFunction.hpp:359-363`), so SPICE `Td` cannot be mapped to a time shift for `SIN` — square/triangle/sawtooth *do* use `time_offset`. Both slots are written by `set_phase`, so the `Phase` argument itself works on every waveform (DR-117, fixed 2026-08-29) |
| `PULSE(V1 V2 Td Tr Tf Pw Per)` | — | nearest match `square`/`sawtooth`/`triangle`, or extend `SourceFunction` with a true pulse (recommended) |
| `DC value` / bare value | constant | `set_constant(value)` |
| `PWL(t1 v1 t2 v2 ...)` | piecewise linear | **no equivalent** — add a PWL `SourceFunction` or load via `userdefined` |
| `EXP(...)`, `SFFM(...)`, `AM(...)` | — | out of scope v1; warn |

Gaps to note honestly: BELFEM `sine` has no offset / delay / damping term; SPICE `SIN`
has all three. v1 should **warn loudly** when a SPICE source uses parameters BELFEM
cannot represent rather than silently dropping them. A clean follow-up is to extend
`SourceFunction` with `Pulse` and `PWL` types (they are generally useful, not just for
the importer).

---

## 5. The two BELFEM-only elements

These cannot be expressed in standard SPICE. There are two options; **RESOLVED 2026-08-27,
Christian → Option A** (comment directives). With linkage decision 6b the `terminal_pair`
directive below is unused — the directives cover `superconductor` and the timed `switch` only.

**Option A — BELFEM extension via comment directives.** ngspice ignores lines starting
with `*`. Encode BELFEM elements in a structured comment the importer recognises:

```
* belfem: terminal_pair Z1 n+=0 n-=1 input_curves=[1,2] output_curves=[3,4] order=1   <- unused under 6b
* belfem: superconductor SC1 n+=0 n-=1 Ic=200 n=20 Ec=1e-4 length=1
* belfem: switch S1 n+=1 n-=2 state=open t_switch=20m
* belfem: order L1 2                                      <- BDF order override (O8, decided 2026-08-27)
```

Real ngspice treats these as comments (so the lumped subset still simulates there);
BELFEM's parser harvests them.

**Option B — non-standard element letters.** Use unused/extended letters (e.g. `Z` for
terminal pair, `P` for superconductor) with a fixed positional syntax. Cleaner to read but
makes the file non-loadable by stock ngspice.

The timed `switch` is a genuine semantic mismatch (SPICE switches are *controlled*, BELFEM
switches are *timed*); represent it with the extension directive rather than pretending a
SPICE `S` card maps over.

---

## 6. Minimising `input.conf` (hybrid file model)

Today `circuit{}` carries the whole topology. Proposed minimal form:

```
circuit
{
    file : magnet.cir ;          // ngspice netlist holds the lumped topology
    output { file : CircuitResults.txt ; currents : ... ; voltages : ... ; }
}
```

Decision needed: **where does the FEM linkage live?** Two coherent choices:

- **6a — linkage in the `.cir`** via `* belfem: terminal_pair ...` directives (Option A).
  Fully self-contained; one file. But mixes FEM-mesh concepts (curve ids) into a circuit
  file.
- **6b — linkage stays in `input.conf`.** The `.cir` holds only lumped elements; the
  `terminal pair` blocks remain in `input.conf` (referencing `topology{curves{}}` as
  today). Cleanest separation of concerns; the `.cir` is pure SPICE.

**RESOLVED 2026-08-27, Christian → 6b.** The deciding realisation: 6a's ngspice-loadability
argument was moot — Option A directives are comments, so *both* variants keep the `.cir` valid
for stock ngspice — leaving 6b's separation of concerns (curve ids stay next to
`topology{curves{}}`, one file to drift instead of two) unopposed. The parser must still
hard-error when `input.conf` names a curve or node that contradicts the netlist. Either way
`number of nodes` is *derived* from the netlist (per §8: distinct non-ground names + 1, not max
index + 1), and the SPICE ground node `0` must be remapped (Section 8).

---

## 7. Standalone circuit runner ~~(extend `electricalCircuit.cpp`)~~

> **RETARGETED 2026-09-04 (DR-39 closing plan, executing Christian's ruling of 2026-08-30).**
> `src/executables/electricalCircuit.cpp` was **deleted** on 2026-09-04 after its four circuits
> were harvested into `tests/circuit/test_ElectricalCircuit.cpp` and run green. The "existing
> FEM-free runner" this section builds on no longer exists, and the shared-`CircuitSolver`
> shape proposed at the end of it is **refuted for good**: `Controller::solve_circuit()` is the
> single circuit Newton loop in the tree, and a helper shared between one production caller and
> a demo buys a seam, not safety. The *need* for a file-driven runner survives, unfunded (Phase 4).
> If it is ever funded, its mechanism is a `circuitrun` that runs the controller's own sequence
> directly — `set_timestep → shift → compute_MNA_matrix → Newton` — exactly as the test helpers
> `take_step`/`solve_attempt` already do; not a third copy of the loop and not a helper extracted
> from the controller. The rest of this section is kept as written for its audited ordering facts.

~~**A FEM-free runner already exists:**~~ `src/executables/electricalCircuit.cpp`. It is the
correct reference loop (not `solve_circuit()`, which is only the Newton inner loop). The
task is to make it **file-driven** instead of hardcoded. Its actual structure
(`electricalCircuit.cpp:96-182`), which the importer must preserve:

> **CORRECTED (audit round 1, 2026-08-27 — both auditors, tree-verified).** The listing below is
> **not** faithful to the executable, and the executable's own reject path is defective; do not
> port either verbatim:
> 1. The real code has **no `continue`** after a rejected step — after `shift_back()` + halving
>    Δt it falls through to `tTime += tDeltaTime; shift()` (`electricalCircuit.cpp:160-180`),
>    i.e. it advances time on a FAILED step without re-solving. The `continue` below is the
>    *correct* behaviour, silently repairing a live defect while claiming fidelity.
> 2. A **first-step rejection** calls `shift_back()` with no prior `shift()` —
>    `ShiftRegister::revert` errors out (`cl_ShiftRegister.hpp:323-325`).
> 3. The reject/grow path below re-stamps via `set_timestep()` + `compute_MNA_matrix()` with no
>    `shift()` in between — but `set_timestep` writes the **BDF1** companion
>    (`cl_Inductor.cpp:81-85`, `mRL = L/Δt`) and only `shift()` writes the true BDF-order stamp
>    (`:55-56`). For order-2 L/C this silently mis-stamps the MNA. Tree-verified.
> 4. **v1 does not touch this loop at all** — `hphirun` keeps the controller's
>    `set_timestep → shift → compute_MNA → solve_circuit` sequence, which stamps correctly for
>    any order. This whole section is Phase-4 material, and the Phase-4 runner must fix (1)–(3),
>    not preserve them.

```cpp
// Phase-4 reference loop — corrected semantics, NOT a faithful copy of electricalCircuit.cpp
ElectricalCircuit * tCircuit = NgspiceCircuitFactory( "magnet.cir" ).circuit();
real tTime = 0.0, tT_end = ... , tDt = ... ;        // from .tran
tCircuit->set_timestep( tDt );
tCircuit->compute_MNA_matrix();                      // once; re-done ONLY when Δt changes
while ( tTime < tT_end )
{
    write_output_row();                              // current solution (prev step)
    real tEps = BELFEM_REAL_MAX, tEps0 = BELFEM_REAL_MAX, tOmega = 1.0 ;
    uint tIt = 0 ;
    while ( tIt++ < tMaxIt && tEps > tTol )          // Newton loop (== solve_circuit body)
    {
        tCircuit->set_omega( tOmega );
        tCircuit->compute_jacobian_and_rhs();
        tCircuit->solve();
        tEps = tCircuit->residual();
        /* adaptive omega, electricalCircuit.cpp:146-155 == cl_FEM_Controller.cpp:986-994 */
        tEps0 = tEps ;
    }
    if ( tIt == tMaxIt )                             // non-convergence: reject + shrink
    {
        tCircuit->shift_back();
        tDt *= 0.5 ; tCircuit->set_timestep( tDt );
        tCircuit->compute_MNA_matrix();              // Δt changed -> restamp L/C
        continue ;
    }
    if ( tIt < tMaxIt/2 && tDt != tDtMax )           // fast convergence: grow Δt
    {
        tDt = min( tDtMax, tDt*1.5 ); tCircuit->set_timestep( tDt );
        tCircuit->compute_MNA_matrix();              // Δt changed -> restamp L/C
    }
    tTime += tDt ;
    tCircuit->shift();                               // advance time + push BDF history
}
```

Key ordering facts (audited):
- `compute_MNA_matrix()` must be re-run **whenever Δt changes** because `shift()` updates
  the L/C companion resistances (`cl_Inductor.cpp:55-56`, `cl_Capacitor.cpp:56-57`) and the
  MNA re-stamps them (`cl_ElectricalCircuit.cpp:577-603`). With **constant Δt** the
  companion resistance is constant, so the existing runner recomputes MNA only on a Δt
  change — the FEM controller instead recomputes it every step (simpler, slightly
  redundant). Either loop strategy is correct; a one-shot MNA build is safe only if Δt is
  fixed for the whole run.
- **Phasing differs between the two existing loops:** `electricalCircuit.cpp` writes
  output, solves at the current state, then `shift()`s at the *end*; the controller
  `shift()`s *before* the solve (`cl_FEM_Controller.cpp:140-142`). Pick one deliberately
  and document it; do not mix them.
- `shift_back()` on non-convergence and an output-file header (`init_output_file`) before
  the first `save_timestep()` are both required.

Notes / constraints:
- **Rank-0 only / serial.** The circuit solver is not distributed
  (`cl_FEM_Controller.cpp:131-161`). `circuitrun` should run on one rank (guard or assert
  `comm_size()==1`).
- A netlist with a `terminal pair` cannot run standalone (no FEM to provide I/V). Either
  reject it with a clear error, or allow a `* belfem: terminal_pair ... as=current_source
  value=...` test stub that substitutes an independent source.
- The Newton/omega control constants (`mAlpha/mBeta/mGamma`, max iterations, tolerance)
  should be read from a small `solver{}`/`.options` block so the runner matches the
  coupled solver's behaviour. The Newton loop currently exists in **two** copies
  (`Controller::solve_circuit()` at `cl_FEM_Controller.cpp:2104-2162` and
  `electricalCircuit.cpp:126-156`); the new runner would be a third — factor it into a
  reusable `CircuitSolver` helper so all share one implementation instead of drifting.

This is independently useful for **unit-testing the circuit module** (e.g. an RLC ring-down
with an analytic solution) without spinning up a mesh.

---

## 8. SPICE semantics we must handle correctly

- [ ] **Ground / node 0.** SPICE: node `0` is ground; ngspice auto-converts `gnd` to `0`
  (unless `no_auto_gnd` is set) and treats node names as strings, so `00` is NOT ground.
  **`GROUND`/`ground` alias dispute RESOLVED 2026-08-27 against the re-fetched manual (master
  tree, "Ground node" subsection): only `0` and `gnd` — Codex was right; v1 treats `ground` as
  an ordinary node name.** The manual also states every circuit must have a ground node, which
  confirms the hard-error-on-groundless-deck rule. Ground aliases are resolved as node *names*
  through the map, never by integer parsing.
  BELFEM: the **last** index is ground. Build a node-name → BELFEM-index map that sends the
  ground alias to `number_of_nodes-1` and packs the rest into `0 … N-2`. (Alternative:
  refactor `ElectricalCircuit` to use index 0 as ground — larger blast radius, deferred.)
- [ ] **Node count, not max index.** `number_of_nodes` = (count of **distinct non-ground**
  node names) **+ 1** for ground — *not* `max integer index + 1` (that allocates holes for
  sparse numbering and is meaningless for alphanumeric names).
- [ ] **Implicit / injected ground.** The parser always reserves a BELFEM ground slot —
  **but (audit round 1, Codex) a deck that never references ground at all is not valid ngspice
  and would produce a singular MNA; hard-error on a groundless deck instead of silently
  injecting one.**
- [ ] **Node names.** SPICE nodes can be alphanumeric (`in`, `out`, `n12`), not just
  integers. Need a string→index map, not `atoi`.
- [ ] **SPICE number scaling suffixes.** `T G MEG K M(=milli) U N P F`, plus `MIL`
  (=25.4e-6) and `A`(=1e-18, atto) found in the local ngspice manual, and engineering
  forms. The classic trap: `M` = milli (1e-3) while `MEG` = 1e6 (`2.5e-6`, `1k`, `100u`,
  `10n`, `1MEG`). Trailing unit letters after the suffix are ignored by SPICE
  (`1kOhm` == `1k`). Write one `spice_number_to_si()` helper with thorough tests; this is
  the most error-prone single piece.
- [ ] **Case-insensitivity** of element letters, keywords, node names, suffixes.
- [ ] **Line continuation** with leading `+`; end-of-line comments (`$`, `;`, and `//`);
  full-line comments (`*`); the mandatory title first line; the required `.end` terminator.
- [ ] **`r=` / `c=` / `l=` keyword vs positional** value forms — these are the ngspice keywords;
  `value=` is **not** standard (audit round 1, Grok) — support it only as a documented BELFEM
  alias, or not at all.
- [ ] **Initial conditions** (`ic=`, `.ic`, `uic`) — would need new state/history setters on
  the components (none exist today); for v1 warn that they are ignored.
- [ ] **Unsupported elements must hard-error, not silently drop.** Controlled sources
  `E`/`F`/`G`/`H`, behavioral `B`, mutual inductance `K`, transmission lines `T`/`O`/`U`,
  and all semiconductor devices beyond the diode have no BELFEM equivalent — reject with a
  clear message naming the offending card.

---

## 9. Proposed code structure

New files under `src/circuit/` (keep it inside the `belfem::electronics` namespace):

```
src/circuit/
├── cl_NetlistParser.{hpp,cpp}      // text  → intermediate representation (IR)
├── cl_NgspiceCircuitFactory.{hpp,cpp}  // IR → ElectricalCircuit via create_* calls
├── fn_spice_number.{hpp,cpp}       // spice_number_to_si( const string & ) -> real
└── doc/ngspice_import_guide.md     // user-facing guide once implemented
```

and a runner:

```
src/executables/circuitrun.cpp      // standalone, FEM-free circuit time loop
```

### 9.1 `NetlistParser` (lexing + structural parse)

Responsibilities (pure text → data; no `ElectricalCircuit` knowledge):

- [ ] Read file; strip title line; strip/normalise comments; join `+` continuations.
- [ ] Tokenise each card; classify by first letter / leading `.`.
- [ ] Emit an **IR**: `Cell<ElementRecord>`, `Cell<ModelRecord>`, `Cell<ControlRecord>`,
  `Cell<BelfemDirective>` (from `* belfem:` lines).
- [ ] `ElementRecord { char tType; string tName; Cell<string> tNodes; Cell<string>
  tValues; Map<string,string> tKwargs; }` (strings at this stage; numeric conversion is
  the factory's job, so the parser stays dumb and testable).
- [ ] Collect the set of node names for the renumbering map.

Follow BELFEM conventions: `Cell<T>` for lists (not `std::vector`), `a/t/m` naming,
`BELFEM_ERROR` for malformed input (always-on), `BELFEM_ASSERT` for internal invariants.

### 9.2 `NgspiceCircuitFactory` (IR → circuit)

- [ ] Build the node map (ground = `0`/`gnd` → last index); set `number_of_nodes`.
- [ ] Construct `ElectricalCircuit`.
- [ ] For each `ElementRecord`, convert values with `spice_number_to_si`, resolve any
  referenced `.model`, and call the matching `create_*`.
- [ ] ~~For each `BelfemDirective`, create terminal pairs (+ push FEM BCs …) or superconductors
  / timed switches.~~ **Superseded by 6b (audit round 1):** directives create ONLY
  superconductors and timed switches. Terminal pairs and their FEM BCs are created by the
  deck-side reader from `input.conf`, resolving node references through the netlist's node map
  (O5) — directives never touch FEM coupling.
- [ ] ~~Expose `circuit()` … as a drop-in replacement so `hphirun` can choose `input.conf` *or*
  `.cir`.~~ **Superseded by 6b (audit round 1): under `circuit { file : ... }`, `hphirun`
  consumes BOTH files** — netlist for lumped topology, deck for terminal pairs / output / FEM
  linkage — feeding one `ElectricalCircuit` instance. "Drop-in" was also underspecified against
  the real entry point: the existing factory ctor takes the input *path*, owns an `InputFile`,
  owns and **deletes** the circuit in its dtor, and `circuit()` returns base `Circuit*`
  (`cl_ElectricalCircuitFactory.hpp:47,58`, `.cpp:23,41`). Ownership/sequencing contract = O6.

Reuse, don't fork: this factory and `ElectricalCircuitFactory` should converge on one set
of `create_*` calls. Consider extracting the per-component construction into shared helpers
so the two front-ends (BELFEM `circuit{}` and SPICE `.cir`) cannot drift apart.

### 9.3 Testing

- [ ] Unit tests for `spice_number_to_si` (the suffix/`M`-vs-`MEG` trap, units glued to
  suffix, scientific notation, signs).
- [ ] Parser tests on small decks (R divider; RLC; V-source + SIN; diode + `.model`).
- [ ] An end-to-end RLC ring-down in `circuitrun` checked against the analytic solution.
- [ ] Cross-check the lumped subset against **stock ngspice** on the same `.cir` (the
  Option-A comment encoding makes this possible).

---

## 10. Phasing / checklist

- [x] **Phase 0 — design sign-off.** Done 2026-08-27 (Christian): **6b** for FEM linkage,
  **Option A** for extension encoding, **no `.subckt` in v1**, ground stays a remap layer
  (plan recommendation, unchallenged — reopen only if the remap turns ugly in Phase 3).
- [x] **Phase 1 — `spice_number_to_si` + tests.** DONE — implemented AND code-audited 2026-08-27
  (`src/circuit/fn_spice_number.{hpp,cpp}`, `tests/circuit/` suite — 89 expectations in 21
  tests, `fast`-labeled, wired into `check` AND `check-fast` foreach lists plus
  `set( LIBLIST circuit )`). Semantics verified against the re-fetched manual (master clone in
  `tmp/ngspice-manuals/`; `ground`-alias dispute resolved → only `0`/`gnd`). **Code-audit
  round complete (Codex + Grok), all findings fixed:** post-scale overflow check ("1e300T"),
  missing LIBLIST (would not have linked — Grok caught what Codex passed), locale-independent
  `std::from_chars` conversion (per `fn_GT_parse.hpp` precedent) with a guarded de_DE
  regression test, magnitude policy pinned toolchain-independently (zero or normal double
  range, before and after scaling — subnormals reject loudly), A-is-atto-not-Ampere warning
  + tests, boundary battery. **Verified by execution: scratchpad probe mirroring the suite,
  88/88 pass** (probe at `tmp/probe_spice/`; negative control demonstrated); both TUs
  syntax-clean. **GATE PASSED, verified by execution: Christian's in-tree `make check`,
  15/15 tests green incl. the new `circuit` suite (fast label 10 tests, 8.9 s), 2026-08-27.**
  Exchange: `ngspice_phase1_code_audit.md`.
- [x] **Phase 2 — `NetlistParser`** (R/L/C/V/I/D + `.model` + `.tran`), IR + tests. DONE.
  Implemented 2026-08-27: `src/circuit/cl_NetlistParser.{hpp,cpp}` (title line, all comment
  forms, `+` continuations, `(`/`)`/`,` normalization, `=` spacing tolerance, case-folding,
  first-appearance node list per O9, `* belfem:` directives incl. the O8 `order` form, hard
  errors per §12.2 with source:line:card context — `S`/`W`/`X` and every non-R/C/L/V/I/D
  letter, every control card except `.model`/`.tran`/`.end`, `.end`-with-trailing-tokens,
  malformed/duplicate kwargs). `.end` optional at EOF and comment-tolerant continuations are
  documented BELFEM extensions (the manual is self-contradictory on the latter — reconciled
  in the exchange). **Code-audit round complete (Codex + Grok), all findings resolved:** the
  real defect was Grok's F1 — directive-carried `n+`/`n-` nodes missed the O9 packing list;
  fixed with a written-order harvest + regression test. Also: `.endc`/`.ends`/`.endif`
  locked by test (P0 blast radius if the matcher is ever "simplified"), `* belfem :`
  colon-space recognized, `r==5`/duplicate-kwarg hard errors via shared `store_kwarg`,
  error-message contract asserted, on-disk `Ascii` ctor smoke test. **Verified by execution:
  real gtest TUs standalone — 47/47 pass (26 parser + 21 Phase-1); GATE PASSED same day:
  Christian's in-tree `test_circuit` binary, 47/47 green.** Exchange:
  `ngspice_phase2_code_audit.md`.
- [x] **Phase 3 — `NgspiceCircuitFactory`** for the lumped subset; node remap. DONE.
  Implemented 2026-08-27: `src/circuit/cl_NgspiceCircuitFactory.{hpp,cpp}` — O9 node map
  ("0"/"gnd" one node, LAST index; others first-appearance-packed, logged at Verbose),
  components created in netlist line order (elements + directives merged by line, freezing
  the unknown-current layout), labels = folded instance names, `node_index()` as the O5
  deck-side lookup, produce-and-hand-over ownership. R/C/L (value XOR keyword, nonzero,
  `ic=` refused), V/I (DC/bare/SIN with vo/td/theta==0 enforced, degrees→radians), diode
  (.model d, is/n only), superconductor/switch/order directives with full key validation.
  **Verified by execution, standalone real-TU build: 59/59, including a SOLVED DC divider
  (v(mid)=5.0 — V-source first-node-is-plus confirmed) and a SOLVED current-source polarity
  circuit (+1000 V — the SPICE n+→through→n− convention CONFIRMED, closing the §12.5
  polarity risk by gate instead of cross-simulation).** **Code-audit round complete
  (Codex + Grok, complementary — each caught a defect the other passed): Codex found the
  parser letting directive nodes jump the O9 queue past a pending card (fixed — directives
  are now statements that flush the card and end continuation chains); Grok found a
  PRE-EXISTING P0 in `ElectricalCircuit::compute_MNA_matrix` (V-source unknown-current
  column miscounted behind a switch — filed as DR-121, reachable from the deck path too;
  the netlist factory REFUSES switch-before-V layouts until it is fixed). Also fixed:
  `unique_ptr` ownership (exception-safe, non-copyable), diode/superconductor positivity
  guards, single-digit order token (kills a uint wraparound), case-folding `node_index()`,
  parser `source()` accessor. Suite 59 → 66 tests, 66/66 on the standalone real-TU build.**
  **GATE PASSED, verified by execution: Christian's in-tree `test_circuit`, 66/66 green
  (SUPERLU solves visible in both analytic circuits), 2026-08-27 — Phases 2 and 3 both
  closed on the same run.** Exchange: `ngspice_phase3_code_audit.md`.
- [x] **Phase 3b — v1 delivery (after: 1–3).** DONE. Wiring implemented 2026-08-27: the `file :`
  branch lives INSIDE `ElectricalCircuitFactory` (`read_circuit_from_netlist`), so `hphirun`
  and `hphiTrun` get it without touching either executable (O6 by construction); the classic
  TERMINALPAIR case and output block were refactored into shared `read_terminal_pair` /
  `read_output` helpers; hybrid node references resolve as netlist node NAMES (O5), currents
  fold, `.tran` logged-and-ignored on rank 0 (O7); guards: `number of nodes` forbidden,
  lumped-in-topology hard error, topology optional. Input Contract updated in BOTH artifacts
  (schema `file` key + §11.1 in the reference). **Audit round complete (Codex + Grok, converged on
  every finding, all fixed same session):** hybrid mode is now case-insensitive THROUGHOUT
  (the one-sided currents fold both vendors caught as P1 — deck pair labels fold at
  creation), cross-source duplicate labels hard-error, `unique_ptr` factory members close
  the ctor-throw leak, node-map logging rank-guarded, schema machine-fields made
  mode-conditional, and a rank-0 decoder line maps voltage names to CSV columns (Grok F6).
  **Verified by execution: 74/74 on the standalone real-TU build (8 hybrid tests incl. the
  solved netlist divider, a mixed-case save_timestep round trip, and the classic-path
  regression divider).** **GATE PASSED, verified by execution: Christian's in-tree `make check` green with the
  74-test circuit suite and the new fem-kernel LIBLIST chain, 2026-08-27.** **DR-115 gate PASSED on the classic deck
  (Christian, np=4 through a timestep rejection) — row struck and archived — and the
  showcase CONVERTED: `examples/circuit/tapestack.cir` + the slim `circuit { file }`
  section. Pre-flight verified by execution against the installed production deck (factory
  probe): 2 nodes with the IDENTICAL index mapping to the old deck (live→0, ground→1), the
  same five components (netlist is/l1/c1 + deck pairs z1/z2), two BCs in deck order — the
  FEM coupling and CSV columns are continuous with the pre-conversion run.** **FINAL GATE PASSED, verified by execution
  2026-08-27: Christian's `./Allrun` on the converted showcase runs — hphirun builds its
  circuit from `tapestack.cir` in production.** (The first attempt aborted on a STALE
  hphirun binary — `make check` rebuilds tests, not solver binaries; one `make` fixed it.)
  Residual: the classic terminal-pair deck test (named coverage gap). Exchange:
  `ngspice_phase3b_code_audit.md`. Wire `circuit { file : ... ; }` into `hphirun`
  (terminal pairs stay in `input.conf` per 6b) with the `* belfem:` directives for
  superconductor / timed switch; convert `examples/circuit`'s lumped part to a `.cir`; update
  `doc/input_file_reference.md` + `doc/input_schema.yaml` in the same session (Input Contract).
  **Sequencing gate: DR-115 must be FIXED first — a waiver is not available, because `Allrun`
  defaults to a multi-rank count and the showcase would abort under its own demo (§12.4.1, which
  also carries the verified root-cause diagnosis). Gates on O5–O8 sign-off (§12.1).**
- [ ] **Phase 4 — `circuitrun`** standalone runner. ~~RLC analytic test. Factor the Newton
  loop out of `Controller::solve_circuit()` into a shared helper — resolving the
  top-vs-bottom `tEpsilon0` drift (header note) deliberately.~~ **Retargeted 2026-09-04 (see the
  §7 note): the shared-helper leg is SUPERSEDED by the DR-39 ruling — the demo executable is
  deleted, `Controller::solve_circuit()` is the single loop, and the `tEpsilon0` placement was
  structural, not numerical. The RLC analytic test already exists
  (`RLCRingAcrossRejectedStep`, plus the four harvested transients beside it). UNFUNDED; a
  runner, if built, drives the controller's sequence directly.**
- [ ] **Phase 5 — remaining extensions:** 6a-style one-file mode (optional), any further
  BELFEM directives beyond superconductor/switch.
- [ ] **Phase 6 — `SourceFunction` extensions** (`Pulse`, `PWL`) to close the SPICE source
  gap and ex-DR-39 item (b); user guide in `src/circuit/doc/`. Today's user-visible behavior is
  the parser's explicit refusal of both by name (`cl_NgspiceCircuitFactory.cpp`, search `PULSE`).
  **Three constraints carried across from the struck DR-39 on 2026-09-04 — tree-checked
  2026-08-30, do not rediscover them:**
  1. **Not circuit-scoped.** `boundary_condition_function_type()` is the deck-string mapper for
     `cl_MaxwellBoundaryConditionFactory.cpp` and `cl_ThermalBoundaryConditionFactory.cpp` as
     well as `cl_ElectricalCircuitFactory.cpp`, so adding the types widens the accepted `type :`
     value set for EVERY boundary condition. That is an Input Contract change in **both**
     `doc/input_file_reference.md` and `doc/input_schema.yaml` (additive, so post-freeze holds).
  2. **The parameter store collides, it is not merely short.** `mValues` is a fixed 7-slot
     `Cell< real >` and `synch()` is `broadcast( mValues )` (`cl_SourceFunction.cpp`); all seven
     slots are already named `BELFEM_BCVAL_AMPLITUDE … FUZZYNESS` (`cl_SourceFunction.hpp`), so
     PULSE's seven parameters collide with that schema — "add slots" is the wrong mental model.
     **PWL does not fit at all:** its variable-length breakpoint list is exactly the payload the
     MPI rule routes through `share`/`receive`, never `broadcast`.
  3. **PWL storage must be `Cell< real >`.** The header forbids `cl_Vector.hpp`/`cl_Matrix.hpp`
     on user-material plugin ABI grounds (`cl_SourceFunction.hpp`).

## 11. Open questions

- Refactor `ElectricalCircuit` so ground is index 0 (SPICE-native) instead of last? Avoids
  a remap layer but touches all the `mNumberOfNodes-1` arithmetic in
  `cl_ElectricalCircuit.cpp`. *(Decided by default 2026-08-27: keep the remap layer for v1;
  revisit only if it bites in Phase 3.)*
- ~~Do we need `.subckt` flattening, or are real magnet decks flat? (Affects Phase 2 scope.)~~
  **RESOLVED 2026-08-27, Christian → no `.subckt` in v1**; hard-error on `X`/`.subckt`.
- Should `circuitrun` share `solver{}` parsing with the FEM path, or read `.options` from
  the netlist? (Leaning: small dedicated block, defaults matching `solve_circuit()`.)
- **Phase 6 split?** Fund `PULSE` (fixed parameter count; fits a widened schema) and gate
  `PWL` on a real deck that needs it. *Claude's recommendation from the DR-39 closing ruling,
  medium confidence — a recommendation, not a decision; it belongs to whoever funds Phase 6.*
- ~~How should the importer report unrepresentable SPICE features — hard error, or warn and
  continue? (Leaning: warn for source-parameter loss, hard error for unsupported elements.)~~
  **RESOLVED 2026-08-27 (audit round 1, adopted at reconciliation): v1 hard-errors on BOTH.**
  Codex's argument: dropping `ic=`/`.ic`/`uic` changes the transient initial state and dropping
  SIN `Vo`/`Td`/`Theta` changes every source sample — value-changing losses are not ignorable
  metadata, and "we ran your deck" after silently altering it is a lie. The showcase deck
  (`SIN(0 500 10)`, no ICs) needs no lossy path anyway. Relaxing to warn is a post-v1 decision.

---

## 12. Audit Round 1 (2026-08-27) — findings, corrections, and new O-items

**Round:** Codex + Grok, blind parallel, read-only, same brief; pre-registration and both full
audits in `tmp/ai_exchange/ngspice_parser_plan_audit.md` (distilled here before sweep).
Claude verified every load-bearing claim against the tree before adoption; verification results
are marked. Belief scoring: pre-registered claims 1 (API current) and 2 (no Pulse/PWL)
CONFIRMED by both; claim 3 (§7 loop) partly refuted (see the §7 correction banner); claim 4
(SPICE checklist) partly refuted (`.tran`, `GROUND`, groundless decks); claim 5 (Phase-0
internally consistent) REFUTED — the decision *list* is coherent but the plan body contradicted
it (now corrected in place); gap G-1 CONFIRMED by both, with the `output{voltages}` sibling.

### 12.1 The three P0s (unanimous) → new O-items

**O5 — node-reference binding across the 6b boundary (G-1).** Today the deck reads
`node +`/`node -` as raw ints checked against node count (`cl_ElectricalCircuitFactory.cpp:76-86`)
and `output{voltages}` as `Vector<id_t>` (`:632-636`). The netlist introduces a name→index remap
(ground → last). An integer deck reference is then ambiguous — SPICE name `"0"` is ground,
BELFEM index 0 is not — and Grok showed the showcase has the *inverted* convention today (live=0,
ground=1), so copying current numbers into a `.cir` silently flips polarity. **Adopted contract
(both auditors: "the only safe rule"; RESOLVED 2026-08-27, Christian — approved):** when `circuit{file:}` is
set, EVERY deck node reference — terminal pairs and `output{voltages}` — is a netlist node
**name** (string, case-folded like the parser map; `"0"`/`"gnd"` resolve through the map). No
post-remap integers, ever. Current outputs already match by label (`:623-628`); the companion
rule is **label = SPICE instance name** (`Is`, `L1`, …; Z1/Z2 stay deck-side), with one declared
case-fold policy — `save_timestep` matches labels with case-sensitive `==`
(`cl_ElectricalCircuit.cpp:915-926`), so the plan must fold at creation time.

**O6 — hybrid construction and ownership.** The existing factory owns and deletes its circuit
(`cl_ElectricalCircuitFactory.cpp:41-48`); `read_circuit` walks all component types with a
silent `default: break`. Adopted requirements (RESOLVED 2026-08-27, Christian — approved): ONE owner; with `file:` set,
`topology{}` may contain **only** `terminal pair` blocks — any lumped element left there is a
hard error, not a duplicate or a silent skip; `number of nodes` is forbidden/derived; the node
map is exported to the deck-side pair reader; all `create_*` calls complete before the first
`compute_MNA_matrix()` (adjacency is built once, `cl_ElectricalCircuit.cpp:538-546`);
**`hphiTrun` gets the same `file:` branch as `hphirun`** or thermal+circuit decks silently stay
on the old path.

**O7 — `.tran` precedence (RESOLVED 2026-08-27, Christian — approved).** Factually settled (see corrected §3 table): SPICE `tstep` is the
print increment, not Δt. Adopted rule: coupled v1 **ignores `.tran` with a log line**;
`solver{timestep{}}` owns time integration. The Phase-4 standalone runner defines its own
`.tran` policy (`tmax`→Δt bound, `tstop`→duration, initial Δt chosen independently).

**O8 — where does L/C BDF `order` live? RESOLVED 2026-08-27, Christian → option (a), the
`* belfem: order <instance> <n>` directive** (per-instance, next to its card, reuses the
Option-A machinery, `.cir` stays vendor-valid; the showcase carries
`* belfem: order L1 2` and `* belfem: order C1 2`). SPICE cards
cannot carry it; the factory default is 1; the showcase sets `order : 2` on both L and C, so
converting as planned would silently change the example's integrator. Note
`solver{timestep{scheme:}}` is NOT this order — circuit L/C carry their own `ode::BDF`
(`cl_Inductor.cpp:20-24`). Options: **(a) `* belfem: order <instance> <n>` directive**
(recommended — per-instance, sits next to the card, reuses the Option-A machinery, `.cir` stays
vendor-valid); (b) a deck-side `circuit{ orders{} }` block keyed by instance name; (c) waiver —
showcase becomes BDF1 (changes shipped-example behaviour; not recommended).

### 12.2 P1 adoptions (folded into the plan; tick when implemented)

- [ ] **O9 — deterministic layout, frozen for restart (RESOLVED 2026-08-27, Christian — approved; box ticks when implemented).** The `x` layout is nodes `0…N-2` then
  unknown currents in `create_voltage_source`/`create_switch` **creation order**
  (`cl_ElectricalCircuit.cpp:297-340`); the restart dump validates dof count ONLY
  (`:984-990`) — no labels, no topology, no L/C history. A netlist factory with unordered name
  packing yields a same-size `x` that `load_state` accepts and misapplies. Adopted: node
  packing = **first appearance in the netlist**, ground last; V/switch emission = netlist card
  order; print the name→index map at setup. The showcase (no V-source, no switch) is
  insensitive; general v1 is not. The missing-history/L-C-order half is pre-existing restart
  debt (see `todo/restart_circuit_verification.md`), not parser scope.
- [ ] **Error-reporting contract (Phase 2 requirement).** `ElementRecord` gains a source line
  number and the raw card text; every parser `BELFEM_ERROR` names file, logical line (after
  `+` joining), and offending token. Hard-error card list extended beyond §8's: `S`/`W`,
  `.control`/`.endc`, `.ac`/`.dc`/`.op`/`.nodeset`/`.options`/`.save`/`.print` (real ngspice
  decks are full of `.control` blocks; ignoring them and claiming "ran your deck" is a lie).
- [ ] **Units seam (Phase 1/3 requirement).** Netlist numbers NEVER pass through `unit_to_si`;
  `spice_number_to_si` is the only converter, and per-front-end validation stays separate from
  the shared `create_*` calls. **The `100F` footgun:** SPICE `F` = femto, so the showcase's
  `100 F` capacitor must be written `C1 <n+> <n-> 100` — `100F` would be 100 fF, a 1e-15
  error. Phase-1 tests must include `100F`→1e-13, longest-match `MEG` before `M`
  (case-folded first), `MIL`, atto-`A` (documented as ngspice extension), signs, `2.5e-6`,
  glued unit letters after a real suffix, and hard-error on old-style `2k5`.
- [ ] **Phase-2 test additions:** a deck whose first line is a real card (title-line
  consumption MUST eat it — dropping the first source is the classic failure), `+`
  continuation across a source function, `$`/`;` end-of-line comments (`//` accepted but
  ngspice-ish — don't detect it inside numeric tokens), `.end` handling documented as
  policy (EOF-tolerant is common SPICE practice; pick and document).
- [ ] **PULSE arity note (§4):** modern ngspice is `PULSE(V1 V2 TD TR TF PW PER [NP])` — the
  7-arg form in §4 is a supported subset, not exhaustive. Diode: area factor on the `D` card
  hard-errors in v1; document the default temperature behind `Vt = 0.026`.

### 12.3 Corrected anchor table (Grok, spot-checked)

| Plan citation (stale) | Today | Searchable token |
|---|---|---|
| `solve_circuit :950-1004` / `:1981-2039` | `:2104-2162` | `Controller::solve_circuit` |
| tEpsilon0 controller `:1999` | `:2122` | `tEpsilon0_Circuit = tEpsilon_Circuit` |
| tEpsilon0 exec `:159` | `:154` | `epsilon0 = epsilon` |
| exec Newton `:131-159` / `:127-157` | `:126-156` | `while ( tIt < tMaxIt` |
| `initialize_timestep :138-156` | circuit block `:237-271` (copy at `:369-401`) | `The circuit MNA solve runs on rank 0 only` |
| omega update exec `:146-155` == ctrl `:986-994` | exec `:145-154`; ctrl `:2140-2148` | `mBeta + mGamma` |
| MNA L/C stamps `:577-603` | `:584-610` | `get_discretized_resistance` |

### 12.4 By-catch (report-only here; filed separately)

1. **DR-115 root cause diagnosed statically (Grok; Claude tree-verified):** `mX` is sized only
   in `compute_MNA_matrix()` (`cl_ElectricalCircuit.cpp:544`), which runs on **rank 0 only**
   (`cl_FEM_Controller.cpp:244-249`), while `reset_timestep` calls `shift_back()` on **every
   rank** (`:2368-2370`) and `update_components()` indexes `mX(j)` (`:252-255`). Non-root ranks
   hold an empty `mX` → bounds abort on first rejection, "three of four ranks", serial clean.
   The register row's "history register depth" hunch was the wrong neighborhood
   (`ShiftRegister::revert` would raise "Cannot revert", not the `cl_BZ_Vector.hpp:349` bounds
   assert). Static diagnosis, high confidence; the fix still gets its own plan+audit round.
   **Consequence for this plan: DR-115 is a FIX gate, not a waive gate** — `examples/scripts`
   `Allrun` defaults to an even multi-rank count, so the showcase would abort under its own demo.
2. **`SourceFunction::set_periodic` time-offset precedence defect (Grok; Claude verified):**
   `set_time_offset((aPhase/2*constant::pi)*aPeriod)` (`cl_SourceFunction.cpp:175`) computes
   `(φ/2)·π·T`; the phase-fraction-of-period is `(φ/(2π))·T` — off by π² ≈ 9.87.
   **FIXED 2026-08-29 (DR-117 struck).** The "latent for sine" reading recorded here was
   WRONG: `set_time_offset` back-computes the phase from the offset, so the same line also
   clobbered `PHASE` to `φ·π²` — and `function_sine` reads `phase()`. Sine with nonzero
   phase was wrong too, including this parser's `SIN` phase argument. Fix was to DELETE the
   line: `set_phase` already writes both slots correctly.
3. **`electricalCircuit.cpp` reject path advances time on failed steps** and errors on a
   first-step rejection (§7 correction banner, items 1–2) — a defect in the shipped standalone
   executable, independent of this plan; matters when Phase 4 touches that file.

### 12.5 Open risks carried (no action yet)

- **Current-source stamp polarity vs SPICE convention** (Grok, medium): BELFEM stamps
  `mRHS(n+) += I; mRHS(n-) -= I` (`cl_ElectricalCircuit.cpp:735-746`); SPICE defines positive
  current as flowing from `n+` through the source to `n-`. Whether the sign conventions agree
  in BELFEM's residual form was NOT proven — settle it with the Phase-3 cross-simulation test
  (stock ngspice vs BELFEM on a trivial R-source deck) before converting the showcase, since
  O5's renumbering is exactly where a polarity slip would hide.
- **DR-75** (ordinal BC ↔ terminal-pair pairing, one bracket group) still applies under 6b;
  the hybrid reader must not reopen the multi-group grammar.
- MPI posture confirmed: parse the `.cir` on every rank exactly like `input.conf` today; no
  `share`/`receive` of the netlist; the circuit solve stays rank 0.

### 12.6 What may proceed

Phase 1 (`spice_number_to_si` + the §12.2 test set) and Phase 2 (lexer/IR + error contract) can
start once Christian approves implementation. Phase 3 must not bind to `hphirun`, and the
showcase must not be converted, until O5–O8 are signed off and DR-115 is fixed.
