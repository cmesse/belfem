# BELFEM Input File Reference {#doc_input_file_reference}
**Date:** 2026-08-10 (last revision; §4 solver keys, §4.2/§4.4 re-baselined)
**Purpose:** Complete reference of every section and key the `input.conf`
parser consumes — syntax, types, units, defaults, aliases, and known
pitfalls. This is the **human** half of the input contract; the machine-readable
half is `doc/input_schema.yaml`.

**Living document — both halves:** every change to the input contract (new key,
removed key, new default, new or removed alias, changed enum, changed
required-ness) must be reflected in **this file AND `doc/input_schema.yaml`** in
the same session that lands the code change. The rule is stated for both AIs in
`CLAUDE.md` §"The Input Contract: Two Artifacts, One Rule", because it is
triggered by editing the C++ — not by editing the docs.

**Maintenance rule:** the parser code is authoritative. When this document and
the code disagree, fix whichever is wrong and note the correction here.

> **`file:line` citations in this document are navigation aids only.** A citation
> can be wrong the day it is written — a line number read from a working copy is
> already stale by the time that copy is committed — and an unrelated edit
> elsewhere in a file shifts every citation into it at once. Some rows also point
> at a consumer rather than a parse call, deliberately. **When a citation and the
> code disagree, believe the code.**
> `doc/input_schema.yaml` anchors each key by a searchable token instead
> (`anchor: '"tolerance"'`), so it does not carry line numbers.

---

## 1. File syntax

Parsed by `src/io/cl_InputFile.cpp` / `cl_Input_Section.cpp`.

```
// comment ( # also truncates the rest of a line; no /* */ )
section [ : label ]
{
    key : value ;          // normal key
    flag ;                 // flag form, stored as the string "true"
    nested { ... }
}
```

- Comments: `//` to end of line (`cl_InputFile.cpp:67`); a `#` also
  truncates the rest of a line via `clean_string`
  (`stringtools.cpp:86-92,160-162`).
- Multiple `key : value ;` statements may share one physical line; blank
  lines are dropped.
- **Section** header is the line before `{`; `type` (before the first `:`)
  and `label` (after it) are both lowercased. Two sections of the same type
  are distinguished by label (`thinshell : tape1 { }`).
- **Sections are stored twice** (`cl_Input_Section.cpp:76-89`): every child
  section is pushed in file order onto an ordered list that preserves
  duplicates, and also inserted into a by-name map keyed by `type`, or
  `type:label` when labeled. A later section with the same map key overwrites
  the earlier map entry. Current consumers of repeatable section sets — the
  Maxwell and thermal BC factories, topology domains, materials, top-level
  `layers`, and circuit components — use the ordered list, so repeated
  same-type sections under those consumers are consumed separately. By-name
  lookup (`section("current")`) sees only the last unlabeled duplicate; a new
  consumer of a repeatable section type must iterate by index.
- **Key names are lowercased** and may contain spaces (`max iterations`).
  **Values are not lowercased** — whether a value is case-sensitive depends
  on the consumer; case-sensitive values are marked below.
- Duplicate keys overwrite each other — except in `layers` blocks, which
  are read line-by-line and allow repeats (see §6).
- Booleans: `true`, `on`, `yes`, `1` (lowercased) are true; anything else
  is false (`stringtools.cpp:300-310`).

### Lists, ranges, groups

- Id lists: comma-separated, brackets optional as separators; `a:b` is an
  **inclusive range in either direction** (`4:9` = `9:4` = 4…9)
  (`cl_Input_Section.cpp:557-628`). Lists are deliberately NOT uniquified
  or sorted (sorting broke periodic node pairing).
- Bracket **groups** (terminal lists): `1,3,5` = three separate groups;
  `[1,2]` = one group of two; `[1],[2]` = two groups
  (`cl_Input_Section.cpp:632-712`). The brackets carry the electrical intent:
  ids inside one `[ ]` are electrically connected and act as a single
  terminal; each separate group is its own independent connection.
  **A range outside a bracket expands the same way**, one group per member:
  `42:48` is seven groups, i.e. seven boundary conditions, not one terminal
  of seven sidesets. Write `[42:48]` for the single terminal. The expansion
  flushes each id as its own group whenever the bracket flag is down, which
  makes this indistinguishable from a hand-written comma list. Since each
  group becomes one current/voltage condition, and each condition consumes
  one cohomology generator, the missing brackets surface much later as
  "more conditions than cohomology generators" from
  `Cohomology::updatekGeneratorsFromHomology`, not as a parse error.
- Intersection form (3-D curves): `a @ b` (`cl_Input_Section.cpp:504-515`).

### Units

`unit_to_si` (`src/core/stringtools.cpp:334-1122`) converts the second word
of a value to SI. Unit checks compare the **dimension**, not the token — any
unit of the right dimension is accepted (`kA` where `A` is expected, `mum`
where `m` is expected). Unknown token → fatal `"Unknown unit : %s"`.

Supported families (each with the usual SI prefixes): length (`nm`…`km`,
`in`, `ft`, `mi`), mass (`mg`…`t`, `lb`, `slug`), time (`ns`…`h`),
frequency (`Hz`…`MHz`), force (`nN`…`MN`, `lbf`), energy (`mJ`…`MJ`,
`cal`), power (`nW`…`GW`), pressure (`muPa`…`GPa`, `bar`, `atm`, `psi`),
current (`mA`, `A`, `kA`), voltage (`muV`…`kV`), flux density (`G`,
`muT`…`MT`, `Oe` family), temperature (`K`, `C`/`°C` add 273.15, `°F`,
`R`/`°R`), amount (`mol`, `kmol`), angle (`rad`, `°`, `deg`), resistance
(`nΩ`…`kΩ`, both `Ohm` and `Ω` accepted, `stringtools.cpp:995`),
capacitance (`nF`…`F`), inductance (`nH`…`H`).

A one-word numeric value is taken as already-SI and dimensionless.

---

## 2. Top-level layout

| Section | Required | Consumer |
|---|---|---|
| `mesh` | yes | `cl_MaxwellFactory.cpp:196,249,292` |
| `solver` | yes | `cl_MaxwellFactory.cpp:560,734`, `cl_FEM_Controller.cpp:2240` |
| `materials` | yes | `cl_MaterialFactory.cpp:36` |
| `layers : <tape>` | one per `thinshell` | `cl_MaxwellFactory.cpp:2463-2534` |
| `homology` | when cohomologies are computed (fresh gmsh mesh or stale `.bfm`) | `cl_MaxwellFactory.cpp:922-935` |
| `topology` | yes | `cl_MaxwellFactory.cpp:326-355, 2433` |
| `boundary conditions` | yes | `cl_MaxwellFactory.cpp:139-151` |
| `initial conditions` | optional | `cl_MaxwellFactory.cpp:154-170` |
| `circuit` | optional | `cl_ElectricalCircuitFactory.cpp:28` |

There is no `output` or `visualization` section; the result file names come
from the executable, not the deck. `belfem` names the Exodus file after the
mesh (`belfem.cpp:210`), except on the segregated coupled path, which writes a
fixed `hphi_results.e-s` (`belfem.cpp:303`). The only `output` subsection
lives under `circuit`. The other outputs have fixed names:
`memdump.hdf5`, written next to the results and read back at startup; see
`timestep → restart` (§4.4).

---

## 3. `mesh`

| Key | Type | Notes | Site |
|---|---|---|---|
| `file` | path | required. `.bfm` loads directly (cohomology skipped); a gmsh file looks for a sibling `.bfm` and reuses it when **both** the base-mesh checksum and the mesh-configuration tag match — see below | `cl_MaxwellFactory.cpp:196` |
| `unit` | length unit | required for non-`.bfm` meshes; scales the mesh to SI. A unit change moves every node coordinate, so it already invalidates the cache through the checksum | `cl_MaxwellFactory.cpp:292-303` |

### Which settings invalidate a cached `.bfm`

A `.bfm` stores the **enriched** mesh: cuts, thin-shell layer blocks, edge-coating walls,
periodicity pairing, hanging entities. Two stamps guard its reuse, and both must match:

- the **checksum**, which is the identity of the base `.msh` (node coordinates and element
  connectivity), and
- the **mesh-configuration tag**, a fingerprint of the settings that decide how that mesh is
  enriched.

Editing any tagged setting rebuilds the mesh, and the run reports which setting changed. The
tag includes all `topology` keys (domain types, thin-shell sidesets, `edge coating` and
width, `periodic` source/target, curves), the `layers` stacks (materials, SI thicknesses
and layer count), `homology { algorithm }`, and current-injection terminal ids (`input
terminals` / `output terminals`, or `input curves` / `output curves`) from both the
`boundary conditions` and `circuit` trees. A deck can declare those ids in either tree;
`examples/tapestack_circuit` declares them only under `circuit`, in a `terminal pair`.

The machine-readable half of this list is the `mesh_config_tag:` field in
`doc/input_schema.yaml`, which records how each section reaches the tag. For example,
`topology` and `layers` are walked `wholesale`, so a new key under them is covered
automatically; `homology` contributes one `included` key, so a new mesh-affecting sibling
would need to be added to the tag builder too.

**Not** in the tag, so you can change them freely and still hit the cache: everything under
`solver`, timestepping and tolerances, output settings, and the *amplitude*, period and
waveform of a boundary condition.

Two consequences worth knowing:

- Thicknesses are compared in SI, so `50 mum` and `0.05 mm` are the same configuration and do
  not trigger a rebuild.
- A change *inside* a user-defined material definition is deliberately not tracked. Layer
  material **labels** are, so swapping `copper` for `silver` rebuilds; redefining a user
  material so that it gains or loses a resistivity does not.

---

## 4. `solver`

Keys read directly in the `solver` section (outside the subsections). Placement matters: the
same key inside `linear magnetic { }` / `linear thermal { }` is read by nobody and silently
ignored (§13). The two deliberate exceptions are `compute conditioning` and
`mumps error analysis`: both are read in either location and documented with the linear
sections in §4.1:

| Key | Type | Default | Site |
|---|---|---|---|

### 4.1 `linear` / `linear magnetic` / `linear thermal`

Magnetic solves prefer `linear magnetic`, thermal solves prefer
`linear thermal`; both fall back to `linear` (`cl_MaxwellFactory.cpp:564-567`,
`cl_ThermalFactory.cpp:165-168`). All parsed by
`SolverParameters` (`src/sparse/cl_SolverParameters.cpp`).

For the `belfem` executable, the thermal solver sections also select the
mode. An unlabeled `linear thermal` or `nonlinear thermal` section requests
the coupled h-ɸ/T problem. Without either section, `belfem` solves the
magnetic problem only. It prints the selected mode at startup. Labeled
sections (`linear thermal : name { }`) do not count. A deck that declares
`boundary conditions → thermal` without one of these sections is refused
at load time as inconsistent.

| Key | Values | Default | Site |
|---|---|---|---|
| `library` | `umfpack`, `superlu`, `mumps`, `strumpack`, `pardiso`, `petsc` (case-insensitive) | best available: STRUMPACK → MUMPS → PARDISO → UMFPACK → SUPERLU | `:80-87` |
| `matrix format` | `csr`, `csc`, `aij` | `csr`; forced `aij` for parallel PETSc unless explicit | `:43-63` |
| `krylov method` | `preonly`, `bcgs`, `cg`, `cgs`, `ibcgs`, `gmres`, `tfqmr`, `auto` | `auto`. Consumed by PETSc and STRUMPACK; MUMPS ignores it (fixed internal refinement). **On STRUMPACK, `auto` resolves to factorization-preconditioned GMRES since 2026-08-16** (uniform semantics with PETSc; strictly more robust than the old Richardson refinement on degraded factors — the refinement floor trace), and **the unmapped methods `cg`, `cgs`, `ibcgs`, `tfqmr` are a hard error at setup since 2026-08-18** — they used to fall through silently to the library's own AUTO, i.e. Richardson-REFINE, a solver nobody asked for. `preonly` = no outer iteration at all — an expert setting for measured cases only: on an ill-conditioned matrix the raw factor delivers ~4 digits and the outer GMRES is what reaches deep residuals. | `:47-50` |
| `preconditioner` | `none`, `asm`, `gamg`, `jacobi`, `bjacobi`, `lu`, `ilu`, `icc`, `hmg`, `spai` | `asm` (PETSc only) — was `gamg` parallel / `jacobi` serial until 2026-08-12 | `:51-54` |
| `reordering scheme` | `natural`, `metis`, `scotch`, `automatic`, `parmetis`, `ptscotch` | `automatic`. **`parmetis` / `ptscotch` (since 2026-09-04)** select the parallel library explicitly. Read the overload carefully: for MUMPS and STRUMPACK, `metis` / `scotch` *already* mean ParMETIS / PT-Scotch at more than one rank, and the new values behave identically there. What they add is the **BELFEM-side parallel nested dissection for PETSc**, where `metis` orders serially on rank 0 before the rows are distributed. A value whose library is not linked (`USE_METIS` supplies ParMETIS, `USE_SCOTCH` PT-Scotch) is a hard error at solver construction, on every rank. If the row count is below the rank count the ParMETIS path announces a fallback to serial METIS on rank 0 and continues. | `:55-58` |
| `metis nodendp` | bool (STRUMPACK only, METIS reordering only) | `true` — selects METIS_NodeNDP, the partitioning entry point STRUMPACK itself declares (undocumented in METIS), which flattens the supernodal etree and avoids the deep-tree stack-overflow risk STRUMPACK warns about. Set `false` to restore plain METIS_NodeND if a METIS build lacks NodeNDP. Ignored by PARMETIS and SCOTCH. | `:75-78` |
| `compression scheme` | `off`, `blr`, `automatic` | `automatic` — **off on every library** (unified 2026-08-16): only an explicit `blr` compresses. Before the unification, MUMPS treated everything except `off` as BLR-on with a silent 1e-8 absolute drop tolerance — six shipped example decks were affected, garber with *negative* headroom against its 1e-10 nonlinear target. BLR is an accuracy-for-memory trade; see `src/sparse/doc/solver_memory_and_compression.md` §4 for when it is sensible. | `:59-62` |
| `compression cutoff` | real, > 0 and finite | `1e-8` — the BLR drop tolerance, consumed **only** with `compression scheme : blr`. One number, two strictnesses: STRUMPACK applies it as a *relative* tolerance, MUMPS as the *absolute* `CNTL(7)` (MUMPS 5.7.3 §5.19). Pair it 2–3 decades below the nonlinear `tolerance`. **Zero-or-negative headroom (cutoff ≥ the nonlinear tolerance, stated or inherited) is a hard error at setup** — the same impossibility the PETSc gate refuses; less than two decades but positive prints a warning (unguarded on rank 0 — even `-v 0` shows it). Zero (MUMPS's lossless-BLR niche) is deliberately unreachable from the deck. | `"compression cutoff"` in `cl_SolverParameters.cpp` |
| `memory budget` | integer, MB per process, > 0 (a bare whole number — no unit token; `1.6` is refused, not rounded) | not stated. **MUMPS only.** A per-process working-memory cap handed to MUMPS as `ICNTL(23)`. **Stated, it applies from the first factorization.** Left unstated, the wrapper measures the machine once when the solver instance is created — available memory (Linux `MemAvailable`, further capped by the tightest cgroup limit above the process when the cgroup filesystem is at its standard mount `/sys/fs/cgroup`; Darwin free + inactive pages), divided by the ranks on the node and by the live MUMPS instances, halved for safety, and reduced to the minimum over all ranks — and applies that measurement **only after** a first out-of-workspace failure (`INFOG(1) = -9` / `-8`): if the measured budget is at least MUMPS's own estimate for the factorization (`INFOG(16)`, `INFOG(36)` under BLR) it becomes `ICNTL(23)` and the factorization is repeated; if it is **below** the estimate the step is handed to the timestep controller at once, with no `ICNTL(14)` climb — the machine cannot hold the estimate and a larger relaxation would only ask for more. Behind a cap that is in place, the `ICNTL(14)` ladder (30 → 60 → 120 → 240 → 480 %) still runs on a recurring `-9`/`-8`, because MUMPS documents that they can recur with the cap set, and on `-17`/`-20` (MPI send / reception buffer too small — sized from `ICNTL(14)` alone, so the ladder is their only remedy, cap or no cap); `-19` (the cap itself cannot be met) hands the step to the controller. Unstated, the *first* factorization runs uncapped and can still be killed by the kernel if the analysis estimate exceeds RAM — state a value on a shared or small-memory node. If the machine could not be measured on any rank (unknown platform, or a cgroup hierarchy declared in `/proc/self/cgroup` but not found under `/sys/fs/cgroup`), the ladder alone runs, as before. | `"memory budget"` in `cl_SolverParameters.cpp` |
| `relative tolerance` | real | `1e-10` — applied **always**: PETSc KSP rtol, and STRUMPACK's outer-GMRES `rel_tol` (since 2026-08-18). **A normal deck does not need this key** — the class default is the production value; state it only as an expert override. History a user does not need but an expert does: until 2026-08-18 STRUMPACK honored the key only when the deck stated it and otherwise kept the library default `1e-6`. That split was a REFINE-era stall fix (2026-08-15: an unmeetable relative target on a numerically zero RHS, uncapped, hung `IterativeRefinementMPI` for ~1 h); the stall protection is now `maxit = 50` plus the always-set `abs_tol = 1e-14`, and `krylov method : auto` has mapped to factorization-preconditioned GMRES since 2026-08-16. The 2026-08-18 tapestack3d A/B (magnetic `1e-8` → `1e-10`, same restart state) exposed the opposite defect: the loose exit test stopped GMRES on a ~4e-9 plateau and that number **was** the printed nonlinear residual (−85 dB); at `1e-10` the same solve finishes its dive (~1e-15), Picard reaches its 1e-11 target in two iterates, and the step-collapse disappears — the cost is a handful of extra Krylov iterations against a ~25 s factorization. The outer loop stays capped at 50 iterations — a cost bound, not a named failure: at the cap STRUMPACK returns the best-effort solution as success (source-verified 2026-08-16), and the nonlinear loop judges the true residual. **Headroom gate (since 2026-08-15): when the field's linear solver is PETSc, a linear relative tolerance looser than that field's nonlinear `tolerance` is refused with a hard error at setup.** STRUMPACK is exempt on purpose: its working default pairing is linear `1e-10` against a 1e-11 nonlinear target, and that relies on GMRES overshooting the exit test — a deck whose linear residual crawls rather than dives should state a tighter number here, not expect the gate to save it. The lossy case — `compression scheme : blr` on MUMPS or STRUMPACK — has its own two-tier gate since 2026-08-16: zero-or-negative headroom errors, thin-but-positive headroom warns (see the `compression cutoff` row). | `"relative tolerance"` in `cl_SolverParameters.cpp` |
| `initial guess` | bool | `false` — PETSc only. Caveat after a soft-failed solve (since 2026-08-15): the diverged iterate is never written back, so a retry under `true` starts from the *previous* increment — a field estimate for Picard, but a correction for a different Jacobian under Newton. The default retries from zero. | `:67-70` |
| `absolute tolerance` | real > 0, finite | `1e-14` for STRUMPACK (always applied); unset for PETSc (library default `1e-50` kept unless the deck states one). Exit test of the ITERATIVE part of a linear solve, measured on the ABSOLUTE residual. Exists because (measured 2026-08-17) STRUMPACK's library default of `1e-10` was never overridden, and on a small right-hand side — a transient's startup, where ‖b‖ ~ 1e-2 — that floor sits exactly at a `1e-11` RELATIVE nonlinear target, so the outer GMRES stopped a factor ~300 short of what Newton needed and step convergence degenerated into an overshoot lottery (measured: stall entry at floor/‖b‖ = −85 dB, iterate tails of 10–26, Δt collapse). The class default `1e-14` restores ≥3 decades of headroom; the `max iterations` cap bounds what it may cost. Loosening this key back toward `1e-10` re-creates that stall in small-RHS regimes — measure before touching it. | `"absolute tolerance"` in `cl_SolverParameters.cpp` |
| `max iterations` | int > 0 | unset — the library default (PETSc: 10000). **Linear** Krylov iteration budget, consumed by PETSc only; distinct from the *nonlinear* key of the same name in §4.2, which the controller parses — this one caps the inner Krylov loop, not Newton/Picard. Negative or zero values are rejected at setup. **Set it from a measurement, not a guess:** a solve that legitimately needs thousands of Krylov iterations (a tight `relative tolerance` on a stiff system will) turns a low cap into a timestep collapse in which the physics is never at fault — measured 2026-08-15, a cap of 500 on a thermal block at rtol 1e-10 failed every solve down to Δt = 12.5 ms. Info level `Verbose` prints the actual per-solve count. A command-line `-ksp_max_it` still overrides it (`KSPSetFromOptions` runs after). Since 2026-08-15 the iteration-class divergences — `DIVERGED_ITS` (budget exhausted), `DTOL`, `NULL`, and the two `BREAKDOWN`s — are soft failures under the controller (timestep cut, retry); `NANORINF`, `PC_FAILED`, and every unlisted reason stay hard aborts, because KSP reuse after them is not established. | `"max iterations"` in `cl_SolverParameters.cpp` |
| `matching` | bool (STRUMPACK MC64) | `true` — load-bearing for h-φ solves, disable only as a verified opt-out | `:71-74` |
| `compute conditioning` | bool | `false` — prints that field's eigenvalue conditioning estimate as a `\|λ\|max/\|λ\|min` row in each timestep footer. **Read per field here, and also in the `solver` section itself, where it is the shorthand for both**; the per-field key wins (`"compute conditioning"` in `cl_FEM_Controller.cpp`; note the key is read from the LINEAR blocks, not the nonlinear ones). ARPACK/PARPACK always computes the estimate, adding one eigen solve per timestep *per field*. Its combined time appears on the footer's Eigen-value-analysis line. The estimate is sampled **once per timestep, at the first iterate**. That iterate is assembled at the converged previous step, so the series can be read as a trend. Once the ratio exceeds about `tolerance / ε_mach`, it cannot resolve the small end. A mixed h-φ system reaches ~1e17; there it reports `n/a`, with a one-line footer notice. The parser warns when either field opts in without MUMPS (the recommended solver for runs that also want the ADD rows below). The ratio equals κ₂ only for a normal matrix, which the h-φ Jacobian is not — hence the honest label. Use it to diagnose a Δt-independent residual floor. If ratio·ε_mach is near the observed floor, conditioning is the likely cause; if it is far below, look to the model or formulation. Nothing in the iteration scheme consumes it. **This key no longer controls the `MUMPS ADD` rows** — see `mumps error analysis` below (split 2026-08-30). See `src/fem/doc/timestepping_strategy.md` §4. | `"compute conditioning"` in `cl_FEM_Controller.cpp` |
| `mumps error analysis` | bool | `false` — prints the **`MUMPS ADD COND1`** and **`MUMPS ADD COND2`** footer rows. They are the Arioli–Demmel–Duff pair from MUMPS's own error analysis (`ICNTL(11)`): a componentwise 1-norm condition estimate for the solved system *with its actual right-hand side*. This is a different quantity from the spectral ratio above, and the two must not be compared. It follows the same placement rules as `compute conditioning`: the `solver` section is shorthand for both fields, while a per-field key in a linear block takes precedence. It is inexpensive: the controller arms the analysis for the **first solve** of each timestep, then disarms it immediately. It uses omega statistics plus a Hager estimator against the existing factorization; it needs no refactorization or extra solve. COND2 is omitted (not `n/a`) when its term cannot contribute — usually because no matrix row fell into its category, and also in the rare case where rows did but omega2 is exactly zero. A **no-op on any library but MUMPS** — the parser warns, which is why the key names the vendor. | `"mumps error analysis"` in `cl_FEM_Controller.cpp` |

### 4.2 `nonlinear` (alias `nonlinear magnetic`, which wins if present)

Section required; `nonlinear magnetic` is preferred over `nonlinear` when both
exist (`cl_FEM_Controller.cpp:3703-3705`).

| Key | Type | Default | Site |
|---|---|---|---|
| `tolerance` (alias `relative tolerance`, consulted only when `tolerance` absent) | real | `1e-6` | `:2911-2918` |
| `absolute tolerance` | real | `0.0` (disabled) — consumed since 2026-08-09: the magnetic loop also terminates once the absolute residual ‖Ax−b‖ drops below this value, matching the thermal loop. Use it per deck when the relative residual sits at a solver noise floor above `tolerance`; otherwise the floor policy in §4.4 cuts, escalates, and can stop at `floor retries`. Dimensional; no universal default. Negative values are rejected at setup. | `:2920-2928` |
| `tolerance switch` | real (Picard→Newton promotion) | `1e-4` | `:2930-2933` |
| `algorithm` | `Newton`, `Newton-Raphson`, `Picard` — **case-sensitive**; guidance §4.5 | `Picard` | `:2935-2958` |
| `max iterations` / `min iterations` | int | `100` / `2`. Since 2026-08-29 (DR-127 certified exit) `min iterations` counts BODY solves: every timestep first measures the committed state's residual from a fresh assembly (the head), takes a solve only while uncertified or under the minimum, and exits on a head-certified state. `min iterations : 0` therefore permits a zero-solve exit when the predictor already meets tolerance — an explicit deck choice. The exit residual printed is always the committed state's. At the timestep floor, §4.4 can temporarily raise `max iterations` up to 4× the deck value. Since 2026-08-28 both keys are validated at setup: `max` must be positive, `min` must not be negative, and `max` must be at least `min` (checked against the other key's default when only one is set). Negative values used to wrap silently into a huge budget. | `:3756-3771` |
| `target iterations` | int (drives Δt adaptation) | `20` — must be positive (rejected at setup since 2026-08-09) | `:2970-2977` |
| `max relaxation` / `min relaxation` | real | `1.0` / `0.001` | `:2979-2987` |
| `stall window` | int, clamped ≥ 2 | `5`. Since 2026-08-28 negative values are rejected at setup — they used to slip past the clamp by wrapping to a huge unsigned value first. `0` and `1` still clamp to `2` silently. | `:2989-2995` |
| `stall tolerance` | real (dB) | `0.2` (raised from `0.001` on 2026-08-09 so the stall band clears the ±0.05–0.2 dB wander of a floored residual; otherwise the Newton→Picard demotion never fires. Genuine convergence moves ≥ O(1) dB per iterate) | `:2997-3000` |
| `anderson depth` | int 0…8 (else fatal) | `0` (off). If this key is absent, `timestep { anderson stabilization : true ; }` fills `3`. | `:3002-3009` |
| `watchdog window` | int, `0` disables | `30`. Since 2026-08-28 negative values are rejected at setup. At the timestep floor, §4.4 can temporarily raise it up to 4× the deck value, like `max iterations`; `0` stays disabled. **Since 2026-08-21 the fire rule has a spare clause:** the controller does not cut the timestep after a no-new-best window if that field's relaxation has grown since the previous watchdog evaluation. This growth indicates that the line search is still recovering from an overshoot; measurements on the coarse tapestack3d mesh showed that the bare rule cut recoveries mid-climb and caused cascading Δt reductions. Omega-pinned stalls fire exactly as before. See `src/fem/kernel/doc/nonlinear_controller_theory.md`. | `"watchdog window"` in `cl_FEM_Controller.cpp` |

#### 4.2.1 Penalty sub-blocks (optional, since 2026-08-24)

Two optional sub-blocks feed the Maxwell IWG's penalty slots. Since 2026-09-01
both are **opt-in**: an absent block means **off**. With no gauge block the
gauge term is not assembled (`chi = 0`); with no ghost block the thin-shell
layers share their interface edges and faces — no duplicate dofs and no ghost
facets exist at all. Between 2026-08-27 and 2026-09-01 both were on by
absence (gauge at `chi = 1e-4`, ghost at `4` / `1e-3 Ohm`); before 2026-08-27
the gauge was off and the ghost on. The reason for the change is in the `chi`
and `eta` rows below. Turning either on is an explicit deck decision.

```text
nonlinear magnetic
{
    ...
    coulomb gauge penalty        // absent = gauging OFF ( opt-in since 2026-09-01 )
    {
        chi : 1e-2 ;             // dimensionless; state a positive value to gauge
    }
    nitsche ghost penalty        // absent = ghost OFF ( shared interface edges )
    {
        eta   : 4 ;              // dimensionless, REQUIRED in the block; 0 = off
        k_reg : 1e-3 Ohm ;       // unit REQUIRED; ignored when eta is 0
    }
}
```

| Key | Type | Default | Site |
|---|---|---|---|
| `chi` | real ≥ 0, dimensionless | block absent → `0` (**off**, opt-in since 2026-09-01). It was on at `1e-4` from 2026-08-27 until then: on the tape decks that value was invisible to the conditioning estimate — κ identical to two digits with and without it — while κ tracked the timestep (κ ≈ 2e15 × Δt[ms]); `chi : 0.01` does register (a factor 2–3). Gauging regularizes the curl–curl operator once the conductor turns resistive, and the Newton tangent has carried the gauge term consistently since 2026-08-27; state a positive `chi` to use it. An empty block is fatal. | `"coulomb gauge penalty"` in `cl_FEM_Controller.cpp` |
| `eta` | real ≥ 0, dimensionless, **required when the block is present** | block absent → ghost **off** (opt-in since 2026-09-01). `eta : 0` is also off; `eta > 0` switches the ghost on with that multiplier (4 was the old default). **Off means the thin-shell layers share their interface edges and faces — no duplicate dofs, no ghost facets are created at all** (the pre-2026-03-18 model; an interface resistance is then a thin resistive layer such as `rint`, not a Nitsche term). One reader serves the thin-shell factory, the controller and the mesh cache tag (`fn_FEM_ghost_switch.hpp`), so the three agree; a flipped switch misses the `.bfm` cache and rebuilds; an explicitly named `.bfm` built the other way is refused (a file from before 2026-09-01 carries no switch line and is treated as ghost ON, which is what the old factory always built); a memdump written since 2026-09-01 of the other layout is refused by its edge / face counts, and any dump whose edge field has a different length is refused by the field-length check. | `"nitsche ghost penalty"` in `fn_FEM_ghost_switch.hpp` |
| `k_reg` | real > 0, **unit mandatory** (Ohm dimension) | `1e-3` Ohm — read only when `eta > 0`; stated next to `eta : 0` it is ignored with a one-line notice | `"k_reg"` in `cl_FEM_Controller.cpp` |

Unit handling, all enforced at setup:

- `k_reg : 1e-3 ;` without a unit is **fatal** ("Invalid unit … expect Ohm").
  `k_reg : 1e-3 mOhm ;` converts to `1e-6` Ohm; `nOhm`/`muOhm`/`kOhm` and `Ω`
  are accepted too. Unit tokens are **case-sensitive** — `ohm` is an
  unknown unit and fails while the deck is read.
- `chi` and `eta` are dimensionless and **reject** any dimensioned value:
  `eta : 4 mOhm ;` is fatal, not silently parsed as `0.004`. The check
  is by SI dimension, so the scale-only angle tokens slip through: `eta : 4 deg ;`
  silently becomes `4·π/180`. Do not write angle units on these keys.
- `k_reg = 0` is rejected: with the unfloored power law (`mRhoMin = 0`) both
  layer stiffnesses of a superconductor–superconductor ghost facet are exactly
  zero, and the regularized harmonic mean would evaluate `0/0`.

Pitfalls and semantics:

- The sub-blocks are read from the **winning** nonlinear section. If both
  `nonlinear` and `nonlinear magnetic` exist, sub-blocks under `nonlinear`
  are silently ignored — same alias rule as every key in this section.
- The blocks must be written without a label. A labeled block
  (`coulomb gauge penalty : foo { }`) is not found, so it is treated as absent.
- Inside the ghost block, `eta` is required and `k_reg` is optional: a block with
  only `k_reg` is fatal (it does not say whether the ghost is wanted), and `k_reg`
  keeps its default `1e-3` Ohm when omitted.
- `0 < chi ≤ ~2.2e-15` (`BELFEM_EPSILON`) parses as "on" but the assembly
  kernels still skip the gauge term — values in that range behave as off.
- A nonzero `chi` (above that threshold) changes the assembled operator, and
  with it the meaning of every residual and tolerance number in this section.
- The penalties are deck state, not solver state: a warm restart takes them
  from the **new** deck; the memdump does not carry them.
- Effectiveness note (gauge campaign, 2026-08): on tetrahedral and thin-shell
  meshes the gauge term does not affect conditioning — the modes it would
  need to lift are invisible to it. It is retained for hex-dominant meshes
  and experiments.
- At setup the controller logs the three resolved slots once (at default
  verbosity or higher; `-v 0`/`-v 1` hide the line):
  `penalty slots : ghost eta …, ghost k_reg … Ohm, gauge chi …`.

### 4.3 `nonlinear thermal` (optional)

Same shape as §4.2; only differences are listed
(`cl_FEM_Controller.cpp:2675-2790`):

| Key | Default | Notes |
|---|---|---|
| `absolute tolerance` | `0.0` (disabled) | Terminates the thermal loop on ‖Ax−b‖; same semantics as the magnetic key (§4.2) since 2026-08-09. Negative values are rejected at setup. |
| `tolerance switch` | `1e-3` | |
| `max iterations` | `100` | The parse and default match §4.2, but the counting depends on the coupling mode. Enforcement in fully-coupled mode began on 2026-08-27. Before that date, the key was silently inert in this mode, and only the magnetic budget, divergence rules, and watchdogs bounded the coupled loop. In fully-coupled mode, the thermal solver runs in lockstep with the magnetic iterates unless the experimental `update gate` freezes it. Frozen iterates do not count. The counter therefore records completed thermal solves, not iterations of a nested thermal loop. A spent budget cuts the timestep only if the magnetic field has reached its relative or absolute tolerance at least once in the attempt, using a per-attempt latch since 2026-08-28. The cut also requires the thermal field to have met neither of its own tolerances. Before the latch, the magnetic field had to be at tolerance on the same iterate. Continuing thermal updates could therefore perturb an already-converged magnetic field, disarm the budget, and leave the loop running until the magnetic ceiling instead (measured: 23–44 min per hopeless step). A converged thermal field waiting for the magnetic field never trips the budget, and the thermal budget never cuts while a magnetic field that has not yet reached tolerance is still converging. Since 2026-08-28 the value must be positive and at least `min iterations`; negative values used to wrap silently and disable the ceiling. A flattening thermal residual defers the cut until the flat streak resolves, allowing at most four additional solves before the flat-stall exit latches. A latched flat-stall is accepted with a warning instead of being cut. In segregated mode, the key directly caps the inner thermal sub-loop, as it always has. At the timestep floor, §4.4 can temporarily raise the budget to as much as 4× the deck value, as it can for the magnetic key. |
| `update gate` | disabled | experimental: skip thermal update while magnetic residual above gate |
| `min relaxation` | `0.1` | |
| `anderson depth` | `0` (off) | If this key is absent, `timestep { anderson stabilization : true ; }` fills `1` (≈ Aitken). |
| `watchdog window` | `30` | Same semantics as §4.2, with one coupled-mode difference since 2026-08-28: the effective window is clamped to `magnetic max iterations − thermal min iterations − 1` (floored at 1), with a setup notice when the clamp engages. With lockstep counting, the last iterate that can run the thermal watchdog is the magnetic ceiling itself, so a wider window could never fire first — measured on the DR-97 record, where a thermal best at iterate 3 put the trip at 33 against a magnetic ceiling of 30. The clamp is sized to an early thermal best; a late best or the spare clauses can still let the ceiling win. `0` stays the disable, and segregated decks are never clamped. |
| `coupling` | `fully coupled` | or `segregated` (case-sensitive); `coupling factor` (int, default 1) read only when segregated. Since 2026-08-28 the factor must be positive: it divides `initial timestep` to seed the thermal substep, so `0` caused division by zero and a negative value reversed the thermal timestep's sign. |

`target iterations`, `stall window`, `stall tolerance` are magnetic-only.

### 4.4 `timestep` (required)

| Key | Type / unit | Default | Site |
|---|---|---|---|
| `initial timestep` | value, s | required (also seeds the thermal substep Δt/coupling factor). Since 2026-08-09 setup validation requires it to be present, positive, and inside the `minimum`/`maximum timestep` window; previously, an absent key silently propagated NaN time. **Since 2026-08-19 it has a second role: it caps the first step after a warm restart**. `load_memdump` keeps a dumped `delta_time` only up to this cap. On tapestack3d, a 50 ms first restart step produced `GMRES it. 0 = 1.1e4` and a rejected step, while a 5 ms first step gave `0.903` in 16 Krylov iterations and converged in two Picard iterates. The run then grows back under the existing BDF growth limiter, about 13 steps on that deck. If restarts feel slow, raise this knob cautiously — **do not jump straight to the working step**: the cap exists because a large first Jacobian at a developed state is exactly what failed, the cold-start argument does not transfer (a cold start sits at t = 0 with negligible current), and no ceiling is enforced, so an over-raised value re-creates the failure the cap mitigates. Prefer moderate raises validated against your own deck's restart behavior. See `src/fem/doc/timestepping_strategy.md` §5. | `:3136-3138`, validation `:3158-3165` |
| `maximum timestep` | value, s | unbounded | `:3142-3146` |
| `minimum timestep` | value, s | `1e-10` s — enforced since 2026-08-09: timestep cuts clamp to this floor and never halve below it. At-floor behavior was revised on 2026-08-10; see the floor policy below and `src/fem/doc/timestepping_strategy.md` §4. | `:3148-3150`, floor policy `:1851-1918` |
| `floor retries` | int, `0` = unlimited | `20` — backstop for the floor policy below. After this many consecutive retries at the `minimum timestep`, the run stops with a diagnosis instead of spinning on a residual that responds to neither the timestep nor a 4× iteration budget; this prevents unattended jobs from burning an allocation. Results up to the last accepted step are already on disk. | `:3284-3289`, backstop `:2040-2054` |
| `restart` | bool | `true` — a rerun resumes from an existing `memdump.hdf5`, restoring the saved state (fields, globals, thermal fields, and circuit state when present) together with the time, and prints a WARM RESTART banner naming the resume point (before 2026-08-10 this resume happened silently). **The dumped Δt is restored only up to the `initial timestep` cap** (since 2026-08-19): whenever the dump carries a larger step — the normal case for a mature run — the first restarted step runs at `initial timestep` and the run grows back under the BDF growth limiter (see the `initial timestep` row). With `adapt timestep : false` there is no growth path, so the cap is **permanent**: such a deck restarts at `initial timestep` and stays there — safe, but not a resumption of the dumped step size. Opt out with `restart : false ;` to ignore the dump and start fresh — it is then overwritten at the first save. Since 2026-08-15 the BDF step-size history is carried over as well (`bdf_h`/`bdf_step_count`/`bdf_last_dt` in the dump's `meta` group), so a warm start resumes at the full BDF order it had earned. **Since 2026-08-30 that history is MANDATORY**: a dump without the triple is refused with a named error instead of silently re-anchoring the order ramp at BDF1, because every dump this binary writes carries it and one that does not was written by an older binary. A coupled run additionally requires the thermal triple (`bdf_h2`/`bdf_step_count2`/`bdf_last_dt2`), so a coupled deck cannot resume a magnetic-only dump. Both checks fire before either BDF integrator history is restored and before the dof manager is initialized and distributed — note that fields, globals and circuit state have already been read on rank 0 by that point, so this is a guarantee about integrator state, not about the whole restore. The escape is `restart : false ;`, which ignores the dump entirely. One further integration detail is deliberately **not** carried over: the coupling factor (returns to its deck value). | `:3278-3280`, gate `:3163-3173` |
| `simulation time` | value, s | **required** — read unconditionally | `:3168-3169` |
| `adapt timestep` | bool | `true` | `:3171-3174` |
| `scheme` (legacy alias `method`; `scheme` wins) | `bdf1`…`bdf5`, `explicit`, `crc`/`crank-nicolson`, `galerkin` | **`bdf1`**: validated baseline (Messe et al. 2023 §4); guidance §4.5. Crank-Nicolson and Galerkin parse, then hard-error downstream with stiffness. The configured scheme labels the per-timestep box in the log (`BDF1`, `BDF2`, …) since 2026-08-10. | `:3178-3225` |
| `anderson stabilization` | bool master switch | **`false`** (opt-in); guidance §4.5. `true` fills default depths (magnetic 3, thermal 1) only where `anderson depth` is absent. Explicit `anderson depth` always wins, including explicit `0`. `false` zeroes both depths and hard-errors with explicit nonzero `anderson depth`. | `:3227-3258` |
| `save every` | value, s | Absent: save every step. | `:3260-3273` |

Floor policy (revised 2026-08-10): a cut requested while Δt is already at the
`minimum timestep` does not abort immediately. The controller doubles the iteration budgets
(`max iterations`, `watchdog window`, both fields), up to 4× the deck values, prints a warning
box, and retries at the minimum. The budgets return to deck values with the next accepted
step; the warning box counts retries against `floor retries`. Before the backstop trips,
diagnose with `compute conditioning` (§4) and accept a genuine noise floor with
`absolute tolerance` (§4.2).

**Not deck-exposed (internal members, `cl_FEM_Controller.hpp`).** The following Δt-controller
members still have no input keys: the step-size mode (`mUseLegacyTimestepControl`, PID by
default since 2026-08-09), PID gains (`mCtrlKp/Ki/Kd` = 0.15 / 0.30 / 0), post-rejection
growth hold (`mPostFailureHoldSteps` = 2), floor-escalation cap (`mFloorEscalationCap` = 4),
and relaxation neutral band (`mOmegaNoiseBand` = 0.05 ≈ 0.2 dB). These remain compiled-in
pending validation; changing them requires a source edit. The rationale for each is in
`src/fem/doc/timestepping_strategy.md`.

### 4.5 Choosing the machinery: Newton, Anderson, BDF order (field guidance)

The parser defaults are `algorithm : Picard`, `scheme : bdf1`, and
Anderson off. `scheme : bdf1` is the published, validated implicit-Euler
baseline (Messe et al. 2023, §4). With `algorithm : Newton`,
each timestep still starts on Picard, promotes at `tolerance switch`, and
uses the published quasi-Newton endgame to ε = 10⁻¹¹ (Messe et al. 2023,
Eq. 13). Treat every knob below as a deliberate departure: useful
on the right deck, harmful on the wrong one. Grep logs for these fingerprints.

**`algorithm : Newton` (terminal quasi-Newton after Picard promotion)**

- *Useful:* smooth decks where Picard's linear contraction stalls before
  a tight tolerance. Newton polishes the last decades cheaply once the
  iterate is inside the basin and the tangent is complete (Messe et al. 2023, Eq. 13).
- *Turn it off (use `Picard`)* on **net-transport-current decks with a
  point bearing**: the Newton tangent carries a near-null φ-gauge mode,
  so promotion can degrade a converging Picard sequence. Since 2026-08-10
  a promoted tangent whose line search rejects all eight trials falls back
  to Picard for the rest of the timestep instead of killing the step (the
  greg5 signature: `Newton n, relax 0.00391` then Δt collapse) — the run
  survives, but each affected step still pays ~9 wasted solves for the
  doomed Newton attempt, so `Picard` remains the right setting for this
  deck class.
- *Log fingerprint:* promotion makes the residual **worse** (e.g. −65 dB →
  −60 dB), then it sits flat while relax sweeps down; or a Newton line
  search rejects every trial — bit-identical residual printed with relax
  collapsing as 2⁻ⁿ — followed by a timestep cut that a pure-Picard retry
  then converges in two iterates. Both mean: this deck wants `Picard`.

**`anderson stabilization : true` (type-II Anderson mixing of the Picard branch)**

- *Useful:* slowly-contracting coupled magnet–thermal iterations, especially
  the quench-adjacent regime where plain Picard creeps for dozens of
  iterations. The mixing extracts missing curvature from the iteration
  history (Walker & Ni 2011). That is the case it was built for.
- *Keep it off* (the default) for strongly nonlinear HTS transport decks
  (power-law n ≳ 25): Anderson has no convergence guarantee on non-smooth
  maps, and extrapolating across the E–J transition can produce non-physical
  iterates. It stays opt-in until validated A/B on such decks.
- *Log fingerprint (historical):* builds before 2026-08-07 reported the
  lagged-operator residual under Anderson — instant "convergence" at
  solver-roundoff level (≈ −120 dB) in exactly `min iterations` iterates,
  with spatial checkerboarding in B/J. Current builds report the honest
  pre-update force residual, so Anderson now iterates realistically; the
  smoothness caveat above still stands.

**`scheme : bdf2` … `bdf5` (higher-order variable-step BDF)**

- *Useful:* smooth transients at fixed or gently varying Δt. Higher order
  can buy accuracy per step. The coefficients are exact for variable steps,
  with a startup ramp.
- *Prefer `bdf1`* (the default and the validated baseline) whenever the
  adaptive controller is working hard: variable-step BDF5 zero-stability
  requires step ratios near 1, while the controller applies ×1.5 growth and
  ×0.5 cuts. That combination is unproven. Flux-front decks that sawtooth
  Δt are the worst case: every cut/regrow cycle stresses exactly the
  property BDF5 does not guarantee.
- *Log fingerprint:* growing or oscillating residuals right after Δt changes,
  on a deck that behaves under `bdf1` with identical tolerances.

Rule of thumb from the greg3 campaign: when a run misbehaves, first return
all three knobs to the baseline (`Picard`, `bdf1`, Anderson off) and
re-diagnose from there. Each knob adds its own failure mode, while the
baseline's failures are the ones the guards and literature actually cover.

**When Δt collapses toward `minimum timestep`.** The controller cuts Δt
whenever an attempt fails to reach `tolerance`; a run that walks to the floor
is therefore reporting a residual that still misses the target. Two causes need
different responses, and `compute conditioning` separates them:

1. *The residual falls with Δt but not far enough* — a genuine timestep-size
   problem. Let the controller work: the cut sequence searches for a usable Δt,
   and the §4.4 floor policy gives slow attempts more iterations at the floor.
2. *The residual sits at the same floor at every Δt* — cutting cannot help.
   The run escalates its iteration budget at the minimum and then stops at
   `floor retries` (§4.4) with a diagnosis. Set
   `solver { compute conditioning : true ; }` (§4) and read the
   `|λ|max/|λ|min` row from the footer as the κ(A) proxy: if κ·ε_mach lands near the observed residual floor, the floor is
   conditioning and the answer is tolerance (accept it with
   `absolute tolerance`, §4.2), not more timestep cuts. If κ·ε_mach is far
   below the floor, the floor is physical or formulation-side — look at the
   model, not the controller.

The relaxation column is the quickest read: a residual pinned within a
fraction of a dB while `relax` decays toward `min relaxation` is case 2.

---

## 5. `materials`

One sub-section per material; the sub-section **type** is the material
label. Four shapes, tested in order:

**Every key and subsection must belong to the shape the material selected.**
Anything else is a fatal error naming the material and the shape that would
have ignored it. This closed a class where a misspelled key, a key belonging
to another shape, or a constant that lost to a `file` was dropped in silence.
Three consequences worth knowing before writing a deck:

- `RRR` is legal only on a builtin whose constructor takes it — the pure
  metals and the alloys. On `ybco`, `hastelloy`/`hastelloyc276` and
  `magnesia`/`mgo`/`buffer` it was read and then discarded, so it is now
  refused rather than ignored.
- On an HTS builtin, `file` and the constants `jc`/`n`/`ec` are exclusive.
  Giving both used to keep `file` and drop the constants without a word.
- A subsection must not carry a label: `defect { }` is read, `defect : mine { }`
  never was. The same holds for `usermat { }`.
- A material name may be defined only once. The later definition used to win
  silently and leak the earlier one.

The one deliberate exemption: a section that selected the **curve** shape may
still carry a valueless `builtin ;`, which some decks use as documentation. It
is inert, and only that exact key is tolerated.

1. **B-H curve ferromagnet** — presence of `curve` selects it:
   `curve : RoxieIron ;` (HDF5 group name), `bhfile : path.hdf5 ;`
   (default `bhdata.hdf5`). Any other key is fatal, except the tolerated
   bare `builtin` described above.
2. **Builtin** — selected by the first matching form: `builtin : <type> ;`
   inside the material section; a material section *labeled* `builtin`
   (`copper : builtin { }`); or a `materials`-level key named exactly like the
   material subsection, with value `builtin`. Types (case-insensitive):
   `aluminum`/`aluminium`/`al`, `chromium`/`cr`, `copper`/`cu`,
   `nickel`/`ni`, `silver`/`ag`, `indium`/`in`, `lead`/`pb`, `tin`/`sn`,
   `hastelloy`/`hastelloyc276`, `ybco`, `magnesia`/`mgo`/`buffer`,
   `iron`/`ferro`/`fe`, or an alloy formula (`Sn40Pb60`). Unknown → fatal.
3. **User material** — `usermat { file : lib.so ; label : name ; }` loads a
   user-material plugin. The name `buffer` is reserved (fatal).
4. None of the three → fatal, naming the material and the three shapes it
   could have selected. A `custom { }` subsection — the name `usermat` carried
   until 2026-08-29 — is refused by name with a message pointing at the rename;
   it is a diagnostic, not an alias, and does not load a material.

**Where a material file is looked for.** Every path-valued key in this section
— `bhfile`, the HTS `file`, `defect { file }`, `heating { file }`, and `usermat { file }` — is
resolved by `material::data_file()` (`fn_material_data_path.cpp`), which
forwards to `search_data_file()` (`filetools.cpp`), in this order:

1. the path as written, relative to the run directory (or absolute);
2. the same relative path below `$BELFEM_DATA/material`;
3. the **file name alone** below `$BELFEM_DATA/material`, so that
   `bhfile : MatData/bhdata.hdf5 ;` still resolves in a run directory that has
   no `MatData` of its own.

The run directory therefore always wins: a local copy overrides the shared
database. If none of the three exists, the path is passed on unchanged, so the
error names the file as written. For the two `.so` keys the fallback never
removes a lookup — an unresolved name is still handed to `dlopen`, which
searches the platform loader path (`$LD_LIBRARY_PATH`, `$DYLD_LIBRARY_PATH` on
macOS) as before.

**The same search applies to a `userdefined` source function's `file` key**
(§9), which before 2026-08-31 was handed to `dlopen` unresolved. One
consequence is worth stating plainly: `dlopen` consults the platform's
dynamic-loader search path (`$LD_LIBRARY_PATH` on Linux, `$DYLD_LIBRARY_PATH`
on macOS) only for a name carrying **no slash**, so for a bare plugin name
`$BELFEM_DATA/material` is now searched *before* it. For a name with a slash —
which is what every shipped example uses — that path was never consulted, and
steps 2 and 3 are new fallbacks that simply did not exist before.

**A caution about `$BELFEM_DATA`.** The search root is not the environment
variable itself but the global `Communicator::set_globals()` fills from it. On
an **installed** tree that global falls back to the compiled-in install data
directory when the variable is unset, so "unset `$BELFEM_DATA`" does *not*
reliably mean "only step 1 applies" — it means that only in a build-tree run
with nothing installed.

Keys inside a material section (the `Applies to` column says which shapes read each one):

| Key | Applies to | Type / unit | Notes | Site |
|---|---|---|---|---|
| `RRR` | pure metals and alloys | real | residual resistivity ratio; legal only when the constructed type is `PureMetal` (pure metals and alloys). On `ybco`, `hastelloy`/`hastelloyc276` and `magnesia`/`mgo`/`buffer` it is **refused at setup** by `check_unused_input`, not discarded | `:100-101` |
| `file` | HTS | HDF5 path | Jc(B,θ,T) and n(B,θ,T) tables; sets both functions. Exclusive with `jc`/`n`/`ec`, which are then not read | `:111-115` |
| `jc` | HTS without `file` | value, A/m² | required together with `n` | `:118-124` |
| `n` | HTS without `file` | real | required | `:126-127` |
| `ec` | HTS without `file` | value, V/m | default `1e-4` | `:129-134` |
| `resistivity type` | HTS / user-material superconductor | `powerlaw`, `power-law`, `piecewise`, `riva` — **case-sensitive** | default power-law | `:147-166, 208-227` |
| `defect { file ; label ; }` | HTS / user-material superconductor | `.so` + symbol | spatial Jc modulation | `:137-144, 198-205` |
| `heating { file ; label ; }` | any `builtin` or `usermat` material | `.so` + symbol | artificial volumetric heat load, W/m³ as a function of x, y, z (m) and t (s); added to the Joule term of the thermal load vector, so it does nothing in a magnetic-only run | `"heating"` |
| `critical temperature` | `usermat` superconductor **only** | value, K | overrides the plugin's own `T_crit`; **fatal on a builtin** | `"critical temperature"` |
| `density correction` | any material | real, dimensionless | default `1.0`; scales **only** the density in the thermal mass matrix | `:275-278` |

**`heating`.** An artificial heat source for quench studies — a heater pulse, a
disturbance, a hot spot — placed on a material rather than on a region: every block
carrying that material evaluates the plugin at each integration point and adds the
returned W/m³ to the Joule heating of the thermal equation. The plugin exports
`extern "C" void <label>_init( Material * )`, which calls `set_user_defined_heating`
with a `real f( x, y, z, t )`; coordinates are in meters whatever the mesh `unit`,
so the function gates its own spot (a tape, a disk) the way a `defect` plugin does.
The load has no temperature dependence and therefore no Newton tangent. On a
thin-shell stack, attach it to the layer that a real heater would warm — the
stabilizer, not the 1–2 µm superconductor — and remember that the heated volume is
the layer thickness under the spot on **every** tape carrying that material unless the
function gates on position. In an adiabatic deck a steady load quenches the conductor
eventually; prescribe the pulse by its energy and choose a period the timestep can
resolve. One `heating` block per material; a second load is refused.

**`critical temperature`.** The temperature above which a superconductor stops
being one, as far as the resistivity laws are concerned. It matters because the
right value belongs to the *jc/n data*, not to the carrier material: a plugin
whose fits are only valid to 90 K needs 90 K, whatever the material it is
attached to would otherwise say.

Four things about it are easy to get wrong.

- **Only the `piecewise` and `riva` laws read it.** They return the normal-state
  resistivity above it and evaluate jc and n below. The default `power-law` does
  not consult it at all, so setting the key changes nothing on a deck that keeps
  the default.
- **It may only be set on a `usermat` material.** On a builtin it is fatal, by
  design rather than by omission. A builtin's critical temperature is part of a
  calibration that its constructor has already consumed — YBCO, for one, samples
  its Callaway thermal-conductivity spline at construction and bakes the value
  into it. Overriding the constant afterwards would move the resistivity gate
  and leave the conductivity built around the old number, which is worse than
  refusing.
- **The comparison is strict.** The superconducting branch is skipped for
  `T > T_crit`; at exactly `T = T_crit` jc and n are still evaluated, so a
  plugin whose fits fall back to a normal-state constant *at* its cutoff will
  still be asked for them there.
- **A unit is required.** `92.5 K` — also `°K`, `C`/`°C`, `°F`, `R`/`°R`, all
  converted to Kelvin on read. A bare number is dimensionless and rejected; so
  is `92.5K` without the space, and a bare `F` is farad, not Fahrenheit.

Set on a material with no jc, or given a non-positive or non-finite value, the
key is fatal. On a `curve` (ferromagnet) section it is fatal too, like every
other key the curve shape does not read; the only exemption there is a
valueless `builtin ;`.

`resistivity type : riva` (2026-08-27) selects the same parallel E-J model as
`powerlaw`: the superconducting power-law channel in parallel with the
normal-state channel (Duron et al. 2004; used by Riva 2021, EPFL thesis 8754,
Eq. 5.4). It is hardened to remain finite for every input that a measured jc/n
table can produce during an iteration. A dead defect (`D = 0`) and spline
under- or overflow fall back to the fully normal branch. When measured n
softens through 1 near T_crit, the law instead uses the finite ohmic closed
form (ec/jc in parallel with the normal channel) rather than producing NaN.
The exponent n is globally floored at 1 in `Material::n_eval` (the ohmic
limit). The `dn_eval_*` derivatives return zero while the floor binds, keeping
the residual and tangent consistent. `riva` works with the `file` jc table or
with constant `jc`/`n`, using the same dependency routing as the other laws.
Since 2026-08-31 selecting `riva` also validates the n source at setup: a
database n table whose stored minimum is at or below 1 is fatal, and a
constant `n` is fatal unless it is finite and greater than 1 — such a
source would run as a plain resistor rather than a superconductor,
silently. The table test covers the whole stored table, including values
above `critical temperature` that riva itself never reads: the build
pipeline floors the entire table at 1.02, so any stored value at or below
1 marks a table that escaped it. The gate applies whichever of the law,
the table, or the constant is set last. Analytic and plugin-supplied fits
carry no cheap bound and are not checked; the stored minimum also does
not bound the interpolant between nodes — the runtime ohmic fallback
remains for both.
It differs from `piecewise` above roughly 1.4 jc: there is no Bézier flux-flow
blend, and the transition to the normal state saturates by about 3 jc instead
of stretching over decades. One deliberate deviation from Riva's published
form is that BELFEM retains its floor semantics for the minimum resistivity
(zero by default) rather than adopting the thesis's additive 1e-17 Ω·m
regularization. The additive form is itself tangent-consistent; the deviation
exists because BELFEM's convention is a floor, and carrying both conventions
in one code invites the value/tangent desynchronization that a positive floor
with an unfloored tangent caused before it was removed on 2026-08-10.

`density correction` exists for the virtual-domain idealization. A thin shell
carries its own mass and heat capacity but occupies no geometry, so the meshed
volume around it spans the full tape pitch and double-counts the shell's own
thickness as filler material. Setting the factor to gap/pitch makes that volume
carry the thermal mass of the thin layer it really contains.

It multiplies the density read by `Calculator`, which reaches exactly one
consumer: `M += Nᵀ ρ cp N dV` in `mt_thermal_h.cpp`. It is deliberately **not**
applied to electrical resistivity or thermal conductivity — those are wrong by
the inverse factor in an over-thick gap, and no single scalar can fix storage
and conduction at once. Correcting the mass while leaving conduction alone is
right only where the temperature field is near-uniform; with a hot spot or a
propagating normal zone the conduction error dominates instead.

---

## 6. `layers : <tapename>` (top level)

One block per `thinshell` label; a missing block is fatal. Read raw
line-by-line
(`cl_MaxwellFactory.cpp:2515-2534`) — repeated material names are allowed
and order is bottom → top:

```
layers : tape
{
    copper : 20 mum ;
    ybco   :  1 mum ;
    copper : 20 mum ;
}
```

The **unit is mandatory** on every thickness (fatal without it); materials
must resolve in `materials`.

Several tape types coexist as separately *labeled* blocks (`layers : tape1`,
`layers : tape2`, …); each `thinshell : <tape>` resolves its stack by that
label (`cl_MaxwellFactory.cpp:2826-2872`). This is label-keyed lookup, distinct
from the unlabeled duplicate-section rule of §1 — two *unlabeled* `layers { }`
blocks would collide on the empty label instead of accumulating.

---

## 7. `homology`

Required only when cohomologies are computed (fatal if missing then).

| Key | Values | Default | Site |
|---|---|---|---|
| `algorithm` | `pellikka`, `ccr`, `belted tree`, `generalized pellikka` (case-insensitive) | `generalized pellikka` | `cl_MaxwellFactory.cpp:927-935` |

`optimize cuts` appears in older decks but is parsed by nothing — inert.

---

## 8. `topology`

Each sub-section is a domain. The sub-section **type** selects the kind
(`en_DomainType.cpp:102-184`, case-insensitive), the section **label** is
the domain name.

| Type strings | Kind | Group key | `material` key |
|---|---|---|---|
| `conductor`, `superconductor` | blocks | `block`/`blocks` | required |
| `buffer` | blocks | ″ | required |
| `ferro`, `iron` | blocks | ″ | required |
| `coil` | blocks | ″ | — |
| `air`, `vacuum`, `void` | blocks | ″ | — |
| `background field`, `air symmetry`, `buffer symmetry`, `ferro symmetry`, `conductor symmetry`, and the four `… antisymmetry` forms | sidesets | `sideset`/`sidesets` | — |
| `thinshell`, `tape`, `shell` | sidesets | `sideset`/`sidesets` | from `layers` |
| `curve`/`curves`, `periodic` | special (below) | | |

**`coil` and `air` blocks are magnetically inert, and that is deliberate.** They
take no `material` key, and if one is supplied it is silently ignored:
`MaxwellFactory::collect_material_labels_from_domains` gives `air` the label
`"air"` and everything it does not recognize — `coil` included — the sentinel
`"inactive"`, which `create_materials` then skips entirely. So a coil block
ends up with **no material object at all**, carries its prescribed transport
current, and contributes nothing but free-space permeability. Do not try to
give a coil a conductor material to "fix" a missing-material complaint: the key
will not be read, and the block is not meant to have one.

The reason is scope. **BELFEM models high-temperature superconductors, not
low-temperature ones.** There is no LTS strand model, so a coil is a current
source rather than a material region: no resistivity, no critical surface, no
quench physics inside it. Conductors that do carry material physics are the
`conductor`/`superconductor`, `buffer` and `ferro` domains, which is why those
three are the ones where `material` is required. A consequence worth knowing
when reading code: any calculator path that dereferences a block's material
must not run for `coil` or `air` blocks, because the pointer is legitimately
null there.

`tape` and `shell` are true aliases of `thinshell`, and singular `curve` of
`curves`. They became so on 2026-08-11: the consumers used to compare the
section-type **string**, so a `tape { }` section yielded the right domain type
and then built no thin shell at all. Consumers now dispatch on the domain-type
enum, which is what makes every accepted spelling behave identically. `cut` was
removed the same day — no parser path had ever read a `topology { cut { } }`
section, so it now fails with `Unknown Domain Type: cut` instead of parsing and
doing nothing. Cuts come from the `homology` section.

`block` xor `blocks` (and `sideset` xor `sidesets`) — exactly one, and a group
key must be present. Both checks are always-active errors since 2026-08-11;
they were assert-level before, so a release build silently read an empty list.
A `label : … ;` **key** inside a domain section is inert — the name comes from
the section label.

Symmetry vs antisymmetry encodes current polarity: preserved → B×n = 0,
inverted → B·n = 0 (natural).

### `thinshell : <tape>` specifics

`sideset`/`sidesets` required (xor, fatal otherwise); the label must match
a `layers` block (`cl_MaxwellFactory.cpp:2433-2582`). Terminal curves on a
shell sideset attach automatically.

**Signed sidesets (since 2026-08-12).** A sideset id may carry a gmsh-style
negative sign:

```
sidesets : -5, -6, 7:20 ;
```

A signed sideset has the orientation of every one of its facets flipped
(master/slave swap plus winding rewrite) after the master normalization and
before the thin-shell pipeline reads the windings for the layer normals
(anchor `"read_signed_sidesets"`, consumer `flip_thin_shell_sidesets`). The
sign exists because the layer stack is laid along the facet normal, and the
master-side normalization is domain-type based — mirror-symmetric, so a tape
stack bounded by air on both faces cannot come out uniform: one outer tape
always flips relative to the rest. Which side the layers face is user intent
that no geometry rule can derive (a corc wrap has no meaningful mean normal),
so the sign carries it per sideset. An unsigned deck behaves exactly as
before.

Rules and pitfalls:

- The sign binds to **single ids only**. `-5:8` and `-5:-8` are refused with
  a fatal error — inside a range the scope of the sign is ambiguous, so the
  flipped ids must be listed individually.
- **Sign all sidesets of one connected sheet together.** A tape meshed as two
  half-surfaces (left/right of the center line) must have both halves signed
  or neither: flipping only one tears the sheet's orientation at the shared
  line, and the node-normal averaging in the thin-shell factory then cancels
  to a zero normal there.
- Downstream consumers see only absolute values — `sidesets()` on the
  protoshell is unsigned, the sign lives in `flipped_sidesets()`.
- Each applied flip is reported at Default level
  (`thin shell <label>: flipped orientation of sideset N ( M facets )`), so a
  run's log states which orientations were touched.
- **The sign selects the jc(θ) lobe (since 2026-08-16).** The solver's
  field-to-normal angle β is unfolded, β ∈ [0, π], and measured jc/n tables
  (`sp-ap.hdf5`, `sst-1.hdf5`) are asymmetric about 90° — up to ~40 % in the
  liquid-nitrogen operating window. β < 90° means the field has a component
  along +n, where +n is the master-side outward facet normal, i.e. the same
  direction the layer stack is laid along. Flipping a sideset therefore swaps
  which measured lobe the tape reads: the deck's normals must be oriented to
  match the tape's physical mounting per the measurement convention (the
  Robinson/SuperCurrent angle is sample-normal-to-field). Analytic Kim-type
  laws are even in β and unaffected; constant jc/n decks are unaffected.

| Key | Values | Default | Site |
|---|---|---|---|
| `edge coating` | bool (`true`, `on`, `yes`, `1` are true; anything else is false) | `off` | `cl_MaxwellFactory.cpp:2760-2769` |
| `edge coating width` | value, m | absent: derive from outer layer thickness when `edge coating` is on | `cl_MaxwellFactory.cpp:2770-2781`; `cl_ThinShellFactory.cpp:478-479` |

`edge coating : on` requests side connector walls (surrounding plating on the
tape slit edges) for this tape. For 3-D meshes, the flag calls
`create_side_connectors`, which builds the walls. With the flag on, the mesh
must be 3-D, the thin shell must be first order, and the tape must have at
least 3 layers. The first and last layers in the `layers` block must be the
matching outer layers: same material and same thickness (relative tolerance
`1e-6`), with positive outer layer thickness (fatal otherwise,
`cl_ThinShellFactory.cpp:437-476`). The wall material is the shared outer
layer material and must be a pure metal (`MaterialType::PureMetal`, fatal
otherwise); a non-copper coating only warns
(`cl_MaxwellFactory.cpp:2380-2397`).

`edge coating width` overrides the wall width. If the key is absent and
`edge coating : on` is set, the width is derived from the outer layer
thickness. Setting the width without `edge coating : on`, or setting it to a
non-positive value, is fatal.

### `curves`

Used only when the mesh carries no curves. Every key is a curve id:
2-D `1 : 5, 6 ;` (from sidesets); 3-D `1 : 12 @ 4 ;` (sideset
intersection) (`cl_MaxwellFactory.cpp:359-387`).

### `periodic`

Used only when the mesh carries no periodicity. `source` and `target` each
take **exactly 3 vertex ids** defining the matched planes, in corresponding
order (fatal otherwise) (`cl_MaxwellFactory.cpp:645-670`). These are **vertex
ids, not node ids**: a vertex is a gmsh type-15 point element, and its id is
that element's *geometry* tag (`cl_Mesh_GmshReader.cpp`, "we want to use the
GeometryTag as Entity ID"), resolved through `mMesh->vertex( id )->node( 0 )`
in `PeriodicityFactory::set_master_plane`. Only points whose node is also used
by another element survive the read. The gmsh `$Periodic` section is not read;
this block is the only way to declare periodicity on a `.msh`. A sideset is
tagged periodic when all of its nodes lie in one of the two planes, so the two
periodic faces must be sidesets of their own. Without this section no periodic
path is exercised — "it ran" does not mean periodicity was active.

---

## 9. `boundary conditions`

If a `maxwell` sub-section exists the Maxwell BCs are built from it, else
from the section itself; thermal BCs come from a `thermal` sub-section
(`cl_MaxwellFactory.cpp:138-149`, `cl_ThermalFactory.cpp:289-301`).

Sub-section type = BC type (case-insensitive): `neumann`, `dirichlet`,
`bearing`, `gauge`, `current`, `voltage`, `background`. `background dirichlet`
parses to a type the creation switch does not handle and is **refused with a
hard error** since 2026-08-11 (§13); any other unrecognised type is refused the
same way, instead of silently reaching the source-function loop.

**Repeated BC sub-sections of the same type are legal.** The factory iterates
sub-sections by index (§1), so a deck may carry, e.g., two `current { }` blocks
with different amplitudes for two coils, and both blocks are consumed. For
`current` and `voltage`, each bracket group in each block becomes a physical
condition and consumes one cohomology generator, so count repeated blocks and
groups together against the available generator budget.

**Mesh globals.** Scalars are written into the Exodus file alongside the
fields through two channels, both on **rank 0 only** — the rank the Exodus
and memdump writers read; workers neither create nor read them, and no value
is replicated over MPI. Both survive memdump warm restarts.

*Terminal conditions* — `current`, `voltage`, and the circuit terminal pairs —
are written by `Controller::save_IV` (since 2026-09-01), **one `I_…`/`U_…`
pair per condition**, under the same names that head the columns of
`iv_results.csv`. The pair carries the imposed value *and* its response: for a
current-driven condition `I` is the imposed current and `U` the terminal
voltage, for a voltage-driven one the reverse. The name is the block's header
label if one is given (`current : coil1 { }` → `I_coil1`, `U_coil1`), else
positional, zero-padded to the width of the generator count (`I_1`, `U_1` below
ten generators; `I_01` … from ten on). A label shared by several conditions —
one block with several bracket groups, or two blocks with one label — gets a
running suffix on **every** occurrence (`I_coil1_1`, `I_coil1_2`), so a bare
name never sits next to a numbered one. Circuit terminal pairs use their
component `label`, else `terminalpair_<n>` where `n` is the component's
1-based position in the circuit topology list (not a terminal-pair counter: a
resistor listed first makes the first unlabeled pair `I_terminalpair_2`).
Generators the deck does not drive — the incidence matrix may carry more cuts
than the deck names conditions — are positional. The pair is refreshed at every
save step, immediately before the Exodus frame is written. Terminal blocks
publish **no** global of their own: before 2026-09-01 an unlabeled `current { }`
also wrote a `current` global carrying the same number as `I_1`.

*Every other value-imposing block* — `background`, `gauge`, `neumann`,
`dirichlet`, and the thermal conditions — publishes **one mesh global per
block** (2026-08-15, narrowed from per-condition to per-block 2026-08-26)
carrying the imposed value, updated every timestep. The base name is the
header label if one is given, else the type; **if the same base occurs more
than once in the deck, every occurrence carries its section ordinal**
(`background_1`, `background_2` — never a bare `background` next to a
`background_2`). One block, one global, whatever its group count: every
condition of a block shares that block's value function, so per-group globals
would be copies of a single number. Remember §1: unbracketed id lists expand to
one group per id, but those groups do not each publish a global. One consequence
to know about: a `userdefined` block builds a separate source function object
per group, each initialized from the same (`file`, `label`) pair. If the
library's init carries per-instance side effects, so that two objects built
from one symbol evaluate differently, the global reports only the first group's
value. The physics still uses each object; only the published scalar collapses.
Thermal conditions prefix `thermal_` to their type when unlabeled, since both
kernels share one mesh. `bearing` imposes a bare node constraint and publishes
nothing. Any duplicate final name is a hard error at setup — label the sections
to disambiguate.

**The temperature publishes as `T_max` or `T_bulk`** (`T_bulk` since 2026-09-05, `T_fixed` from 2026-09-01,
`temperature` before), created at the first save and refreshed at every one.
In a coupled run `T_max` is the **live maximum of the nodal field `T`** — the
same number the console box prints. In a magnetic-only run `T_bulk` is the
bulk temperature `gTbulk` at which every material law is evaluated, taken from
`initial conditions { temperature : … }` or the executable's 77 K default; the
name says that it is held, not solved. Exactly one of the two exists in a given
file. Both names are reserved: a boundary condition labeled either aborts the
run rather than being overwritten at the first save.

**Warm restarts and renamed globals:** memdump load re-adds any dumped label
that the current run has not created, and the ID-mismatch check does not look
at the label *set*. A dump written before a naming change therefore resurrects
the retired names alongside the new ones — a pre-2026-09-01 dump brings back
`current` and `temperature`. Delete `memdump.hdf5` after any change to how
globals are named. One restart caveat: the globals' IDs follow deck order, so a
deck whose BC *order* changed across a memdump will refuse the warm restart
with an ID-mismatch error.

Domain keys (`cl_MaxwellBoundaryConditionFactory.cpp:46-179`):

| BC type | Key | Notes |
|---|---|---|
| `bearing`, `gauge` | `nodes` | id list, required |
| `current`, `voltage` | `input terminal(s)` or `input curve(s)` | required; `curve` forms mark thin-shell terminals; bracket groups = separate BCs |
| ″ | `output terminal(s)` / `output curve(s)` | optional; absent = 2-D (input list reused); group counts must match |
| `voltage` (2-D) | `length` | value, m; required, scales by 1/length |
| `background`, `neumann`, `dirichlet` | `sideset` xor `sidesets` | required |
| `background` | `direction` | 3 comma-separated reals. **Pitfall:** components are read through `std::stoi` (`cl_Input_Section.cpp:731`), so fractional components are truncated — `0, 0.5, 0` becomes `0,0,0`. Use integer component ratios until fixed |

A `bearing` fixes the magnetic potential at its nodes; in practice it removes
the null space, the free additive constant in the potential, from the linear
system. **As of 2026-08-28 the bearing is OPTIONAL:** when no bearing section
is present, gauge pinning is automatic — one pin per connected φ component
(the air region and each enclosed buffer patch), chosen as the node farthest
from any sideset, computed once on the fresh mesh and persisted in the `.bfm`
(dataset `pinned`), so reloaded meshes keep their pins. A deck-written
bearing takes precedence and **disables all automatic pins** — including the
buffer-patch anchors, which then float. **Writing a bearing by hand is an
expert function and is not recommended.** The automatic pinning covers every
φ component; a hand-written bearing names one node and silently gives up the
pins for every component the deck does not mention, so the potential floats
wherever you did not look. Reach for it only when you need a specific node
held and can say why. A reloaded `.bfm` written before the autopin feature
carries no pins: without a bearing the run warns that the potential is
unpinned (delete the `.bfm` to recompute, or set a bearing).

The imposition machinery (bearing, thermal, `gauge`) is periodicity-aware
since the same date: a pin landing on a periodically condensed node is
rerouted onto the node's partner (the one independent dof of the pair; the
named node follows automatically), reported once on the console. Any other
constrained target — condensed onto several sources, or with a non-unit
weight — is a fatal error asking for an unconstrained node, as is a target
that names a nonexistent point or a node without dofs (previously all silent
no-ops). For `current`/`voltage`, omitting the output
list selects the 2-D form:
the input groups are reused for the return path
(`cl_MaxwellBoundaryConditionFactory.cpp:123-128`).

**Sign of a declared current.** In 2-D, a positive `amplitude` drives the
current along **+z**, out of the plane, so that the field circulates counter-
clockwise around the conductor by the right-hand rule. In 3-D the direction is
set by the terminals instead: a positive amplitude flows from the `input`
terminal to the `output` terminal, and the amplitude's sign reverses that.
Mixed-polarity models — the two halves of a pancake cross section, a return
leg, an antisymmetric pair — set the polarity by the sign of `amplitude` per
condition, since each condition owns its own cut.

A `voltage` condition drives the same cut, so it inherits the convention: a
positive 2-D `voltage` amplitude drives a current along +z, the same direction
a positive `current` amplitude gives. Current- and voltage-driven models of the
same geometry therefore agree on polarity.

Both statements are as of 2026-08-28, when the 2-D convention was corrected:
until then a declared positive 2-D current ran along −z. The inversion was
uniform across bulk and thin-shell conductors, so sign-insensitive magnitudes
(|B|, |J|, J/Jc magnitude, losses, temperature) were unaffected and relative
polarity within a model was preserved — but every signed current and field
direction an older 2-D run reported was reversed, which matters wherever
absolute direction does. Decks that carried compensating signs to work around
the old behavior need those signs removed.

Source-function keys. These apply to every BC type **except `bearing`**, which
reads only `nodes` and ignores a source function entirely (true in both the
Maxwell and the thermal factory). Amplitude dimension: `current` → A,
`voltage` → V, `background` → A/m **or tesla** (see below), other **Maxwell**
types dimensionless. The thermal factory uses its own table — see the thermal
note below.

**A `background` amplitude may be written as a flux density.** The natural way
to describe a background field is "1 T", not "795774.7 A/m", so the Maxwell
factory accepts either. A Tesla-family unit — `T`, `mT`, `muT`, `kT`, `MT`, `G` —
is read as *B* and converted to the field strength the weak form imposes,
`H = B/µ₀`, using `constant::nu0`. Anything else is still required to carry the
A/m dimension. Four things to know:

- **Write the unit with a space.** `1 T` parses; `1T` is read as the bare number
  `1` with no unit and then fails against A/m, because the value/unit split is
  by whitespace.
- **`µ₀` is the CODATA value** (`1.25663706212e-6`), so `1 T` is `795774.71503`
  A/m. The Oersted spelling `10 kOe` gives `795774.71546` A/m — the Oersted is
  defined against the retired exact `4π×10⁻⁷`. The two differ by 5·10⁻¹⁰
  relative, far below any solver tolerance, but they are not the same number.
- **The conversion assumes the boundary is not inside a magnetic material.**
  `B = µ₀H` holds outside one; a far-field boundary bordering a `ferro` block
  normally draws a warning at setup, and the run continues with `H = B/µ₀` as
  imposed. State the amplitude in A/m if you want to reason about the local `B`
  there. Treat the warning as a courtesy, not a guarantee — it is printed at the
  default log level (so `-v 0` hides it), and it stays silent on a sideset whose
  facets or block links it cannot read, which a `.bfm` may present.
- For a `userdefined` source the dimension comes from its `units` key, not from
  `amplitude`, so `units : T ;` is what converts.

| Key | Applies to | Notes |
|---|---|---|
| `type` | all but `bearing` | **case-sensitive**: `ramp`, `sigmoid`, `sine`, `square`, `triangle`, `sawtooth`, `constant`, `userdefined`; unknown = fatal. Absent leaves the `SourceFunction` in its UNDEFINED state — no usable source is selected, so treat `type` as required in practice |
| `amplitude` | all but `bearing`, and except userdefined | required, dimension as above |
| `period` / `frequency` | ramp, sigmoid (period required); periodic family (one of the two required) | `period = 1/frequency`. **If both are given, `period` wins** — `frequency` is only consulted when `period` is absent |
| `offset` | ramp, sigmoid | required, s |
| `phase` | periodic family | rad (write `deg`, it converts); default 0. **A positive phase advances the waveform**, following SPICE: the shape is evaluated at `t + phase/omega`, so a quarter-period phase gives the zero-phase wave sampled a quarter period later. This holds for all four periodic types. Do not confuse it with `offset`, which is a start *delay* and belongs to `ramp`/`sigmoid` only |
| `fuzzyness` | sigmoid | default `0.01`. The only knob on the logistic rate, which `function_sigmoid` builds from `fuzzyness` and `period` |
| ~~`expk`~~ | — | **REMOVED 2026-08-29 and now refused.** It was parsed and stored for years without ever reaching the waveform. **Delete the key to keep the results you get today.** Do not translate it: a deck carrying `expk` has always run at whatever `fuzzyness` said ( `0.01` when absent ), so writing `fuzzyness = 1/(1+expk)` would CHANGE the waveform rather than preserve it. That formula is the equivalence `expk` was *meant* to express — `expk = 99` is `fuzzyness = 0.01` — and is given only for someone who intended the slope it names |
| `file` / `label` / `units` | userdefined | `.so` path, symbol, and a **unit token** whose dimension must match; all required. The path is resolved exactly like a material plugin's — run directory, then `$BELFEM_DATA/material`, then the name alone there, then handed to `dlopen` unchanged (§5). Before 2026-08-31 it went to `dlopen` unresolved |

`boundary conditions → thermal`: `bearing`/`gauge` (`nodes`) and
`dirichlet` (`sideset(s)`) with the same source-function keys; `neumann` is
recognized but hard-errors (not implemented). **Both `dirichlet` and `gauge`
take their amplitude in K** — a thermal gauge imposes a Dirichlet temperature
on its nodes, so it carries the same dimension as a thermal Dirichlet. (Before
2026-08-11 the units switch skipped `gauge`, leaving its amplitude checked
against an empty required unit.) This is why the "others dimensionless" rule
above is a *Maxwell* rule: amplitude dimension depends on the consuming factory,
not on the BC type alone.

---

## 10. `initial conditions`

| Key | Type | Behavior | Site |
|---|---|---|---|
| `t` / `temp` / `temperature` | value, K (`C`, `°F`, `R` convert) | Sets the bulk temperature `gTbulk`; a magnetic-only run publishes it in the Exodus output as the global `T_bulk` (§9). If no temperature is set — the section absent, or present with none of the three recognized keys — `belfem` substitutes 77 K in either mode and prints a console note. A deck with thermal boundary conditions but no thermal solver section is rejected as inconsistent. (The retired `hphiTrun` aborted at thermal setup instead.) | parse `cl_MaxwellFactory.cpp:154-170`; 77 K fallback `belfem.cpp:152-161`; rejection `belfem.cpp:123-136` |

Other keys in this section are silently ignored.

---

## 11. `circuit` (optional)

`number of nodes` is required. `circuit → topology` holds one sub-section
per component; the **type** selects it (case-insensitive): `resistor`,
`capacitor`, `inductor`, `switch`, `voltage source`, `current source`,
`diode`, `superconductor`, `terminal pair`
(`cl_ElectricalCircuitFactory.cpp`). Do not label component section headers
(`resistor : R1 { }` fails the type lookup) — use the `label` key.

Every component: `label` (optional), `node +`, `node -` (ints
< number of nodes, required).

| Component | Keys (dimension) | Defaults |
|---|---|---|
| `resistor` | `value` (Ω / Ohm) | required |
| `capacitor` | `value` (F), `order` | order 1 |
| `inductor` | `value` (H), `order` | order 1 |
| `voltage source` / `current source` | same source-function keys as §9, amplitude V / A; `type` required | |
| `switch` | `initial state` (`closed`/`open`, case-sensitive), `switch time` (s) | both required |
| `diode` | `Is` (A), `Vt` (V) | `0.1 mA`, `26 mV` |
| `superconductor` | `Ic` (A), `n`, `Ec` (V/m), `length` (m) | all required |
| `terminal pair` | `input`/`output terminal(s)`/`curve(s)`, **exactly one bracketed group** (since 2026-08-15) — `[1,2]`, not `[1],[2]` and not the unbracketed `1,2`, which §1 expands to one group per id. More than one group is a hard error at setup: each terminal pair is one lumped component, the controller pairs conditions with components by position, and extra groups would silently overrun that pairing. For several pairs, define several components. | |

`circuit → output`: `file` (required), `currents` (component **labels**,
comma-separated), `voltages` (node id list)
(`cl_ElectricalCircuitFactory.cpp`, anchor `set_output_voltages`).

### 11.1 Hybrid mode: `circuit → file` (since 2026-08-27)

Setting `file : <netlist>.cir` switches the section into **hybrid mode**
(anchor `read_circuit_from_netlist`): the lumped topology comes from an
ngspice netlist instead of the deck, and the deck contains only what SPICE
cannot express.

```
circuit
{
    file : tapestack.cir ;

    topology
    {
        terminal pair
        {
            label : Z1 ;
            node + : tape ;          // netlist node NAME, not an index
            node - : 0 ;             // "0" and "gnd" are the same node
            input curves  : [1,2] ;
            output curves : [3,4] ;
        }
    }

    output
    {
        file : CircuitResults.txt ;
        currents : Is, L1, Z1 ;      // netlist instance names + pair labels
        voltages : tape, 0 ;         // node NAMES here too
    }
}
```

The following violations are hard errors:

- `number of nodes` must **not** be set — the count is derived from the
  netlist.
- `topology` may contain **only** `terminal pair` sub-sections; every lumped
  element (R, L, C, sources, diode, superconductor, timed switch) belongs
  in the netlist. `topology` may be omitted entirely for a purely lumped
  circuit.
- Every node reference — terminal-pair `node +` / `node -` and
  `output → voltages` — is a netlist node **name** (case-insensitive;
  ground is `0` or `gnd`), resolved through the netlist's node map. Raw
  indices are not accepted: BELFEM assigns ground the *last* internal
  index while SPICE calls it node `0`, so an integer would be ambiguous.
- Hybrid mode is **case-insensitive throughout**: netlist instance names,
  terminal-pair labels, and the `output → currents` list are all folded to
  lower case, and the folded form is what appears in output headers and
  exodus globals. A terminal-pair label that collides with a netlist
  instance name (after folding) is a hard error — it would make the current
  output ambiguous.
- The circuit CSV voltage columns are labeled with the resolved internal
  index (`V0`, `V1`, ...); rank 0 prints the name→column mapping at
  startup so each column can be identified.

A `.tran` card in the netlist is ignored with a log message — the deck's
`solver → timestep` controls time integration (SPICE's `tstep` is a print
increment, not a solver step). The netlist grammar itself — supported
cards, the `* belfem:` directives for superconductor, timed switch and
per-instance BDF `order`, and the refusal list — is documented in
`src/circuit/doc/` and in the parser headers (`cl_NetlistParser.hpp`,
`cl_NgspiceCircuitFactory.hpp`). Two classic traps worth repeating: SPICE
`F` is *femto* (a 100 F capacitor is written `100`, never `100F`) and `A`
is *atto*, never Ampere.

---

## 12. Aliases and precedence

| Preferred | Alias | Rule |
|---|---|---|
| `scheme` | `method` | `scheme` wins when both present |
| `tolerance` | `relative tolerance` | alias consulted only when `tolerance` absent |
| `block` / `sideset` | `blocks` / `sidesets` | exactly one (xor) |
| `input terminal` | `input terminals`, `input curve`, `input curves` | first found wins; `curve` forms = thin shell |
| `period` | `frequency` | `period` wins when both present; otherwise period = 1/frequency |
| `linear` | `linear magnetic`, `linear thermal` | specific section wins, `linear` is the fallback — so `linear` is required unless the specific section exists |
| `nonlinear` | `nonlinear magnetic` | specific section wins |

---

## 13. Known dead keys and pitfalls

- `nonlinear [magnetic] → absolute tolerance` — dead until 2026-08-09, now
  live (see §4.2); kept here as a historical warning so old advice
  ("magnetic termination is relative-only") is not re-propagated.
- **A rerun in a directory that already holds `memdump.hdf5` continues the
  previous run** — `timestep → restart` defaults to `true` (§4.4). Since
  2026-08-10, a WARM RESTART banner names the resume time; for a fresh run,
  set `restart : false ;` or delete the dump. Before that date the resume was
  silent, which is how the behavior went unnoticed.
- `homology → optimize cuts` — parsed by nothing; inert legacy key.
  Removed from all in-tree decks 2026-08-07; external decks may still
  carry it (harmless).
- `label : … ;` key inside `topology` sub-sections — inert (name comes
  from the section label). Removed from in-tree decks 2026-08-07. NOT to
  be confused with the live `label` keys in `materials → usermat`/`defect`,
  circuit components, and userdefined BC sources.
- `background dirichlet` — parses to an enum value that the Maxwell creation
  switch does not handle. **Refused since 2026-08-11:** it now hard-errors
  naming the type. Until then it created no BC *and* let the following loop
  write a source function into a neighboring condition (or underflow the
  index when it was the first), so it silently corrupted an unrelated BC.
  `impose_bc` still carries a BackgroundDirichlet branch, so the enum and its
  dead code survive; do not read that as work in progress.

  **The background field is imposed weakly by design, and the strong variant
  was deliberately dropped (2026-08-31).** `background` adds its own stiffness
  matrix rather than writing nodal values, which is why it composes with any
  gauge treatment and works unchanged on a periodic model. Imposing the same
  field strongly would mean fixing φ to the spatial value −B·r on every
  boundary node, and that costs two things. First, the automatic gauge pins
  would have to be removed: with φ prescribed over the whole boundary the
  potential is already determined, so an autopin fixing an interior node to
  0.0 becomes an inconsistent extra constraint rather than a redundant one.
  Second, the source value is spatial, so on a periodic model each condensed
  target needs the value at its *source's* position under the periodic
  transform — a gradient derivation over the φ domain, not a copied constant.
  The whole return on that is a handful of prescribed dofs. It is not worth
  the complexity, and the weak form is the supported path.

  A corollary for anyone tempted to repair the dead branch cheaply: it must
  **not** be routed through `fem::pin_dirichlet_dof`. That helper pins a
  condensed dof's source to a single constant, which is correct for a
  constant-valued Dirichlet sideset and wrong here — it would write the
  target's value onto the source and so impose a wrong constraint on the
  opposite face.
- `direction` components pass through `std::stoi` — fractional values
  truncate to integers (`cl_Input_Section.cpp:731`).
- A key written without `:` (`order 1 ;`) becomes a flag named `"order 1"`
  — a silent no-op, not an error. (In-tree instances removed 2026-08-07;
  they sat in `terminal pair` sections, where `order` is not consumed
  anyway.)
- ~~The `block`-xor-`blocks` (and `sideset`-xor-`sidesets`) rule is
  assert-level, so a release build silently uses the singular when both
  keys are present.~~ **Fixed 2026-08-11:** both checks are always-active
  errors, so a missing group key and both-keys-present each fail the same
  way in debug and release (`cl_FEM_Domain.cpp:84-90`).
- Value case-sensitivity is inconsistent by consumer: `algorithm`,
  `coupling`, BC/source `type`, `resistivity type`, `initial state` are
  case-sensitive; `library`, `scheme`, domain types, builtin material
  names are not. When in doubt, match the tables above exactly.

---

## 14. Extending this document

When adding an input feature:
1. Add the key to the right section table with type, unit dimension,
   default, and the parse-site `file:line`.
2. State the behavior when the key is absent (default vs fatal).
3. If the key supersedes another, add the pair to §12 with the precedence
   rule, and keep the old form parsing as an alias.
4. If a key is retired, move it to §13 rather than deleting it silently —
  decks in the wild will still carry it.
