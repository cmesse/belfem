# AC Loss Postprocessing — Energy, Decomposition and Honesty on Top of the Power Global

> **CLOSED 2026-09-03** (todo/ currentness sweep, round 3): purpose served — `dotQ` is collected every timestep (`Controller::collect_dotQ`, `cl_FEM_Controller.cpp:2333`) and written to the exodus output as a global; the energy/decomposition refinements were judged not worth a live plan. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-31 (rewritten 2026-09-01 against `7323d528`)
**Purpose:** Turn the instantaneous dissipated-power global into a reportable AC loss: integrate it
over a cycle, decompose it per block, disclose where the resistivity clamp corrupted it, and
document it. The volumetric summation and the MPI reduction are **already done**.
**Module:** `src/fem/maxwell`, `src/fem/kernel` (controller), documentation
**AIs involved:** Claude (survey + plan), peer session (read `7323d528` from the Mac)
**Status:** OPEN — scope reduced. The original plan assumed nothing computed loss; that was true at
`bc6fb67d` and false at `7323d528`, which lands the `dotQ` global. R1/R2 of the first draft are **done in
that commit** and struck below. What remains is smaller, and is the part that makes the number
publishable rather than merely present.

> **Scope guards:**
> - Postprocessing and reporting. No change to the assembled physics.
> - Out of scope: hysteretic/coupling loss decomposition by mechanism.
> - **Read `devlog/dl20260901_heatloss_and_iv_globals.md` (in `7323d528`) first** — it records that
>   commit's own jury findings, including a `save_IV` second element pass left open by decision and
>   a two-rank run still owed.

---

## 1. What `7323d528` already provides

`save_dotQ( Calculator *, const real )` (`mt_maxwell_h.hpp:33-36`) accumulates into a mesh
global:

```cpp
aCalc->mesh()->global_variable( "dotQ" )->value() += aDotQ ;
```

fed from each kernel's integration-point loop by `dotQ += rho * dot( j, j ) * wdV`.

| Level | Handled? | Mechanism |
|---|---|---|
| integration points | ✅ | local `dotQ` accumulator |
| elements | ✅ | `+=` into the global |
| MPI ranks | ✅ | `Controller::collect_dotQ()`, `cl_FEM_Controller.cpp:2329-2344`, guarded `mCommSize > 1` |
| **timesteps** | ❌ | `reset_dotQ()` zeroes it before every `compute_jacobian_and_rhs()` |
| **per block** | ❌ | one scalar for the whole mesh (`cl_MaxwellFactory.cpp:1039`) |

Kernel coverage: `h_picard`, `h_newton_mu0`, `h_newton_mu`, `h_side_connector`,
`h_side_connector_newton`. **Not** `h_ghost` — correctly, as that is a Nitsche interface penalty
whose rho enters a stabilization coefficient rather than a volumetric dissipation term. Nothing in
`mt_maxwell_phi.cpp` accumulates, correctly, the phi domain being non-conducting.

### 1.2 The name `dotQ` is reused, deliberately

`dotQ` was already taken: `cl_IWG_StaticHeatConduction.cpp:63` sets `mFluxFields = { "dotQ" }` and
`:181` reads it via `node_data( "dotQ" )` — verified here, and it is the only other `"dotQ"` in the
tree. That one is a **nodal field on a static heat conduction problem**; the new one is a **mesh
global on a magnetic problem**. Not a collision — `Mesh` keeps `field_exists()` and
`global_variable_exists()` in separate namespaces, Exodus writes `EX_NODAL` and `EX_GLOBAL`
separately, and in a coupled run they live on different meshes. Recorded so that a reader meeting
`dotQ` twice knows it was examined rather than overlooked.

**Bottom line:** the stored value is the instantaneous **power** `∫ρ|j|²dV` in **watts**, at the last
assembly of each timestep, summed over ranks. That is the right quantity, and it is not yet AC loss.

### 1.1 Units — verified, because everything below depends on it

`mJ = C( aIndex ) * q()` (`cl_FEM_Calculator.hpp:2586`): the curl operator against the H-field edge
dofs, so `j` is A/m². Then `ρ·|j|²·dV` = Ω·m · (A/m²)² · m³ = **W**. Independently corroborated by
the thermal kernels, which push the same product into `f()` of a thermal problem
(`mt_thermal_h.cpp:49,124`) — an RHS that must already be in watts. Confidence: high.

## 2. Ordered steps

- [x] ~~**R1** — element-level reduction~~ — done in `7323d528`
- [x] ~~**R2** — MPI reduction over owned elements~~ — done, `collect_dotQ()`
- [ ] **R3 — Time integration.** AC loss is energy per cycle [J]; the global is power [W]. Nothing
      in the tree integrates it. Accumulate `∫P dt` over the period with the active timestepping
      scheme's weights, not a naive rectangle rule, so the reported energy is consistent with the
      integrator that produced the power samples.
- [ ] **R4 — Per-block decomposition.** One number for the whole mesh cannot separate tape from
      stabilizer from former, which is the first question anyone asks of a loss figure. Requires a
      per-block accumulator rather than a single global.
- [ ] **R5 — Clamp contamination flag** (LOW priority — see O2; the clamp is unreachable with stock
      materials, so this is defensive cover for plugin materials, not a correctness fix).
- [◐] **R6 — Documentation and schema.** **Done upstream, unverified here.** A peer session reports
      §9.3 "Ohmic Dissipation Global (dotQ)" now exists in
      `src/fem/maxwell/doc/maxwell_usage_guide.md`, with a pointer from the module README, covering
      units (watts), the mesh-wide scope, the lifecycle, the last-assembly-not-average consequence,
      kernel coverage, and the clamp caveat naming `dotQ` a lower bound. Confirm after pulling.
      Still open: a schema row **if** a deck key is ever added to control this — per the two-artifact
      rule in `CLAUDE.md`, `doc/input_file_reference.md` and `doc/input_schema.yaml` change together.
- [ ] **R7 — Validation.** A Norris/Brandt strip solution on `2D_Tapestack`, the cheapest gate deck.
      `7323d528` is a plumbing commit; there is no evidence the number has been checked against
      anything. Until it is, this is an unvalidated output.
- [ ] **R8 — Two-rank gate.** The commit's own devlog records that a two-rank run is still owed, and
      `collect_dotQ()` is guarded by `mCommSize > 1` — so the reduction path is exactly the part
      no one has executed. `make check` runs multi-rank only for `tests/comm`, so a green suite is
      not this gate.

## 3. Open design questions

- [x] ~~**O1 — which resistivity, and does it work isothermally?**~~ **RESOLVED by `7323d528`:**
      the magnetic side, via `mx->compute_rho( k )` in `mt_maxwell_h.cpp`. It therefore works in an
      isothermal magnetic-only run, which is how AC loss is usually measured. Good outcome.

- [ ] **O2 — The clamp is unreachable with stock materials and stock bounds. Keep the flag anyway,
      as cheap insurance for the one case that can reach it.**

      **Two earlier framings in this plan were wrong and are retracted. Both are kept, per the
      convention that negative results prevent re-audits.**

      - ~~"the clamp under-reports precisely the quench regime"~~ — **FALSE POSITIVE (retracted
        2026-09-01).** A quenched HTS sits at ~1e-6…1e-2 Ω·m; the cap is 1e10, insulator territory.
      - ~~"an E-J power law driven deep over-critical exceeds 1e10 at a modest multiple of Jc"~~ —
        **FALSE POSITIVE (retracted 2026-09-01).** True of `rhoPL`, and irrelevant to the returned
        value, which is the point I missed.

      **Why it cannot fire — verified in this checkout at `bc6fb67d`.** All three HTS law families
      put the power-law channel *in parallel* with the normal-state channel, so the power law can
      only ever pull the result **down**:

      - `rho_powerlaw`, all eight overloads (`powerlaws.hpp:193, 214, 242, 269, 290, 309, 331, 348`):
        `return 1.0/((1.0/rhon) + (1.0/rhoPL)) ;`
      - `rho_riva`, with the author's own comment *"parallel combination, branch-stable against a
        huge ρPL"*: `rhoPL > rhon ? rhon/(1.0 + rhon/rhoPL) : rhoPL/(1.0 + rhoPL/rhon)`
      - `rho_piecewise` returns the unbounded `rhoPL` branch only **below** its knee, where `rhoPL`
        is bounded by construction; above it returns `rhoFF`, then `rhon`.

      Algebraically `1/(1/ρn + 1/ρPL) = ρn·ρPL/(ρn+ρPL) ≤ min(ρn, ρPL) ≤ ρn`, so **ρ ≤ ρn always**,
      however far over-critical an iterate drives `|j|`.

      **The author states this in the source, in both prose and code.** Cited by greppable sentence
      rather than line number, per the anchor convention in `CLAUDE.md` — line numbers go stale, and
      each of these strings is unique in the tree:

      - the design note, in `powerlaws.hpp`: *"past a cap the parallel combination is ρn to machine
        precision, so the residual returns ρn and the tangents return the matching normal-branch
        values instead of the raw inf/inf"*
      - the same statement made executable, at the cap in `riva_rho_pl`: *"past this cap 1/ρPL
        vanishes to machine precision against any physical ρn — the caller takes the fully-normal
        branch"*, guarding `if ( ! ( lg <= 250.0 ) ) return false ;`
      - and `cl_Communicator.cpp`: *"The defaults make the clamp a no-op"*

      **Two mechanisms, one destination.** `rho_riva` early-outs — past the cap `riva_rho_pl`
      returns false and `rho_riva` returns `ρn` outright. `rho_powerlaw` has no early-out and needs
      none: an overflowing `ρPL` sends `1/(1/ρn + 1/ρPL)` to `ρn` by arithmetic alone.

      The guard is also deliberately **NaN-aware** — the comparison is negated (`! ( lg <= 250.0 )`)
      so that a NaN `lg`, from a NaN `n` in a measured table, a negative `ec`, or `inf·0` at
      `J == jc` with infinite `n`, also lands on the fully-normal branch, since NaN fails every
      ordinary comparison. Worth knowing independently of the clamp question: it means a corrupt
      jc/n table degrades to normal-state rather than poisoning the assembly with NaN.

      **The two ways it can still fire**, neither reachable from a deck:
      1. a **user plugin material** whose normal-state `rho(T)` exceeds 1e10 Ω·m. No shipped
         material comes near it — verified, and there is no insulator material class — but a plugin
         is arbitrary code (and plugins are already known to be ABI-fragile).
      2. an **executable narrowing the window** via `Communicator::set_globals()`. None in this tree
         does — verified, zero hits in `src/executables/`.

      The low side is not "never" either: `gRhoMin = 0` trips on a **negative** resistivity, which
      would be a material defect worth surfacing rather than a guard doing its job.

      **Consequence for R5:** the flag is no longer a correctness fix, it is a cheap assertion that
      the no-op stays a no-op. Priority drops accordingly — but it does not drop to zero, because
      the plugin case is real and is precisely the case where nobody would otherwise find out. One
      accumulator settles it for any deck. **What was wrong in both retracted framings was
      asserting a frequency; the flag measures one instead.**

- [ ] **O3 — Thin shells: likely already covered, confirm after pulling.**
      Thin-shell tape **layers** are ordinary conductor blocks dispatched through
      `h_picard`/`h_newton`, so they should accumulate; coating connectors are covered by the two
      `h_side_connector` kernels. If so, O3 of the first draft is resolved and the geometries
      BELFEM targets are included. Confidence: medium — this rests on a peer's static read of a
      tree not present here, plus the layer-blocks-are-Blocks plumbing. **Verify directly against
      `mt_maxwell_h.cpp` after pulling `7323d528`.**

- [ ] **O4 — Where does the energy live?** A second mesh global, a per-block vector, or a
      time-history column. Interacts with R4; decide them together.

## 4. Definition of done

- [ ] Cycle energy [J] reported, not just instantaneous power [W]
- [ ] Per-block breakdown
- [ ] Clamp contamination reported, not hidden (O2)
- [◐] `dotQ` documented with its **units** stated (done upstream in the maxwell usage guide,
      §9.3 — verify after pulling); schema row only if a deck key is added (R6)
- [ ] Validated against a published strip solution (R7) — not "it produced a number"
- [ ] Two-rank gate run (R8)
- [ ] `src/fem/maxwell/doc/README.md` postprocessed-field list updated — **losses** was removed from
      it on 2026-08-31 on the evidence then available, which `7323d528` has since overtaken

## 5. Provenance of this plan

Everything about `7323d528` here came from a peer session reading Christian's Mac; the commit is
**not on this filesystem** and is not pushed to `origin` or `backup`. Independently verified in this
checkout at `bc6fb67d`: the clamp mechanism and its tangent-zeroing precedent, `rho_clamped()`'s
existence, `gRhoMin`/`gRhoMax` not being deck-settable, and the units derivation in §1.1. Everything
else — the `save_dotQ` body, the call sites, `collect_dotQ()`, and the kernel coverage list
— is **unverified here** and must be re-read after pulling. Reviewed, not verified; nothing built
or run.
