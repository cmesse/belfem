# Todo Directory Sweep: Active Set 47 → 33

**Date:** 2026-08-26
**Purpose:** Record the currentness sweep over `todo/*.md` that moved finished plans to
`closed/` and parked ones to `deferred/`, and the deletion questions it surfaced.
**Module:** cross-cutting (todo/)

## Method, stated because it bounds the result

Classified from **each file's own Status line and its checkbox state**, cross-referenced against
the live and archived register rows. Two mechanical signals did most of the work:

1. **Checkbox census** — open `- [ ]` against ticked and struck, per file. Ten files carry zero
   open boxes; those are reports and references, not plans.
2. **Register cross-reference** — which plans are named by a *live* (unstruck) row versus only by
   archived rows. A plan cited by a live row still has forward work by definition.

Neither signal is sufficient alone, and both were checked against the prose. `falsification_tooling.md`
is the clearest example of why: **zero open checkboxes, but its prose says items (3) and (4) are
still open** and Christian's review is owed. It stayed active. A checkbox count is not a status.

## Moved (14) — nothing deleted

**→ `closed/` (7).** Work finished, nothing forward remaining:

| file | why |
|---|---|
| `dr92_restart_step_ramp.md` | COMPLETE 2026-08-19, verified by execution, all four acceptance criteria, five-step round finished |
| `thermal_picard_freeze_dr08.md` | its own Status says CLOSED; DR-08 struck |
| `bdf_nonlinear_mass_verification.md` | Status says CLOSED; DR-33 and DR-64 both struck; its one open step V1 is itself BLOCKED — the driver it names is gone |
| `tapestack3d_J_asymmetry.md` | investigation complete, verdict delivered, nothing to do |
| `coulomb_gauge_penalty_report.md` | a closing report by its own title; ends with a recommendation, by scope |
| `instruction_doc_currency.md` | R1–R4 all done, DR-60 struck; only the open *question* O3 remains |
| `nonlinear_iteration_strategy_near_quench.md` | its own Status: "THEORY REFERENCE — keep it for its literature argument, not as a work plan" |

**→ `deferred/` (5).** Sound, parked, nothing scheduled:

`2d_thinshell_todo.md` (PAUSED; the sprint calendar expired and the file says to read it as an
ordering, not a schedule) · `ngspice_parser_plan.md` (a pure proposal — no parser, no
`circuitrun`, nothing written) · `jacobian_init_remaining_bottlenecks.md` ("Open, not scheduled …
only worth attacking if that 23 s starts to matter" — the textbook deferred shape) ·
`handoff_scls_strumpack_openmp.md` (an SCLS **build-system** issue, explicitly "not a BELFEM one")
· `gauging_context_from_the_quench_run.md` (an evidence pack handed forward to gauging work that
nobody has scheduled).

## Two calls I did NOT make on my own

- **`iterate_refactor_plan.md`** reads "ACTIVE, but **substantially overtaken**". Deferring it is
  defensible, but the file declares itself active and a sweep should not overrule a plan's own
  status line on a judgement call. Left active, flagged here.
- **`nonlinear_iteration_strategy_near_quench.md`** went to `closed/` because `todo/` is for
  forward work, but `doc/` may be its better home — it is a literature argument, not a record of
  finished work. Worth a ruling.

## Link hygiene

Every reference inside `todo/` was rewritten to the new paths. Three **pre-existing** broken links
in `todo/README.md`, left by *earlier* moves that never updated the index, were repaired:
`anderson_picard_acceleration_plan.md`, `2d_thinshell_gap_analysis.md` and
`coreduce_complexPellikkaGeneralized_performance_findings.md`. **`todo/` now resolves zero broken
markdown links** — checked mechanically, not by eye.

The archived register (`debt_register_closed.md`) had its paths rewritten too. That mutates a
historical record, which is worth naming: it repairs a *pointer* so it still resolves, and changes
no claim.

## The deletion question, which is really a repo-boundary question

**`fvm_implementation_plan.md` and `fvm_module_next_steps.md` describe `src/fvm`, and `src/fvm`
no longer exists.** Verified on disk: the directory is absent, `nonfree/fvm` is present, and DR-41
records Christian's 2026-08-14 ruling that FVM is out of the open 1.0 release.

So these are two open-source planning documents for a module that is now proprietary. By the same
rule that sends nonfree devlogs to `nonfree/devlog`, they belong in **`nonfree/todo/`** — which
already exists and already holds two FVM files (`fvm_mpfa_o_vs_fem_prompt.md`,
`fvm_mpfa_o_vs_fem_recommendation.md`, both about the Manta spacecraft-thermal application rather
than the HTS quench shell, so there is no duplication).

**EXECUTED on Christian's ruling ("yes, move the fvm plans to nonfree/todo").** Both files were
copied to `nonfree/todo/`, byte-compared against the originals, then `git rm`'d from the
open-source tree. Each gained a relocation banner, and their internal paths were rewritten:
`src/fvm` → `nonfree/fvm` (7 and 4 occurrences), and the cross-references between the two plans.

**The move repaired a reference that was already broken.** Both plans cite
`fvm_mpfa_o_vs_fem_recommendation.md`, and that file has always lived in `nonfree/todo/` — so the
citation never resolved from the open-source tree. It resolves now. That is independent evidence
these documents belonged on the other side of the boundary.

One claim was deliberately **not** rewritten: both plans assert things about `src/CMakeLists.txt:12`
carrying a commented-out FVM line. Verified — that line is gone, so the claim is stale. The
relocation banner says so rather than silently editing the plan's technical content.

They are untracked in the nonfree repository; staging and committing there is Christian's.

No other file was judged worth deleting. `ngspice_parser_plan.md` comes closest (a proposal never
started, and DR-39 recommends the `.subckt` half be re-scoped or dropped) but it is a real design
document; deferred is the honest home and deletion is Christian's call, not a sweep's.

## The three outside references — fixed, minimally

All three were repointed rather than deleted, and the reasoning differs per site:

- **`src/fem/kernel/cl_FEM_Controller.cpp:4346`** — `// todo/dr92_restart_step_ramp.md` →
  `todo/closed/…`. **Path only.** The surrounding comment states DR-92's mechanism as a
  "current working explanation, not proof", and rewording it is exactly what the **live** row
  **DR-95** asks for — with its suggested resolution explicitly *not applied* pending approval.
  Doing half of an open row's work while fixing a path would have pre-empted it.
- **`src/circuit/doc/README.md`** and **`circuit_usage_guide.md`** — repointed to
  `todo/deferred/ngspice_parser_plan.md`, and each now says the work is **deferred and nothing is
  written**. Neither cites the todo as the *source of a claim* — they are "related / future work"
  pointers, which the inline-facts rule does not forbid — but a reader following them deserves to
  know they lead to a parked proposal rather than an imminent feature.

## DR-95's record hygiene — fixed

The row flagged that `dr92_restart_step_ramp.md` O3 still concluded **"T3 and T4 hold; a ramp is
sufficient"** while T3 and T4 are struck **REFUTED** 150 lines above it.

The two claims are about different objects, which is why the sentence survived. **T3 claimed
in-process healing** — that a sick process cures itself after a few small steps — and measurement
killed it: the sick 50 ms restore rejected step 222 and retried at 25 ms with `it.0 = 9888`,
essentially unchanged from 10851. **What the ladder actually measured was a NEW process *born* at
5 ms** climbing cleanly to 50 ms. The ramp conclusion is sound; the attribution to T3/T4 was not,
and it invited T3 to be quoted as a survivor.

O3 now reads "a ramp **started small** is sufficient to reach full step" and carries an inline note
stating that T3 and T4 remain refuted, what the ladder measured instead, and the measurement that
killed T3. The file was swept for other instances — the only other mention is the header's own
"T3 and T4 REFUTED" line, which is correct.

**DR-95 stays live on its other half:** the `cl_FEM_Controller.hpp:53-56` comment demotion is still
not applied, and that half needs an approved code step.

## Path audit, by-catch of the same fix

DR-95's own row cited `` `dr92_restart_step_ramp.md:240-243` `` — a bare filename with a line-number
suffix, which is exactly the shape that escaped the sweep's earlier path rewrite. A mechanical audit
of both registers for `todo/…md` references that no longer resolve found **six** such stragglers:
one in the live register (DR-92/DR-95's plan) and five in the archive
(`bfm_stale_cache_detection.md`, `example_deck_and_material_db_repair.md`,
`ferro_cut_flux_island.md`). All repointed. **Both registers now resolve every `todo/` path they
cite — checked mechanically, not by eye.**

This is the `file:line` rot the input-contract rule warns about, in its cheapest form: the citation
carried a line range, so a naive filename rewrite missed it.

## Register integrity check, and the one real find

Auditing the registers after the path repairs turned up an ID gap: **DR-30 was in neither file.**
It had not been lost — a parallel session **struck** it (gate discharged on a bulk-superconductor
substitute deck at 1 and 2 ranks, since the 2-rank thermal-coupled corc deck the row named does not
exist) but left the struck row sitting in the **live** file instead of filing it into the archive.

That is the register's own two-actions rule showing up in its mirror image: *closing the status
cell and striking the ID are two separate actions* — and here, striking the ID and **filing the row
into the archive** are also two separate actions. The row was moved to `debt_register_closed.md`.

Both registers now check out: **43 live + 67 archived = 110 IDs, no gaps, no overlap, no struck
rows left in the live file, and every live row tagged.**

The autopilot list drops back to **two** rows (DR-87, DR-107) — `run_gate_batches.md` updated, with
Batch C kept as an explicit "DISCHARGED, do not queue" entry so nobody re-runs it. Its method is
worth keeping: an **unreachable** gate was discharged by substituting a deck that exercises the same
code path, rather than by waiting for the deck the row named. That is the counter-move to the
evidence decay in the blocked table.

## Owed

- `doc/lessons_learned_evidence.md` cites several moved files, but it is a dated evidence record
  and sweep-excluded — left as written.
- **DR-95's other half:** the `cl_FEM_Controller.hpp:53-56` comment demotion. Its record-hygiene
  half is discharged; the comment wording still needs an approved code step.

Nothing was compiled or run. Every change is a file move, a path rewrite, or index prose.
