# DR-134 Closed: §20 of the DofManager Guide Rewritten Against the Live Thin-Shell Code

**Date:** 2026-08-28
**Purpose:** Close DR-134 (the stale-doc twin split off DR-85 at its strike the same day) by
rewriting the thin-shell sections of `dof_manager_usage_guide.md` against the current tree, and
fix the two sibling sites carrying the same defect.
**Module:** src/fem/kernel/doc, src/fem/maxwell/doc

---

## What was wrong

DR-134 was filed as documentation debt, and it was: §20 of the DofManager guide walked the reader
through a function, `h_ts_metal`, that no longer exists. It had a second, sharper property. The
walkthrough's listing showed the thin shell's normal field as

```cpp
bn = tCalc->Bm(0) * phi_m + tCalc->Bs(0) * phi_s;
bn *= -0.5 * constant::mu0;
```

which is the **unprojected** master/slave average. That is exactly the defect DR-85 found in the
postprocessor and Christian fixed on 2026-08-16. The fix landed in the code; the prose kept the
bug. A reader routed to §20 by the module README would have met a deleted kernel presented as the
current implementation, carrying a defect the tree no longer has, and would have had no way to
tell.

## What the code actually does

Two layers, and the guide had them collapsed into one:

| Layer | Anchor | Decides |
|---|---|---|
| Kernel selection | `DomainType::ThinShell` in `IWG_Maxwell::link_to_group` | solver algorithm, and whether mu is constant |
| Per-point field/material math | `tIsThinShell` in the `calculator::MaxwellData` ctor | edge vs. node h, the normal term, the rho family, the T source |

`Conductor` and `ThinShell` share one `case`, and the choice inside it never looks at the
material's shape. `mt_maxwell_h.hpp` exports six kernels total: `h_picard`, `h_newton_mu0`,
`h_newton_mu`, `h_ghost`, `h_side_connector`, `h_side_connector_newton`. Everything that makes a
block a shell is a function pointer bound at construction.

The normal field comes from `compute_hn` (a free function in `cl_FEM_Calculator.hpp`), and its
last line is the one the old doc was missing:

```cpp
hn = 0.5 * ( hm + hs ) ;
n  = tCalc->normal( k );
hn = dot( hn, n ) * n ;      // <- this
```

## What changed

`src/fem/kernel/doc/dof_manager_usage_guide.md`:

- **§19.7** the `mFunMKF` example replaced with the real dispatch (one `case` for Conductor and
  ThinShell, branch on `SolverAlgorithm::NewtonRaphson` and `is_constant( MaterialProperty::mu )`).
  Its `cl_IWG_Maxwell.cpp:165-350` line citation replaced by a greppable anchor.
- **§20.3** retitled "Where the Thin-Shell Assembly Lives". Now documents the two-layer split with
  the real `h_picard` body, and states outright that there is no thin-shell kernel to call.
- **§20.4** new, "The Normal Field". Quotes the live `compute_hn` with the projection, and calls
  out the four load-bearing properties: the projection itself, the `k == 0` cache and why it is
  exact rather than approximate, `get_normal_calculator` as the sanctioned bridge, and the fact
  that the normal's **sign** is a modeling input (master to slave, flipped by a minus prefix on
  the `thinshell` sidesets key).
- **§20.5** rewritten around `bn_angle`. The old text's physics was **backwards**: it claimed
  `beta = 0` (B perpendicular to the tape) gives maximum Jc. For a REBCO tape it is the minimum;
  maximum is B in the tape plane. Also records that the angle is deliberately unfolded onto
  `[0, pi]` and must not be re-folded, and separates `bn_angle` (tape normal, HTS) from `bj_angle`
  (Kohler, metal) with the reason `mFundRhodBeta` stays `return_zero` on HTS blocks.
- **§20.6** the kernel-variant table (which listed six deleted functions) replaced by the deck's
  `resistivity type` enum mapped onto the `compute_rho_{powerlaw,piecewise,riva}_ts[_defect]`
  bindings, with the `jc`/`n` vs `file` alternative and the `ec` default from the schema.
- **§20.7** rewritten onto `compute_T_fem` / `compute_T_const`, including the `block_exists()`
  guard against adopting the EmptyBlock's calculator as thermal peer, and the clamp contract (dT
  derivatives forced to zero while clamped; a converged solution at a clamp is a modeling error).
- **§21.6** literature mapping de-staled.
- Revision history gets a **v1.2** entry naming the three substantive corrections.

Two sibling sites, same defect class, fixed in the same pass:

- `src/fem/maxwell/doc/maxwell_usage_guide.md` §8 carried the same deleted `h_ts_metal` listing
  with the same unprojected average, and it routes the reader to §20 of the DofManager guide. It
  now states the no-thin-shell-kernel fact and quotes the projection.
- `src/fem/maxwell/doc/contact_impedance_theory.md:133` named `h_ts_metal` / `h_ts_hts`; one line,
  renamed to the live kernels.

After the pass, `grep -rn 'h_ts_metal\|h_ts_hts' src/ doc/` returns nothing. The remaining
tree-wide hits are all in `devlog/`, `todo/` and `doc/lessons_learned_evidence.md`, which are
records of what was true when written and stay as written.

## By-catch

**`compute_h_ts_node` was dead code, and Christian ruled it out the same session** ("if
compute_h_ts_node is dead code, let's remove it"). See the removal section below.

**A register row was wrong about its own evidence.** DR-134 asserted that
`src/fem/maxwell/matrices/` carries no `mt_maxwell_h.*`, "grep-verified 2026-08-28". The files are
there and are live. What no longer exists is the per-material `h_ts_*` **function family** inside
them. The row's anchor (`h_ts_metal`) was the correct half and it drove the fix, but the file-level
claim was overstated; the correction is recorded in the struck row. Cheap lesson: "the family is
gone" and "the file is gone" are different claims, and a grep for a function name settles only the
first.

## Status

Register: DR-134 struck, `[P]` count 25 → 24, recounted mechanically. Doc-only, so under the
2026-08-28 ruling (comments and docs free, code takes the jury) it landed without an audit round.
The Codex language sweep required for user-facing documentation was run over the rewritten
sections. No source file was touched. Uncommitted.

---

## Follow-on the same session: the dead function removed

Christian's call on the by-catch. Because this touches code that compiles into the binary, it took
the full standing round (plan + audit, code + audit, Codex **and** Grok) rather than the doc-only
exemption above.

### Evidence it was dead

| Check | Result |
|---|---|
| Access | `private` member of `MaxwellData` (`cl_FEM_Calculator.hpp:313` governs) |
| Friends | none declared |
| Textual refs, whole tree incl. `nonfree/`, `tests/`, `todo/` | exactly two: the declaration and the definition |
| `mFunH` bindings | four sites, all other functions |
| Git history, all branches | `-S "mFunH = & MaxwellData::compute_h_ts_node"` returns nothing: never bound in any commit |
| Emitted symbol | `nm -D libbelfem.so.0.9.0` gives 0 hits, while `compute_h_ts_edge`, `compute_h_bulk_node` and `compute_h_bulk_edge` each emit one in `libbelfem_kernel.a` |

It entered in `82a1efc1` (2026-07-21, "more work on matrix collapse") and was never wired.

### What landed

- `cl_FEM_Calculator.hpp`: declaration and inline definition deleted, nothing else.
- `cl_FEM_Calculator.cpp`: the Buffer/`F4` comment no longer names the deleted function (it now
  warns against *any* thin-shell h on a volume block, which is the more general and more useful
  form); the `F4` marker is preserved. A new comment on the `tIsThinShell` branch records why no
  edge/node test is needed there.
- `dof_manager_usage_guide.md` §20.3: the sentence written earlier this session describing the
  function as "declared and defined" was deleted. **Both auditors independently required this**,
  and they were right: leaving it would have re-created, within hours, exactly the stale-name
  class DR-134 had just closed.

### What the jury contributed

**Grok supplied the argument that actually settles "leftover or placeholder".** The one open plan
that might have wanted a nodal thin shell,
`todo/deferred/thinshell_hphi_formulation.md` §3.1, defines its φ layer as nodal `phi` with a
Laplace form. Recovery there is `h = -grad phi`, which is `compute_h_bulk_node`. The deleted body
computed `hn + (-B q)`, adding the air-side normal on top of a gradient that a φ element already
carries in full: it would **double-count** `hn`. So the function was not merely unbound, it was the
wrong formula for its only plausible future use. If insulator-φ layers are ever built they need a
purpose-designed bind, not a resurrection of this one.

**Grok also caught a near-miss in the comment I first wrote.** It cited
`cl_Maxwell_FieldList.cpp` for "ThinShell shares the Conductor dof table (edge_h)". That file has
*two* ThinShell-related tables: `collect_block_dofs` routes ThinShell **blocks** through
`Conductor`, but the `FieldList::ThinShell` member is the **sideset** table and does carry `phi`
(pushed inside the Air loop). A reader grepping `ThinShell` in that file would land on the one that
appears to contradict the comment. Both the comment and the doc paragraph now name
`collect_block_dofs` and warn about the near-miss. Verified on disk before keeping it.

**Codex independently confirmed the load-bearing claim** that no supported path builds a nodal
ThinShell block: the assembly IWG pushes fields only under `Formulation::HPhi`, and the other
enum values (`L2PhiH`, `L2PhiB`, `L2EdgeH`) belong to `IWG_MaxwellPostproc`, which never gets a
`MaxwellData`. I checked the same independently: `MaxwellFactory::create_equation` is private and
only ever called as `create_equation( mFormulation )` with the member defaulted to `HPhi`.

Verdicts: plan **APPROVE WITH CHANGES** from both; code **APPROVE** from both (Codex high, Grok
high ~90%).

### Gates

`g++ -fsyntax-only -std=gnu++17 -Wall` with the module's own `flags.make`, on
`cl_FEM_Calculator.cpp` (baseline green **before** the change, so the result is attributable),
plus `cl_MaxwellPostprocessor.cpp`, `mt_maxwell_h.cpp` and `cl_IWG_Maxwell.cpp`: all exit 0,
stable across re-runs. `nm` shows the function emitted no symbol even before removal, so there is
no codegen and no supported ABI change.

**Reviewed, not verified.** No full library build and no runtime gate was run. `make check-fast`
is owed at the next rebuild, and is the honest remaining gate.

### Tooling and tree hazards hit on the way

Two collisions with concurrent work in this shared checkout, both worth recording:

1. **`ask_grok.sh` was edited at 22:24:39, during the Grok code-audit run**, leaving it
   transiently unparseable ("break: only meaningful in a loop", then a syntax error near `else`).
   The run died at exit 2 and Grok's response was lost. This looked exactly like the known
   read-only-breach signature and was not: `bash -n` passes on the file now. Re-run against a
   private snapshot of the script succeeded. Lesson: when a vendor wrapper fails with a shell
   syntax error, check the wrapper's mtime before blaming the vendor or the quota.
2. **A parallel session is implementing DR-135 in this same working tree** (`cl_FEM_SideSet.{cpp,hpp}`,
   `cl_FEM_Bearing.cpp`, `cl_FEM_DofManager*.hpp`, `cl_MaxwellFactory.cpp`, newest write 22:22:47).
   No overlap with this work (zero `MaxwellData`/`compute_h` hits in their diff), and the audit
   diff was scoped to the two files touched here so nothing leaked into the review. **Nothing was
   committed**, deliberately: a commit from this session would have swept up their in-progress
   work.
