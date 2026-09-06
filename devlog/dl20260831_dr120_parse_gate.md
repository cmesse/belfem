# Devlog 2026-08-31 — DR-120 negative-deck parse gate RUN

**Date:** 2026-08-31
**Purpose:** Execute the one run gate DR-120 still owed, refresh the row's stale
citations, and record what the gate does and does not settle
**Module:** fem/kernel

## What this session was, and what it was not

The session opened on the premise that DR-120 "would be an easy fix". It is not a
fix at all: **DR-120's code landed 2026-08-28 and has been in the tree ever since.**
Both `build/bin/belfem` and `cmake-build-debug/bin/belfem` (Aug 30) were confirmed to
carry the error strings before anything was run, so the Aug-29 suite ran against a
binary that already had the checks compiled in.

What the row actually owed was its **run gate** — a deck with a negative
`max iterations` must abort at setup naming the key. That gate has now RUN.

## The gate

`Controller::set_params` is called from `MaxwellFactory::create_controller()`
(`cl_MaxwellFactory.cpp:1047`), i.e. *after* mesh load and kernel construction. The
gate therefore needs a mesh a current deck accepts, and no `.msh` in the tree
qualified: `cmake-build-debug/tapestack3d.msh` is a Jul-2 two-volume mesh, while
`examples/3D_tapestack/input.conf` now wants blocks `1` and `2:8`. The deck chosen was
`examples/2D_Tapestack` — 2D, self-contained, no user-material plugin — meshed by
gmsh 4.13.1 into the scratchpad (62865 nodes, 130495 elements). Release binary,
serial. Nothing was written into the repository.

| case | deck change | result |
|---|---|---|
| control | none | accepted; 12 completed timesteps |
| A | magnetic `max iterations : -1` | abort `:4211`, "nonlinear: max iterations must be positive ( is -1 )" |
| B | magnetic `max 30` / `min 200` | abort `:4226`, cross-check names both values |
| C | `stall window : -1` | abort `:4256`, "stall window must not be negative ( is -1 )" |
| D | magnetic `watchdog window : -1` | abort `:4282`, "( is -1, 0 disables )" |
| E | thermal `max iterations : -1` | abort `:4412` — **the row's literal gate** |
| F | segregated `coupling factor : 0` | abort `:4483`, "must be positive ( is 0 )" |

Every abort names `belfem::fem::Controller::set_params()` and occurs after the mesh
read, cuts, edges and connectivities — so these are genuine setup aborts, not a
parse-time short circuit that never reached the controller.

## What the gate discriminates, and what it does not

**It discriminates on the live code.** The control deck ships `watchdog window : 0`,
the documented disable sentinel, and it passes; `-1` on the same key aborts. So the
checks separate legal from illegal values rather than rejecting the key outright, and
they do not fire on a legal deck (12 timesteps).

**It does not re-measure the pre-fix behaviour.** That a negative used to wrap into a
huge `uint` budget is a static argument — `Section::get_int` is `round( get_real )`
with no range check (`cl_Input_Section.cpp:349-352`), and the members are `uint`
(`cl_FEM_Controller.hpp:77-78, 100, 109-110, 121-122`). Reproducing the silent wrap
needs a rebuild with the checks reverted, which was judged disproportionate for a P3.
The row and this entry both say so rather than letting "gate ran green" imply more.

## Citations were stale again, by ~210 lines

The row cited the magnetic block at `:3997-4018` and the thermal at `:4198-4219`; the
actual sites are `:4208-4227` and `:4408-4428`. The 2026-08-30 currentness sweep had
already re-anchored two DR-120 citations, and they had drifted again inside a day —
`cl_FEM_Controller.cpp` is under active uncommitted modification. This is the exact
failure mode `CLAUDE.md` names for the input contract: *a line number read from a
working copy is already wrong by the time it is committed.* The register currently
anchors by line number; `doc/input_schema.yaml` anchors by searchable token for
precisely this reason. Whether the register should adopt the same discipline is a
question for Christian, not a change this session made.

## By-catch (observation only, no row filed)

`coupling factor` is read with `get_int` — i.e. `round( get_real )` — into a `real`
member `mCouplingFactor`. A deck writing `coupling factor : 2.5` therefore gets 3
with no notice. Harmless for a step count, and filing a row for it would be the
register inflation the extend-don't-branch policy exists to prevent.

## Status

**DR-120 STRUCK and archived 2026-08-31 on Christian's ruling,** in the same session
its gate ran. The row was re-anchored and its status cell rewritten first, then lifted
out of `todo/debt_register.md` into `todo/debt_register_closed.md` with ID and
description struck and the status column left intact — that column is where the
closure evidence lives. The `[P]` count recounted 9 → 8 (8 `[P]` + 3 `[W]` = 11 live
rows, 138 closed); `check_doc_claims.py` 37/37 both before and after.

Two things were deliberately carried into the archived row rather than dropped with
the strike: the `coupling factor` rounding by-catch, and the citation-drift
observation. **Struck is not verified** — and this row is the mild case of that rule:
the post-fix behaviour is verified by execution, the pre-fix claim is reasoned.

No code changed this session, so no audit round was owed.

## By-catch: two malformed rows in the closed register

Checking my own archived row's cell count turned up two rows carrying **seven cells in
a six-column table** — content shifts one column right, so the "blocking-1.0?" column
renders the wrong text. Both predate this session (confirmed against a pre-edit
backup; they only moved by the one line my insertion added) and both came from other
sessions' archival work. Fixed on Christian's instruction, delimiters only, no wording
changed:

- **`~~DR-111~~`** — a literal shell pipe inside a code span,
  `` `strings <binary> | grep -c 'Negative Jacobian determinant'` ``. GFM splits a
  table row on `|` **even inside backticks**; escaped to `\|`, which renders back as
  `|` so the command stays copy-pastable. This is the trap to remember: a code span
  does not protect a pipe in a table.
- **`~~DR-129~~`** — a stray ` | ` between the STRUCK closure paragraph and the
  original "Fix shape:" status text the archiving session prepended to. Replaced with
  an em-dash, the register's own idiom for that transition.

Both registers now pass a mechanical cell-count check on every row, and the closed
register's row inventory is unchanged at 138.
