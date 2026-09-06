# Fujikura `sch04` Tables Removed from the Shipped Material Data

**Date:** 2026-08-31
**Purpose:** Record the removal of the two brochure-derived critical-current tables from
`share/material/`, the reasoning, and what was checked before deleting them.

## What happened

Christian asked what distinguished `fesc-sch04.hdf5` / `fysc-sch04.hdf5` from `fesc.hdf5` /
`fysc.hdf5`, then ruled that the `sch04` pair should not ship. Both were deleted (`git rm`) and
`share/material/README.md` was updated to describe six tables instead of eight.

## Why they are different, and why that made them a problem

They are not another processing of the same measurements — they are a different wire and a different
kind of data.

| | `fesc` / `fysc` | `fesc-sch04` / `fysc-sch04` |
|---|---|---|
| Conductor | Fujikura FESC-S12 / FYSC-S12, 12 mm (Robinson FK088 / FK094) | Fujikura FESC-SCH04 / FYSC-SCH04, 4 mm |
| Origin | measured `Ic(T,B,θ)` scans, figshare, CC BY 4.0 | digitised from the Fujikura brochure Rev JUL2024, pp. 12–14 |
| Angular structure | measured | **a model** — the brochure gives `B‖c` only for FESC (replicated over θ) and two directions for FYSC (elliptical interpolant) |
| `t_eff` | 1.0 µm convention | 2.4 µm / 1.9 µm (published layer thickness) |
| Licence | CC BY 4.0 | © 2024 Fujikura Ltd., digitised interpolant, no manufacturer guarantee |

Two independent reasons to drop them, and the copyright one is decisive with the release imminent:

1. **Copyright.** `meta/license` in each file says as much in its own words. The encumbrance is not
   only in the interpolant: `source/points` carries the ~308 digitised marker centres, `python/`
   carries the generator, and `meta/pdf_sha256` fingerprints the brochure PDF. Shipping the file
   redistributes the digitised plot data itself.
2. **They are a model wearing a measurement's clothes.** Every other REBCO table in the directory
   has a measured angular dependence, and `doc/input_file_reference.md:756-757` tells decks that the
   asymmetry about 90° is measured and must not be folded away. The `sch04` tables are flat in θ
   (FESC) or a fitted ellipse (FYSC). A user picking one by name would get an angular response with
   no measurement behind it, and nothing in the deck would say so.

## What was checked before deleting

Tree-wide `grep` for `sch04` outside `.git` and `tmp/`:

- **No** reference from any source file, example deck, test, `doc/input_schema.yaml`, or
  `doc/input_file_reference.md`. The only citations of specific tables in the input docs are
  `sp-ap.hdf5` and `sst-1.hdf5`.
- `CMakeLists.txt:364` installs `share/` as a whole directory, so no CMake edit was needed —
  deleting the files removes them from `make install` automatically.
- `scripts/check_doc_claims.py` makes no claim about the material directory; run anyway, 37/37.
- The generator copies under `tmp/bscco/**` are untracked (`.gitignore:17`) and were left alone;
  their disposal is R10 of the todo plan.

## Files touched

- `share/material/fesc-sch04.hdf5`, `share/material/fysc-sch04.hdf5` — deleted (≈18 MB).
- `share/material/README.md` — eight tables → six; "seven REBCO tables … except for the two `sch04`
  files" → "five REBCO tables"; the `t_eff` comparison paragraph no longer needs to except the
  brochure files; the licence paragraph drops its `sch04` sentence; date bumped.
- `todo/sch04_history_purge.md` — new.

## The part that is not done

`git rm` does not remove the blobs from history. They are reachable from `568a0fff` and `b99ba890`
(both 2026-08-29) on `main`, `devel` and `claude`, and all three branches on `origin`
(belfem.lbl.gov GitLab) still carry them. The GitHub `backup` mirror does **not**, as of the last
fetch. Christian decided to purge history; the plan is `todo/sch04_history_purge.md` and **none of it
has been run** — it rewrites shared history and needs collaborator coordination first. Its O1 asks
the question that decides whether the rewrite is needed at all: if the public release is cut as a
tarball or a squashed repository, the blobs never reach the public.

Nothing was built or run this session beyond `check_doc_claims.py`, so this is *reviewed*, not
verified. The README revision is a set of deletions and count corrections to an already-swept
document; a Codex language sweep over it is owed but not blocking.
