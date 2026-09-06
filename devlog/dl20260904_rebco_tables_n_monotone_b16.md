# REBCO jc/n tables: n non-increasing in B, field axis to 39.8 T

**Date:** 2026-09-04
**Purpose:** regenerate the five REBCO critical-current tables with a monotonicity contract on `n` along `log10 B` and a wider field axis
**Files:** `share/material/{superox,sp-ap,sst-1,fesc,fysc}.hdf5`, `share/material/README.md`

## Question and finding

Christian asked whether `d log10(jc) / d log10(B) <= 0` and the same for `n` would be a
physical contract to add beside the existing `d/dT <= 0` and the angular continuity. Checked
against the shipped tables:

- `jc` already carried it: the B-spline control net is isotonic in T and `log10 B` for `jc`
  (`meta/monotonicity`), and the nodal values have zero positive B-steps in all six files.
- `n` did not. The shipped `n` had no B contract and no T contract; only the clip to
  `[n_floor, n_ceiling]`. In the measured band (T <= 78 K, B <= 8 T) the violations are
  measurement scatter (superox at 20 K: at most 0.04 dex above a non-increasing profile).
  Outside it they are not:

  | table | where | n goes from -> to |
  |---|---|---|
  | sp-ap | 25 K, 6.3 -> 7.9 T, theta 80-87 deg | 13 -> 58 |
  | sst-1 | 40 K, above 8 T | 12.7 -> 17.9 |
  | fysc | 78 K, above 6 T | 4 -> 10 |
  | superox, fesc, fysc | 79-85 K, above 2 T | 4 -> 12-14 |

  The sp-ap case sits inside the measured range and comes from one bad n-value fit at one
  angle scan (the raw rows at 75 and 90 deg give 12 and 22 at 8 T); the 79-85 K cases are the
  log-log tail fitted on floor-limited Ic rows. Both are one-sided upward excursions.

Physics: in the creep picture `n` scales with the pinning energy over kT, which falls with
field; empirically `n` tracks `jc` in these datasets; near the irreversibility field `n -> 1`.
The only mechanism for `n` to rise with B at fixed T and angle is a peak effect, which would
appear in `jc` too, and none of the five sources shows one. Confidence medium-high; rests on
the raw data and standard creep phenomenology, not on a paper checked this session.

## What was done

Work directory `tmp/rebco_b16/` (gitignored, ephemeral). Kit unpacked from the `python` group
of the shipped `superox.hdf5` (identical across the five files).

1. **Reproducibility gate.** `build.py` (Grok's generator) on the unchanged grid reproduces
   the pre-projection tables in `tmp/bscco/grok2/orig/` to round-off (max 2.7e-15 dex on
   `n`, exactly 0 on `jc`, sp-ap and superox). One patch was needed to run it at all: the
   Ic(T) sweep of superox, fesc and fysc contains repeated temperatures (2, 5, 5 rows) and
   scipy's PCHIP refuses a non-strict abscissa; repeated T are now averaged. The round-off
   agreement on superox shows this matches what produced the shipped tables.
2. **Field axis.** `jc/points` and `n/points` set to `89 x 37 x 181` and the tables rebuilt
   from `source/`. Christian asked for `10^1.5 T`; an order-2 axis needs an odd node count,
   so the axis ends at `log10 B = 1.6` (39.81 T), which contains 1.5. Below 1.4 the
   pre-projection `jc` changes only where the Bernstein isotonic passes now pool the two new
   tail nodes: T >= 75 K, `log10 B >= 0.1`, at most 0.23 dex (fesc) and 0.014 dex inside the
   measured band. `n` is unchanged there to round-off.
3. **n gate.** `project_table.py`: for the periodic (REBCO) mode the `n` control net is
   replaced by its **running minimum along `log10 B`** (low to high field) before the clip.
   A running minimum of a control net is a control net, so the table stays a B-spline and
   C1; the clip is order preserving so it cannot undo it. Running minimum rather than the
   L2 isotonic fit used for `jc`, because the violations are one-sided: pooling the 81 K
   trough (n = 4) with its rebound (n = 14) would lift the physical value to about 11, the
   running minimum holds 4 and cuts the rebound. Interior monotonicity of `n` along B is
   now checked element-wise like `jc` (0 bad elements, all five).
4. **Metadata.** `tidy_meta.py` rewrites the storage-format block, the `extrapolation`
   ceiling and the construction list from the grid the file carries (the shipped docs still
   said `45 x 31 x 181`, 2 K, 25.12 T), states the `n` contract in `meta/monotonicity` and
   the documentation, and embeds the patched `build.py` in the kit. Outputs repacked
   (`repack.py`, object-by-object copy) because HDF5 does not reclaim deleted-dataset space:
   19.5 MB -> 10.4 MB per file.

## Gates (executable, all ran)

- `verify_tables.py` through the shipped reader `howtoread.py`, all five: C1 slope mismatch
  <= 2e-12 on every axis, `d jc/dT`, `d jc/dB` non-positive to round-off over 400 000 random
  points, `n` in `[1.02, 90]`, theta seam matched to 1e-14. **All pass.**
- `report_accuracy.py`, `log10 jc` vs the measured rows, 0.01-8 T: unchanged at the fourth
  decimal in every band for all five tables (the widened axis and the `n` gate do not touch
  `jc` in the measured band). fesc 70-85 K worst case 0.848 -> 0.892; that point is in the
  near-Tc collapse where the table is already 0.85 dex off.
- Cost of the `n` gate, new vs shipped, measured band (T <= 78 K, B <= 8 T):

  | table | p95 | max | where the max is |
  |---|---|---|---|
  | superox | 0.013 dex | 0.14 | 75-78 K collapse |
  | sp-ap | 0.016 dex | 0.63 | the 25 K, 7.9 T bad fit (58 -> 9.9) |
  | sst-1 | 0.029 dex | 0.15 | low-B plateau scatter |
  | fesc | 0.023 dex | 0.13 | 75-78 K collapse |
  | fysc | 0.011 dex | 0.45 | 78 K, 6-8 T rebound (10 -> 4) |

  The running minimum sits on the lower envelope of the scatter, so on the low-B plateau it
  biases `n` down by the scatter amplitude (2-7 %). Accepted.

## Not done, and what remains

- **bscco-2223 untouched** by decision: mirrored build, already monotone in T and B for both
  tables, and its field axis stays at 25 T.
- No C++ changed. `JcFunction_Database` reads the grid from the file, so the wider axis
  needs no code; `make check` was not run and has nothing to test here.
- **Jury round run** (Codex gpt-5.6-terra/high, Grok grok-4.6/high, blind;
  `tmp/ai_exchange/review_rebco_table_stage.md`). No P0. Both auditors independently
  confirmed the control-net running minimum (axis, B-spline space, C1, theta periodicity,
  clip order) with two different proofs. Confirmed defects, all in the **embedded metadata**
  and none in the values: `tidy_meta.py` rewrites failed open, so fesc/fysc (different doc
  layout) still said `89 x 35 x 181` / `25.12 T`, sp-ap/sst-1 still said "2 K tensor table"
  and carried the old step-8 wording without the n contract (Codex found the grid text by
  reading the files; Grok found the step-8 path by reading the code); `repack.py` was not in
  the embedded kit or the provenance hash; `build.py`'s docstring said "2 K T-mesh";
  `verify_tables.py` never sampled `dn/dB`. My pre-registration claim that all five docs had
  been inspected was retracted: only superox's storage block had been read.
  **Fixed on approval, one metadata pass, no re-projection:** rewrites now count their
  matches and abort on a miss, cover both layouts and the "Self-describing" purpose variant,
  rewrite an existing L2 step in place; `repack.py` embedded (11 kit files); docstring fixed;
  `verify_tables.py` now gates `d n/dB <= 0` (max +4e-11 against |d n/dB| up to 2e3) and the
  n theta seam (1e-14). Values byte-identical before and after the pass; the installed files
  re-measured: B-spline in-space residual <= 1.2e-13, face jump <= 7.1e-13, all six tables.
  The three single-raiser P2 items were then settled on approval: the sweep dedupe keeps
  the arithmetic mean (it is what reproduces the 2026-08-29 tables to round-off; the
  geometric mean would move them by 0.009 %) and now says so in a comment; the control-net
  periodicity is printed before and after the gate (8.7e-14 / 8.7e-14 for jc, 1.1e-14 /
  1.1e-14 for n on superox); the kit is read from the script's own directory, proven by
  running the re-tidy from `/`. Values byte-identical through both metadata passes;
  installed files re-measured a third time: residual <= 1.2e-13, face jump <= 7.1e-13.
- No `dn/dT` contract for the REBCO tables, **by decision, not omission.** The measured
  `n(T)` for B in the tape plane is not monotone: at 3 T it falls to a minimum near 50-60 K
  and rises again by 0.2-0.26 dex to about 70 K before the near-Tc collapse, in all five
  tapes (superox 10.8 -> 14.4, sp-ap 8.0 -> 14.6, sst-1 7.4 -> 11.4, fesc 9.4 -> 13.8,
  fysc 9.1 -> 14.1). Five tapes, two manufacturers, with and without artificial pinning:
  a measured feature, not scatter. For B parallel to c the profiles are monotone apart from
  a small 40-55 K plateau. A running minimum along T would clamp in-plane `n` at its 55 K
  value up to the collapse: measured-band p95 0.010-0.087 dex, max 0.15-0.28 dex, and the
  large moves are the measured rise. Mechanism reading (low confidence): a crossover
  between pinning regimes for in-plane field. `bscco-2223` keeps its T contract; it was
  digitised from figures and carries no such structure.
- The consumer clamps B at 39.8 T now instead of 25.1 T; decks that relied on the clamp at
  25 T for jc above that field get the tail model instead.

## Second pass, same day: n smoothing (Christian's n.png)

Christian's ParaView slice of superox `n` at 20 K showed jagged contours in angle and
horizontal streaks that `jc` does not have. Measured: in rms second-difference terms `n` is
no rougher than `jc` along either axis (about 3e-3 dex); the look comes from `n` being flat
in B over the low-field plateau, so a 0.01 dex angular wiggle moves a contour far in B. But
the table followed the raw `n` to a median 0.004 dex while the raw scatter between
neighbouring 5-degree angle scans is 0.017 dex, so it was partly tracking one-IV-fit noise;
the streaks are single angle scans off their neighbours, carried to 1 degree by PCHIP.

Two levers tested on the pre-projection stage (before the gates and the projection), on
superox and fysc, scored by roughness, fit to the raw `n` through the reader, ab-peak height
and the 81 K collapse profile:

- **Lever 1, adopted:** `n` gets its own periodic Whittaker in theta, lambda = 300 (jc keeps
  20) with the ab-peak weight dip deepened to w_min = 0.01 so the in-plane peak survives.
  Theta roughness at 20 K halves (superox 3.2e-3 -> 1.7e-3, fysc 3.9e-3 -> 2.3e-3); fit to
  raw `n` unchanged within noise (superox median 0.0043 -> 0.0049, residual rms 0.0092 ->
  0.0100 against 0.017 raw); fysc in-plane peak 0.138 -> 0.123 dex; 81 K collapse
  unchanged (8.0 / 4.2 / 3.0 / 3.0 at 1 / 2.5 / 4 / 10 T). Streaks gone in all five.
- **Lever 2, rejected:** a Whittaker along log10 B for `n` (lambda 5-20). Halves the B
  roughness but blurs the near-Tc collapse, which is measured (superox scans run to 95 K):
  n at 81 K, 4-10 T rises from 3 to 6. Fading it out above the last measured scan does not
  help for that reason. `n` is already no rougher than `jc` along B.

Both levers are in `build.py` behind meta keys `whittaker_lambda_theta_n`,
`whittaker_peak_wmin_n`, `whittaker_lambda_logB_n` (the last unset = off), with environment
overrides for sweeps; `tidy_meta.py` states the n penalty in `meta/smoothing` and in the
construction list. Full rebuild of the five tables: `jc` values byte-identical to the
previous install, `n` changed by mean 0.003 dex (max 0.19-0.35 at the removed streaks);
verifier passes; installed files byte-compared to the outputs and re-measured in a fresh
process: B-spline residual <= 1.2e-13, face jump <= 7.1e-13, all six tables.

Incident, own script: a closing-gate read of fysc reported residual 0.5 and then an HDF5
"bad object header" error. Cause: the install copy overwrote a file that the same process
still held open in h5py, and HDF5 shares open handles by path, so the gate read through a
stale handle. Byte comparison and a fresh process showed the installed file intact. Lesson:
never overwrite an HDF5 file a live handle points at; close, copy, reopen.
