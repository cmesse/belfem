# Material data files

**Date:** 2026-08-31
**Purpose:** What each shipped data file is, and what a deck can assume about it

`material::data_file()` looks for these files in this order: the path as written,
then `$BELFEM_DATA/material/<path>`, then `$BELFEM_DATA/material/<name>`.
A copy in the run directory therefore always takes precedence over the copy
installed here.

## Superconductor tables

There are six critical-current tables. Each uses an order-2 (triquadratic)
tensor mesh for `log10 jc` and `log10 n`. The mesh axes are temperature,
`log10 |B|`, and the angle from the tape normal. A deck selects a table with
the `file` key in its `jc` block.

| File | Conductor | T range | `t_eff` | Source |
|---|---|---|---|---|
| `bscco-2223.hdf5` | Bi-2223 hermetic tape, AMSC, stainless laminate | 4 – 108 K | 310 um | 10.1063/1.2900382 |
| `sp-ap.hdf5` | SuperPower Advanced Pinning, 12 mm, SCS12050-AP | 3 – 91 K | 1.0 um | 10.6084/m9.figshare.4256624 |
| `sst-1.hdf5` | Shanghai Superconductor High Field Low Temperature, 4 mm | 3 – 91 K | 1.0 um | 10.6084/m9.figshare.5331145 |
| `superox.hdf5` | SuperOx YBCO, 4 mm | 3 – 91 K | 1.0 um | 10.6084/m9.figshare.13708690 |
| `fesc.hdf5` | Fujikura FESC-S12, 12 mm, artificial pinning | 3 – 91 K | 1.0 um | 10.6084/m9.figshare.15095703 |
| `fysc.hdf5` | Fujikura FYSC-S12, 12 mm, no artificial pinning | 3 – 91 K | 1.0 um | 10.6084/m9.figshare.3759321 |

The five REBCO tables span **0.01 – 39.8 T** on 89 x 37 x 181 nodes.
`bscco-2223` spans **0.01 – 25 T** on 105 x 35 x 181 nodes. Every table has a
temperature resolution of **1 K** and a field resolution of 0.1 in `log10 B`.
Above the last measured field, 8 T for the REBCO tapes, the tables are a
log-log tail fit and not measurement; `meta/extrapolation` in each file says
so.

The five REBCO tables come from measured `Ic(T, B, angle)` scans.
`bscco-2223` was digitised from figures. It is mirrored about 90 degrees
because the published data covers only 0 – 90 degrees.

`jc` is the current density in the superconducting **layer**:
`(Ic / width) / t_eff`. Therefore, `jc * t_eff` recovers the measured sheet
current `Ic/w` in A/m. Read `meta/t_eff_m` before comparing tables. All five
REBCO files are Robinson-derived and use the same 1 um convention, so they can
be compared directly. `bscco-2223` uses 310 um.

The 310 um value is neither a convention nor an engineering current density
based on the whole tape. A Bi-2223 conductor has a superconducting layer of
that order. The value excludes the silver coating and stainless laminate; the
finished tape is about 340 um thick. This thickness explains why
`bscco-2223` reports `jc` near 1e8 A/m2 while the REBCO tables report
1e10 to 1e11. Most of the difference comes from the thickness ratio, not from
conductor quality. Compare sheet currents, or compare `jc` only between tables
that use the same `t_eff`.

### What a deck may assume

- `jc` is non-increasing in T and in `|B|` throughout the interpolant, not only
  at the nodes. `n` is non-increasing in `|B|` in every table, and also in T
  for `bscco-2223`; there is no `dn/dT` contract for the REBCO tables. `n`
  stays within the bounds given by `meta/n_floor` and `meta/n_ceiling`. All of
  these hold everywhere, not only at the nodes, because they are imposed on
  the B-spline control net.
- First derivatives are continuous. The nodal values sample a quadratic
  B-spline, so `d jc/dT`, `d jc/dB`, and `d jc/dtheta` do not jump between
  elements. This continuity helps Newton iterations on these tables converge
  cleanly.
- The angle axis covers `[0, pi]`, and the consumer wraps it with period `pi`.
  For `bscco-2223`, the two halves are mirrors, and the angular gradient
  vanishes at 0 and 90 degrees. For the REBCO tables, the asymmetry about
  90 degrees is measured and must not be folded away.
- Outside the T and B windows, the consumer clamps the values and sets the
  corresponding tangents to zero.
- The tables reproduce the measurements used to build them. When evaluated at
  the measured points with the shipped reader, the mean ratio of table values
  to measured `Ic/w` is 0.997 – 1.003 for all five Robinson tapes at 77.5 K
  over 0.01 – 8 T. The 95th percentile of the pointwise error is 5 – 11 %.
  This error is the cost of the B-spline projection. The projection provides
  the continuous derivatives described above but uses about half the degrees
  of freedom per axis, so it smooths structure finer than an element. The bias
  is small, but the scatter is not zero.

### What is in the file besides the table

Every table has a self-documenting `meta` group. Start with
`meta/documentation`. Read `meta/accuracy`, `meta/extrapolation`, and
`meta/continuity` before trusting a value near the edge of the range. Every
table also carries its inputs in `source` and its generator in `python`, so you
can rebuild and check it without BELFEM. For the five REBCO tables `source`
holds the raw measurement rows. For `bscco-2223` it holds the digitised figure
extractions instead, under `source/fig5` and `source/figs7and8.csv`. `python/howtoread.py` reads any of these files with the same
interpolation as the solver.

The same reader is installed beside this directory at
`share/python/database/howtoread.py`:

```python
import sys; sys.path.insert(0, '../python/database')
from howtoread import JcFunctionDatabase
jc = JcFunctionDatabase('sp-ap.hdf5', 'jc')
jc.eval(5.0, 0.0, 77.0)          # |B| [T], angle from the normal [rad], T [K]
```

Licences are recorded in `meta/license` in each file. The five
figshare-derived tables are CC BY 4.0. `bscco-2223` lists its sources in
`meta/reference` instead.

## Other files

| File | Contents |
|---|---|
| `bhdata.hdf5` | B-H curves for the ferromagnetic materials |
