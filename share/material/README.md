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
- Outside the T window and above the top of the B window, the consumer
  clamps the values and sets the corresponding tangents to zero.
- Below `Bmin` = 10 mT, the tables have no data. BELFEM continues `n` and
  `jc` for the five REBCO tables down to B = 0. The value and its field slope
  are continuous at `Bmin`, and the value does not increase with `|B|`; `jc`
  also does not increase with T. `n` follows a parabola with zero slope at
  B = 0, so `n(0)` is slightly above `n(Bmin)`. `jc` increases toward the
  self-field level measured at B = 0, which the file stores in `source`. If
  that level is too low for a monotone join, the reader raises it, most
  strongly below the lowest measured temperature and near Tc. `bscco-2223`
  has no self-field measurement, so its `jc` remains clamped below `Bmin`.
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
interpolation as the solver for B >= `Bmin`. Below `Bmin`, this reader clamps
the values, while the solver continues them as described above.

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

### Why the tables are monotone and smooth

The shape contracts listed above are modelling choices made for the solver, not
laws of physics: `jc` is non-increasing in T and `|B|`, `n` is non-increasing
in `|B|` and bounded, and the interpolant has continuous first derivatives
across element faces.

**Monotonicity suppresses noise-driven features.** The two tables are treated
differently:

- `jc` gets an isotonic fit in T and `log10 B`. Inside the measured range,
  where the raw scans are close to monotone, the fit moves the values by less
  than the measurement scatter. It therefore removes noise without changing
  the measured trend.
- `n` gets a running minimum along `log10 B`, from low to high field. This
  removes one-sided upward excursions, such as a bad fit in one angle scan or
  the log-log tail fitted to rows where `Ic` is at the measurement floor.
  These excursions can be much larger than the scatter: in `sp-ap`, one bad
  fit sent `n` from 13 to 58. `n` has no T contract for the REBCO tables
  because the measured in-plane `n(T)` rises again between about 55 and 70 K.

The shape of `jc(|B|)` depends on the pinning landscape. With columnar
BaZrO3 defects along the c-axis, `jc` falls much more slowly when the field is
parallel to c. The defect density that works best in a given field range
follows the matching field, where the vortex density equals the defect density
(Wee et al. 2013, DOI 10.1038/srep02310).

Do not assume that `jc` falls steeply with field. The tables also cannot
represent a conductor whose `jc` rises over part of the field range. None of
the five REBCO tapes behind the shipped tables shows such a rise. If the source
data contained one, two stages of the build would flatten it: a monotone
Bernstein projection before the B-spline stage, and the constraint on the
B-spline control net. More sophisticated regularizations exist, including
machine-learning approaches. BELFEM uses the simplest one that gives the
solver what it needs.

**Smoothness helps Newton converge.** The tangent of the power-law
resistivity contains `d jc/dB` and `d jc/dT` (`powerlaws.hpp`). If the table
values are continuous but their derivatives jump between elements, the
tangent jumps too. The table resolution is 1 K, so during a quench the
temperature crosses many table elements within a few time steps. In our
experience, derivative jumps of this kind are what stall the Newton iteration
in that phase. This has not been isolated in a reproducer.

The B-spline projection removes the jumps. Well-behaved material functions
and the projection together are a significant part of the solver's
robustness on HTS problems. Keep both when you change a table.

The C1 property belongs to the table interpolant. The functions BELFEM serves
are not C1 everywhere: below `Bmin`, the `jc` continuation can have
discontinuous T and angle tangents, and the `n` floor is C0 only
(`resistivity_laws.md`).

### Using an older or modified table

When BELFEM opens a `jc` or `n` table, `JcFunctionDatabase::open_checked`
(`cl_JcFunction_Database.hpp`) checks the number of axes, the angle coverage,
and one property of the data.

On the lowest-field face, `d log10 value / d log10 B` must be certified
non-positive by a Bernstein bound over the whole face. A table that fails is
refused with

    database jc in sst-1.hdf5 cannot be certified not to rise with |B| at its
    lowest field ( bound on d log10 value / d log10 B: 0.0198602 ), so it
    cannot be continued monotonically below Bmin

For `n`, the reader also requires `n >= 1.01` on the same face. The bound is
conservative: a face that falls everywhere but is strongly concave can still
fail it.

The check protects the continuation below `Bmin`, which cannot be
non-increasing if the table already rises at `Bmin`. It applies to every
table, including a `jc` table without a measured self-field curve, such as
`bscco-2223`, which is clamped below `Bmin` rather than continued.

One cause of the refusal is a table built before the current versions. This
can be an older copy of a shipped file left in the run directory, which takes
precedence over the installed copy (see the top of this page), or a
hand-edited file. Replace it with the current file before doing anything else.

The check covers only the lowest-field face. A rise elsewhere in the table,
or derivative jumps from values edited without re-projection, load without
complaint and may show up later as a Newton iteration that does not converge.

If you edit a table, re-run the B-spline projection and verify the result. The
tools are in the file's `python` group:

```python
import h5py, pathlib
with h5py.File('my-table.hdf5', 'r') as f:
    for k, v in f['python'].items():
        pathlib.Path(k).write_bytes(v[()])
```

```sh
python3 project_table.py my-table.hdf5 -o my-table-c1.hdf5            # REBCO
python3 project_table.py my-table.hdf5 -o my-table-c1.hdf5 --mirror   # bscco-2223
python3 verify_tables.py my-table-c1.hdf5     # C1, monotonicity, n bounds, theta seam
```

`--mirror` keeps the Bi-2223 table symmetric about 90 degrees. Without it,
the projection treats the angle as periodic and breaks that symmetry. Note
that `howtoread.py` in the same group clamps values below `Bmin`, so it does
not show the continuation the solver uses.

The projection re-imposes `jc` non-increasing in T and `|B|`, and `n`
non-increasing in `|B|`, so it flattens any such rise you intended to keep.
It does not touch `n(T)`. `project_table.py` has no option to switch these
constraints off.

Keeping a genuine `jc` rise therefore means editing the script, and also the
monotone pass in `build.py` if you rebuild from `source`. What else is needed
depends on where the rise is:

- **A rise inside the field range** passes the load check. Only the scripts
  need to change.
- **A rise at the lowest field** is refused at load. Keeping it also requires
  disabling the check in `open_checked`. Below `Bmin`, the continuation is
  then no longer guaranteed to be non-increasing.

Neither configuration is tested. Treat both as a last resort, and prefer the
current shipped table unless the rise is the physics you are studying.

## Other files

| File | Contents |
|---|---|
| `bhdata.hdf5` | B-H curves for the ferromagnetic materials |
