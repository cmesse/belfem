# Building a jc / n Table for an HTS Conductor

**Date:** 2026-10-26
**Purpose:** How to turn measured or fitted critical-current data for an HTS tape or cable into an
HDF5 file that `JcFunctionDatabase` accepts, and how to make sure its B-spline mapping holds: the
table is C¹, monotone and periodic, and stays that way
**Module:** `src/physics/materials` (reader), `src/physics/database` (interpolant and projector)

---

## 1. What you are building

A deck points an HTS material at a table with the `file` key:

```
materials
{
    tape
    {
        builtin          : ybco ;
        file             : my-cable.hdf5 ;
        resistivity type : riva ;
    }
}
```

`MaterialFactory` then opens the file twice: once for the `jc` group and once for the `n` group
(`cl_MaterialFactory.cpp`, the `MaterialType::HTS` branch). Both are read by `JcFunctionDatabase`
(`cl_JcFunction_Database.hpp`). With `file` set, the deck's `jc`, `n` and `ec` keys are not read, so
the electric-field criterion is fixed at the default `ec = 1e-4 V/m` (1 µV/cm). Section 7 explains
what that means for your data.

Each group is a `Database` table (`src/physics/database/doc/database_usage_guide.md`): a structured
tensor grid with **order-2 (HEX27) Lagrange interpolation**. The table is only as good as its node
values. If those values are raw samples, the interpolant is C⁰: values are continuous, but
`d jc/dT`, `d jc/dB` and `d jc/dθ` jump at every element face. Those derivatives go into the Newton
tangent of the power-law resistivity (`powerlaws.hpp`). A jumping tangent can stall the Newton
iteration. This is most visible during a quench, when T crosses many 1 K elements within a few time
steps.

The node values must therefore be **samples of one quadratic B-spline** on the same mesh. This makes
neighbouring HEX27 patches reproduce the same polynomial pieces and match their first derivatives
across every face. This guide calls that property *the B-spline mapping*. Most of the work below is
about establishing it and not breaking it afterwards.

## 2. The file format

### 2.1 Required datasets

Two groups, `jc` and `n`, on **identical grids**. Each holds five datasets:

| Dataset | Type | Content |
|---|---|---|
| `order` | `uint32` scalar | `2` |
| `origin` | `float64[3]` | first node on each axis: `{ Tmin, log10 Bmin, 0 }` |
| `points` | `uint32[3]` | number of **nodes** per axis (not elements) |
| `step` | `float64[3]` | **node** spacing per axis |
| `values` | `float64[nT·nB·nA]` | node values, flattened in C order: axis 0 (T) slowest, axis 2 (θ) fastest |

The axes and the stored quantity are fixed by the reader:

| Axis | Quantity | Scale |
|---|---|---|
| 0 | temperature T [K] | linear |
| 1 | field magnitude | **log₁₀**(B / 1 T) |
| 2 | angle θ [rad] between B and the **tape normal** (`angleNxB`); θ = 0 is B ∥ c | linear |
| values | `jc` [A/m²] and `n` [–] | **log₁₀** (not ln: the metal ρ tables use ln, these do not) |

The C-order layout matches `TensorMeshConfig::node_index`, which is `k + nA·( j + nB·i )`. In numpy,
build the array with shape `( nT, nB, nA )` and `indexing="ij"`, then write `v.ravel( order="C" )`.
A Fortran-order write loads without error and silently transposes the table.

### 2.2 Grid rules for the B-spline mapping

The order-2 B-spline layer (`bspline.py`, §4) uses `Ne + 2` control points for `Ne` elements.
Element `e` uses controls `e, e+1, e+2`. The resulting grid rules are:

| Rule | Why |
|---|---|
| **Odd** node count on every axis (`points = 2·Ne + 1`) | an order-2 element spans three nodes; an even count leaves a half element. `bspline.node_matrix` asserts this |
| θ covers exactly `[ 0, π ]`: `origin[2] = 0`, `step[2] · ( nA − 1 ) = π` | `open_checked` refuses a table that does not reach both 0 and π (within 1e-6). `make_extension.py` asserts the period |
| `nA − 1` even (with `nA` odd this holds automatically) | the θ period has to be a whole number of elements, so the periodic extension (§4.2) tiles cleanly |
| last θ column equals the first: `v[:, :, -1] == v[:, :, 0]` | the reader wraps θ with period π. A mismatched seam is a jump the projection cannot repair |
| `origin[1]` is the lowest field you trust, typically `−2` (10 mT) | below `Bmin` the table holds no data, and the reader continues it analytically (§6) |

The shipped REBCO tables use 89 × 37 × 181 nodes (T 3–91 K in 1 K steps, log₁₀ B from −2 to 1.6 in
0.1 steps, θ in 1° steps). That is a good default. A coarser grid works too, but the projection
smooths out structure finer than about one element (§5).

### 2.3 The `meta` group

The solver uses only the self-field calibration in `meta` (§6). The projection and verification
scripts need two entries. The remaining entries document the file for future users:

| Entry | Needed by | Content |
|---|---|---|
| `meta/n_floor`, `meta/n_ceiling` | `project_table.py`, `verify_tables.py` | bounds imposed on `n`. The shipped tables use 1.02 and 90. The reader requires `n ≥ 1.01` on the lowest-field face |
| `meta/t_eff_m` | the self-field calibration, and every human comparing tables | thickness [m] that turns `Ic/w` into `jc` (§7) |
| `meta/Icw_77p5K_sf_A_per_m` | self-field calibration | measured `Ic/w` [A/m] at 77.5 K and B = 0 (§6) |
| `meta/angle_folding` | `verify_tables.py` | present only for a table mirrored about π/2 (§4.2) |
| `meta/ec_criterion`, `meta/reference`, `meta/doi`, `meta/license`, `meta/documentation` | people | the Ic criterion, the data source and its licence, and a free-text description |

Copy the `meta` layout of a shipped table (`share/material/sp-ap.hdf5`) rather than inventing one.
`share/material/README.md` explains what each shipped entry means.

## 3. Step 1: fill the raw node values

Build the raw values from the data available to you. Typical sources are:

- **A parametric fit** (Kim, Hilton, a Fleiter-type fit for REBCO, a vendor formula): evaluate it
  directly at the nodes.
- **Measured scans** `Ic( T, B, θ )`: interpolate onto the nodes in log space. PCHIP keeps
  interpolants free of overshoot. Extrapolate beyond the measured T and B ranges with a
  physically motivated law, and record in `meta/extrapolation` where measurement ends and model
  begins. The `build.py` shipped inside every REBCO table is a worked example. It is written for
  the HTSDB/Robinson export layout, so treat it as a template, not a general tool.

Before projection, apply these checks to the raw values:

1. **No NaN or Inf anywhere.** Add a small floor before taking the logarithm
   (the shipped tables add `1e4 A/m²` to `jc`) so that the region near Tc stays finite.
2. **`jc` non-increasing in T and in log₁₀ B**, `n` non-increasing in log₁₀ B and inside
   `[ n_floor, n_ceiling ]`. The projection re-imposes these on the control net (§4.3), but it
   moves values less, and the result stays closer to your data, when the input already roughly
   satisfies them. An isotonic (PAVA) pass along each axis is the simplest fix for measurement
   noise.
3. **Close the θ seam**: `v[:, :, -1] = v[:, :, 0]`.
4. **Smooth noisy θ scans if needed.** Do not fold an asymmetric angular dependence about 90°:
   if the asymmetry was measured, it is physics.

A minimal writer, which creates the skeleton every later step expects:

```python
import h5py, numpy as np

T0, dT, nT = 4.0, 1.0, 87                # 4 .. 90 K, odd node count
L0, dL, nB = -2.0, 0.1, 37               # B = 0.01 .. 39.8 T, odd
nA = 181; dA = np.pi / (nA - 1)          # theta = 0 .. pi in 1 deg steps

T  = T0 + dT * np.arange(nT)
lB = L0 + dL * np.arange(nB)
th = dA * np.arange(nA)
Tg, Lg, Ag = np.meshgrid(T, lB, th, indexing="ij")   # shape ( nT, nB, nA )

log10_jc = ...   # your model or data at ( Tg, 10**Lg, Ag ), log10 of A/m^2
log10_n  = ...

with h5py.File("my-cable.hdf5", "w") as f:
    for g, v in (("jc", log10_jc), ("n", log10_n)):
        v[:, :, -1] = v[:, :, 0]                     # close the theta seam
        G = f.create_group(g)
        G["order"]  = np.uint32(2)
        G["origin"] = np.array([T0, L0, 0.0])
        G["points"] = np.array([nT, nB, nA], np.uint32)
        G["step"]   = np.array([dT, dL, dA])
        G["values"] = np.asarray(v, np.float64).ravel(order="C")
    m = f.create_group("meta")
    m["n_floor"] = 1.02 ; m["n_ceiling"] = 90.0 ; m["t_eff_m"] = 1.0e-6
```

At this point the file loads, but it is only C⁰. `verify_tables.py` (§5) normally reports failures
on the three C1 axes and on the θ seam slope.

## 4. Step 2: the B-spline projection

### 4.1 Get the tools

Every shipped table carries its generator pipeline in the group `python`. Unpack it next to your
file:

```python
import h5py, pathlib
with h5py.File("share/material/sp-ap.hdf5", "r") as f:
    for k, v in f["python"].items():
        pathlib.Path(k).write_bytes(v[()])
```

| Script | Role |
|---|---|
| `bspline.py` | numpy reference of the B-spline layer: collocation matrix, L2 projection, control-net round trip, isotonic regression, face-jump diagnostics |
| `make_extension.py` | writes the θ-extended copy that BELFEM's projector needs (§4.2) |
| `project_table.py` | projects, folds the extension back, re-imposes the contracts on the control net, writes the C¹ table |
| `verify_tables.py` | checks the finished file through the reference reader (§5) |
| `helpers.py` | element-local Bernstein tools used by `build.py` and `project_table.py` |
| `howtoread.py` | the reference reader, also installed as `share/python/database/howtoread.py` |

`verify_tables.py` looks for `howtoread.py` beside itself first, then under
`$BELFEM_DATA/python/database`.

### 4.2 Why θ is extended, and which extension to use

`database::Projector` imposes **no boundary condition**. The projected spline is therefore free at θ
= 0 and θ = π, and its slope need not match across the seam. The condition has to come from the
domain: θ is extended before projection, and only the middle period is kept:

| Table | Extension | What it gives at θ = 0, π |
|---|---|---|
| π-periodic, asymmetry about 90° measured (all REBCO tapes) | three **translated** copies on `[ −π, 2π ]` (default) | value **and** slope match across the seam |
| only `[ 0, π/2 ]` measured, mirrored (`bscco-2223`) | **even reflection** at both planes (`--mirror`) | zero angular gradient |

These two cases are not interchangeable. `--mirror` on an asymmetric table destroys the asymmetry.
The default on a mirrored table breaks the symmetry. A mirrored table must also carry
`meta/angle_folding` so that `verify_tables.py` holds it to the mirror contracts.

### 4.3 Run the projection

Two routes produce the same table to round-off.

**numpy only** (no BELFEM build needed):

```sh
python3 project_table.py my-cable.hdf5 -o my-cable-c1.hdf5            # periodic
python3 project_table.py my-cable.hdf5 -o my-cable-c1.hdf5 --mirror   # folded about pi/2
```

**BELFEM's own projector, with numpy as the cross-check.** `mattest` is a development executable
built with the materials module (`src/physics/materials/mattest.cpp`):

```sh
python3 make_extension.py my-cable.hdf5 my-cable-ext.hdf5   # [--mirror]
mattest my-cable-ext.hdf5                                   # writes values_bspline into jc/ and n/
python3 project_table.py my-cable.hdf5 my-cable-ext.hdf5 -o my-cable-c1.hdf5   # [--mirror]
rm my-cable-ext.hdf5
```

With the second argument, `project_table.py` compares BELFEM's projected values with the numpy
reference. It aborts with `FATAL: the two projections disagree` when they differ by more than 1e-9.

What `project_table.py` does, in order:

1. Extends θ (§4.2) and L2-projects onto the quadratic B-spline space. The projection uses the
   same element matrices as `database::Projector`: `A = Tᵀ M T`, `Y = Tᵀ M f`
   (`cl_DatabaseProjector.cpp`).
2. Keeps the middle period and reports the *end effect*, which should have decayed to round-off,
   and the *in-space residual*, which should be near 1e-13 and confirms that the values lie in
   the spline space.
3. Converts to the control net and **re-imposes the contracts there, not on the node values**:
   - `jc`: isotonic (non-increasing) in T and log₁₀ B; θ is left free.
   - `n`: running minimum along log₁₀ B (low to high field), then clipped to
     `[ n_floor, n_ceiling ]`. It uses a running minimum rather than an isotonic fit because the
     violations in `n` are one-sided upward excursions from bad IV fits.
   - with `--mirror`: both tables are made monotone in T and log₁₀ B, and θ is driven
     non-decreasing on `[ 0, π/2 ]` with the mirror imposed on the control net.
4. Evaluates the node values from the constrained control net, checks interior monotonicity on
   every element, and writes `jc/values` and `n/values`.

Step 3 preserves the mapping because the quadratic B-spline basis is non-negative and sums to one. A
control net that is non-increasing along an axis gives an interpolant that is non-increasing along
that axis **everywhere**, and a bounded net gives a bounded interpolant. Both changes leave the
field a B-spline. The same contracts imposed on the node values would move the values off the spline
and bring back the C⁰ kinks.

## 5. Step 3: verify

```sh
python3 verify_tables.py my-cable-c1.hdf5
```

The script reads the file through `howtoread.py`, the consumer path, and checks:

| Gate | Pass criterion |
|---|---|
| C1, per axis | one-sided `d/dx` limits agree at every interior face to 1e-8 of the reference slope |
| monotone | `d jc/dT ≤ 0`, `d jc/dB ≤ 0` at 400 000 random interior points; periodic tables also `d n/dB ≤ 0` |
| n bounds | `n` inside `[ n_floor, n_ceiling ]` |
| θ seam (periodic) | value and slope match across θ = 0 / π, for `jc` and `n` |
| θ mirror (folded) | `jc( θ ) = jc( π − θ )`, zero edge slope, `jc` and `n` rising on `[ 0, π/2 ]`, `d n/dT ≤ 0` |

It ends with `all tables pass every gate` or `FAILED: { … }` and a non-zero exit code. A table that
has not been projected fails the C1 gates. That is a quick way to tell whether a file you were
handed has been through §4.

Also compare the projected table with your data. The projection costs accuracy. On the shipped REBCO
tables the mean ratio of table to measurement is 0.997–1.003 at 77.5 K, with a pointwise
95th-percentile error of 5–11 %, because the spline has about half the degrees of freedom per axis
and smooths structure finer than one element (`share/material/README.md`). If that is too coarse for
your conductor, refine the grid. Do not skip the projection. Record the measured error in
`meta/accuracy`.

### What BELFEM checks at load

`JcFunctionDatabase::open_checked` refuses a file, on every rank together, when:

- the table does not have exactly 3 axes;
- θ does not cover `[ 0, π ]` (a folded `[ 0, π/2 ]` or `[ −π/2, π/2 ]` export);
- the table cannot be certified **not to rise with |B| on its lowest-field face**. This is a
  Bernstein bound on `∂ log₁₀ value / ∂ log₁₀ B` over the whole face `log₁₀ B = origin[1]`, with
  tolerance 1e-10. The bound is conservative, so a face that falls everywhere but is strongly
  concave can still fail;
- for `n`: `n < 1.01` anywhere on that face.

Loading does not check C¹ continuity or interior monotonicity. A table that skipped §4, or whose
values were edited afterwards, loads without complaint. The problem shows up later as a Newton
iteration that does not converge. **Run `verify_tables.py` on every file before it goes into a
deck.**

## 6. Below Bmin: the optional self-field anchor

The log₁₀ B axis cannot reach B = 0. Below `Bmin` the reader continues each table C¹ and
non-increasing in |B| from the values and slopes on the lowest-field face
(`cl_JcFunction_Database.hpp`, "below Bmin"):

- `n` always: a parabola with zero slope at B = 0.
- `jc` only if the file carries a measured self-field level: a monotone cubic Hermite from that
  anchor to the table. Without the level, `jc` is **clamped** at its `Bmin` value.

The `jc` anchor is read from three places, and all three must be present:

| Item | Content |
|---|---|
| `source/columns` | string array naming the columns; it must contain `temperature`, `field` and `icw` |
| `source/points` | `( N × ncols )` float64 rows; rows with `field == 0` exactly are the self-field rows. NaN and non-positive values are skipped |
| `meta/Icw_77p5K_sf_A_per_m` and `meta/t_eff_m` | the curve is **calibrated**: scaled so that its 77.5 K value equals `Icw_77p5K_sf / t_eff`. The unit of the `icw` column therefore does not matter. The calibration holds the raw curve at its end value outside the measured temperatures, so data that stops short of 77.5 K is calibrated against the wrong temperature. In evaluation, the anchor is instead continued linearly from its end knots |

Rows are binned by temperature (0.25 K tolerance) and averaged over angle. At least two temperature
bins are needed, otherwise the clamp is kept. If any of the items is missing, the reader falls back
to the clamp without an error.

## 7. Cable-specific points

The table describes the conductor **at the scale resolved by the mesh**. Before building a cable
table, decide which case applies.

**Tape-resolved cable** (Rutherford, CORC or Roebel geometry with the tapes as thin shells, as in
`examples/rutherford`). The table is a **tape** table. Use the tape data, not the cable's `Ic`, and
keep θ measured from the **local tape normal**, which is what the solver supplies at every
quadrature point. The same tape table serves every tape in the cable.

**Thickness must match.** `jc` is a current density in the superconducting layer,
`( Ic / w ) / t_eff`. The solver multiplies it by the layer thickness the deck assigns:

```
layers : tape
{
    bscco : 310 mum ;     // equals meta/t_eff_m of bscco-2223.hdf5
}
```

If the `layers` thickness differs from `t_eff_m`, every model `Ic` is scaled by that ratio, and
BELFEM does not warn you. Build the table with the `t_eff` you mesh, or mesh the `t_eff` you built.

**The Ic criterion is fixed at 1 µV/cm.** With `file` set, `ec` is not read and is 1e-4 V/m. If your
`Ic` used another criterion `E₁`, convert before building the table:
`Ic( E₀ ) = Ic( E₁ ) · ( E₀ / E₁ )^( 1 / n )`. Store the criterion in `meta/ec_criterion`.

**Do not anchor on a cable-level self-field measurement** (medium confidence). The self-field `Ic`
of an assembled cable is set largely by the cable's own field distribution, which the FEM model
computes itself. Using it as the B = 0 anchor of a tape table would count that field twice. Leave
out `source` and accept the clamp below `Bmin`, or anchor on a short-sample measurement of the
single tape.

**A homogenized cable table**, with one `jc` for the whole cable cross section, is not covered by
this guide. The θ axis would no longer be an angle to a tape normal, and the anisotropy of the tapes
inside the cable does not survive averaging in a way this reader models. Treat it as untested.

## 8. Pitfalls

| Pitfall | Symptom | Avoid by |
|---|---|---|
| Raw samples stored without projection | loads fine; Newton stalls in transients, especially quenches | §4, then `verify_tables.py` |
| Contracts (monotone, bounds, clipping) applied to node values **after** projection | C1 gate fails again; derivative kinks come back | impose them on the control net, which `project_table.py` does; never hand-edit `values` |
| Even node count on an axis | `assert` in `bspline.node_matrix` | `points = 2 · Ne + 1` |
| θ range `[ 0, π/2 ]` or `[ −π/2, π/2 ]` | refused at load: "unfolded bn_angle contract" | export the full `[ 0, π ]`; mirror if only half was measured |
| Values in ln instead of log₁₀ | `jc` off by orders of magnitude | log₁₀ for `jc` and `n`; ln is only for metal ρ |
| Fortran-order `values` | loads, silently transposed | `ravel( order="C" )` on a `( nT, nB, nA )` array |
| Wrong extension mode (`--mirror` vs. periodic) | broken asymmetry, or a θ seam failure | §4.2; set `meta/angle_folding` for mirrored tables |
| Table rises with \|B\| at `Bmin` | refused at load: "cannot be certified not to rise with \|B\|" | start the B axis where the data is monotone, or let the projection's isotonic pass flatten it |
| `layers` thickness ≠ `t_eff_m` | every `Ic` scaled by the ratio, no warning | §7 |
| `Ic` at a criterion other than 1 µV/cm | `jc` biased by `( E₀/E₁ )^(1/n)` | convert before building |
| An old copy of the file in the run directory | old table used instead of the new one | a file named as written, in the run directory, wins over `$BELFEM_DATA/material` (`share/material/README.md`) |

## 9. Checklist

- [ ] grid: order 2, odd node counts, θ = `[ 0, π ]`, last θ column equals the first
- [ ] values: log₁₀ of A/m² and of n, C order, no NaN
- [ ] `meta`: `n_floor`, `n_ceiling`, `t_eff_m`, `ec_criterion`, `reference`, `license`, `accuracy`
- [ ] projected with `project_table.py` (and cross-checked against `mattest` if a build is available)
- [ ] `verify_tables.py` prints `all tables pass every gate`
- [ ] table compared with the source data, error recorded in `meta/accuracy`
- [ ] `layers` thickness in the deck equals `meta/t_eff_m`
- [ ] the file loads in BELFEM without an `open_checked` refusal

## Related documentation

- `src/physics/database/doc/database_usage_guide.md`: the interpolant, the projector and the
  clamping conventions
- `share/material/README.md`: the shipped tables, their contracts, and how to modify them
- `resistivity_laws.md`: what `powerlaw`, `piecewise` and `riva` do with `jc` and `n`
- `share/python/database/howtoread.py`: the reference reader

This guide has not been tested on a table built from cable-level data. The end-to-end recipe in
§3–§5 was run on a synthetic anisotropic Kim-model table through the numpy route. The `mattest`
route and the BELFEM load of the result were not run for this guide.
