# Fluid Property Tables

**Date:** 2026-08-05
**Purpose:** Build the fluid property tables BELFEM ships.

Two tools:

- **`build_tables.py`** produces `thermo.inp` and `trans.inp`. It extracts the species listed
  in `species.txt` from the NASA CEA database and merges a generated low temperature interval
  into each record where one can be fitted. The output is a complete table per property:
  there is no overlay file, no second table, and no load order contract.
- **`nasa9_lowT.py`** holds the fitting machinery and can also be run alone to produce just
  the low temperature records.

Extracting rather than shipping the vendor file is deliberate. The `thermo.inp` in
circulation is not purely a US government work: later versions carry propellant components
contributed by third parties and one block from a Russian database, none licensed. All of
those additions are solid and propellant materials, so a table restricted to the species
BELFEM needs contains none of them. That removes a redistribution question rather than
working around it, and drops 1.3 MB to 51 kB while doing so.

## Why, and where the handover belongs

The CEA thermodynamic tables start at 200 K or 300 K. The 2021 revision raised many of the
lower bounds to 300 K because the fits had not been validated below that, keeping the
coefficients unchanged.

Measured against the reference equations of state, those coefficients are good well below
the declared bound. Relative error in cp:

| | 100 K | 150 K | 200 K | 250 K | 273 K | 300 K |
|---|---|---|---|---|---|---|
| Ar, He | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % | 0.00 % |
| N2 | 7.15 % | 0.56 % | 0.00 % | 0.01 % | 0.00 % | 0.00 % |
| O2 | 15.44 % | 1.46 % | 0.00 % | 0.02 % | 0.01 % | 0.01 % |
| CH4 | 56.45 % | 4.53 % | 0.02 % | 0.00 % | 0.04 % | 0.05 % |
| H2 | 17.06 % | 2.82 % | 0.61 % | 0.12 % | 0.04 % | 0.01 % |

So the 300 K bound is caution about validation rather than a statement about the data, and
extrapolation below 200 K collapses quickly.

The handover therefore defaults to **250 K**. That is comfortably inside the useful range of
the CEA coefficients whichever lower bound the file declares, and it leaves 273.15 K and
everything above on pure CEA. Use `--junction` to move it, or `--junction 0` to hand over at
the start of the CEA interval instead.

## Method

Ideal gas heat capacity is sampled from CoolProp, which implements the same reference
equations of state the standard property tables are generated from. Two points matter:

- **Sample the ideal gas heat capacity, not the real gas one.** The NASA-9 form represents
  the ideal gas at 1 bar. Fitting real gas cp at a finite pressure would fold the departure
  into the polynomial, and the equation of state would then add it a second time.
- **Both target forms are linear in their coefficients**, and continuity of value and slope
  at the junction is a pair of linear equality constraints. The fit is therefore an equality
  constrained linear least squares problem with a closed form solution through the
  Karush-Kuhn-Tucker system. There is no optimiser, no starting guess and no convergence
  criterion, so the result is reproducible to the last digit.

The two integration constants are then solved exactly, so enthalpy and entropy are
continuous with the CEA data at the junction. This is not cosmetic. Those constants carry
the integral below the start of the fit, and the equilibrium constants are built from
absolute entropy on the third law scale, so a step there would move every one of them.

Transport coefficients are fitted the same way to the CEA form,

    ln( y / scale ) = A ln T + B/T + C/T^2 + D

against the dilute gas limit, which is what the CEA transport data describes. This matters:
at 1 bar the density contribution reaches 6 % in the conductivity of hydrogen at 25 K and
about 2 % for methane at 115 K, and folding that into a form meant for zero density would
be wrong. Four parameters cannot span a very wide temperature range, so intervals are added
until the tolerance is met. All intervals are solved **simultaneously**, with continuity between
neighbors as constraints. Fitting them one after another, each constrained to the one
above, lets the error at a shared edge propagate downward and can be worse than a single
interval: for hydrogen conductivity that approach gave 18.9 % where the simultaneous fit
gives 2.8 %.

## Usage

To build the tables BELFEM ships, which is what this directory is for:

```bash
python build_tables.py \
    --vendor-thermo <thermo.inp> \
    --vendor-trans  <trans.inp> \
    --species       species.txt \
    --out-dir       ../../share/fluid \
    --nist-transport
```

To produce only the generated low temperature records, without extracting a table:

```bash
python nasa9_lowT.py \
    --thermo  <thermo.inp> \
    --trans   <trans.inp> \
    --species Ar CH4 H2 He N2 O2 Ne Xe \
    --out-thermo lowT_thermo.inp \
    --out-trans  lowT_trans.inp
```

Requires `numpy` and `CoolProp`.

Useful options: `--margin` sets how far above the triple point to start, `--tolerance` the
transport residual that triggers another interval, `--samples` the number of sample points,
and `--dilute-density` the density used for the dilute gas limit.

The vendor files are inputs, not outputs. They are not redistributed, and nothing in BELFEM
reads them at runtime: the build is run once and its result is what ships.

## File conventions

Both tables are **fixed column and never exceed 80 characters**, including comments. The
readers work by byte offset, and while they ignore everything past the last parsed field, a
data file that quietly grows a 109 character line is one an editor or a downstream tool can
still mangle. The build checks the width itself and warns rather than trusting that it held.

Anything that is not the original vendor data is therefore marked with a **footnote reference**
on the record line, with the text of the note in the file header:

- `@1` — a low temperature interval was generated for this species
- `@2` — *(`trans.inp` only)* the record is derived from a lighter isotopologue by mass scaling

The notes carry the per species detail: fitted range, interval count and residual. This keeps
the record lines short and puts the explanation in one place instead of repeating a fragment of
it on every line. The marker separator narrows from two spaces to one where the citation is long,
so the 80 column rule holds without truncating anything.

**`thermo.inp` holds caloric properties only; `trans.inp` holds viscosity and thermal
conductivity only.** Neither file comments on the other's contents — a note about transport in
the caloric table would be in the wrong place, however true it is. Each header says where the
other property lives.

## Sources and authority

CoolProp is the sampling engine because it is scriptable, reproducible and can be pinned to
a version, none of which is true of a live web service. **NIST is the authority where the
two disagree.**

`--validate-nist` fetches a small table per species from the NIST Chemistry WebBook and
reports the largest relative difference. For most fluids the two implement the same
reference correlations and agree to the printed digits. Methane is the exception found so
far, where CoolProp differs by 2.18 % in viscosity and 0.77 % in conductivity, so its
transport fit should be treated as the weaker one.

`--nist-transport` covers species CoolProp has no transport model for, neon, xenon and
deuterium among them, by sampling NIST at low pressure instead. Without it those species get a
heat capacity fit only.

Both options issue one request per species. Neither is a bulk download.

## Output quality

Measured against CoolProp over the fitted range:

| species | range [K] | cp | viscosity | conductivity | junction jump |
|---|---|---|---|---|---|
| Ar | 84.3 - 250 | 0.00 % | 0.16 % | 0.09 % | 2e-16 |
| CH4 | 91.2 - 250 | 0.01 % | 1.43 % | 1.36 % (2) | 8e-16 |
| H2 | 14.5 - 250 | 0.51 % | 1.67 % | 0.97 % (3) | 8e-16 |
| He | 2.7 - 250 | 0.00 % | 1.25 % (2) | 1.34 % (2) | 3e-16 |
| N2 | 63.7 - 250 | 0.01 % | 0.34 % | 1.93 % (2) | 4e-16 |
| O2 | 54.9 - 250 | 0.03 % | 0.75 % | 0.63 % | 2e-16 |
| Ne | 25.1 - 250 | 0.00 % | 0.24 % | 0.59 % | 2e-16 |
| Xe | 161.9 - 250 | 0.00 % | 0.75 % | 1.99 % (2) | 4e-16 |
| D2 | 19.2 - 250 | 0.08 % (2) | 1.14 % | 5.85 % (3) | 8e-09 |

Neon, xenon and deuterium transport comes from NIST, everything else from CoolProp.

The bracketed number is how many intervals the fit needed, for cp in the first column and for
transport in the others. The junction jump is the largest discontinuity in cp, its slope, h
and s, scaled by a physical magnitude; all of them are at round off level. Deuterium is a
digit or two larger only because it is measured across two generated edges rather than one.

Deuterium conductivity is the one entry that misses the tolerance, and the residual is not a
fitting failure: **CEA and NIST disagree by 4.6 % at the junction itself**, and the fit is
constrained to meet CEA there while following NIST below. The cause is visible in the vendor
file. CEA's deuterium viscosity is exactly sqrt(2) times its hydrogen viscosity, to within
0.05 % at every temperature, so that record is mass scaled rather than measured; the
conductivity was derived the same way with a partial correction for the internal modes. NIST
carries measurements. Continuity is kept regardless, since a step in conductivity at 250 K
would be worse than a fit that absorbs a disagreement already present in the sources.

The generated records are checked by parsing them back with the same reader. The heat
capacities recover the right physical limits: 5/2 R for the monatomic species, 7/2 R for
nitrogen and oxygen where the rotational modes are active and the vibrational ones are
frozen, and 5/2 R for hydrogen at 14 K where even rotation is frozen out.

Residuals should be read against the uncertainty of the source rather than driven to zero.
The published uncertainties of the underlying correlations are roughly 0.5 to 2 % for
viscosity and 1 to 5 % for thermal conductivity in the dilute gas region, and about 1 % for
the heat capacity of hydrogen. Fitting tighter than that tracks one particular correlation
more closely without being more true.

The heat capacity escalates to more intervals the same way transport does, and for the same
reason, but the escalation almost never fires: one interval already holds every species here
to better than 0.6 %. It exists for the case where cp is **not monotonic** over the fitted
range, which a single seven term polynomial cannot follow.

Normal deuterium is that case. Its rotational heat capacity overshoots the classical 7/2 R,
peaking at 3.62 R near 98 K before settling back to 3.51 R, because the ortho and para forms
of a boson pair weight the rotational levels differently than in hydrogen. One interval gives
4.66 %; two, split near 69 K, give 0.08 %.

Hydrogen is the instructive contrast. Its heat capacity climbs from 5/2 R to 7/2 R across the
same range and never turns over, so a single interval reaches 0.51 % and splitting it does
not help — an earlier attempt at a fixed 100 K breakpoint made it worse. Shape, not range,
decides whether another interval is worth anything.

## Deuterium

Fusion fuel handling needs the hydrogen isotopologues, so the table carries `D2`, `D`, `HD`,
`OD`, `D2O` and the deuterium ions. CEA has thermodynamic data for all of them, and transport
for `D2` and `D2O` only.

`D`, `HD` and `OD` therefore get their transport **derived from the corresponding hydrogen
species** by mass scaling. In the dilute gas the collision integral depends on the electronic
potential, which isotopologues share, so only the mass prefactor of Chapman-Enskog changes:

    eta_child    / eta_parent    = sqrt( m_child / m_parent )
    lambda_child / lambda_parent = sqrt( m_parent / m_child )

Since the CEA form is `ln( y/scale ) = A ln T + B/T + C/T^2 + D`, a constant factor moves only
`D`. The other three coefficients are the parent's, unchanged, which makes every derived record
checkable by eye against its parent.

This is the same construction NASA used: the ratio of their `D2` to `H2` viscosity is 1.4135
to 1.4137 against sqrt(2) = 1.41421 over 300 to 4000 K. Derived records say so in their header
rather than inheriting the parent's citation, which would credit a measurement of hydrogen for
data that is not one.

The conductivity relation is exact only for the translational part; internal modes do not scale
with mass. So it is exact for `D` from `H`, which is an atom with no internal modes at all, and
an approximation for `OD` and `HD`. Its size is visible in CEA's own numbers, where the D2/H2
conductivity ratio drifts from 0.707 to 0.738 against the pure 1/sqrt(2) = 0.7071. For `OD` the
mass changes by 6 %, so the drift is correspondingly smaller.

Two limits are worth stating. **`D2O` gets no generated low temperature interval**: heavy water
freezes at 276.97 K, above the junction, so there is no gas phase to fit. And CoolProp will
refuse `Cp0molar` below the melting line even though `Cp0` is a function of temperature alone;
since the pressure is immaterial to the answer, the sampler retries at a lower one, which is
what makes deuterium reachable down to 19.2 K at all.

Tritium is not included. Neither CEA nor CoolProp carries it.

## Liquid transport

The correlations fitted here are the **dilute gas** curves, which is what the CEA transport
form represents. Liquid transport reaches BELFEM through the Lucas and Stiel-Thodos residual
corrections instead, and their accuracy there is worth knowing, because the two HTS coolants
sit at opposite ends of it.

For saturated liquid, the residual relative to the dilute gas base:

| | T [K] | rho_r | viscosity | conductivity |
|---|---|---|---|---|
| helium | 4.2 | 1.80 | 3.0 x | 2.3 x |
| nitrogen | 77.35 | 2.57 | 30 x | 21 x |

For liquid helium the residual is a genuine correction and the correlations are on firm
ground. For liquid nitrogen the residual **is** the answer, so the accuracy is that of the
correlation alone, with no help from the well characterised dilute term. Expect roughly
10 to 20 % there against 1 to 2 % in the gas.

CoolProp is more accurate in the liquid because its transport models are fluid specific fits
over the whole surface rather than corresponding states correlations. For nitrogen, from
Lemmon and Jacobsen 2004: viscosity as a Chapman-Enskog dilute term plus a higher order
polynomial in reduced temperature and density; conductivity as a dilute term, a residual in
temperature and density, and a critical enhancement. Stiel-Thodos has no analogue of that
last term, so conductivity near the critical point is a regime BELFEM cannot represent.

Fitting liquid curves with the machinery here would work along the saturation line, where
transport is single valued in temperature, but not off it: the CEA form is a function of
temperature only and liquid transport is a surface in temperature and density. That needs a
different table type, not a different fit.

**Helium has a hard floor at 2.1768 K**, the lambda transition. Below it the fluid is
superfluid and thermal conductivity is not a transport coefficient at all, since heat moves
by counterflow. No correlation applies and none should be extrapolated there.

## Limitations

- **Ozone is not in CoolProp or the WebBook fluid pages.** It needs statistical
  thermodynamics from spectroscopic constants, which is a separate piece of work.
- **Methane transport carries a known 2.2 % disagreement** between CoolProp and NIST. Until
  that is resolved, prefer `--nist-transport` for it or treat the fit as provisional.
- The species mapping from CEA label to CoolProp fluid is an explicit table at the top of
  the script. A label that is not in it is reported and skipped rather than guessed.
