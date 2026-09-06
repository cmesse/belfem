# Reviving the 2022 gantry-dipole deck on the modern input format

**Date:** 2026-08-12
**Purpose:** Port `cmake-build-debug/gantry/` from the 2022 input format (464 hand-written cuts) to the current `input.conf` contract with automatic cohomology cuts and an SAE1010 yoke
**Module:** fem/maxwell (deck only — no source change)

## What the model is

`gantry.msh` is a 2-D upper-half model of a gradient bending dipole: a C-shaped
iron yoke whose shaped pole runs from `y = 25 mm` over `x = 59.7 … 91.65 mm` up
to `y = 52 mm` at `x = 203.3 mm` (so the full gap opens from 50 mm to 104 mm
across the aperture), and a return leg touching the midplane over
`x = 381.13 … 686.67 mm`. Two coil packs of four pancakes each sit under the
yoke bottom at `x ≈ 31.5 … 55.1 mm` and `x ≈ 336.4 … 360.0 mm`.

The winding is BSCCO tape, 5 mm wide and 0.396 mm thick, on a 0.414 mm pitch —
96 % fill. Each turn is meshed as a **line**, not an area: the tape is a
thin shell whose 5 mm width is the sideset and whose 0.396 mm thickness is the
layer laid along the facet normal. The rectangles between consecutive turns are
air, which is why the 2022 deck put *every* block except the yoke in `air`.

8 pancakes × 58 turns = **464 tapes**, all in series.

## The mesh carries no physical groups

`gantry.geo` defines no `Physical Surface`/`Physical Curve`, so gmsh 4.1 writes
every geometric entity and BELFEM uses the entity tags directly as block and
sideset ids. Classification, derived from `gantry.msh` itself rather than from
the old deck:

| ids | meaning | how identified |
|---|---|---|
| block 1 | iron yoke | boundary is `{14,16..32}` = `Curve Loop(1)` in the geo |
| blocks 2:932 | air (including all inter-turn gaps) | everything else; all 932 surfaces carry elements |
| sidesets 1, 2 | far-field arcs | the two `Circle` entities |
| sidesets 3:13, 15 | `y = 0` midplane in air | both endpoints at `y = 0` |
| sideset 14 | `y = 0` midplane inside the yoke | the one axis curve in block 1's boundary |
| 928 sidesets | the two 2.5 mm halves of each tape | vertical curves of length `w/2 = 2.5 mm` |

The tape sidesets are eight contiguous runs of 116:
`207:322, 497:612, 790:905, 1080:1195, 1370:1485, 1660:1775, 1953:2068, 2243:2358`.
They start at 207 — the same id the 2022 deck's `tape` list starts at, so gmsh
entity numbering has not drifted since 2022 and the old deck can still be read
as a cross-check.

Mesh size: 42 414 nodes, 84 591 TRI3, 13 229 LINE2, 1 427 point elements
(every geometric point becomes a vertex). 10 facets per tape half, i.e. 20
elements across the 5 mm tape width.

## What changed against the 2022 deck

**Cuts.** The 2022 deck listed 464 `cut tapeNNN { … }` entries, and the geo
still carries the geometry that supported them: eight `Circle` arcs from the
midplane points 23–26 / 28–31 up to the first turn of each pancake, plus the
mid-height horizontal line that splits every inter-turn gap. Those are now
inert — `homology { algorithm : generalized pellikka ; }` generates the cuts.
The arcs survive as internal air–air interfaces, which
`Topology::sideset_type( Air, Air )` classifies `Inactive`, so no geo change was
needed.

**Currents.** The 2022 amplitudes ran 20137.6 A down to 347.2 A in 58 steps of
347.2 A, repeated for each of the eight pancakes. Those are **cumulative cut**
values, not transport currents: 20137.6 = 58 × 347.2. The physical current is
347.2 A in every turn. In the new deck each tape is one `curve` built from its
two half-sidesets, and

```
current
{
    input curves : 1:464 ;
    amplitude    : 347.2 A ;
}
```

relies on `Section::get_id_groups` treating an **unbracketed** id list as one
group per id, so the range expands to 464 separate boundary conditions rather
than one 464-curve terminal.

**Midplane.** The 2022 deck declared the whole axis `antisymmetry`. That is the
B·n = 0 magnetic wall, which is also what an undeclared one-sided air face gets
by default (`Topology::sideset_type`), so the old deck effectively declared the
default. For this magnet the lower half is the mirror image with the current
**preserved**, which gives B×n = 0, so the new deck declares
`air symmetry : 3:13, 15` and `ferro symmetry : 14`. `MaxwellFactory`
implements those as `impose_dirichlet( 0.0 )` — legitimate here because φ is odd
about the midplane. Sideset 14 has to be declared explicitly: left undeclared,
a one-sided ferro face is classified `FerroAntiSymmetry` and deactivated, which
would force B·n = 0 inside the return leg.

Declaring the midplane Dirichlet raised the obvious worry — a loop that runs
through the upper half from one midplane point to another and closes along the
midplane encircles a nonzero ampere-turn count, so it must cross a cut, and if
the cut *terminates* on the midplane both sides of its landing node look pinned
to φ = 0. **That worry is unfounded, but not for the reason first assumed.**
`CutFactory::unflag_symmetry_sidesets` does not steer cuts away from the
symmetry plane; it removes the symmetry simplices from the complex, which if
anything makes the midplane an attractive landing (and the 2022 geo's `Circle`
arcs landed there deliberately). What saves it is that the cut jump is **not**
a difference of two nodal values: `CutSet::create_duplicates` gives each
duplicate node the source list `{abstract nodes, original}` with unit weights,
so φ_dup = I + φ_org, and `IWG_Maxwell::set_currents` fixes the abstract dof to
the terminal current. `SideSet::impose_dirichlet` pins only the nodes carried by
the midplane facets, which keep their originals. Landing on a Dirichlet boundary
therefore gives original = 0, duplicate = I, jump = I — consistent. Confidence
medium-high; the honest gate is a run that writes the cut meshes and checks
∮H·dl on a midplane-closing loop.

**Yoke material.** `RoxieIron` → `SAE1010`, both groups of `bhdata.hdf5`
(copied into the run directory, which wins over `$BELFEM_DATA/material`).

**Tape material.** BSCCO is a multifilamentary composite, so it stays one
homogenized layer spanning the full 0.396 mm, as in 2022. The 2022 deck used a
constant `rho = 1e-10 Ohm*m`; the new deck uses the E-J power law with an
engineering `jc = 2.5e8 A/m²` over 5 mm × 0.396 mm → Ic = 495 A, 1.43× the
operating current, and `n = 15`. `builtin : ybco` is only the power-law carrier
— there is no BSCCO builtin. **These two numbers are placeholders** and should be
replaced with the measured Ic(B,T) of the actual tape.

## Files

- `cmake-build-debug/gantry/input.conf` — the new deck (generated; the 464 curve
  definitions carry coil / turn / coordinate comments)
- `cmake-build-debug/gantry/input.txt` — the 2022 deck, kept for reference
- `cmake-build-debug/gantry/input.conf.draft.bak` — the partially-ported stub
  that was in the directory beforehand

## Open

- Not yet run. The load-bearing unknown is whether `generalized pellikka` plus
  the Smith normal form finishes at 464 generators on 84 591 triangles at
  acceptable cost. Both auditors named this, not the boundary condition, as the
  thing most likely to kill the revival.
- **Ampere-turn sanity.** 232 turns × 347.2 A across the 25 mm half-gap at the
  pole tip works out to B ≈ 4 T, well past SAE1010 saturation, which would make
  the shaped pole pointless. Either the design current per turn is lower than
  the 2022 arithmetic suggests, or the magnet is deliberately a saturated
  superferric. Worth checking against the design value before trusting any
  field number this deck produces.
- Timestepping is copied from greg5/costheta (10 s ramp, 15 s simulated), not
  from 2022 (5 ms sigmoid over 20 ms) — a different Faraday regime. Nonlinear
  `tolerance : 1e-7` is looser than the ε < 10⁻¹¹ that Messe et al. 2023 asks for
  after the Picard stage to avoid HTS checkerboarding.
- `initial conditions { temperature : 77 K ; }` is the `hphirun` default written
  out explicitly, not a statement about the operating point.
- Reviewed only — see `tmp/ai_exchange/gantry_revival.md` for the Codex and Grok
  audit round.
