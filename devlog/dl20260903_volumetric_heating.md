# Artificial volumetric heat load: `heating { }` plugin on a material

**Date:** 2026-09-03
**Purpose:** Session record — a heater-pulse trigger for quench runs as an alternative to
the Ic defect: literature check, what the thermal side lacked, the wiring, and the
tapestack3d deck that uses it.
**Module:** `src/physics/materials`, `src/fem/kernel`, `src/fem/thermal`

## Question

`cmake-build-debug/tapestack3d_heating` quenches an eight-tape stack by removing 90 % of jc
on a 1.5 mm disk of the top tape (`src/defect.cpp`). Christian asked whether depositing an
artificial heat load at the same spot is a realistic alternative, whether BELFEM can do it,
and which function and W/m³ to prescribe.

## Literature

Read-only agent sweep of `./literature/` (extractions, `rg`). No FEM quench paper in the
library triggers with a heater: Wozniak et al. 2025 use an Ic defect, Badel et al. 2019
random Ic scatter (and argue in §1 that MQE/NZPV framing does not transfer to HTS), Riva et al.
2023 overcurrent. The heater lives in the books: Iwasa 2009 Eq. 6.1 carries the disturbance
term `g_d(t)` beside the Joule term (§6.2.5 taxonomy; Table 6.4 energy margin ≈ 8 J/cm³ per
5 K at 70 K; Disc. 6.8 MPZ radius 0.1–10 mm at 77 K; §8.5 heater-pulse-triggered NZP on
YBCO at 77 K, Table 8.5), and Russenschuck 2010 §18.7.1 models protection heaters as
volumetric sources with an exponential decay (Eq. 18.23, 20 W/m per cable, 30 ms delay).
Verdict: realistic, but energy-defined, and physically different from the defect — jc stays
intact, so the spot can recover (the classical MQE experiment).

## What was missing

The thermal load vector carried only `Nᵀ ρ |j|² dV` (`mt_thermal_h.cpp`); the thermal BC
factory accepts `dirichlet`/`gauge`/`bearing` and aborts on `neumann`; `Material` had a
defect hook but no heat hook. Christian added `HeatFunc`, `set_user_defined_heating`,
`read_heating`, `have_heating`, `volumetric_heatload` on `Material`, a
`compute_volumetric_heatload` stub on `MaxwellData`, and duplicated the two thermal kernels
as `T_h_*_heating`.

## Audit of the draft

- `MaxwellData::compute_volumetric_heatload` was unfinished (`return mMaterial->have_heating()`
  without `;`, never called the function).
- `~Material` did not `dlclose( mHeatHandle )`.
- `read_heating` carried the defect comment verbatim; the header had a doubled `@brief`.
- The duplicated kernels needed a material lookup in `link_to_group`, which the collapsed
  kernel deliberately does not do, and copied 120 lines that every future fix would have to
  mirror.

## Wiring (this session, source edits approved by the request)

- `MaxwellData` gets `mFunHeat`, bound once in the constructor:
  `compute_heatload_user` (calls `compute_x`, then the plugin with `mX, mY, mZ, mTime`) when
  the material has a heating plugin, `return_zero` otherwise — the same dispatch the ρ
  tangents use. The duplicated kernels were dropped; `T_h_picard` and `T_h_newton` add
  `mx->compute_volumetric_heatload( k )` to the Joule term. No tangent: the load is
  prescribed in (x, t).
- `MaterialFactory` reads `heating { file ; label ; }` on any builtin or usermat material
  after `density correction`; `heating` is on the unknown-input allow-lists of both shapes
  and the subsection is checked for `file`/`label` only (DR-146 discipline).
- `~Material` closes the heating handle; comments fixed.
- Input contract: row and paragraph in `doc/input_file_reference.md` §materials, `heating:`
  entry in `doc/input_schema.yaml` beside `defect:`.
- Deck: `src/heating.cpp` + `ramps::heating_pulse` (derivative of the deck's own logistic,
  unit integral), `CMakeLists.txt` compiles it into `userdefect.so`, `copper` gets the
  `heating` block, the ybco `defect` block is commented out.

Numbers behind the deck values: disk area π(R² + e²/2) = 7.21 mm² (the erf-smoothed disk
integrates to that exactly); both 20 µm copper layers → 2.88e-10 m³; adiabatic enthalpy of
the disk through the full 95.6 µm tape 77 K → T_cs(700 A) ≈ 13 mJ, → Tc ≈ 20 mJ, consistent
with Iwasa Table 6.4. Chosen: E = 20 mJ over a 50 ms period at t = 3.0 s → mean 1.4e9 W/m³
in the copper (peak 4.8e9), against 7.5e8 W/m³ Joule heating of a fully normal disk at
87.5 A per tape. Sweep 10 / 20 / 50 mJ.

## Evidence

Compiled the four touched sources with the build's own `-Werror` flag set
(`c++ <flags.make> -c … -o /dev/null`) and the three plugin sources with the plugin's
flags: clean. **Reviewed, not verified** — no `make`, no run. Owed: build + deck run;
check the deposited energy (∫ q̇ dV dt over the pulse should give 20 mJ — confirms that a
thin-shell layer block's `dV` is the layer volume, which was assumed, medium confidence);
Codex language sweep of the new `input_file_reference.md` paragraph; an
`example_user_heating.cpp` template next to `example_user_defect.cpp`.
