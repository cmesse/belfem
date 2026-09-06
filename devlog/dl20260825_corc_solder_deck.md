# CORC Solder Deck: Air Block 2 Becomes Solder, Quench Drive, Density Correction Tooling

**Date:** 2026-08-25
**Topic:** `cmake-build-debug/corc_solder` deck conversion + mesh-generator tooling
**Scope:** deck + python tool only — no framework source touched

## What changed

### Deck (`cmake-build-debug/corc_solder/input.conf`)

- **Block 2 (winding annulus, r 2.149–2.449 mm) converted from air to solder**:
  `conductor : solder { blocks : 2 ; material : solder ; }`, air now `1, 3`.
  Material: `builtin : Pb38Sn62 ; RRR : 1.5` (alloy-formula path,
  `cl_MaterialFactory.cpp:260-265`, RRR consumed by the `Alloy` ctor).
- **Current BC collapsed from six per-tape conditions to one bracketed group**
  (`[1,3,5,7,9,11]` / `[2,4,6,8,10,12]`): the solder shorts the six tapes into one
  connected conductor, so the six per-tape air generators no longer exist. Six
  separate conditions would die later as "more conditions than cohomology
  generators".
- **Second bearing node**: the solder annulus disconnects inner air (block 1)
  from outer air (block 3); each φ-component needs its own gauge node. Chosen:
  vertex **5** (fourth front quarter point, outer rim) — stable type-15 id that
  survives mesh regeneration and stays out of the periodic triples. An earlier
  in-session choice of a raw node id (21664) was replaced for exactly that
  regeneration-fragility reason.
- **Quench drive**: `sigmoid`, amplitude 600 A, period 10 s, offset 0,
  fuzzyness 1e-4 (tapestack3d's rationale: zero value *and* slope at t = 0).
  At constant jc = 1e10 A/m² the perpendicular tape width (3.25/3.10 mm at the
  35.6° lay) gives cable Ic ≈ 305 A, so 600 A ≈ 2×Ic with a 5 s plateau inside
  the 15 s simulation. Added `minimum timestep : 0.1 ns` (tapestack3d's guard).
- **`density correction : 0.735`** — see below; the first in-session value
  (0.469) was **wrong** and is corrected here.

### Density-correction convention (the error worth recording)

The annulus double-counts tape mass that the thin shells already carry. The
first derivation put *both* bounding layers' stacks inside the annulus
(→ 0.469). That contradicts the generator's radial-slot construction
(`r += tapeThickness + solderThickness` per layer): each layer's stack lies
radially **outward** of its shell, so the annulus between layers l and l+1
carries **one** stack — layer l's — and the outermost stack extends into outer
air, where there is no solder mass to correct. Cross-check: tapestack3d's
0.044 is exactly `(100 − 95.6)/100` per slab, the same one-stack convention.

    1 − 3 · 4 mm · 95.6 µm / (π (2.4486² − 2.1486²) mm²) = 0.735

Confidence high: analytic value and a direct mesh integration (tet volumes vs.
sideset areas) agree to 0.001.

### Mesh generator (`corc_solder/python/`)

- `Winding.solder_fractions(stackThickness)` — computes the factor per
  inter-layer annulus (cable length cancels; exact for straight centerlines).
- `Cable.stackThickness` (µm, physical BELFEM layers sum = 95.6, distinct from
  the 100 µm geometric slot), `Cable.density_corrections()`, reported by
  `C.print()` and emitted as a ready-to-paste comment by `write_topology`
  (maps annulus l → `Volume_{l+2}` per the assembler's region ordering).

### Realistic-geometry attempt: negative result, documented in `main.py`

Three configurations toward a realistic 45° lay / thinner solder bed all ran
gmsh's volume pass past the 1 h assembler timeout:

| pitch | solder bed | tapeRes | outcome |
|---|---|---|---|
| 2.15 (45°) | 50 µm | 0.15 | 198k nodes at timeout, insertion decaying |
| 2.15 | 100 µm | 0.15 | 163k nodes at timeout, worst-tet-radius spikes to 246 |
| 2.15 | 200 µm | 0.2 (proven regime) | 131k+ nodes vs. 74k for pitch 3 — **pitch itself is the driver** |

The row count along the cable is pitch-invariant (delta fixes it), so the
third experiment isolates the lay angle: at 45° the counter-wound layers cross
at 90° (vs. 71° at pitch 3), and gmsh's Delaunay refinement diverges between
the crossing constrained strips. `main.py` reverted to pitch 3 / 200 µm /
0.2 mm (the existing `corc.msh` is exactly that geometry — no regeneration
needed) with the findings in comments. Steepening the lay or thinning the bed
needs a dedicated meshing campaign (candidate levers: gmsh `Mesh.Algorithm3D`,
resolution grading near the crossings, longer assembler timeout).

## Findings flagged, not acted on

- **Modified Kim rule is not selectable from `input.conf`.**
  `cl_JcFunction_ModifiedKim.hpp` and the factory overload
  `create_jc_function(jc0, B0, k2, alpha)` (`cl_MaterialFactory.cpp:337`)
  exist, but the HTS parse branch knows only `file` or constant `jc`+`n`.
  `set_jc_function` composes with a constant n, so the wiring is a small third
  parse branch + the Input Contract updates. Awaits Christian's go and the
  plan+audit round.
- **Curves-vs-terminals for the current BC** (medium confidence): tapestack3d's
  deck comment argues a soldered stack must use bulk *terminal* sidesets —
  the thin-shell curves branch of `suggest_Homology` builds its generator from
  per-tape chains that miss the solder's current share, and all twelve corc
  terminal curves now lie inside the conductor. `corc.msh` has no solder-only
  end-face sidesets (11/12 span the whole plane), so switching needs either new
  physical surfaces from the generator or verification that the linked current
  is the total transport current. During a quench this is exactly the current
  that matters.
