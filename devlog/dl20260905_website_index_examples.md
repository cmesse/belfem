# Website: index overhaul and examples prose check

**Date:** 2026-09-05
**Purpose:** Record the rewrite of the public site's `index.html`, the fact-check of `examples.html` against the gallery images and the shipped decks, and the Codex language sweep over both.
**Scope:** `~/html/` (the belfem.lbl.gov working copy, not part of this repository). `scls.html` untouched by instruction.

---

## index.html

The old page was three paragraphs from before the release: it still said BELFEM "is not
available to the public" and promised a future BSD release, named VIPER cables (nothing in the
tree models one), and used the pre-redesign header (no brand, no Examples link, Roboto font).
Rewritten on the design vocabulary `scls.html` and `examples.html` already use (`eyebrow`,
`section-band`, `section-head`, `lead`, `footer-text`):

- hero: what BELFEM is, and the range from one tape to a full cable or magnet cross-section
- "What It Does": h-ɸ, thin shells with resolved layers, automatic cohomology cuts; coupled
  heat conduction, lumped circuits (inline RLC and SPICE), nonlinear materials, periodic
  translation + twist
- "Built for HPC": C++17 + Fortran, MPI, the solver list, Armadillo/Blaze, gmsh in, Exodus out,
  SCLS as the dependency stack
- "Applications": the gallery and the eighteen shipped decks; BCMT + Polytechnique Montréal
- "Availability": 0.9.0 released 2026-09-04 under BSD-3-Clause-LBNL (from `CITATION.cff` and
  `LICENSE`), GitLab link, four-line build, citation line
- footer carries the Regents copyright line and the CORC® trademark attribution

Facts were taken from `README.md`, `CITATION.cff`, `LICENSE`, `examples/README.md` and the
deck READMEs. Open point for Christian: the GitLab URL
(`https://belfem.lbl.gov/gitlab/codes/belfem`) is the origin remote; whether it is publicly
readable was not checked.

## examples.html

Each of the six cards was checked against its image (frame extracted from the GIF) and the
matching deck README. Changes:

| card | was | problem | now |
|---|---|---|---|
| Undulator | "Periodic magnet array; 2D magnetostatic field" | `undulator2d` is a 1200 s transient h-ɸ run with REBCO thin shells at 15 K and a user current source; not magnetostatic, not a permanent-magnet array | HTS undulator, two rows of REBCO tape coils on an iron yoke, flux density in yoke and beam gap |
| Cos-theta | "used to validate field level and harmonic content" | the deck does no harmonic analysis; it demonstrates symmetry/antisymmetry planes with a B-H iron and an 8498 A ramp | 68 conductor blocks, nonlinear yoke, symmetry planes, 8.5 kA ramp |
| Gantry | generic | fine, but unspecific | upper half, 464 Bi-2223 thin-shell tapes, `j/jc` plus `|B|` (both legends are in the image) |
| CORC | "resolves inter-tape current sharing and AC loss" | neither is modeled or shown; the decks model thin-shell tapes under a 50 Hz transport current with a periodic pitch | thin-shell tapes, `j/jc` and direction, translated-or-twisted periodic pitch |
| Tape stack | "a local defect" | the image shows two defects | two defects; matches no shipped deck (`tapestack3d` has 8 tapes and no defect, `tape_quench_usermat` has one tape and one defect) so the prose describes the picture only |
| Racetrack | fine | image shows one coil with two leads, `j/jc` and B arrows; `racetrack_usermat` has three coils, so no deck count is quoted | unchanged in substance |

Lead paragraph corrected: the old one claimed every case uses meshes for thin structures; the
cos-theta deck is bulk conductors. Now: every case is h-ɸ; tapes are thin shells; bulk, iron
and air are volumes.

**CORC®:** per Christian's standing rule, every occurrence on the site carries `&reg;` (prose,
meta description, alt text), and the trademark attribution sentence from the deck READMEs is
in the examples footer and the index footer.

## Nav consistency

`publications.html` and `contact.html` still had the old header (no brand, no Examples link)
and the Roboto font links. Header block and font links replaced to match the other pages;
three journal citation typos fixed in `publications.html` (Riva et al. "Volume 3" → 33, a stray
`)` and a missing space). Nothing else on those pages touched.

## Language sweep

Codex luna/medium, slug `website_prose`, prose only, technical claims and CORC® frozen.
23 proposals, no suspected errors. 21 applied; rejected: "magnetothermal" (the README spells it
magneto-thermal), "Simulation Examples" for the eyebrow (the index links to the "example
gallery"), and a caption that would repeat the "2D Examples" heading.

## Status

All four touched pages parse with balanced tags. Not rendered in a browser. Reviewed, not
verified.

## Addendum: gallery merged into the index

Christian's call, same session: the six field plots should be the first thing a visitor sees.
`examples.html` is retired and its two gallery bands now sit directly under a shortened hero on
`index.html`; the "Applications" band went away (it was the card list in sentence form), its
collaboration sentence moved into "What It Does" and its eighteen-decks sentence into
"Availability". Page order: hero, 2D examples, 3D examples, What It Does, Built for HPC,
Availability. The Examples nav item is removed from `index`, `publications` and `contact`
(`scls.html` still carries it; not touched by instruction, so its Examples link now points at
a deleted page until Christian edits that file). CORC® attribution stays in the index footer.
Tags balanced on the three edited pages; not rendered.

## Addendum 2: plugins and a Materials section

Christian's three requests: the hero says "HTS tapes" instead of "REBCO and Bi-2223 tapes";
"What It Does" gets a **User plugins** item (sources and boundary conditions, defect functions,
heat loads, complete materials, each a shared library named in the deck); and a new
**Materials** band between "What It Does" and "Built for HPC", listing what the database
ships. Roster taken from `cl_MaterialFactory.cpp` (nine pure metals, Hastelloy C-276,
magnesia, YBCO), formula alloys from `alloy_homogenization.md` (mass-percent formulas over
Al, Cr, Fe, Ni, Cu, Ag, In, Sn, Pb — the "mix solder materials" line), the three B-H groups
from `bhdata.hdf5` (RoxieIron, SAE1010, SAE310@4K), the six conductor tables from
`share/material/README.md`, and the three resistivity laws (`powerlaw`, `piecewise`, `riva`)
from `doc/input_file_reference.md`. Section rules re-alternated. `style.css` gained a scoped
`table.material-table` rule (full width, narrow label column); nothing else in the stylesheet
changed, so `scls.html` renders as before. Tags balanced; not rendered.

## Addendum 3: Materials section detail

Right column of the table turned into bulleted lists (pure metals, solders and alloys, other
materials, B-H curves, HTS conductors, superconductor models). Added on request: RRR is a free
input for pure metals and alloys; Kohler tables give T- and B-dependent transport properties;
the REBCO tables come from the Robinson Research Institute HTS database (the README calls the
five figshare tables "Robinson-derived"); the Riva model links the EPFL thesis
(doi:10.5075/epfl-thesis-8754); closing paragraph on B-spline storage giving continuous first
derivatives for Newton-Raphson (`SplineLookupTable`, `BhSplineCurve`, and the quadratic
tensor-mesh conductor tables). Two scoped CSS rules for `p` and `ul` inside the table cells.

## Addendum 4: three list entries sharpened

Magnesia named as MgO representing the buffer layer; YBCO entry states the normal-state
resistivity above the critical current and above T_c plus cp, λ (Callaway) and α over the whole
range (all present in `cl_Material_YBCO.hpp`); ROXIE iron identified as the default B-H curve
of CERN's ROXIE code, with a link to roxie.docs.cern.ch.

Addendum 5: the "fields up to 8 T" clause removed from the HTS conductors entry on Christian's request; the tables extend past the measured range and the number misleads.

Addendum 6: nav entry "Documentation" -> https://belfem.lbl.gov/doc/ (the Doxygen site) on index, publications and contact, placed after BELFEM; the Availability paragraph links it too. scls.html nav left to Christian.
