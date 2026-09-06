# Devlog 2026-08-15 — tapestack3d Outer-Copper |J| Asymmetry (Overnight, Read-Only)

**Date:** 2026-08-15
**Topic:** Why one outer copper layer of the tapestack3d thin-shell stack renders as
jagged z-periodic |J| fingers while the other renders as smooth bands
**AIs involved:** Claude (primary + data forensics), Codex (formulation-side audit),
Grok (postprocessor-side audit) — parallel blind dispatch, reconciled
**Claude Confidence:** high (electrical topology, solution-side localization);
medium (~60 %) on the mechanism (as corrected — the overnight ~70 % tolerance-floor
claim was retracted next morning, see Summary)
**Codex Audit Confidence:** medium | **Grok Audit Confidence:** high (mapping), medium (mechanism)
**Literature References:** Messe et al. 2023 §4 (checkerboarding / 1e-11 criterion);
Burman & Zunino 2006, Rivière-Wheeler-Girault 2001 (ghost penalty context)
**Verification:** read-only forensics on the run's own exodus/log output (numeric
probes, level 4) + static source trace (level 5); overnight part read-only, "reviewed"
not "verified"; the morning `_Volumes` fix is syntax-checked only, executable gate =
Christian's rebuild + `make check-fast`

## Summary

The premise "the two outer copper layers should look alike" is wrong for this stack:
the magnesia buffer has no `rho`, so `ThinShellFactory::create_buffers` turns it into a
`DomainType::Buffer` φ-region and `create_ghost_facets` skips its interfaces — each tape
is two separate conducting halves (Cu/Ag/YBCO vs hastelloy/Ag/Cu) joined only at the
coated edges and, across tapes, through the solder. The jagged face is the YBCO-side
copper; the noise is a **standing, slowly growing solution mode** (inter-save
correlation 0.996–0.998) concentrated in the layers that current-share with the
saturating YBCO, worst on the outermost tapes, absent on the buffered top face. [Corrected next
morning: the overnight draft tied it to an unmet 1e-11 tolerance — wrong, a dB
misreading; the iteration exits at −108…−117 dB = 1.6e-11…2e-12, i.e. AT tolerance
under the code's 10·log10 convention. The mode lives in the converged solution;
candidate controls are the ghost-penalty calibration, the edge z-resolution, and the
collapsed-Δt BDF5 dynamics. Newton's bit-identical-residual stalls near the floor
remain an efficiency observation only.] No defect in EF_PENTA6TS (E/C hand-verified
analytically), the orientation/tie machinery, the T-matrix path (inactive on linear
meshes), or the postprocessor's layer handling. Full report with tape×plane noise
matrix, evidence chain, and morning verification plan:
**`todo/tapestack3d_J_asymmetry.md`**.

## Key Findings

- Noise matrix (exodus s.00119, t = 5.95 s): planes 0/1 (YBCO-side copper) 3.8–4.7 %
  relative roughness on tapes 1/8, plane 7 (top copper outer face) 0.5–0.7 % on all
  tapes; tape 2 anomalously clean (open question).
- Every tape's layer elements are owned by a single MPI rank (exodus `ElementOwner`) —
  all MPI/aura/partition postprocessing hypotheses dead for the tape faces.
- **By-catch defect (real, inert here):** `Kernel::compute_element_volumes` thin-shell
  branch (`cl_FEM_Kernel.cpp:944-950`) indexes facet areas with an owned-only counter;
  wrong SPR weights on any mesh where a rank's owned shell elements are not a prefix.
  Inert in this run only because tapes are single-rank, contiguous, and identically
  meshed (proven by reconstructing the corrupted lookup from `ElementOwner`). **Fixed
  next morning on Christian's approval** (counter now advances per element; new
  size-precondition assert), report §6 A1 ticked.
- ~~By-catch observation: magnetic loop exits at residual ≈ 3·10⁻⁶ despite
  `tolerance : 1e-11`.~~ **RETRACTED next morning** — the dB convention is 10·log10,
  the exits are at ≈1e-11; report §3b/§5 carry the correction and the lesson.
- Known interface-node SPR leak confirmed at work: nodal J on the buffer/silver-L1
  faces shows the YBCO's ~2.5·10¹⁰ A/m² (insulator block rendering at 1.3·10¹⁰);
  inner-plane nodal J is unusable by design, outer planes are unaffected.
- Doc drift flagged by Grok (recovery-theory line ranges, node-sharing doc's factory
  path, side-connector "WIP stop" vs built walls) — refresh queued (report §6 A5).

## Changes Made / Proposed

- Overnight: no source, input, or run-state modifications (read-only mandate).
- Angle-convention unification applied (Christian's approval, afternoon):
  `compute_superconductor_ts` now computes β via `Calculator::bn_angle` (the
  assembly's folded |n·b| convention + small-field guard) instead of its own
  signed `acos` — display-only, syntax-checked; closes the proven inconsistency
  from the round-3 audit; the ~+30 % JJC bias question remains A7 step 1.
- Morning follow-up (Christian's approval): `_Volumes` indexing fix applied in
  `cl_FEM_Kernel.cpp` (counter advances for every shell element; new size-precondition
  assert) — **both auditors approved**; plus A1b, a second fix found by Codex's audit
  and verified by reading (rank 0 sized remote not-my-volume buffers with its own
  count; now `tProcIDs.length()`). Both syntax-checked with build-tree flags;
  rebuild/rerun is Christian's.
- **Cause of the oscillation closed, three-AI verified** (report §3b-ii/iii): the
  copper is an E-field microscope of mesh-scale overcriticality texture in the YBCO.
  Round-3 corrections folded in: true Ic ≈ 1.5–1.8 kA (sp-ap jc ≈ 3·10¹⁰, deck
  comment stale) so the run sits at ≈ 1.0–1.2 · Ic; local n_eff ≈ 17–22 (database),
  true J/Jc ≈ 1.10 with ~0.2 % texture, ×19 amplification → the 3.8 % fingers;
  signature triple measured (JJC +0.59 / |J|_ybco −0.19 / |B| −0.04). Tape-2 anomaly
  resolved (uniform J/Jc). Literature: closest analog Messe 2023 §5 staircase;
  mitigation = graded tip refinement; NOT §4 checkerboarding. By-catch: assembly
  |n·b| vs postproc signed-acos angle convention; written JJC biased ≈ +30 %.
- Tolerance enforcement measured over the full log: 568/569 step-final magnetic
  residuals ≤ 1e-11.
- Created `todo/tapestack3d_J_asymmetry.md` (report + action/verification checklists).
- Remaining proposals: A3 ghost-penalty calibration check (Ag|YBCO, deep saturation —
  physics call), A4 extend the circulation battery to the TS element family, A5 doc
  refresh. (A1 applied; A2 obsolete after the retraction.)

## Open Questions

- Tape 2's cleanliness (not ownership, not geometry — cut/terminal topology?).
- ~~Origin of the −110 dB residual floor.~~ Resolved: −110 dB IS the 1e-11 tolerance
  (10·log10 convention).
- Whether the electrically-severing rho-less buffer is intended modeling; if yes it
  belongs in `doc/input_file_reference.md`.

## Appendix: forensics recipe (scipy netcdf, read-only)

Per-plane extraction used throughout (tape k = 1..8, layer block `tape_7X_...`,
`which` = 'low'/'high' node plane of that thin layer volume):

```python
from scipy.io import netcdf_file
import numpy as np
f = netcdf_file('hphi_results.e-s.00119','r',mmap=True)
names  = [b''.join(n).decode().strip() for n in f.variables['name_nod_var'][:]]
bnames = [b''.join(n).decode().strip() for n in f.variables['eb_names'][:]]
x,y,z = (f.variables['coord'+c][:].copy() for c in 'xyz')
g = lambda s: f.variables['vals_nod_var%d'%(names.index(s)+1)][0].copy()
J = np.sqrt(g('Jx')**2+g('Jy')**2+g('Jz')**2)
def plane(bn, tk, which):
    yc0 = (-0.35+0.1*(tk-1))*1e-3
    conn = f.variables['connect%d'%(bnames.index(bn)+1)][:]
    nid = np.unique(conn)-1
    idx = nid[np.abs(y[nid]-yc0)<0.05e-3]; yy = y[idx]
    ym = 0.5*(yy.min()+yy.max())
    return idx[yy<ym] if which=='low' else idx[yy>ym]
# roughness metric: mean |second difference| of J along z within 0.05 mm x-bins;
# render: matplotlib tripcolor of (x,z,J[sel]) with vmin/vmax 6.7e4/1.8e5.
# ElementOwner reconstruction of the _Volumes bug: for each rank r, positions
# idx = where(owner==r) in block order, area_used[idx] = area[arange(len(idx))].
```

## Files Updated

- src/fem/kernel/cl_FEM_Kernel.cpp (the `_Volumes` fix, morning)
- todo/tapestack3d_J_asymmetry.md (new; corrected next morning)
- todo/README.md (registration)
- devlog/README.md (this entry)
- tmp/ai_exchange/tapestack3d_j_asymmetry.md (ephemeral thread, distilled here)
