# Plain QUAD4 unblocked in the Calculator; the rectangular-only rule for QUAD/HEX documented

**Date:** 2026-08-27
**Purpose:** Two dispatch cases let quad-meshed 2D blocks allocate Calculators for physics
that never touch nedelec data, gated by a QUAD4 twin of the DofSeeding fixture; and the
known edge-element distortion hazard is now written down with its literature trail
**Register:** no new row — by-catch of the DR-91 close-out, fixed same-session
**Thread:** `tmp/ai_exchange/quad4_nedelec_dispatch.md`

## The defect, and how it grew from one site to two

Found while gating DR-91's H-C fix: `Calculator::allocate_memory()`'s nedelec-data switch
runs unconditionally for every block/sideset Calculator and had no plain-QUAD4 case, so a
quad-meshed 2D block aborted with "Unsupported Element Type" under *any* physics —
including thermal/Poisson, which never use nedelec data. The pre-registered plan proposed
one line; both plan audits independently found a second rejection site behind it: the 2D
`mFundV` volume dispatch accepts only TRI geometry and `QUAD4TS`, so plain QUAD4 died
again at the "Higher order thin shells are not implemented!" error. The plan's stop rule
fired, the scope was amended in the open, and both edits landed:

1. `case ElementType::QUAD4 :` joined the `PENTA6TS/QUAD4TS/HEX8TS/HEX8TB/HEX8` group in
   the nedelec switch. The bound function is a pure dof-gather (`nedelec_data_linear`
   over `edge_h`), only ever invoked by h-φ machinery; both audits rejected the `nullptr`
   alternative because the wrapper calls the pointer unguarded.
2. The 2D volume branch binds `QUAD4 → dV_hex`. Grok's plan audit had explicitly warned
   against the tempting `dV_tri6_tet10` (its `mIsCurved` shortcut yields a constant
   `det J` — silently wrong on a general bilinear quad); `dV_hex`'s fallback caches
   `det( J(aIndex) )` per integration point, which is the correct isoparametric weight,
   and its edge-function branch takes the null path because the `EdgeFunctionFactory`
   cannot build a plain QUAD4.

Gate: the `DofSeeding` fixture is now parameterized by element type (TRI3 default keeps
the two existing tests behavior-identical; a guard fails any other type loudly), and
`DofSeeding.SeedModeOnQuad4Block` drives two CCW quads through `initialize( true )` with
the full seeding assertions. Red evidence exists for both sites. The gate proves dispatch
acceptance and seeding — not the assembled Jacobian; the debug `adV >= 0` assert backstops
inverted elements at the first real assembly.

Process: light two-vendor round both stages, no blocker anywhere; Grok's comment-accuracy
advisories applied (cache-key wording; anti-fallthrough guard). Known leftovers, recorded
not fixed: the `mesh::number_of_nedelec_dofs(QUAD4)==4` vs `fem::num_nedelec_dofs==0`
table split, and the future footgun that adding QUAD4 to the edge factory would silently
route `dV_hex` onto its hex `update_nabla` branch.

## The rectangular-only rule (Christian's instruction, literature-verified)

Quadrilateral and hexahedral elements in a Maxwell problem must be **perfectly
rectangular** — distorted (trapezoidal/sheared) quads and hexes break the edge-element
convergence theory, silently. The literature check confirmed and sharpened the claim from
the local library itself:

- Monk 2003 builds the hexahedral Nédélec theory on "parallelepipeds with edges parallel
  to the coordinate axes" (§6.1) and warns of "non-optimal convergence rates (or even
  non-convergence)" on non-parallelepiped hexes, even under a plain trilinear map
  (§8.2–8.3).
- Boffi et al. 2013, Remark 2.5.5: on general quad meshes the lowest-order H(div) family's
  divergence does not converge at all.
- The definitive papers: Arnold/Boffi/Falk(/Gastaldi) 2001, 2002, 2005 and
  Falk/Gatto/Monk 2011 (the 3D H(curl) case). All four DOIs verified against the
  publishers this session.

Written down in three places: `src/fem/interpolation/doc/nedelec.md` §6.6 (full mechanism
and meshing rule; simplex elements are exempt — any non-degenerate TET is affine),
`maxwell_usage_guide.md` §1.7 (new Common Pitfalls entry), and a new subsection with the
four citations in `doc/literature_references.md`. Codex language sweep applied to all
three passages; its by-catch fixed too — the nedelec.md factory-dispatch table was missing
the `HEX8`/`HEX8TS`/`HEX8TB` rows the factory actually serves.

## Status

Syntax gates clean (target `flags.make`, all TUs). **Reviewed, not verified** — runtime
gate owed: `make check`, expecting all four `DofSeeding` tests green including
`SeedModeOnQuad4Block`. `tests/fem/test_DofSeeding.cpp` is untracked and needs `git add`
at commit time.
