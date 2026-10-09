# Pancake Inner Lead-Out: From the Hard Corner to a Bent Lead

**Date:** 2026-09-24
**Purpose:** Replace the inner lead of the single-tape pancake generator
(`examples/pancake_tape/python/pancake/`), currently an L-shaped hard-way corner in the bore wall, by a
lead that first separates from the inner turn on a smaller circle arc whose center is off the
axis and then bends gradually until the tape direction runs along z. The mechanism is a
piecewise-analytic centerline with a moving (t, n, b) frame, meshed as a structured strip, with
the solder bar offset along n.
**Module:** `examples/pancake_tape/python/pancake` (mesh generator, Python); consumer `src/fem/maxwell`
thin-shell pipeline and `src/homology` cut factory (read-only constraints)
**AIs involved:** Claude (exploration + plan), Codex (audit), Grok (third voice); Christian
(physics ruling on the lead-out shape)
**Status:** PARKED 2026-10-09 — the `pancake_tape` deck was removed from devel on Christian's
decision (it never produced a solution and does not ship in 0.9.3); the generator this plan
rewrites is recoverable from git (`git show 00b7972:examples/pancake_tape/python/pancake/`).
The brief stands as written for whenever the deck returns; no code changed

> **Scope guards:**
> - The outer lead (easy-way bend in the plane, solder bar inside) is OUT of scope.
> - The block/sideset rules of the cut factory stay as they are: every tape facet separates two
>   different blocks, air is only ever on the outer side, half-turn pieces
>   (`doc`: see `todo/pancake_inner_leadout.md` §1 and the generator's `validate()`).
> - The tape stays a zero-thickness thin shell; bending limits of real REBCO are physics
>   context, not a model constraint.

---

## 1. Current Behaviour and Why It Is Wrong

Frame nomenclature (Christian): **t** tape direction, **n** normal to the tape face,
**b** binormal = across the width, b = t × n.

In the winding, n is radial and b = e_z. The inner lead today is the sheet continued on the
bore wall r = ri: in the unrolled (r_i θ, z) plane the sheet is an L, the arm (width b_w = 4.11
mm, columns 0..nL at θ ∈ [θ_L, 0]) runs straight up to the cap
(`sheet.py`, docstring and `_build_columns`; `pack.py::_build_sheet`). That is a **hard-way
bend of radius zero**: t turns by 90° about n at a reentrant corner (θ = 0, z = b/2).

| Failure | Mechanism | Evidence |
|---|---|---|
| Unphysical joint | a 90° corner about n has infinite geodesic curvature; no tape does this | geometry, `sheet.py` L-construction |
| Current crowding | sheet current turning a reentrant corner has a singular density at the corner | classic L-shaped conductor; medium confidence, not simulated |
| Lead glued to the first turn | the arm sits on r = ri with the sliver (block 1/2) outside it for |z| < b/2, so the lead-out shares the winding's surface up to the cap | `pack.py` sliver cells k < 0, `ARM_UPPER` |

**Bottom line:** the lead must leave the cylinder surface of the winding before it turns, and it
must turn gradually.

## 2. Candidate Strategies

**S1 — off-center arc, then gradual hard-way bend (Christian's proposal as read by Claude).**
At the inner end the centerline leaves the spiral on a circle arc of radius r_1 < r_i whose
center C_1 is displaced from the origin (tangent-continuous with the spiral at the departure
point). The tape stays wrapped on the cylinder of radius r_1 about the axis through C_1
(n radial to that cylinder, b = e_z). Then, still on that cylinder, the centerline bends in the
unrolled plane with radius R_h from horizontal to vertical: t → e_z, b → −e_θ' (tangential to
the small cylinder), n stays radial to it. Then straight up to the cap. Surface: a strip of a
cylinder, developable, exactly the `Bend(hard=True)` of `examples/pancake/python/pancake/lead.py`
but on an offset cylinder. Meshing: structured strip in the unrolled (s, z) plane of the small
cylinder, mapped onto it; the solder bar is an offset along n (toward the axis C_1) with
columns copied at radius r_1 − t_bar. Strain of a real tape: hard-way bending needs
R_h ≫ b_w (order 1 m for REBCO), so the model shape is not a manufacturable tape; the joint
region is really the solder/copper block and the current transfer is through the block.

**S2 — off-center arc, twist, easy-way bend (the manufacturable lead).**
Same departure arc; then a twist of 90° about t over a length L_t (n: radial → e_z,
b: e_z → radial), then an easy-way bend of 90° about the new b (radius R_e ≥ 5 mm, as the
outer lead), so t → e_z; straight up to the cap with the tape face horizontal... i.e. the tape
rises with its flat face in a vertical plane containing the radial direction. Surface: the
twist piece is a ruled, non-developable helicoid strip; the rest developable. Meshing: columns
along the centerline with cross-section direction b(s); the bar offset along n(s) rotates with
the frame; the quads of the twist are warped but the Kuhn split tolerates it. The terminal on
the cap is a rectangle b_w × (thin-shell edge) plus the bar's face.

**S3 — helical exit without turning t to z.** Leave on the arc, then run a helix of constant
lay angle on the small cylinder to the cap; t never reaches e_z, the cap cuts the tape
obliquely. Cheapest, but the terminal edge is slanted and the tape must still change its lay
angle from 90° (hard-way, gradual) — S1 in disguise.

## 3. Open Design Questions (for the jury)

- **O1** Which of S1/S2 is the right lead-out for the h-φ thin-shell model given that the tape
  is a current sheet with layers stacked along n? Does the twist of S2 break any assumption of
  the thin-shell factory (layer normals from facet windings, node-normal averaging)?
- **O2** Departure arc geometry: r_1, the center offset, and the clearance from the winding
  (the sliver is 50 µm; the bar sits inside the arc). Is a single tangent arc enough or should
  the departure be a clothoid to avoid a curvature jump against the spiral?
- **O3** Hard-way bend radius R_h for S1 (model only): what governs it — element aspect
  ratio, current crowding, or nothing physical? Proposal: R_h = 2 b_w.
- **O4** Meshing of the transition: is the unrolled-plane strip on the offset cylinder (S1)
  or the moving-frame extrusion (S2) safer with respect to the cut-factory rules (the bar
  block must be inside, air outside, the sideset piece boundaries must stay on columns)?
- **O5** Where should the terminal be: the cap (z = H) as now, or a flat patch on a
  cylinder r = r_1 at some height, so the lead needs no full rise?

## 4. Plan Steps (after the jury and Christian's ruling)

- [ ] R1 Ruling on S1/S2 and O1–O5 (Christian).
- [ ] R2 `sheet.py`: centerline pieces (departure arc, transition, straight) with a frame
      (t, n, b) per column; keep the sliver and the half-turn pieces.
- [ ] R3 `pack.py`: lead columns carry a cross-section direction b(s) instead of the fixed
      e_z; the arm bar becomes an offset along n(s).
- [ ] R4 Boundary: the cap hole follows the new bar face; the tube unchanged.
- [ ] R5 Validation: existing rules plus a smoothness check of the centerline (no curvature
      jump above a threshold) and a positive-Jacobian check of the strip.
- [ ] R6 20-turn smoke run in both lead directions; 200-turn build; deck update.

## 5. Exchange Threads

- `tmp/ai_exchange/review_pancake_leadout.md` — jury round on this brief (2026-09-24).
