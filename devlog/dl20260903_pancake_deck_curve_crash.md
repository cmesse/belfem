# Pancake deck: two-chain curve crash, tape sign, 20 K current

**Date:** 2026-09-03
**Purpose:** Review of the pancake wedge deck (`cmake-build-debug/pancake`, python generator in its
`python/` subdirectory); root cause of the `sort_end_nodes` abort; tape orientation sign; operating
current against the SuperOx table.
**Module:** mesh (`cl_CurveFactory`), fem/maxwell (deck contract)

## Deck review

Every key of the deck resolves against `doc/input_schema.yaml`: mesh unit, the MUMPS block,
the nonlinear keys, `method : bdf5`, Anderson, the builtin names including the `Pb38Sn62` alloy
formula, `mum`, layer labels against the materials block, `generalized pellikka`, bracketed
terminal groups (one condition with seven input and seven output faces), `sigmoid` with
`fuzzyness`, `superox.hdf5` under `$BELFEM_DATA/material`. All 34 deck ids match the mesh's
physical names (the generator writes geometry tag = physical tag). The mesh's `$Periodic`
section is not read by `cl_Mesh_GmshReader`; the deck's `periodic { source ; target }` block is
what builds the periodicity, and its ids are **vertex** ids (`PeriodicityFactory::
set_master_plane` calls `mMesh->vertex( A )->node( 0 )`). The generator writes vertex elements
1..6 on those nodes and orders the nodes first, so the two id spaces coincide.

Minor: `initial conditions` appears twice with the same value (last wins); solder
`density correction : 0.044` is copied from tapestack3d and is inert in a magnetic-only run —
for this stack (155.2 µm of layers at 250 µm pitch) gap/pitch is 0.38; the alloy rho cache is
named `%s_RRR%u.hdf5`, so RRR 1.5 lands in `Pb38Sn62_RRR1.hdf5` — the stored RRR is compared on
reload (hard error, not silent reuse), so changing RRR within [1,2) needs the cache deleted.

## The crash

```
Row index 4294967295 out of bounds, which must be smaller than 8.
  CurveFactory::sort_end_nodes  cl_CurveFactory.cpp:1144
  CurveFactory::intersect       cl_CurveFactory.cpp:597   (aID = 1)
```

Curve 1 was `18 @ 9`: the outer cylinder (the wedge's terminal face, both leads end on it)
against tape 1. Tape 1 borders the air, so the generator used the domain face for both of its
ends — and on the wedge both end edges lie on that one face. Measured in the mesh: the
intersection has 8 shared nodes in **two disconnected chains**; every other pair in the deck has
4 nodes in one chain. `sort_end_nodes` counts one- and two-segment nodes (4 + 4 = 8, the
"branch" check passes), walks from the lowest-id tip along the first chain, leaves the other
four nodes at `gNoIndex`, and line 1144 indexes the adjacency with the sentinel. Curve 9 (outer
end) was the identical pair.

**Fix (deck + generator, applied):** tape 1's curves use the end face of the solder slab it
shares with tape 2 — `1 : 21 @ 9`, `9 : 28 @ 9` (one chain each, measured). Equivalent
semantics: `MaxwellFactory` attaches every deck curve touching a tape sideset as that tape's
terminal curve regardless of the other sideset's type (`cl_MaxwellFactory.cpp:3327-3355`).
`fem.py::topology_section` now takes tape k's faces from slab `max(k-1, 0)`; the box domain had
the same defect (both stack ends are holes in one top face).

**Diagnostic gap, not fixed (reported):** a two-chain intersection is a deck error, but the
factory reports it as an out-of-bounds index in a debug build and would run off the sentinel in
release. `sort_end_nodes` could count connected components after the branch check (a chain has
exactly two one-segment nodes; `tOne == 2 || tOne == 0` is the cheap form) and name the curve id
and the two sidesets. The cut factory's later "expect 2 tips" guard would also have caught it,
but the abort comes first.

## Tape orientation sign

`fix_facet_masters` puts the master of an air-adjacent tape on its solder side; inner tapes keep
the as-loaded master, the neighbour with the smaller element ids (devlog
`dl20260812_tapestack3d_setup_and_terminal_guard.md`, Addendum 2). The pancake mesh numbers air
first, then slab 1 to slab 7 (element ids 7–312446, 312447–321662, …), so every inner tape has
its master on the slab above it; tape 8's solder is above it and agrees, tape 1's solder is below
it and comes out mirrored. Same outcome as tapestack3d (`-5,-6`), same reason. Deck and
generator now write `sidesets : -9, 10:16 ;`.

**Jury round (`tmp/ai_exchange/pancake_tape_sign.md`), unanimous:** Codex terra/high and Grok
4.6/high (second attempt; the first hit its 30-turn budget and reversed itself three times,
discarded) both confirm the chain with high confidence, and the inner-tape rule is now derived,
not borrowed: the worklist is a FIFO (`cl_Queue.hpp:73-76`), the reader orders facets by
sideset (`cl_Mesh_GmshReader.cpp:1063-1070`), `check_facet_orientation` is the manifold
shared-edge test (`fn_check_facet_orientation.hpp:52-110`), and a same-type facet is oriented
once by the first seed that reaches it (`cl_MaxwellFactory.cpp:1553-1567`). The wall of slab
k−1 (sideset 33+k) precedes the wall of slab k, so tape k takes the outward sense of the slab
above; the facet winding is the outward face of the master (TET4 table = Exodus outward table,
checked on the unit tet), so master = slab above = as-loaded. Refuted in Codex: solder end faces
as seeds (one tet each, no slave). **Gate owed (Grok's falsifier):** master block of an interior
tape-2 facet after `fix_facet_masters`, expected block 2; and the first frame showing eight
stacks facing the same way. Open physics: `9, -10, …, -16` is the other uniform choice with the
opposite sense; `bn_angle` maps theta → π−theta under it and the SuperOx table is not folded
about 90°, so the sense is Christian's physical input, not a code question.

## Operating current

`share/material/superox.hdf5`: `jc` is log10 of the layer current density at t_eff = 1 µm,
array order (T, log10 B, theta), theta from the tape normal. Read at the nodes:

| T | B | Ic ⊥ (4 mm) | Ic ∥ |
|---|---|---|---|
| 20 K | 0.01 T | 2386 A | 2388 A |
| 20 K | 0.5 T | 1715 A | 2299 A |
| 20 K | 1.0 T | 1249 A | 2209 A |
| 20 K | 2.0 T | 877 A | 2080 A |

Check value: 77 K, 0.01 T, ∥ → 52.4 kA/m against `meta/Icw_77p5K_sf_A_per_m` = 51.2 kA/m at
77.5 K. (The `meta/documentation` text says 45×31×181 nodes at 2 K; the dataset's own
`points`/`step` say 89×35×181 at 1 K — the dataset wins.)

Self-field of the 4×2 mm stack, equivalent round conductor of 12 mm perimeter: ~0.105 T/kA,
perpendicular at the tape edges. Self-consistent at 20 K: 6 kA → I/Ic ≈ 0.48, 7 kA → 0.63,
10 kA → 1.00; a flat stack concentrates more at the edges than the round model, so 10 kA is at
or above Ic. The neighbouring turn adds ~0.07 T/kA in-plane, which the table barely penalises.
Recommended 6 kA; **Christian kept 10 kA** ("if it breaks, it breaks") — the run is a transition
study on the way up the 10 s sigmoid, not an operating point.

## Input contract

`doc/input_schema.yaml` said `periodic { source ; target }` take node ids. They are **vertex**
ids: `create_periodic` hands them to `PeriodicityFactory::set_master_plane`, which resolves
`mMesh->vertex( id )->node( 0 )`; a vertex is a gmsh type-15 point element whose id is its
geometry tag (`cl_Mesh_GmshReader.cpp`, "we want to use the GeometryTag as Entity ID"), kept only
if its node is used by another element. Both artifacts corrected on Christian's approval (schema
`id_kind` + `id_kind_note`; reference §`periodic`, which also now states that the gmsh
`$Periodic` section is not read and that the periodic faces must be sidesets of their own,
per `tag_periodic_sidesets`). `check_doc_claims.py` 37/37 after the edit. Codex language sweep
of the touched reference paragraph owed.

## Status

Reviewed, not verified: no BELFEM binary run this session. Deck edits verified against the mesh
by a node-chain count only. Regeneration of the mesh with the fixed generator and the first run
are Christian's.
