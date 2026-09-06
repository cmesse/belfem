# Inert coil blocks kept out of the magnetic equation (costheta repaired)

**Date:** 2026-08-27
**Topic:** commit `4adf3f9e` — the equation-side block selection now skips
`DomainType::Coil`, which is what made `examples/costheta` runnable from its own mesh.

## The defect

Coil blocks are current sources, not material regions (no LTS model in BELFEM), so they
deliberately carry no material — and, as Christian pointed out mid-investigation, no dofs
either: `Cell< string > Coil` in the Maxwell FieldList has never been populated in the
tree's history, so the coil doftable was always empty.

Two block-selection paths disagreed about what to do with that. The kernel-parameter path
(`create_magnetic_kernel`) already excluded coils; the equation-side path
(`set_block_types_in_magnetic_equation`) handed every mesh block to the magnetic IWG — and
the IWG's selected list is the one `BlockData::create_blocks` and the Group constructor
actually build FEM blocks and Calculators from. A coil block therefore got a Calculator,
and `MaxwellData` dereferenced its null material. On `examples/costheta` — the only shipped
deck with a `topology { coil }` section — that was a segfault (later a named error, once
Christian's null-material guard landed in the working tree), and the reason the deck could
not run from its mesh.

## The fix and its fences

One filter in `set_block_types_in_magnetic_equation`, skipping `DomainType::Coil` only.
Both auditors independently insisted on the "only": the two selection lists are
deliberately not identical — coating blocks belong on the equation side despite being
absent from the kernel-parameter filter, with their own dof mapping and side-connector
dispatch. Current imposition is untouched by construction: terminal ids reach the
cohomology as mesh-block lookups (`Homology::suggest_Homology` → `mMesh->block(id)`) and
the current is fixed onto abstract cut dofs (`IWG_Maxwell::set_currents`), never resolved
through the selected block list.

Two of my own claims were corrected in audit: the equation list is built on rank 0 and
broadcast by `IWG_Maxwell::initialize`, not built per rank; and the comment was reworded to
stop implying the two lists should agree.

## Gates

- costheta from its `.msh`: runs clean to −84 dB, and the imposed 169.96 A produces
  |B| = 0.34 T at the coils and 83 mT in the aperture — the current arrives, not silently
  dropped. Christian confirmed the full run on his side.
- gantry (coil-free): residual sequence bit-identical to the same-mesh pre-filter
  reference — the filter is a proven no-op where it must be one.
- `make check-fast` 9/9.

## Also documented this session

`doc/input_file_reference.md` §8 now states that coil and air blocks are magnetically
inert by design and why (commit `132cbf24`), after the missing-material complaint on a
coil block briefly sent this investigation toward "fix the deck" — the deck was right.

## Residues

- Dead code, deliberately not removed here: `Cell< string > Coil` and the
  `case DomainType::Coil` in `collect_block_dofs` — future hygiene bundle.
- `examples/costheta` still lacks `bhdata.hdf5` next to its deck; the canonical copy is
  `share/material/bhdata.hdf5`. Provisioning is an Allrun/packaging question, not physics.
- costheta being healthy again is input to the DR-65 fix-or-drop decision on the freeze
  list: the shipped-examples story is one deck stronger.

Round record: `tmp/ai_exchange/coil_block_exclusion.md`.
