# The first-facet sideset campaign leaves devel: e30bd98 and 411e69a reverted on Christian's ruling, cfaf7e8 stays, the Topology commit archived

**Date:** 2026-10-09
**Purpose:** Christian's ruling of 2026-10-09 (evening), in his words: "we do not include this
in 0.9.3 because it tries to fix something that is not broken, and I would even remove these
fixes from devel and archive them." By quiz he set the scope: the unpushed Topology commit
and the two devel commits e30bd98 (admission and set key from all duplicated facets, the
refusals for a multi-block private-copy sideset, a sideset whose facets do not all touch a
coil, and a sideset mixing interface facets with boundary or conductor|conductor facets) and
411e69a (the four ordered passes, the three reasons, the `masters|slaves` key in the docs)
come out of devel; cfaf7e8 (the mixed-master-type refusal and the id-ordered interface copies)
stays. This commit reverts the two on top of the tip, as new history; nothing is rewritten.
**Module:** `src/homology` (`cl_InterfaceProcessor`), `tests/homology`, `src/homology/doc`,
`doc/input_file_reference.md`, `CHANGELOG.md`, `todo/`
**AIs:** Claude Fable 5.1 (DrClean session); the Server session (ruling relay, the
Gmsh-reader argument); Codex `gpt-6-astra` and Grok `grok-4.7` (the juries of the reverted
commits, recorded in the archive)
**Status:** COMPLETE — reverted on devel; `src/homology/cl_InterfaceProcessor.{cpp,hpp}` and
`tests/homology/test_InterfaceProcessor.cpp` are byte-identical to their cfaf7e8 versions
(`git diff cfaf7e8 <this commit> -- <the three files>` is empty), so that code needs no new
audit: it is blind-audited, pipeline-green (pipeline 61) and deck-compared. Gate: private debug
tree, `make check` 17/17. The campaign's work is on the branch
`archive/first-facet-sidesets` (tip 93d8318, on top of 411e69a), never on main
**Archive:** `archive/first-facet-sidesets` — the Topology commit as it stood (classification
from all facets, the boundary rule, the symmetric Ferro|Coil arm, full-pipeline tests), the
round-2 plan, both jury threads, the consumer survey, the sideset-type probe, the per-deck
evidence table and the rule (a) script, under `devlog/archive_first_facet_sidesets/`

## Summary

- **Why it was not broken.** `GmshReader::create_sidesets` builds one sideset per geometric
  entity tag (`cl_Mesh_GmshReader.cpp:878-902`). A Gmsh boundary surface lies on one volume
  and an interior surface between exactly two, so a native Gmsh mesh cannot contain a
  mixed-master, partly-interface, multi-block or mixed-boundary sideset; only
  generator-written or hand-built meshes can. The first-facet classification that the
  campaign set out to make order-independent therefore has nothing to act on in any deck
  BELFEM reads from Gmsh.
- **What the gates showed.** The strict Topology commit (every disagreement refused) ran all
  21 shipped decks through setup with no `.bfm`: 21 of 21 pass, and a probe printing every
  sideset's type after `Topology::run` gives, for every deck, types identical to 6cb9685
  (8 to 2364 sidesets per deck); the first dof banner is identical on all 21. No shipped deck
  has a sideset whose classification depends on facet order.
- **What the consumer survey showed, and why the Topology commit could not ship.** An
  auto-typed antisymmetry sideset that became `Inactive` is seen by three readers:
  `Topology::select_sidesets` (`cl_Topology.cpp:480-491`) puts Air/Buffer/FerroAntiSymmetry
  on the φ-boundary list that the cut processor uses for its pruning exemptions
  (`cl_CutProcessor.cpp:51-63`); `DofData::count_node_dofs` (`cl_FEM_DofMgr_DofData.cpp:946-968`)
  with the field list's per-type dof tables puts the air dofs on every node of such a sideset;
  `CutFactory::compute_poisson_problem` (`cl_CutFactory.cpp:2669-2675`) flags the nodes of
  air (anti)symmetry sidesets. The test_fem `ThinShellAlignment` fixture, whose box faces
  span the tape stack and the air, proved the first: 12 tests died with "Thin cut 0 is empty"
  when its lateral boundary sideset lost the AirAntiSymmetry type its first facet used to give
  it. The readers that do not see the difference: the factory's deactivation of every
  auto-typed antisymmetry group (`cl_MaxwellFactory.cpp:851-859`), the IWG's type switch,
  `FEM_Domain`, the autopin selection and contact set, periodic tagging (runs before
  `Topology::run`), deck re-typing (before, and typed sidesets are skipped), the `.bfm` key.
- **What devel keeps (cfaf7e8).** The mixed-master-domain-type refusal over every facet
  with a slave of a sideset that has a duplicated facet (so conductor|air beside air|air is
  refused as a master-type mix); interface copies allocated set by set in ascending original
  id and appended in ascending copy id; set order the sorted hex of the union bitset of the
  first facet's block pair; the classification of a sideset from its first facet, as in
  0.9.2. The docs say that and nothing more, keeping the corrections that are true of
  cfaf7e8: the periodic backup order (within-set pairs set by set, then the cross-set pass),
  the check running after `CutProcessor` on the sidesets that existed before it, both facets
  named in the message, "a reused `.bfm`", and the numbering scope as it stands until the
  cut-copy commit (interface copies deterministic for a given incoming node order; cut copies
  get deterministic ids but are appended in hash order).
- **What was taken out.** From the processor: the all-facet admission and the
  `masters|slaves` key, the coil-contact, interface-facets-only and one-master-block
  refusals, the four ordered passes and the three reasons, the handed master block, the
  value-member bitsets. From the tests: every test of those rules; the file is cfaf7e8's. From
  the docs: every sentence about them.

## Key Findings

- The campaign's juries (plan and code, two rounds each for e30bd98 and 411e69a, one for the
  Topology commit) found real defects in each step and are preserved in the archive; the
  decisive finding — that the all-facet typing is not inert — came from `make check`, not
  from a jury.
- Residue for after 0.9.3, if Christian wants it: the archive branch, with option (a) (a
  boundary-only sideset typed by its most φ-like block) ready and through the debug gates.

## Changes Made

- `src/homology/cl_InterfaceProcessor.cpp`, `src/homology/cl_InterfaceProcessor.hpp`,
  `tests/homology/test_InterfaceProcessor.cpp` — cfaf7e8 content
- `src/homology/doc/homology_usage_guide.md`, `doc/input_file_reference.md`, `CHANGELOG.md`
  — describe cfaf7e8
- `todo/closed/interface_sideset_order_independence_2.md` (closed by ruling),
  `todo/closed/interface_sideset_first_facet_gaps.md` and
  `todo/closed/interface_copy_refusal_and_ordering.md` (reverted notes), `todo/README.md`,
  this devlog, `devlog/README.md`
