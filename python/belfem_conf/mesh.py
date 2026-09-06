"""Reading just enough of a mesh to check a deck against it.

Deliberately partial. The point is not to load a mesh — it is to answer the
handful of questions a deck raises: does this file exist, is it 2-D or 3-D, and
do the block and sideset ids the deck names actually occur in it.

Two facts shape everything here.

**Ids are ELEMENTARY (geometry) tags, not physical tags.** `GmshReader` sets
both on an element but groups on `geometry_tag()` everywhere — `create_group_ids`,
`create_blocks`, `create_sidesets` — and `read_physical_tag` only advances the
buffer cursor, so `$PhysicalNames` never reaches a block or sideset id. Reading
physical tags here would agree with a deck by luck and disagree by design.

**The mesh is usually not there.** The examples ship `.geo` sources and generate
the `.msh` with gmsh, so five of the six carry no mesh at all. A missing mesh is
therefore normal and must never be an error — it makes a class of check
unavailable, and the validator says so rather than passing in silence.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# gmsh element type -> topological dimension. `element_type_from_gmsh` is a
# bare static_cast, so this must cover every id in Mesh_Enums.hpp that gmsh can
# write — omitting one silently drops its elementary tag from the id pools and
# reports a valid block as missing.
ELEMENT_DIM = {
    1: 1, 8: 1, 26: 1, 27: 1, 28: 1,                       # lines
    2: 2, 3: 2, 9: 2, 10: 2, 16: 2, 21: 2, 23: 2, 25: 2,   # triangles / quads
    32: 2,
    4: 3, 5: 3, 6: 3, 7: 3, 11: 3, 12: 3, 13: 3, 14: 3,    # tets / hexes /
    17: 3, 18: 3, 19: 3, 29: 3, 30: 3, 92: 3,              # prisms / pyramids
    15: 0,                                                 # points (vertices)
}

ELEMENT_NODES = {
    1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6,
    10: 9, 11: 10, 15: 1, 16: 8, 17: 20,
}


@dataclass
class MeshInfo:
    path: Path
    kind: str                        # "gmsh" | "bfm" | "missing" | "unreadable"
    detail: str = ""
    version: str = ""
    dimension: int | None = None
    block_ids: set[int] = field(default_factory=set)
    sideset_ids: set[int] = field(default_factory=set)
    vertex_ids: set[int] = field(default_factory=set)
    sibling_bfm: Path | None = None

    @property
    def usable(self) -> bool:
        return self.kind == "gmsh" and self.dimension is not None


def resolve(deck_path: Path, filename: str) -> Path:
    """Mesh paths in a deck are relative to the deck, as the solver reads them."""
    return (deck_path.parent / filename.strip()).resolve()


def inspect(deck_path: Path, filename: str) -> MeshInfo:
    path = resolve(deck_path, filename)

    # a gmsh file may be superseded by a sibling .bfm that the factory prefers
    sibling = path.with_suffix(".bfm")
    sibling = sibling if sibling.is_file() and sibling != path else None

    if not path.is_file():
        return MeshInfo(path, "missing", sibling_bfm=sibling)

    if path.suffix.lower() == ".bfm":
        # HDF5. Reading it would need h5py, which this package deliberately
        # does not depend on: the checks must run in a hook and on a machine
        # that has never built BELFEM. Recorded as known-but-unread.
        return MeshInfo(path, "bfm",
                        detail="enriched mesh; not read (HDF5 needs h5py)",
                        sibling_bfm=sibling)

    try:
        return _read_gmsh(path, sibling)
    except (OSError, ValueError, IndexError) as exc:
        # IndexError included: a truncated section indexes past the line list,
        # and that must surface as "unreadable", not as a traceback
        return MeshInfo(path, "unreadable", detail=str(exc), sibling_bfm=sibling)


def _read_gmsh(path: Path, sibling: Path | None) -> MeshInfo:
    lines = path.read_text(errors="replace").splitlines()

    version, filetype = "", ""
    for i, line in enumerate(lines):
        if line.strip() == "$MeshFormat" and i + 1 < len(lines):
            words = lines[i + 1].split()
            version = words[0] if words else ""
            filetype = words[1] if len(words) > 1 else ""
            break

    info = MeshInfo(path, "gmsh", version=version, sibling_bfm=sibling)

    if filetype == "1":
        # binary payload after an ASCII header; parsing the text lines would
        # produce garbage coordinates and ids that LOOK checked
        info.detail = f"binary .msh not parsed (version {version})"
        return info

    # GmshReader accepts exactly 2.2 and 4.1 and hard-errors on anything else,
    # so a wider net here would happily validate a mesh the solver refuses.
    if version == "2.2":
        _read_nodes(lines, info)
        _read_elements(lines, info)
        return info

    if version == "4.1":
        # Not an optional extra: gmsh has defaulted to 4.1 since 4.0, so every
        # freshly generated mesh is 4.1 and a 2.2-only reader checks nothing.
        _read_nodes_41(lines, info)
        _read_elements_41(lines, info)
        return info

    info.detail = f"format {version} not parsed (the solver accepts only 2.2 and 4.1)"
    return info


def _read_nodes(lines: list[str], info: MeshInfo) -> None:
    """Dimension from the z extent, which is how GmshReader decides it.

    Not from topology: a curved 2-D surface embedded in 3-D reads as 3-D there,
    so matching the rule matters more than being right in the abstract.
    """
    try:
        start = lines.index("$Nodes")
    except ValueError:
        return
    count = int(lines[start + 1])
    zs = []
    for k in range(count):
        parts = lines[start + 2 + k].split()
        if len(parts) >= 4:
            zs.append(float(parts[3]))
    if zs:
        # exact comparison, because GmshReader tests min(tZ) == max(tZ); a
        # tolerance here would call 2-D what the solver treats as 3-D
        info.dimension = 2 if min(zs) == max(zs) else 3


def _read_elements(lines: list[str], info: MeshInfo) -> None:
    try:
        start = lines.index("$Elements")
    except ValueError:
        return
    count = int(lines[start + 1])

    by_dim: dict[int, set[int]] = {0: set(), 1: set(), 2: set(), 3: set()}
    for k in range(count):
        parts = lines[start + 2 + k].split()
        if len(parts) < 3:
            continue
        etype, ntags = int(parts[1]), int(parts[2])
        if ntags < 2:
            continue
        dim = ELEMENT_DIM.get(etype)
        if dim is None:
            continue
        # tags are [physical, elementary, ...]; the SECOND is what groups
        by_dim[dim].add(int(parts[4]))

    info.vertex_ids = by_dim[0]
    if info.dimension == 2:
        info.block_ids, info.sideset_ids = by_dim[2], by_dim[1]
    else:
        info.block_ids, info.sideset_ids = by_dim[3], by_dim[2]


def _read_nodes_41(lines: list[str], info: MeshInfo) -> None:
    """4.1 splits nodes into per-entity blocks, tags first then coordinates.

        numEntityBlocks numNodes minTag maxTag
        entityDim entityTag parametric numNodesInBlock
          <numNodesInBlock tag lines>
          <numNodesInBlock coordinate lines>
    """
    try:
        start = lines.index("$Nodes")
    except ValueError:
        return

    header = lines[start + 1].split()
    if len(header) < 1:
        return
    num_blocks = int(header[0])

    zs: list[float] = []
    cursor = start + 2
    for _ in range(num_blocks):
        block = lines[cursor].split()
        if len(block) < 4:
            return
        in_block = int(block[3])
        cursor += 1 + in_block                 # skip the tag lines
        for k in range(in_block):
            parts = lines[cursor + k].split()
            if len(parts) >= 3:
                zs.append(float(parts[2]))
        cursor += in_block

    if zs:
        # exact, matching GmshReader — see _read_nodes
        info.dimension = 2 if min(zs) == max(zs) else 3


def _read_elements_41(lines: list[str], info: MeshInfo) -> None:
    """In 4.1 the elementary tag is the block's entityTag.

        numEntityBlocks numElements minTag maxTag
        entityDim entityTag elementType numElementsInBlock
          <numElementsInBlock element lines>

    `GmshReader` calls `set_geometry_tag( entityTag )` on this path and never
    `set_physical_tag`, so entityTag is exactly the id a deck names.
    """
    try:
        start = lines.index("$Elements")
    except ValueError:
        return

    header = lines[start + 1].split()
    if len(header) < 1:
        return
    num_blocks = int(header[0])

    by_dim: dict[int, set[int]] = {0: set(), 1: set(), 2: set(), 3: set()}
    cursor = start + 2
    for _ in range(num_blocks):
        block = lines[cursor].split()
        if len(block) < 4:
            return
        entity_dim, entity_tag, _etype, in_block = (
            int(block[0]), int(block[1]), int(block[2]), int(block[3]))
        if entity_dim in by_dim:
            by_dim[entity_dim].add(entity_tag)
        cursor += 1 + in_block

    info.vertex_ids = by_dim[0]
    if info.dimension == 2:
        info.block_ids, info.sideset_ids = by_dim[2], by_dim[1]
    else:
        info.block_ids, info.sideset_ids = by_dim[3], by_dim[2]


def expand_ids(raw: str) -> list[int]:
    """Deck id list to integers: commas, whitespace, and `a:b` ranges.

    A range always expands ASCENDING, matching the reader: `9:4` and `4:9` both
    give 4..9. Bracket group markers are separators here, since this is used
    for existence checks and not for grouping.
    """
    out: list[int] = []
    for token in re.split(r"[,\s\[\]]+", raw.strip()):
        if not token:
            continue
        if ":" in token:
            a, _, b = token.partition(":")
            try:
                lo, hi = int(a), int(b)
            except ValueError:
                continue
            out.extend(range(min(lo, hi), max(lo, hi) + 1))
        else:
            try:
                out.append(int(token))
            except ValueError:
                continue
    return out
