"""Locating the repository and its sources.

Everything here is deliberately path-based rather than build-based: the drift
check must run without a configured build tree, in a git hook, or on a machine
that has never compiled BELFEM.
"""

from __future__ import annotations

import os
import re
import sys
from functools import lru_cache
from pathlib import Path

# Files that read `input.conf` through input::Section. The code -> schema check
# scans exactly these: widening it to all of src/ would sweep up XML, gas-table
# and HDF5 readers, whose string keys are not input-deck keys.
#
# This list is transcribed, which is the failure mode the tool exists to kill —
# a new factory that reads the deck is invisible until someone adds it here.
# There is no mechanical way to derive it (any file may take an
# input::Section*), so keep it reviewed: candidates are grep hits for
# `input::Section` in src/. cl_FEM_Kernel.cpp and cl_ThinShellFactory.cpp were
# removed 2026-08-13 — their key_exists calls are on Maps, not deck sections.
CONSUMERS = (
    "src/io/cl_Input_Section.cpp",
    "src/fem/maxwell/cl_MaxwellFactory.cpp",
    "src/fem/maxwell/cl_MaxwellBoundaryConditionFactory.cpp",
    "src/fem/maxwell/fn_mesh_config_tag.hpp",
    "src/fem/kernel/cl_FEM_Controller.cpp",
    "src/fem/kernel/cl_FEM_Domain.cpp",
    "src/fem/thermal/cl_ThermalFactory.cpp",
    "src/fem/thermal/cl_ThermalBoundaryConditionFactory.cpp",
    "src/physics/materials/cl_MaterialFactory.cpp",
    "src/sparse/cl_SolverParameters.cpp",
    "src/circuit/cl_ElectricalCircuitFactory.cpp",
    "src/circuit/cl_Component.cpp",
    "src/circuit/cl_Resistor.cpp",
    "src/circuit/cl_Capacitor.cpp",
    "src/circuit/cl_Inductor.cpp",
    "src/circuit/cl_CurrentSource.cpp",
    "src/circuit/cl_VoltageSource.cpp",
    "src/homology/cl_Topology.cpp",
)

SCHEMA = "doc/input_schema.yaml"
REFERENCE = "doc/input_file_reference.md"


def find_root(start: Path | None = None) -> Path:
    """Walk up until we see the schema.

    The working directory is tried BEFORE the tool's own location: with two
    checkouts on disk, running A's copy while standing in B must validate B,
    not silently fall back to A. Exit status 2, the documented bad-usage code.
    """
    starts = [start] if start is not None else [Path.cwd(), Path(__file__)]
    for origin in starts:
        here = origin.resolve()
        for candidate in [here, *here.parents]:
            if (candidate / SCHEMA).is_file():
                return candidate
    print(
        f"belfem-conf: not inside a BELFEM checkout (no {SCHEMA} above "
        f"{' or '.join(str(s) for s in starts)})", file=sys.stderr)
    raise SystemExit(2)


class Sources:
    """A lazily-built index of the C++ sources, with whitespace-blind search.

    Whitespace normalisation is not a convenience: call sites genuinely vary
    between `key_exists("nodes")` and `key_exists( "nodes" )`, so a literal
    search reports false misses. Comment lines are dropped, because a key named
    only in a comment is not a parse site.
    """

    def __init__(self, root: Path):
        self.root = root
        self._by_name: dict[str, Path] = {}

    @lru_cache(maxsize=None)
    def lines(self, path: Path) -> tuple[tuple[int, str], ...]:
        out = []
        for n, line in enumerate(
            path.read_text(errors="replace").splitlines(), start=1
        ):
            stripped = line.lstrip()
            if stripped.startswith("//") or stripped.startswith("*"):
                continue
            out.append((n, line))
        return tuple(out)

    def resolve(self, name: str) -> Path | None:
        """Find a source file by bare name, e.g. 'cl_FEM_Domain.cpp'."""
        if name in self._by_name:
            return self._by_name[name]
        for base, dirs, files in os.walk(self.root / "src"):
            dirs[:] = [d for d in dirs if d != ".git"]
            if name in files:
                self._by_name[name] = Path(base) / name
                return self._by_name[name]
        self._by_name[name] = None
        return None

    def find(self, token: str, in_file: str | None = None) -> list[tuple[str, int]]:
        """Whitespace-blind search. Returns [(basename, line), ...]."""
        needle = re.sub(r"\s+", "", token)
        if not needle:
            return []

        if in_file:
            path = self.resolve(in_file)
            paths = [path] if path else []
        else:
            paths = [p for p in (self.resolve(Path(c).name) for c in CONSUMERS) if p]

        hits = []
        for path in paths:
            for n, line in self.lines(path):
                if needle in re.sub(r"\s+", "", line):
                    hits.append((path.name, n))
        return hits

    def consumer_paths(self) -> list[Path]:
        out = []
        for rel in CONSUMERS:
            path = self.root / rel
            if path.is_file():
                out.append(path)
        return out
