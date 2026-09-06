#!/usr/bin/env python3
"""
Sweep the C++ tree for direct third-party calls that bypass their dedicated wrapper.

L-21 (`doc/lessons_learned.md`) is hard policy: third-party libraries are never called
directly; every access goes through the wrapper layer that carries this tree's contracts —
`comm_type<T>`'s 32/64-bit `int_t` resolution, `comm_abort`'s MPI lifecycle guards, the
chunking and tag discipline in `commtools.hpp`.  A direct vendor call starts from zero of
that knowledge, every time.

The recipe the card originally shipped could not enforce it:

    grep -rn 'MPI_' src/ --include='*.cpp' --include='*.hpp' | grep -v src/comm/

On 2026-08-31 that produced fourteen hits, **none** of them a violating call — comments, a
serial-build `typedef int MPI_Comm`, `MPI_Comm` parameters and members, a Blaze configuration
macro, a PETSc error enumerator — while being structurally blind to the thirty-three `MPI_`
lines in `src/sparse/*.f90`.  A check that reports only noise trains its reader to ignore it.

This script matches **calls** — a vendor-prefixed identifier immediately followed by an open
parenthesis, after comments and string literals are removed — and it names every exclusion it
applies rather than hiding one inside a `grep -v`.

Scope note: it finds direct *calls*.  A vendor *constant* (`MPI_SUM`, `H5T_NATIVE_DOUBLE`)
passed across a module boundary is also a leak of the token the wrapper exists to hide, but it
is not syntactically distinguishable from a legitimate mention, so it is out of scope here and
belongs to the review round.

    python3 scripts/check_wrapper_policy.py [--verbose] [--root DIR]

Exit code 0 if no violation, 1 otherwise.
"""

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# The vendor families, their call syntax, and the module that owns each.
#
# Owners are directories relative to the scan root.  A call inside its owner is
# the wrapper doing its job; a call outside is the finding.
# ---------------------------------------------------------------------------

FAMILIES = [
    {
        "name": "MPI",
        # MPI_Init(, MPI_Comm_rank(, MPI_ALLREDUCE( ... but never the MPI_Comm
        # type, the MPI_COMM_WORLD constant, or BLAZE_MPI_PARALLEL_MODE
        "pattern": re.compile(r"\bMPI_[A-Za-z][A-Za-z0-9_]*\s*\("),
        "owners": ["comm"],
        "wrapper": "commtools.hpp / cl_Communicator.hpp",
    },
    {
        "name": "HDF5",
        "pattern": re.compile(r"\bH5[A-Z][A-Za-z0-9_]*\s*\("),
        "owners": ["io", "physics/materials"],
        "wrapper": "the io module",
        "owner_notes": {
            "physics/materials": "ruled by design (Christian, 2026-08-31): the database currency "
                                 "probes open a file to test for a SIBLING group of the one the "
                                 "wrapper opened, and the wrapper exposes no root handle to probe "
                                 "with -- see cl_JcFunction_Database.hpp and "
                                 "fn_rho_database_is_current.hpp, each of which says so in place",
        },
    },
    {
        "name": "PETSc",
        "pattern": re.compile(r"\b(?:KSP|PC|Mat|Vec|PetscOptions|PetscInitialize|PetscFinalize)"
                              r"[A-Z][A-Za-z0-9_]*\s*\("),
        "owners": ["sparse", "comm"],
        "wrapper": "the Solver* classes",
        "owner_notes": {
            "comm": "ruled by design (Christian, 2026-08-31): PETSc's requirements are unusual "
                    "enough that its scope reaches into the communicator -- PetscInitialize is "
                    "paired with MPI_Init in Communicator::init and cannot live behind Solver*",
        },
    },
    {
        "name": "MUMPS",
        "pattern": re.compile(r"\b[sdcz]mumps_c\s*\("),
        "owners": ["sparse"],
        "wrapper": "the Solver* classes",
        "expect_empty": "MUMPS is driven from mumpstools.f90; no direct C++ call is expected",
    },
    {
        "name": "STRUMPACK",
        # a C++ namespace API, not a C prefix: strumpack::StrumpackSparseSolver< ... >( ),
        # strumpack::ReturnCode::SUCCESS ( no parenthesis, correctly not a call )
        "pattern": re.compile(r"\bstrumpack::[A-Za-z_][A-Za-z0-9_]*\s*(?:<[^;{}]*>)?\s*\("),
        "owners": ["sparse"],
        "wrapper": "the Solver* classes",
    },
    {
        "name": "SuperLU",
        "pattern": re.compile(r"\b(?:[sdcz]gssvx?|Destroy_(?:SuperMatrix|CompCol|SuperNode)"
                              r"[A-Za-z0-9_]*|set_default_options|StatInit|StatFree)\s*\("),
        "owners": ["sparse"],
        "wrapper": "the Solver* classes",
    },
    {
        "name": "LAPACK/BLAS",
        # Fortran-linkage symbols: dgesvd_(, zgetrf_(, dgemm_( ...
        "pattern": re.compile(r"\b[sdcz][a-z][a-z0-9]{1,5}_\s*\("),
        "owners": ["linalg"],
        "wrapper": "the fn_* LAPACK interface",
    },
    {
        "name": "ARPACK",
        "pattern": re.compile(r"\b(?:p?[sdcz](?:na|sa|ne|se)upd)_?\s*\("),
        "owners": ["sparse"],
        "wrapper": "the eigensolver classes",
        "expect_empty": "ARPACK is driven from arpacktools.f90 / parpacktools.f90; "
                        "no direct C++ call is expected",
    },
]

CXX_SUFFIXES = {".cpp", ".hpp", ".h", ".c", ".cc", ".inl", ".tpp"}
FORTRAN_SUFFIXES = {".f90", ".f", ".F90", ".F"}

# ---------------------------------------------------------------------------
# Named exclusions.  Every one prints; none is silent.
# ---------------------------------------------------------------------------

EXCLUSIONS = [
    {
        "name": "Fortran drivers",
        "test": lambda rel, path: path.suffix in FORTRAN_SUFFIXES,
        "reason": "Christian, 2026-08-31: the wrapper requirement is specific to C++.  The "
                  "Fortran drivers for MUMPS, PARDISO and ARPACK ARE the wrapper the C++ tree "
                  "uses, and their communication stays inside the vendor package's own "
                  "collectives (MUMPS's internal MPI traffic is confined to MUMPS).  The "
                  "contracts differ too: MPI_INTEGER, not comm_type<T>.",
    },
    {
        "name": "third-party sources vendored into the tree",
        "test": lambda rel, path: rel.parts[0] in {"extern", "thirdparty"},
        "reason": "not our code",
    },
]


def strip_cxx(text):
    """Blank out comments and string literals, preserving line structure."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                out.append(" ")
                i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            while i < n and not (text[i] == "*" and i + 1 < n and text[i + 1] == "/"):
                out.append("\n" if text[i] == "\n" else " ")
                i += 1
            out.append("  ")
            i = min(i + 2, n)
        elif c in "\"'":
            quote = c
            out.append(" ")
            i += 1
            while i < n and text[i] != quote:
                if text[i] == "\\":
                    out.append(" ")
                    i += 1
                    if i < n:
                        out.append(" ")
                        i += 1
                    continue
                out.append("\n" if text[i] == "\n" else " ")
                i += 1
            out.append(" ")
            i += 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def owned_by(rel, owners):
    """True when rel sits under one of the owner directories.

    Owners may be nested ("physics/materials"), so this is a path-prefix test, not a
    first-component test.
    """
    posix = rel.as_posix()
    return any(posix == owner or posix.startswith(owner + "/") for owner in owners)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--root", default=str(REPO / "src"),
                    help="directory to scan (default: src/)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"no such directory: {root}", file=sys.stderr)
        return 2

    findings = []
    skipped = {e["name"]: [] for e in EXCLUSIONS}
    scanned = 0
    owner_hits = {f["name"]: 0 for f in FAMILIES}

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root)

        excluded = next((e for e in EXCLUSIONS if e["test"](rel, path)), None)
        if excluded is not None:
            if path.suffix in CXX_SUFFIXES or path.suffix in FORTRAN_SUFFIXES:
                skipped[excluded["name"]].append(rel)
            continue

        if path.suffix not in CXX_SUFFIXES:
            continue

        scanned += 1
        try:
            raw = path.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            print(f"cannot read {rel}: {exc}", file=sys.stderr)
            return 2
        code = strip_cxx(raw)

        for lineno, (line, rawline) in enumerate(zip(code.splitlines(),
                                                     raw.splitlines()), start=1):
            for fam in FAMILIES:
                # every match on the line, not just the first: two violating
                # calls on one line are two work items
                for m in fam["pattern"].finditer(line):
                    if owned_by(rel, fam["owners"]):
                        owner_hits[fam["name"]] += 1
                        continue
                    findings.append((rel, lineno, fam, m.group(0).rstrip("( \t"),
                                     rawline.strip()))

    # ---- report ----------------------------------------------------------
    print(f"wrapper-policy sweep (L-21) over {root}")
    print(f"  {scanned} C/C++ translation units and headers scanned\n")

    print("  excluded, by name:")
    for exc in EXCLUSIONS:
        files = skipped[exc["name"]]
        if not files:
            print(f"    - {exc['name']}: no files matched")
            continue
        print(f"    - {exc['name']}: {len(files)} file(s)")
        for rel in files:
            print(f"        {rel}")
        print(f"      reason: {exc['reason']}")
    print()

    # A pattern that matches nothing anywhere is a vacuous check that reads as a clean
    # tree.  Say so, unless the family is known to be driven from Fortran.
    print("  in-wrapper calls (the wrapper doing its job, not findings):")
    for fam in FAMILIES:
        hits = owner_hits[fam["name"]]
        found = any(f[2] is fam for f in findings)
        line = (f"    - {fam['name']:<12} {hits:>4} call(s) "
                f"inside {', '.join(o + '/' for o in fam['owners'])}")
        if hits == 0 and not found:
            note = fam.get("expect_empty")
            line += f"  — {note}" if note else "  — WARNING: this pattern matched nothing " \
                                               "anywhere; verify it before trusting a green run"
        print(line)
        for owner, note in fam.get("owner_notes", {}).items():
            print(f"        {owner}/ is an owner by ruling, not by module: {note}")
    print()

    if not findings:
        print("no direct third-party calls outside their wrapper module.")
        return 0

    print(f"{len(findings)} direct call(s) outside the owning wrapper module:\n")
    for rel, lineno, fam, token, text in findings:
        print(f"  src/{rel}:{lineno}: {token} — {fam['name']} belongs behind "
              f"{fam['wrapper']}")
        if args.verbose:
            print(f"      {text}")
    print("\nL-21: extend the wrapper, then call it. See doc/lessons_learned.md.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
