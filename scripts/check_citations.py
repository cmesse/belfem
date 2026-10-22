#!/usr/bin/env python3
"""Check that every author-year citation in the tree resolves to a known work.

Two tables are authoritative:

  doc/literature_references.md                             citation key index (general literature)
  src/physics/materials/doc/material_property_sources.md   bibliography (property data provenance)

The index is delimited so it survives edits to the surrounding prose:

    <!-- BEGIN CITATION INDEX -->
    | Key | Cite as | Reference | DOI / ISBN |
    ...
    <!-- END CITATION INDEX -->

Scanned: src/ sources and module docs, doc/. Not scanned: devlog/ and todo/ (records, kept as
written), tmp/, literature/, nonfree/, archive/.

Usage:
    scripts/check_citations.py [--verbose] [--list-keys]

Exit status 1 if an unresolved citation is found.
"""

import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INDEX_FILE = "doc/literature_references.md"
MATERIALS_FILE = "src/physics/materials/doc/material_property_sources.md"

INDEX_BEGIN = "<!-- BEGIN CITATION INDEX -->"
INDEX_END = "<!-- END CITATION INDEX -->"

SCAN_DIRS = ("src", "doc")
SCAN_EXT = (".hpp", ".cpp", ".h", ".c", ".f90", ".md")

SKIP_FILES = (INDEX_FILE, MATERIALS_FILE)

SKIP_DIRS = {"devlog", "todo", "tmp", "literature", "nonfree", "archive", ".git",
             "cmake-build-debug", "build", "bin", "lib"}

# A citation token: one or two capitalized surnames, optional "et al.", then a four-digit year.
# Surnames may carry accents, hyphens and apostrophes (Luethi, Strom-Olsen, O'Neal, Agelas).
SURNAME = r"[A-ZÀ-Þ][A-Za-zÀ-ÿ'’.-]{2,}"
CITATION = re.compile(
    r"\b(?P<first>" + SURNAME + r")"
    r"(?:\s*(?:&|and)\s*(?P<second>" + SURNAME + r"))?"
    r"(?P<etal>\s+et\s+al\.?)?"
    r",?\s+\(?(?P<year>1[89]\d\d|20\d\d)(?P<letter>[ab])?\b(?!-\d)"
)

DOI = re.compile(r"\b10\.\d{4,9}/[^\s,;)\"'\]]+")

# Capitalized words that precede a year without being a citation. Extend only for genuine
# false positives -- a real work belongs in the index, not here.
NOT_A_SURNAME = set("""
    Page Figure Table Section Chapter Equation Version Release Appendix Volume Vol No Ed Edition
    Series Ser Part Chap Note Draft Rev Commit Branch Tag Issue Item Step Tier Phase Round
    January February March April May June July August September October November December
    Jan Feb Mar Apr Jun Jul Aug Sep Sept Oct Nov Dec
    Monday Tuesday Wednesday Thursday Friday Saturday Sunday
    Copyright ISBN ISSN DOI HAL LCCN
    The This That These Those Since Until From Before After When While Their There Here Between
    Standard Update Updated Date See Added Removed Fixed Changed Committed Merged Found Corrected
    CORRECTED Addendum Devlog Session Status Summary Result Verdict Outcome Decision Ruling
    Springer Wiley Elsevier Pearson Dover Prentice Cambridge Oxford Academic CRC SIAM IEEE ACM
    Plenum Pergamon Woodhead Butterworth Heinemann Addison Wesley McGraw Hill
    Press Publishing Publishers Verlag Books Ltd Inc GmbH
    Addison-Wesley Butterworth-Heinemann McGraw-Hill Wiley-VCH Prentice-Hall
    Comput Compute Phys Appl Supercond Numer Math Sci Technol Trans Proc Acta Conf Cryogenics
    Physica JETP Pramana Mater Chem Ref Geom Magn Suppl Lett Eng
    CLAUDE CODEX Claude Codex Grok GPT CLion CppCon Gmsh COMSOL GetDP MFEM PETSc MUMPS STRUMPACK
    SuperLU PARDISO METIS HDF5 Exodus ParaView Valgrind GitLab GitHub Linux
    EPFL LBNL LBL MIT KIT NIST NBS ORNL FNAL TPRC HTS LTS REBCO YBCO BSCCO CORC
    Boulder Berkeley Montreal Torino Karlsruhe Watertown Zurich Zürich
    SODA SoCG SUST TASC BIT TOMS FVCA ICCS NBSIR CBMS
    Temperatures Properties Matter Handbook Procedures Equations Elements Methods
    Christian Gregory Sirous Nicolo Emily
    III IV VI VII VIII IX XI XII
""".split())

def _bibliography_aliases(key):
    """The materials bibliography keys composite works as "Kamm-Alers 1964".

    Citations in the tree write those as "Kamm & Alers 1964" or just "Kamm 1964", so both
    spellings are registered alongside the table's own form.
    """
    match = re.match(r"^(.+?)\s+((?:1[89]|20)\d\d)$", key)
    if not match:
        return []
    names, year = match.group(1), match.group(2)
    parts = [part for part in names.split("-") if part]
    aliases = ["%s %s" % (part, year) for part in parts]
    for left, right in zip(parts, parts[1:]):
        aliases.append("%s & %s %s" % (left, right, year))
    if len(parts) > 1:
        aliases.append("%s & %s %s" % (parts[0], parts[-1], year))
    return aliases


def _cells(line):
    parts = [c.strip() for c in line.strip().strip("|").split("|")]
    return parts


def load_index(verbose=False):
    """Return {key: reference} from both authoritative tables."""
    keys = {}

    path = os.path.join(ROOT, INDEX_FILE)
    if not os.path.exists(path):
        sys.stderr.write("error: %s not found\n" % INDEX_FILE)
        sys.exit(2)

    with open(path, encoding="utf-8") as fh:
        text = fh.read()

    if INDEX_BEGIN not in text or INDEX_END not in text:
        sys.stderr.write(
            "error: %s carries no citation index (missing %s / %s)\n"
            % (INDEX_FILE, INDEX_BEGIN, INDEX_END))
        sys.exit(2)

    block = text.split(INDEX_BEGIN, 1)[1].split(INDEX_END, 1)[0]
    for line in block.splitlines():
        if not line.strip().startswith("|"):
            continue
        cells = _cells(line)
        if len(cells) < 3 or cells[0] in ("Key", "") or set(cells[0]) <= set("-: "):
            continue
        # The key cell may carry several accepted spellings, separated by ";"
        for key in cells[0].split(";"):
            key = key.strip().strip("`")
            if key and key != "-":
                keys[key] = cells[2]

    mat = os.path.join(ROOT, MATERIALS_FILE)
    if os.path.exists(mat):
        with open(mat, encoding="utf-8") as fh:
            in_bib = False
            for line in fh:
                if line.startswith("## "):
                    in_bib = line.strip().lower().startswith("## bibliography")
                    continue
                if not in_bib or not line.strip().startswith("|"):
                    continue
                cells = _cells(line)
                if len(cells) < 2 or cells[0] in ("Key", "") or set(cells[0]) <= set("-: "):
                    continue
                raw = cells[0].strip("`")
                keys.setdefault(raw, MATERIALS_FILE)
                for alias in _bibliography_aliases(raw):
                    keys.setdefault(alias, MATERIALS_FILE)
                # the Authors cell names co-authors the key leaves out
                year = re.search(r"((?:1[89]|20)\d\d)$", raw)
                if year and len(cells) > 1:
                    for surname in re.findall(SURNAME, cells[1]):
                        if len(surname.rstrip(".")) > 2 and "." not in surname:
                            keys.setdefault("%s %s" % (surname, year.group(1)),
                                            MATERIALS_FILE)
    elif verbose:
        print("note: %s not found, property-data keys unavailable" % MATERIALS_FILE)

    return keys


def _clean(name):
    """Strip a possessive and trailing punctuation from a matched surname."""
    name = name.rstrip(".")
    if name.endswith("'s") or name.endswith("\u2019s"):
        name = name[:-2]
    return name


def _is_surname(token):
    """Initials such as "J.-C" are not surnames: require three alphabetic characters."""
    return sum(ch.isalpha() for ch in token) >= 3


def normalize(match):
    """Reduce a matched citation to the key form used in the index: 'Surname YYYY'."""
    first = _clean(match.group("first"))
    second = match.group("second")
    if second:
        second = _clean(second)
    year = match.group("year")
    letter = match.group("letter") or ""
    if second and not match.group("etal"):
        return "%s & %s %s%s" % (first, second, year, letter)
    return "%s %s%s" % (first, year, letter)


def scan_files():
    for top in SCAN_DIRS:
        base = os.path.join(ROOT, top)
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            for name in sorted(filenames):
                if not name.endswith(SCAN_EXT):
                    continue
                path = os.path.join(dirpath, name)
                if os.path.relpath(path, ROOT) in SKIP_FILES:
                    continue
                yield path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true", help="report resolved citations too")
    ap.add_argument("--list-keys", action="store_true", help="print the known keys and exit")
    args = ap.parse_args()

    keys = load_index(args.verbose)

    if args.list_keys:
        for key in sorted(keys):
            print(key)
        return 0

    unresolved = {}
    resolved = 0
    bad_doi = []

    for path in scan_files():
        rel = os.path.relpath(path, ROOT)
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                lines = fh.readlines()
        except OSError:
            continue

        in_fence = False
        is_markdown = rel.endswith(".md")

        for lineno, line in enumerate(lines, 1):
            # A fenced block in a document quotes a data file, a deck or program output.
            # Whatever looks like a citation in there belongs to the quoted thing, not to us.
            if is_markdown and line.lstrip().startswith("```"):
                in_fence = not in_fence
                continue
            if in_fence:
                continue

            for match in CITATION.finditer(line):
                if _clean(match.group("first")) in NOT_A_SURNAME:
                    continue
                if not _is_surname(match.group("first")):
                    continue
                if match.group("second") and match.group("second").rstrip(".") in NOT_A_SURNAME:
                    continue
                prefix = line[max(0, match.start() - 24):match.start()]
                if "http" in prefix or "10." in prefix or "doi" in prefix.lower():
                    continue
                key = normalize(match)
                if key in keys:
                    resolved += 1
                    continue
                # A lettered citation may be indexed without the letter
                if key[-1] in "ab" and key[:-1] in keys:
                    resolved += 1
                    continue
                # A letterless citation where only lettered works exist is the a/b trap,
                # and is reported rather than silently resolved.
                hint = ""
                if (key + "a") in keys or (key + "b") in keys:
                    hint = "  <- ambiguous, needs the a/b letter"
                unresolved.setdefault(key + hint, []).append("%s:%d" % (rel, lineno))

            for doi in DOI.finditer(line):
                value = doi.group(0).rstrip(".,);")
                if len(value) < 10 or "/" not in value:
                    bad_doi.append("%s:%d  %s" % (rel, lineno, value))

    print("citations resolved : %d" % resolved)
    print("distinct unresolved: %d" % len(unresolved))

    if bad_doi:
        print("\nmalformed DOI:")
        for entry in bad_doi:
            print("  " + entry)

    if unresolved:
        print("\nunresolved citations (add to the citation key index, or to the materials")
        print("bibliography if the work is property data):\n")
        for key in sorted(unresolved):
            sites = unresolved[key]
            shown = ", ".join(sites[:4])
            more = "" if len(sites) <= 4 else "  (+%d more)" % (len(sites) - 4)
            print("  %-34s %s%s" % (key, shown, more))
        return 1

    print("\nall citations resolve")
    return 0


if __name__ == "__main__":
    sys.exit(main())
