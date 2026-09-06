"""Re-deriving C++ enum value lists, so the schema's copies stay honest.

Several deck values are closed enums resolved like this:

    string tString = string_to_lower( aString );
    for ( uint k = 0; k < static_cast< uint >( Enum::UNDEFINED ); ++k )
        if ( string_to_lower( to_string( static_cast< Enum >( k ) ) ) == tString )
            return static_cast< Enum >( k );
    BELFEM_ERROR( false, "unknown ..." );

so the accepted spellings are exactly the `to_string` outputs, lowercased. That
makes the `to_string` switch the authority, and it is machine-readable — which
is the whole reason this exists. Hand-transcribing those lists into the schema
would create one more copy to rot, and this repository has already measured how
that ends.

A schema key opts in by naming its enum:

    krylov method:
      values: [preonly, bcgs, ...]
      enum_ref: KrylovMethod
"""

from __future__ import annotations

import re
from pathlib import Path

# Sentinels that exist to terminate the search loop, not as deck values.
NOT_A_VALUE = {"undefined", "unknown"}

_TO_STRING = re.compile(
    r"to_string\(\s*const\s+(\w+)\s+\w+\s*\)\s*\{(.*?)\n    \}", re.S
)
_RETURNS = re.compile(r'return\s+"([^"]+)"')


def enum_tables(root: Path) -> dict[str, list[str]]:
    """Every `to_string( SomeEnum )` in src/, as {EnumName: [values]}.

    Values are lowercased because every lookup lowercases both sides; SolverType
    spells its enumerators `STRUMPACK` while decks write `strumpack`.
    """
    tables: dict[str, list[str]] = {}
    for path in sorted((root / "src").rglob("en_*.cpp")):
        text = path.read_text(errors="replace")
        for match in _TO_STRING.finditer(text):
            name, body = match.group(1), match.group(2)
            values = [
                v.lower() for v in _RETURNS.findall(body)
                if v.lower() not in NOT_A_VALUE
            ]
            if values:
                # a duplicate spelling is legal (PETSc maps ASM and GASM both
                # to "asm"); keep first-seen order and drop the repeat
                tables.setdefault(name, list(dict.fromkeys(values)))
    return tables


def collect_refs(node, path: str = "") -> list[tuple[str, str, list[str]]]:
    """Find schema entries carrying `enum_ref`.

    Returns [(yaml_path, enum_name, declared_values), ...].
    """
    found: list[tuple[str, str, list[str]]] = []
    if isinstance(node, dict):
        ref = node.get("enum_ref")
        if isinstance(ref, str):
            declared = node.get("values")
            found.append(
                (path, ref, [str(v).lower() for v in declared]
                 if isinstance(declared, list) else [])
            )
        for key, value in node.items():
            found.extend(collect_refs(value, f"{path}.{key}" if path else str(key)))
    elif isinstance(node, list):
        for i, item in enumerate(node):
            found.extend(collect_refs(item, f"{path}[{i}]"))
    return found


def check(schema_data: dict, root: Path) -> list[tuple[str, str]]:
    """Compare declared value lists against the C++. Returns [(where, problem)]."""
    tables = enum_tables(root)
    problems: list[tuple[str, str]] = []

    for where, enum_name, declared in collect_refs(schema_data):
        actual = tables.get(enum_name)
        if actual is None:
            problems.append(
                (where, f"enum_ref {enum_name!r} matches no to_string() in src/")
            )
            continue
        if not declared:
            problems.append(
                (where, f"enum_ref {enum_name!r} declared but no `values:` to check")
            )
            continue

        missing = [v for v in actual if v not in declared]
        extra = [v for v in declared if v not in actual]
        if missing:
            problems.append(
                (where, f"{enum_name} accepts {missing} which the schema omits")
            )
        if extra:
            problems.append(
                (where, f"schema lists {extra} which {enum_name} does not accept")
            )
    return problems
