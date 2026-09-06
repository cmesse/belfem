"""Reading `doc/input_schema.yaml`.

The schema is walked as a parsed data structure, never with line regexes. That
is not a style preference — the first prototype of the anchor check matched
`anchor:` with a regex anchored at line start, silently skipped the 60 anchors
written inside inline flow mappings, and reported "all 53 resolve" for a file
that held 115. It passed a check it had never run.

Walking the parsed tree also means findings are reported by YAML path rather
than by line number, which suits a file whose whole design premise is that line
numbers rot.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

# Fields that carry a searchable token into the sources. `anchor` is the main
# one; the rest are the same idea under names that read better in context.
ANCHOR_FIELDS = (
    "anchor",
    "gate_anchor",
    "evidence_anchor",
    "accepts_all_three",
    "but_gates_on_one",
)

# Where key NAMES live, for the code -> schema direction.
KEYNAME_FIELDS = ("keys", "domain_keys", "keys_by_component_type")
KEYLIST_FIELDS = (
    "aliases", "value_aliases", "inventory", "required", "optional",
    # a fallback chain names real section spellings: `linear magnetic` exists
    # only here, not as a `sections:` entry of its own
    "resolution_chains",
)

# Section names are part of the contract too: the code looks them up with
# `section( "..." )` / `section_exists( "..." )`, so they appear as parse-call
# literals exactly like keys do. Collecting only `keys:` reported all thirteen
# section names as undocumented on the first run.
SECTION_FIELDS = ("sections", "thermal_subsection", "subsections")


@dataclass
class Anchor:
    token: str
    path: str          # YAML path, e.g. sections.solver.sections.nonlinear.keys.tolerance
    consumer: str | None


@dataclass
class Schema:
    data: dict
    anchors: list[Anchor] = field(default_factory=list)
    keys: set[str] = field(default_factory=set)
    section_names: set[str] = field(default_factory=set)
    anchor_fields_seen: int = 0

    @property
    def sections(self) -> dict:
        return self.data.get("sections", {})


def load(root: Path) -> Schema:
    path = root / "doc/input_schema.yaml"
    data = yaml.safe_load(path.read_text())
    schema = Schema(data=data)
    _walk(data, "", None, schema)
    return schema


def _walk(node, path: str, consumer: str | None, schema: Schema) -> None:
    """Recursive descent, inheriting the nearest enclosing `consumer:`."""
    if isinstance(node, dict):
        # a consumer declared here applies to everything below it
        local = node.get("consumer", consumer)
        if isinstance(local, str):
            local = local.strip().strip("'\"").split("/")[-1]
        else:
            local = consumer

        for key, value in node.items():
            child = f"{path}.{key}" if path else str(key)

            if key in ANCHOR_FIELDS:
                # count DECLARED per token, exactly as the checked side counts,
                # or the declared-vs-checked self-report can mask a skip: a
                # list-valued anchor counted as one declaration makes checked
                # exceed declared, and mixed shapes cancel out. A field whose
                # value yields no token still counts one, so it shows up as a
                # declared-but-unchecked warning instead of vanishing.
                tokens = _as_tokens(value)
                schema.anchor_fields_seen += len(tokens) if tokens else 1
                for token in tokens:
                    schema.anchors.append(Anchor(token, path or "(root)", local))
                continue

            if key == "anchors" and isinstance(value, list):
                for token in value:
                    schema.anchor_fields_seen += 1
                    if isinstance(token, str):
                        schema.anchors.append(Anchor(token, path or "(root)", local))
                continue

            if key in KEYNAME_FIELDS and isinstance(value, dict):
                schema.keys.update(
                    k.lower() for k in value.keys() if isinstance(k, str)
                )

            if key in SECTION_FIELDS and isinstance(value, dict):
                schema.keys.update(
                    k.lower() for k in value.keys() if isinstance(k, str)
                )
                schema.section_names.update(
                    k.lower() for k in value.keys() if isinstance(k, str)
                )

            if key in KEYLIST_FIELDS:
                schema.keys.update(_as_keynames(value))

            _walk(value, child, local, schema)

    elif isinstance(node, list):
        for i, item in enumerate(node):
            _walk(item, f"{path}[{i}]", consumer, schema)


def _as_tokens(value) -> list[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return [v for v in value if isinstance(v, str)]
    return []


def _as_keynames(value) -> set[str]:
    """Key names hide in several shapes: a string, a list, a dict, or a dict of
    lists (`resolution_chains: {magnetic: [linear magnetic, linear]}`), so this
    recurses rather than handling one level."""
    out: set[str] = set()
    if isinstance(value, str):
        out.add(value.lower())
    elif isinstance(value, list):
        for item in value:
            out |= _as_keynames(item)
    elif isinstance(value, dict):
        out.update(k.lower() for k in value.keys() if isinstance(k, str))
        for item in value.values():
            out |= _as_keynames(item)
    return out


def known_keys(schema: Schema) -> set[str]:
    """Every key name the schema is aware of, lowercased.

    Anchors are included because the file's own convention is that a key anchor
    IS the bare quoted key literal, so the anchor set and the key set largely
    coincide. Lowercasing matters: `Section::get_string` and `key_exists` both
    call `string_to_lower`, so `get_reals( "Direction" )` reads the key stored
    as `direction`. Comparing case-sensitively would report it as unknown.
    """
    keys = set(schema.keys)
    for anchor in schema.anchors:
        token = anchor.token.strip()
        if token.startswith('"') and token.endswith('"') and len(token) > 2:
            keys.add(token[1:-1].lower())
    return {k for k in keys if k}
