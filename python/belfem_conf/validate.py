"""Deck validation against `doc/input_schema.yaml`.

What this can and cannot check is bounded by the schema, not by ambition, and
the report says so explicitly at the end: `-v` lists every class of check that
was skipped, and why. Conditional requiredness (`audit_todo` T3) is enforced;
what remains out of reach needs data the deck does not carry — mesh-state
gates, bracket-group arity, and the material-shape prose (T4).

A validator that silently skips a class of check is worse than one that says it
skipped it — the reader cannot tell "passed" from "never looked".
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import units
from .document import Document
from .mesh import expand_ids
from .parser import Section, Statement

ERROR, WARNING, NOTE = "error", "warning", "note"

# Sections whose OWN statements are data rather than schema-named keys: the
# `materials` body holds `<name> : builtin ;` selectors, a `layers` stack holds
# material names, a `curves` block holds curve ids. An unknown-key check is
# meaningless there, so it is skipped and counted rather than fabricated.
# Deliberately NOT listed: boundary-condition and topology subsections — their
# NAMES are data but their KEYS are schema'd (the merged domain/source-function
# maps), so a typo'd `amplitde` there must be reported.
DATA_KEYED = {"materials", "layers", "curves"}


@dataclass
class Diagnostic:
    severity: str
    message: str
    line: int
    where: str
    hint: str = ""


@dataclass
class Result:
    diagnostics: list[Diagnostic] = field(default_factory=list)
    checked: dict[str, int] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)

    def add(self, sev, msg, line, where, hint=""):
        self.diagnostics.append(Diagnostic(sev, msg, line, where, hint))

    def tick(self, what: str, n: int = 1):
        self.checked[what] = self.checked.get(what, 0) + n

    @property
    def errors(self) -> int:
        return sum(1 for d in self.diagnostics if d.severity == ERROR)


def _line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def validate(doc: Document, schema: dict, root: Path) -> Result:
    res = Result()
    sections = schema.get("sections", {})
    rootstr = str(root)

    units_ok = bool(units.table(rootstr))
    if not units_ok:
        res.skipped.append(
            "unit dimensions — unit_to_si could not be parsed out of "
            "src/core/stringtools.cpp (moved, or reformatted beyond the "
            "extractor's regex)")

    _check_top_level(doc, sections, res)
    _check_domain_types(doc, sections, res)
    _check_bc_types(doc, sections, res)
    _check_circuit(doc, sections, res)

    bc_keys = _bc_keys(sections)
    topo_keys = _topology_keys(sections)
    mat_keys = (sections.get("materials") or {}).get("keys")
    select_first = (sections.get("materials") or {}).get("select_first")
    families = _families(sections)
    bc_root = sections.get("boundary conditions", {})
    no_source = set((bc_root.get("source_function_keys") or {})
                    .get("does_not_apply_to_bc_type") or [])
    source_names = set((bc_root.get("source_function_keys") or {}).get("keys") or {})

    # `one_of` sits on the section that DESCRIBES the subsections, so the rules
    # for a material body live on `materials`, and the rules for a BC source
    # function on `source_function_keys`. Map each to the deck sections it
    # governs rather than to the section it is written on.
    one_of_by_owner = {
        "materials": (sections.get("materials") or {}).get("one_of") or [],
        "boundary conditions": (bc_root.get("source_function_keys") or {}).get("one_of") or [],
    }

    for sec in doc.walk():
        spec = _spec_for(sec, sections)
        keys = _keys_for(sec, spec, bc_keys, topo_keys, mat_keys, sections)
        _check_keys(doc, sec, keys, res, rootstr,
                    select_first if sec.parent is not None and sec.parent.type == "materials" else None,
                    units_ok)
        _check_duplicates(doc, sec, res)
        if spec is not None and isinstance(spec.get("keys"), dict):
            # only where the schema owns this section's keys outright — the
            # merged BC/topology maps carry per-type applicability instead,
            # and a blanket required check there would be wrong for most types
            _check_required_keys(doc, sec, spec["keys"], res)
        if isinstance(keys, dict):
            _check_declared_constraints(doc, sec, keys, res, families)

        # a material body is a direct child of `materials`; a BC source
        # function lives on the BC-type section, wherever under
        # `boundary conditions` it sits
        if sec.parent is not None and sec.parent.type == "materials":
            _check_one_of(doc, sec, one_of_by_owner["materials"], res, families,
                          overrides={"builtin": _builtin_of(sec)})
        elif _under_boundary_conditions(sec):
            _check_one_of(doc, sec, one_of_by_owner["boundary conditions"],
                          res, families)

        if _under_boundary_conditions(sec):
            _check_applicability(doc, sec, bc_keys, res, families)
            # a BC type that reads no source function at all
            if sec.type in no_source:
                for stmt in sec.statements:
                    if stmt.name in source_names:
                        res.tick("key applicability")
                        res.add(WARNING,
                                f"`{sec.type}` reads no source function, so "
                                f"`{stmt.key}` is ignored",
                                _line_of(doc.text, stmt.key_offset), sec.path,
                                "its second switch only re-reads `nodes`")

    _check_cross_references(doc, res, rootstr, units_ok)
    _check_landmines(doc, res)

    # T3 is complete: conditional requiredness, one-of groups, section
    # fallbacks and the ordered material-shape resolution are all enforced.
    # What is left needs data the deck does not contain.
    unavailable = _check_against_mesh(doc, res)
    if unavailable:
        res.skipped.append(unavailable)
    res.skipped.append(
        "the mesh-STATE gates — a `curves` or `periodic` block is silently "
        "ignored when the mesh already carries them. Only a `.bfm` can, and "
        "reading one needs h5py, which this package avoids depending on")
    res.skipped.append(
        "bracket-group arity — that input and output terminal group counts "
        "match. `id_groups` semantics are not modelled yet")
    res.skipped.append(
        "material SHAPE applicability — `applies_to_material_shape` values are "
        "still prose ('pure metals', 'HTS without file'), so a `RRR` on an HTS "
        "is not flagged")
    return res


def _under_boundary_conditions(sec: Section) -> bool:
    # outermost ancestor, for the same reason as _keys_for: a nested section
    # may share a name with a top-level one
    outermost = None
    node = sec.parent
    while node is not None and node.header_span != (0, 0):
        outermost = node.type
        node = node.parent
    return outermost == "boundary conditions"


# ---------------------------------------------------------------- structure

def _check_top_level(doc: Document, sections: dict, res: Result) -> None:
    _check_required_sections(doc, sections, res)

    known = set(sections) | {"layers"}
    for sec in doc.root.sections:
        if sec.type not in known:
            res.add(ERROR, f"unknown top-level section `{sec.type}`",
                    _line_of(doc.text, sec.header_offset), sec.type,
                    "not described in doc/input_schema.yaml")
        res.tick("top-level sections")

    # Text the parser discards, exactly as the C++ does. Harmless, but it is
    # inert input, and inert input that nobody mentions is how a typo survives
    # for years. Two kinds: text above a section header, and text before a
    # closing brace with no `;` — the latter is usually a statement that lost
    # its semicolon, and the C++ drops the whole statement in silence.
    for sec in [doc.root, *doc.walk()]:
        if sec is not doc.root and sec.header_discarded:
            res.tick("stray text")
            res.add(NOTE,
                    f"ignored text before the `{sec.type}` header: "
                    f"{sec.header_discarded!r}",
                    _line_of(doc.text, sec.header_offset), sec.path,
                    "the parser keeps only the line immediately before `{`")
        if sec.trailing_discarded:
            res.tick("stray text")
            res.add(WARNING,
                    f"text with no `;` is dropped: {sec.trailing_discarded!r}",
                    _line_of(doc.text, sec.trailing_offset),
                    sec.path or "(deck)",
                    "a statement without its semicolon is silently ignored")


def _check_required_sections(doc: Document, sections: dict, res: Result,
                             prefix: str = "", parent=None) -> None:
    """`required: true` and `required_unless_section`, at every nesting level.

    Plain `required` used to be enforced for top-level sections only, so a deck
    with `solver` but no `timestep` validated clean and died in the C++. A
    nested requirement only applies when its holder actually exists — absent
    `solver` is one finding, not one per missing child.

    `required_unless_section` (required unless a named sibling exists) needs
    its own field because the predicate is the existence of a SIBLING SECTION,
    not a key value. The real case is `solver { linear }`, which the magnetic
    and thermal paths both fall back to — a deck with neither dies at
    `section( "linear" )` with no useful message.
    """
    holder = parent if parent is not None else doc.root
    siblings = {s.type for s in holder.sections}
    line = (_line_of(doc.text, holder.header_offset)
            if holder is not doc.root else 1)
    where = prefix.rstrip("/") or "(deck)"

    for name, spec in sections.items():
        if not isinstance(spec, dict):
            continue

        spellings = {name} | _section_aliases(spec)

        if spec.get("required") is True:
            res.tick("required sections")
            if not (spellings & siblings):
                other = sorted(spellings - {name})
                res.add(ERROR, f"required section `{prefix}{name}` is missing",
                        line, where,
                        ("or one of: " + ", ".join(other)) if other else "")

        alternates = spec.get("required_unless_section")
        if isinstance(alternates, list):
            res.tick("section fallbacks")
            if name not in siblings and not (siblings & set(alternates)):
                res.add(ERROR,
                        f"`{prefix}{name}` is missing and so is every fallback",
                        line, where,
                        "one of: " + ", ".join([name] + [str(a) for a in alternates]))

        # Descend through whichever spelling the deck used. `nonlinear` and
        # `nonlinear magnetic` share a schema entry, so its subsections
        # (`coulomb gauge penalty`, `nitsche ghost penalty`) hang off the alias
        # too -- and the C++ reads them from the WINNING section only.
        nested = spec.get("sections")
        here = None
        for spelling in _spellings_in_precedence(name, spec):
            here = holder.section(spelling)
            if here is not None:
                break
        if isinstance(nested, dict) and here is not None:
            _check_required_sections(doc, nested, res, f"{prefix}{name}/", here)


def _check_domain_types(doc: Document, sections: dict, res: Result) -> None:
    topo = doc.root.section("topology")
    catalogue = sections.get("topology", {}).get("domain_types")
    if topo is None or not catalogue:
        return

    accepted, rejected = {}, {}
    for entry in catalogue.get("accepted", []):
        for spelling in entry["spellings"]:
            accepted[spelling] = entry
    for entry in catalogue.get("rejected", []):
        for spelling in entry["spellings"]:
            rejected[spelling] = entry

    for sub in topo.sections:
        res.tick("topology domain types")
        line = _line_of(doc.text, sub.header_offset)
        if sub.type in rejected:
            entry = rejected[sub.type]
            res.add(ERROR, f"`{sub.type}` is no longer a valid domain type", line,
                    f"topology/{sub.type}",
                    f"use instead: {entry.get('use_instead', '—')}")
        elif sub.type not in accepted:
            res.add(ERROR, f"unknown domain type `{sub.type}`", line,
                    f"topology/{sub.type}",
                    "domain_type() would abort with 'Unknown Domain Type'")
        else:
            entry = accepted[sub.type]
            if entry.get("material") == "required" and not sub.key_exists("material"):
                res.add(ERROR,
                        f"`{sub.type}` requires a `material` key", line,
                        f"topology/{sub.type}",
                        "the Domain constructor hard-errors without it")
            _check_group_keys(doc, sub, entry, res)


def _check_group_keys(doc, sub, entry, res) -> None:
    kind = entry.get("groups")
    if kind not in ("blocks", "sidesets"):
        return
    singular, plural = ("block", "blocks") if kind == "blocks" else ("sideset", "sidesets")
    line = _line_of(doc.text, sub.header_offset)
    res.tick("group keys")
    # Only the "neither is present" half lives here. Which of the two pairs a
    # domain type needs comes from `domain_types.groups`, which the schema can
    # express; that BOTH must not appear is `xor_with`, which it also expresses
    # and which _check_declared_constraints now reads. Keeping both halves here
    # reported the same fault twice.
    if not sub.key_exists(singular) and not sub.key_exists(plural):
        res.add(ERROR, f"`{sub.type}` defines no {plural}", line,
                f"topology/{sub.type}", f"needs `{singular}` or `{plural}`")


def _check_bc_types(doc: Document, sections: dict, res: Result) -> None:
    """A boundary-condition subsection is NAMED for its BC type, and the name
    is resolved through `boundary_condition_type()` — an unknown one is fatal,
    and `background dirichlet` is refused outright. The names are data to the
    unknown-key machinery, so without this check a typo'd type sailed through.
    """
    catalogue = (sections.get("boundary conditions") or {}).get("bc_types")
    if not isinstance(catalogue, dict):
        return
    legal = {str(v).lower() for v in catalogue.get("values") or []}
    if not legal:
        return

    for bc_root in doc.root.find_all("boundary conditions"):
        for sub in bc_root.sections:
            res.tick("bc types")
            line = _line_of(doc.text, sub.header_offset)
            if sub.type == "background dirichlet":
                res.add(ERROR,
                        "`background dirichlet` is refused at parse", line,
                        sub.path,
                        "it used to corrupt an unrelated BC — see the schema's "
                        "bc_types.background_dirichlet history")
            elif sub.type not in legal:
                res.add(ERROR, f"unknown boundary-condition type `{sub.type}`",
                        line, sub.path,
                        "one of: " + ", ".join(sorted(legal)))


def _check_circuit(doc: Document, sections: dict, res: Result) -> None:
    """The circuit contract from the schema fields that were never read.

    `required_subsections`, `component_section_types`,
    `universal_component_keys` and `keys_by_component_type` describe what
    `ElectricalCircuitFactory::read_circuit` hard-errors on: a missing
    `topology`, an unknown component type, a labelled component header, and a
    component without `node +` / `node -`. Per-type required keys are enforced
    only when the schema entry names a real key from `component_keys.inventory`
    — some entries ("source function family") describe a key GROUP in prose,
    and inventing a key name for them would report ghosts.
    """
    spec = sections.get("circuit")
    circuit = doc.root.section("circuit")
    if not isinstance(spec, dict) or circuit is None:
        return
    line = _line_of(doc.text, circuit.header_offset)

    for name in spec.get("required_subsections") or []:
        res.tick("circuit structure")
        if circuit.section(str(name)) is None:
            res.add(ERROR, f"`circuit` has no `{name}` subsection", line,
                    circuit.path, "the factory hard-errors without it")

    types_spec = spec.get("component_section_types") or {}
    legal = {str(v).lower() for v in types_spec.get("values") or []}
    inventory = _as_keynames_lower(
        (spec.get("component_keys") or {}).get("inventory"))
    # sources mirror the BC source-function keys (the schema's component_keys
    # note names them), so those spellings and their aliases are legal too
    sfk = ((sections.get("boundary conditions") or {})
           .get("source_function_keys") or {}).get("keys") or {}
    for kname, kspec in sfk.items():
        inventory.add(str(kname).lower())
        if isinstance(kspec, dict):
            inventory.update(str(a).lower() for a in _as_list(kspec.get("aliases")))
    inventory.update({"file", "units"})
    universal = (spec.get("universal_component_keys") or {}).get("required") or []
    by_type = spec.get("keys_by_component_type") or {}

    topo = circuit.section("topology")
    if topo is None:
        return
    for comp in topo.sections:
        cline = _line_of(doc.text, comp.header_offset)
        res.tick("circuit components")

        if legal and comp.type not in legal:
            res.add(ERROR, f"unknown circuit component `{comp.type}`", cline,
                    comp.path, "one of: " + ", ".join(sorted(legal)))
            continue
        if comp.label and types_spec.get("forbid_label_in_header"):
            res.add(ERROR,
                    f"component header `{comp.type} : {comp.label}` must be "
                    f"the bare type", cline, comp.path,
                    "a labelled header does not resolve — use `label : ...` "
                    "inside the section")

        for key in universal:
            if not comp.key_exists(str(key)):
                res.add(ERROR, f"`{comp.type}` has no `{key}`", cline,
                        comp.path, "every component needs node + and node -")

        wanted = (by_type.get(comp.type) or {}).get("required") or []
        for key in wanted:
            key = str(key)
            if key.lower() not in inventory:
                continue            # prose descriptor, not a literal key
            if not comp.key_exists(key):
                res.add(ERROR,
                        f"a `{comp.type}` requires `{key}`", cline, comp.path)

        if inventory:
            for stmt in comp.statements:
                if stmt.name not in inventory:
                    res.add(WARNING, f"unknown key `{stmt.key}`",
                            _line_of(doc.text, stmt.key_offset), comp.path,
                            f"not in the component-key inventory for `{comp.type}`")


def _as_keynames_lower(value) -> set[str]:
    if not isinstance(value, list):
        return set()
    return {str(v).lower() for v in value}


# ---------------------------------------------------------------- key level

def _section_aliases(spec: dict) -> set[str]:
    """Every deck spelling that resolves to this schema section.

    Two schema fields name the same thing in different shapes. `aliases` is the
    flat form (`nonlinear` <- `nonlinear magnetic`). `resolution_chains` is the
    per-field form: `linear` is reached as `linear magnetic` or `linear
    thermal`, with the bare name last in each chain. Reading only the first
    made `solver/nonlinear` "missing" on every deck that spells it
    `nonlinear magnetic`; reading neither left `linear magnetic` matched to no
    schema entry at all, so NOTHING inside it was checked.
    """
    if not isinstance(spec, dict):
        return set()
    names = {str(a) for a in _as_list(spec.get("aliases"))}
    chains = spec.get("resolution_chains")
    if isinstance(chains, dict):
        for chain in chains.values():
            names.update(str(c) for c in _as_list(chain))
    return names


def _spellings_in_precedence(name: str, spec: dict) -> list[str]:
    """Every spelling of this section, most-specific first.

    Both schema fields put the winner first, and both agree with the C++.
    `resolution_chains` is already ordered (`[linear magnetic, linear]`), and
    `alias_precedence` says in prose that `nonlinear magnetic` beats
    `nonlinear`. So: chain/alias spellings first, the bare name last. It
    matters when a deck carries BOTH — the controller then reads the specific
    one and ignores the other, and so must anything that descends into it.
    """
    ordered = []
    chains = spec.get("resolution_chains") if isinstance(spec, dict) else None
    if isinstance(chains, dict):
        for chain in chains.values():
            for step in _as_list(chain):
                if str(step) not in ordered:
                    ordered.append(str(step))
    for alias in _as_list(spec.get("aliases") if isinstance(spec, dict) else None):
        if str(alias) not in ordered:
            ordered.append(str(alias))
    if name in ordered:
        ordered.remove(name)
    return ordered + [name]


def _lookup_section(name: str, table: dict) -> dict | None:
    """The schema entry governing a deck section called `name`.

    An exact entry always wins over an alias, so `nonlinear thermal` -- which
    has its own entry AND differs from `nonlinear` in three keys -- is never
    silently validated against the magnetic one.
    """
    if not isinstance(table, dict):
        return None
    if name in table and isinstance(table[name], dict):
        return table[name]
    for other, spec in table.items():
        if other != name and name in _section_aliases(spec):
            return spec
    return None


def _spec_for(sec: Section, sections: dict) -> dict | None:
    """Find the schema entry for a deck section, by walking the same path."""
    chain = []
    node = sec
    while node is not None and node.header_span != (0, 0):
        chain.append(node.type)
        node = node.parent
    chain.reverse()

    spec = _lookup_section(chain[0], sections) if chain else None
    for step in chain[1:]:
        if not isinstance(spec, dict):
            return None
        nxt = spec.get("sections", {})
        if not isinstance(nxt, dict):
            return None
        hit = _lookup_section(step, nxt)
        if hit is not None:
            spec = hit
        elif _names_its_subsections_by_data(spec):
            # A user-named level: a material section is `materials/<its own
            # name>`, so the name is data and the schema has no entry for it.
            # Stay on the current spec and keep descending, otherwise a real
            # subsection below it — `materials/ybco/defect` — cannot resolve
            # and every one of its keys is reported as unknown.
            continue
        else:
            return None
    return spec if isinstance(spec, dict) else None


def _names_its_subsections_by_data(spec: dict) -> bool:
    """True when this section's subsections are named by the user, not the schema."""
    return bool(spec.get("subsection_naming")
                or spec.get("subsection_type_is_bc_type"))


def _bc_keys(sections: dict) -> dict:
    """Key specs for a boundary-condition subsection.

    BC types do not each get a `keys:` map — the schema factors them into
    `domain_keys` (nodes, sidesets, terminals, length, direction) and
    `source_function_keys.keys` (type, amplitude, period, ...), which apply
    across types. Without joining the two, the section with the most keys in a
    deck is skipped in silence.
    """
    spec = sections.get("boundary conditions", {})
    merged = {}
    for name, kspec in (spec.get("domain_keys") or {}).items():
        if isinstance(kspec, dict):
            merged[name] = kspec
    sfk = (spec.get("source_function_keys") or {}).get("keys") or {}
    for name, kspec in sfk.items():
        if isinstance(kspec, dict):
            merged.setdefault(name, kspec)
    return merged


def _topology_keys(sections: dict) -> dict:
    """Key specs for a topology domain subsection.

    Like the BC section, topology factors its keys out into `domain_keys` plus
    the thin-shell extras, rather than giving each domain type a `keys:` map.
    Until these are joined, a topology subsection is checked only for its group
    keys and nothing else.
    """
    spec = sections.get("topology", {})
    merged = {}
    for source in ("domain_keys", "thinshell_keys"):
        for name, kspec in (spec.get(source) or {}).items():
            if isinstance(kspec, dict):
                merged[name] = kspec
    return merged


def _families(sections: dict) -> dict[str, list[str]]:
    sfk = (sections.get("boundary conditions", {})
           .get("source_function_keys") or {})
    fam = sfk.get("function_families") or {}
    return {k: list(v) for k, v in fam.items() if isinstance(v, list)}


def _check_required_keys(doc, sec, keys, res) -> None:
    """Plain `required: true` on a key, where the schema owns the section's
    keys outright. Was never read: a deck with no `mesh { file }` or no
    `simulation time` validated clean and died in setup code. Only the literal
    boolean counts — `required: conditional` and prose values are handled by
    the predicate fields or not at all, and are reported by `-v` as skipped."""
    line = _line_of(doc.text, sec.header_offset)
    for name, kspec in keys.items():
        if not isinstance(kspec, dict) or kspec.get("required") is not True:
            continue
        res.tick("required keys")
        present = sec.key_exists(name) or any(
            sec.key_exists(a) for a in _as_list(kspec.get("aliases")))
        if not present:
            res.add(ERROR, f"`{sec.type}` requires a `{name}` key", line,
                    sec.path)


def _builtin_of(sec) -> str | None:
    """The builtin material name this section selects, or None.

    Mirrors the THREE paths of MaterialFactory (tested after `curve`, which
    wins outright): (1) `builtin : <type> ;` in the section; (2) the section is
    labeled builtin — `ybco : builtin { }`; (3) a `<name> : builtin ;`
    statement at the `materials` level. The schema's HTS rule is conditioned on
    the selected NAME, so evaluating only path (1) let `ybco : builtin { }`
    skip the file-or-jc+n check that the factory then hard-errors on.
    """
    if sec.key_exists("curve"):
        return None
    value = (sec.get("builtin") or "").strip()
    if value and value != "true":              # "true" is a bare flag
        return value
    if sec.label == "builtin":
        return sec.type
    parent = sec.parent
    if parent is not None and (parent.get(sec.type) or "").strip() == "builtin":
        return sec.type
    return None


def _check_duplicates(doc, sec, res) -> None:
    """Duplicate keys and duplicate subsections, which the C++ collapses
    without a word: keys land in a map (last one wins) and named-section
    lookup does the same, while index-based iteration still sees every
    duplicate. Data-keyed sections are exempt — a `layers` stack legitimately
    repeats a material name, and its order IS the physical stack."""
    if sec.type in DATA_KEYED:
        return

    seen: dict[str, int] = {}
    for stmt in sec.statements:
        seen[stmt.name] = seen.get(stmt.name, 0) + 1
    for name, count in seen.items():
        if count > 1:
            res.tick("duplicate keys")
            res.add(WARNING,
                    f"`{name}` appears {count} times — only the last is used",
                    _line_of(doc.text, sec.header_offset), sec.path or "(deck)",
                    "duplicate keys overwrite in the C++ key map")

    kids: dict[str, int] = {}
    for child in sec.sections:
        key = f"{child.type}:{child.label}"
        kids[key] = kids.get(key, 0) + 1
    for key, count in kids.items():
        if count > 1:
            res.tick("duplicate sections")
            res.add(NOTE,
                    f"section `{key.rstrip(':')}` appears {count} times",
                    _line_of(doc.text, sec.header_offset), sec.path or "(deck)",
                    "named lookup uses the last; per-index consumers "
                    "(materials, circuit) process every one")


def _check_declared_constraints(doc, sec, keys, res, families=None) -> None:
    """Enforce the constraint fields the schema already carries.

    `xor_with`, `required_for` and `pairs_with` were declared but never read —
    `validate.py` re-implemented the block/sideset xor by hand instead, which
    is the same "schema knowledge living in Python" pattern the drift checker
    exists to prevent. Reading them here removes the duplicate and makes the
    fields worth writing.

    `required_for` is conditional requiredness in its narrowest useful form:
    required when this section's own `type` value is in the list.
    """
    line = _line_of(doc.text, sec.header_offset)
    function = (sec.get("type") or "").strip()
    seen_xor: set[frozenset] = set()

    for name, kspec in keys.items():
        if not isinstance(kspec, dict):
            continue
        present = sec.key_exists(name) or any(
            sec.key_exists(a) for a in _as_list(kspec.get("aliases")))

        # -- xor_with: at most one of the pair. (The "at least one" half is
        #    not universal — a BC may use `nodes` instead of either sideset
        #    spelling — so only domain types enforce it, via their group rule.)
        partner = kspec.get("xor_with")
        if isinstance(partner, str):
            pair = frozenset((name, partner))
            if pair not in seen_xor and sec.key_exists(partner) and present:
                seen_xor.add(pair)
                res.tick("xor constraints")
                res.add(ERROR,
                        f"`{name}` and `{partner}` are mutually exclusive",
                        line, sec.path, "give exactly one")

        # -- required_for: keyed on this section's `type`
        needed = kspec.get("required_for")
        if isinstance(needed, list) and function and function in needed:
            res.tick("conditional requiredness")
            if not present:
                res.add(ERROR,
                        f"`{name}` is required for a `{function}` source function",
                        line, sec.path)

        # -- required_unless: a sibling key's value decides
        rule = kspec.get("required_unless")
        if isinstance(rule, dict):
            res.tick("conditional requiredness")
            if not present and not _predicate_holds(sec, rule, families):
                res.add(ERROR,
                        f"`{name}` is required unless {_describe_rule(rule)}",
                        line, sec.path)

        # -- pairs_with: both present, and the same number of entries
        twin = kspec.get("pairs_with")
        if isinstance(twin, str) and present and sec.key_exists(twin):
            a_raw, b_raw = sec.get(name) or "", sec.get(twin) or ""
            if "[" in a_raw or "[" in b_raw:
                # bracket groups compare by GROUP, and id_groups semantics are
                # not modelled — counting ids here would invent errors
                continue
            res.tick("paired lists")
            a = len(expand_ids(a_raw))
            b = len(expand_ids(b_raw))
            if a != b:
                res.add(ERROR,
                        f"`{name}` and `{twin}` must have the same number of "
                        f"entries, have {a} and {b}", line, sec.path,
                        "they are matched element by element; ranges like 4:9 "
                        "count expanded")


def _predicate_holds(sec, rule: dict, families: dict | None = None,
                     overrides: dict | None = None) -> bool:
    """Evaluate a sibling-key predicate. Unknown shapes are treated as unmet.

    `in` may name a family (`periodic family`), which is expanded through
    `function_families` — so a condition can be written once against the family
    rather than restated as four source types that could drift apart.

    `overrides` supplies a RESOLVED value for a sibling that the deck does not
    spell as a plain key — the material `builtin` selection has three spellings
    in the factory, and the caller resolves them to one name.
    """
    sibling = rule.get("sibling")
    if not sibling:
        return False
    if overrides and sibling in overrides:
        value = overrides[sibling]
    else:
        value = sec.get(sibling)
    if value is None:
        return False
    value = value.strip()

    if "ends_with" in rule:
        return value.lower().endswith(str(rule["ends_with"]).lower())
    if "equals" in rule:
        return value.lower() == str(rule["equals"]).lower()
    if "in" in rule and isinstance(rule["in"], list):
        expanded: set[str] = set()
        for entry in rule["in"]:
            expanded.update((families or {}).get(str(entry), [str(entry)]))
        return value.lower() in {str(v).lower() for v in expanded}
    return False


def _check_one_of(doc, sec, constraints, res, families, overrides=None) -> None:
    """Exactly one alternative must be satisfied, when the condition holds.

    Generalises `xor_with`, which links two keys. Neither real rule is a pair:
    an HTS takes `file` OR (`jc` AND `n`), and the periodic family takes
    `period` OR `frequency`.
    """
    line = _line_of(doc.text, sec.header_offset)
    for rule in constraints:
        if not isinstance(rule, dict):
            continue
        alternatives = rule.get("alternatives")
        if not isinstance(alternatives, list) or not alternatives:
            continue

        when = rule.get("when")
        if isinstance(when, dict) and not _predicate_holds(sec, when, families,
                                                           overrides):
            continue                      # condition not met; rule is dormant

        res.tick("one-of constraints")
        satisfied = [
            alt for alt in alternatives
            if all(sec.key_exists(k) for k in _as_list(alt))
        ]
        shown = " or ".join(
            "(" + " + ".join(map(str, _as_list(a))) + ")" for a in alternatives)

        if not satisfied:
            res.add(ERROR, f"none of {shown} is given", line, sec.path,
                    rule.get("note", "").strip().split(".")[0] or None)
        elif len(satisfied) > 1:
            res.add(WARNING, f"more than one of {shown} is given", line,
                    sec.path, "only the first is used")


def _describe_rule(rule: dict) -> str:
    sibling = rule.get("sibling", "?")
    for op, word in (("ends_with", "ends with"), ("equals", "is"), ("in", "is one of")):
        if op in rule:
            target = rule[op]
            shown = ", ".join(map(str, target)) if isinstance(target, list) else target
            return f"`{sibling}` {word} {shown!r}"
    return f"`{sibling}` satisfies its condition"


def _check_applicability(doc, sec, keys, res, families) -> None:
    """Enforce the two axes that used to share the name `applies_to`.

    A boundary-condition subsection is named for its BC type, and it may carry
    a source-function `type`. A key can be constrained on either axis, and a
    key constrained on an axis it does not match is inert: the factory never
    reads it, so it silently does nothing.
    """
    bc_type = sec.type
    function = (sec.get("type") or "").strip()

    for stmt in sec.statements:
        kspec = keys.get(stmt.name)
        if not isinstance(kspec, dict):
            continue
        line = _line_of(doc.text, stmt.key_offset)

        allowed_bc = kspec.get("applies_to_bc_type")
        if isinstance(allowed_bc, list) and bc_type not in allowed_bc:
            res.tick("key applicability")
            res.add(WARNING,
                    f"`{stmt.key}` does not apply to a `{bc_type}` boundary "
                    f"condition", line, sec.path,
                    "applies to: " + ", ".join(map(str, allowed_bc))
                    + " — here it is read by nobody")

        allowed_fn = kspec.get("applies_to_function")
        if isinstance(allowed_fn, list) and function:
            expanded = set()
            for entry in allowed_fn:
                expanded.update(families.get(str(entry), [str(entry)]))
            if function not in expanded:
                res.tick("key applicability")
                res.add(WARNING,
                        f"`{stmt.key}` does not apply to a `{function}` source "
                        f"function", line, sec.path,
                        "applies to: " + ", ".join(sorted(expanded))
                        + " — here it is read by nobody")


def _per_field_keys(sec, spec, sections) -> dict:
    """Parent-level keys a FIELD-SPECIFIC spelling of this section may repeat.

    `compute conditioning` and `mumps error analysis` are declared on `solver`,
    and the controller reads each from `linear magnetic` / `linear thermal`
    as a per-field override (cl_FEM_Controller.cpp, anchor
    '"compute conditioning"'). It does NOT read them from a bare `linear`, so
    the merge is restricted to the chain spellings — otherwise the checker
    would bless a key the C++ never looks at there.
    """
    override = spec.get("per_field_override")
    if not isinstance(override, dict) or sec.parent is None:
        return {}
    if not isinstance(sections, dict):
        return {}
    parent = _spec_for(sec.parent, sections)
    if not isinstance(parent, dict):
        return {}
    siblings = parent.get("sections")
    parent_keys = parent.get("keys")
    if not isinstance(siblings, dict) or not isinstance(parent_keys, dict):
        return {}

    # `resolution_chains` lists the BARE name too, so an alias test alone would
    # also bless `compute conditioning` inside a plain `linear` — a section the
    # controller never reads it from. Identify the canonical spelling by
    # identity against the schema table and exclude exactly that one.
    canonical = next((n for n, v in siblings.items() if v is spec), None)
    if canonical is None or sec.type == canonical:
        return {}

    return {name: kspec for name, kspec in parent_keys.items()
            if name in override}


def _keys_for(sec, spec, bc_keys, topo_keys, mat_keys=None,
              sections=None) -> dict | None:
    """The key specs that govern this section.

    Most sections carry their own `keys:`. Boundary conditions and topology
    factor theirs out into shared maps instead, because the keys apply across
    types rather than per type — so those two are resolved by ancestry.
    """
    if spec is not None and isinstance(spec.get("keys"), dict):
        keys = spec["keys"]
        extra = _per_field_keys(sec, spec, sections)
        return {**keys, **extra} if extra else keys

    # The OUTERMOST ancestor decides, not any ancestor. `circuit` has its own
    # nested `topology` subsection, so matching the name anywhere in the chain
    # validated circuit components against mesh-topology keys and reported
    # every one of them as unknown.
    outermost = None
    node = sec.parent
    while node is not None and node.header_span != (0, 0):
        outermost = node.type
        node = node.parent

    if outermost == "boundary conditions" and bc_keys:
        return bc_keys
    if outermost == "topology" and topo_keys:
        return topo_keys
    if outermost == "materials" and mat_keys:
        return mat_keys
    return None


def _check_keys(doc, sec, keys, res, root: str, select_first=None,
                units_ok: bool = True) -> None:
    if not isinstance(keys, dict):
        return
    inert = frozenset(_inert_by_shape(sec, select_first))

    lookup = {}
    for name, kspec in keys.items():
        if not isinstance(kspec, dict):
            continue
        lookup[name.lower()] = (name, kspec)
        for alias in _as_list(kspec.get("aliases")):
            lookup[str(alias).lower()] = (name, kspec)

    data_keyed = sec.type in DATA_KEYED

    for stmt in sec.statements:
        line = _line_of(doc.text, stmt.key_offset)
        hit = lookup.get(stmt.name)
        if hit is None:
            if not data_keyed:
                res.tick("keys")
                res.add(WARNING, f"unknown key `{stmt.key}`", line, sec.path,
                        f"not described for section `{sec.type}`")
            continue

        res.tick("keys")
        _name, kspec = hit
        _check_value(doc, stmt, kspec, res, sec, line, root, inert, units_ok)


def _inert_by_shape(sec, select_first) -> set[str]:
    """Selectors that belong to a shape this section did not select.

    `select_first` says the first selector present decides, and the losers are
    never read. `examples/costheta` has an inert `builtin ;` beside its `curve`;
    enum-checking it reported a legal deck as broken.
    """
    if not isinstance(select_first, list):
        return set()
    chosen = None
    for entry in select_first:
        if isinstance(entry, dict) and sec.key_exists(entry.get("selector", "")):
            chosen = entry.get("selector")
            break
    if chosen is None:
        return set()
    return {e["selector"] for e in select_first
            if isinstance(e, dict) and e.get("selector") != chosen}


def _also_summary(kspec: dict):
    """The one-line description of what a key accepts BESIDES its enum."""
    also = kspec.get("also_accepts")
    if isinstance(also, dict):
        return also.get("summary")
    if isinstance(also, str):
        return also
    return None


def _check_also_accepts(stmt, kspec, res, sec, line, value) -> bool:
    """Apply a structured `also_accepts` rule; True when the value is legal.

    Only `alloy_formula` is modelled, because it is the only one the schema
    describes structurally. A prose-only `also_accepts` (a bare string) cannot
    be checked, so it PASSES rather than fails: the enum is then known to be an
    incomplete statement of the contract, and reporting a legal deck as broken
    is the worse error. That is exactly what happened to `Pb38Sn62`.
    """
    also = kspec.get("also_accepts")
    if isinstance(also, str):
        return True                     # prose only — cannot check, must not fail
    if not isinstance(also, dict):
        return False
    if also.get("id") != "alloy_formula":
        return True                     # a rule this checker does not model yet

    res.tick("alloy formulas")
    pattern = str(also.get("token_pattern") or r"([A-Za-z]+)([0-9]+)")
    tokens = re.findall(pattern, value)
    if not tokens:
        res.add(ERROR,
                f"`{stmt.key} : {value}` is neither a builtin nor an alloy "
                f"formula", line, sec.path,
                "to_pair finds no <symbol><percentage> token, so "
                "create_material hard-errors")
        return True                     # reported here; do not also report the enum

    rrr = str(also.get("rrr_token") or "RRR").lower()
    pure = {str(m).lower() for m in also.get("pure_metals") or []}
    barred = {str(m).lower() for m in also.get("not_pure_metals") or []}

    total = 0.0
    for symbol, percent in tokens:
        low = symbol.lower()
        if low == rrr:
            continue
        total += float(percent)
        if low in barred:
            res.add(ERROR,
                    f"`{symbol}` cannot be an alloy component", line, sec.path,
                    "Alloy::create_component requires MaterialType::PureMetal "
                    f"— `{symbol}` is not one")
        elif pure and low not in pure:
            res.add(ERROR,
                    f"unknown alloy component `{symbol}`", line, sec.path,
                    "one of: " + ", ".join(sorted(pure)))

    # the last fraction is overwritten with the balance, so a total that misses
    # 100 is legal and merely suspicious — never an error
    if also.get("fractions_are_percent") and tokens and abs(total - 100.0) > 0.5:
        res.add(WARNING,
                f"alloy fractions in `{value}` sum to {total:g}%, not 100%",
                line, sec.path,
                "legal — set_components overwrites the LAST fraction with the "
                "balance — but the stated number is then not what is used")
    return True


def _as_number(value: str):
    """The leading word of a deck value as a float, or None if it is not one."""
    words = value.split()
    if not words:
        return None
    try:
        return float(words[0])
    except ValueError:
        return None


def _round_half_away(number: float) -> int:
    """C `round()`: half away from zero, unlike Python's banker's rounding."""
    return int(math.floor(number + 0.5) if number >= 0
               else math.ceil(number - 0.5))


def _check_range(stmt, kspec, res, sec, line, number) -> None:
    """Enforce a declared `range: [lo, hi]`.

    Only a closed numeric interval is modelled. The sibling field `constraint:`
    stays prose ("> 0, and >= min iterations") and is deliberately NOT parsed —
    guessing at English would give exactly the confident-and-wrong findings this
    checker exists to avoid.

    `out_of_range: fatal` is the schema's own word for what the C++ does, so it
    selects the severity rather than being assumed: the two `anderson depth`
    keys are both `BELFEM_ERROR( tDepth >= 0 && tDepth <= 8, … )`
    (cl_FEM_Controller.cpp:4273 and :4445).
    """
    bounds = kspec.get("range")
    if not isinstance(bounds, list) or len(bounds) != 2:
        return
    try:
        lo, hi = float(bounds[0]), float(bounds[1])
    except (TypeError, ValueError):
        return

    res.tick("declared ranges")
    if lo <= number <= hi:
        return

    fatal = str(kspec.get("out_of_range", "")).lower() == "fatal"
    shown = _round_half_away(number) if float(number).is_integer() else number
    res.add(ERROR if fatal else WARNING,
            f"`{stmt.key} : {shown:g}` is outside {bounds[0]}..{bounds[1]}",
            line, sec.path,
            "rejected at setup" if fatal
            else "the schema declares this range but not what breaches it")


def _check_value(doc, stmt: Statement, kspec: dict, res, sec, line, root: str,
                 inert: frozenset = frozenset(), units_ok: bool = True) -> None:
    value = stmt.value
    kind = kspec.get("type")

    if stmt.name in inert:
        res.tick("inert-by-shape keys")
        res.add(NOTE,
                f"`{stmt.key}` is ignored here — another shape selector wins",
                line, sec.path,
                "the first selector present decides; the rest are never read")
        return

    values = kspec.get("values")
    if isinstance(values, list) and values:
        res.tick("enum values")
        legal = [str(v) for v in values]
        probe = value.lower() if kspec.get("value_normalize") == "lower" else value
        pool = [v.lower() for v in legal] if kspec.get("value_normalize") == "lower" else legal
        if probe not in pool:
            # An enum is not always the whole contract. `builtin` falls through
            # to a formula parser when the name matches no member, so a closed
            # list alone rejected `Pb38Sn62` -- a legal solder -- on two
            # shipped decks. The schema says so under `also_accepts`.
            if _check_also_accepts(stmt, kspec, res, sec, line, value):
                return
            sensitive = kspec.get("case_sensitive") or kspec.get("value_normalize") == "exact"
            res.add(ERROR, f"`{stmt.key} : {value}` is not a legal value", line,
                    sec.path,
                    "one of: " + ", ".join(legal)
                    + (", or " + str(_also_summary(kspec))
                       if _also_summary(kspec) else "")
                    + ("  (case-sensitive)" if sensitive else ""))
            return

        # a value can be legal to the PARSER and still refused downstream —
        # crc / crank-nicolson / galerkin parse into EulerMethod and then hit
        # a hard error in cl_IWG_Timestep. The schema marks those, and says in
        # as many words that a validator must reject them.
        status = kspec.get("value_runtime_status")
        if isinstance(status, dict):
            aliases = kspec.get("value_aliases") or {}
            canon = str(aliases.get(probe, probe)).lower()
            verdict = status.get(probe, status.get(canon))
            if verdict in ("accepted_then_hard_error", "rejected"):
                res.tick("runtime-rejected values")
                hint = str(kspec.get("runtime_status_note", "")).strip()
                res.add(ERROR,
                        f"`{stmt.key} : {value}` parses but is refused at "
                        f"runtime", line, sec.path,
                        hint.split(". ")[0] if hint else "marked "
                        f"{verdict} in the schema")
                return

    if kind in ("int", "uint"):
        res.tick("numeric values")
        number = _as_number(value)
        if number is None:
            res.add(ERROR, f"`{stmt.key}` expects an integer, got {value!r}",
                    line, sec.path)
        else:
            # `Section::get_int` is `round( get_real( … ) )`
            # (cl_Input_Section.cpp:349), so a fractional value is ROUNDED, not
            # refused. Reporting it as an error would be the same mistake as
            # rejecting `Pb38Sn62`: stricter than the code, on a deck that runs.
            rounded = _round_half_away(number)
            if rounded != number:
                res.add(WARNING,
                        f"`{stmt.key} : {value}` is not a whole number", line,
                        sec.path,
                        f"get_int rounds it to {rounded} — state whole numbers")
            _check_range(stmt, kspec, res, sec, line, rounded)
    elif kind == "real":
        res.tick("numeric values")
        number = _as_number(value)
        if number is None:
            res.add(ERROR, f"`{stmt.key}` expects a number, got {value!r}",
                    line, sec.path)
        else:
            _check_range(stmt, kspec, res, sec, line, number)
    elif kind == "bool":
        res.tick("boolean values")
        if value.lower() not in ("true", "false", "on", "off", "yes", "no", "1", "0"):
            res.add(WARNING,
                    f"`{stmt.key} : {value}` is not a recognised boolean",
                    line, sec.path,
                    "true/on/yes/1 are true; ANYTHING else is silently false")
    elif kind == "value" and units_ok:
        _check_dimension(stmt, kspec, res, sec, line, root, value)
    elif kind == "unit_token" and units_ok:
        # the VALUE is a bare unit (mesh.unit). The C++ runs it through
        # unit_to_si + check_unit with no key_exists guard, so both an unknown
        # token and a wrong dimension are fatal there.
        res.tick("unit dimensions")
        actual = units.dimension_of(value, root)
        if actual is None:
            res.add(ERROR, f"unknown unit `{value}`", line, sec.path,
                    "unit_to_si would abort with 'Unknown unit'")
        else:
            _compare_dimension(stmt, kspec, res, sec, line, value, value, actual)


def _check_dimension(stmt, kspec, res, sec, line, root, value) -> None:
    """Mirror of how `create_key` + `get_value` treat a `value` key.

    The C++ stores a (scaled number, dimension) pair ONLY for `<number>` and
    `<number> <unit>`. A bare number is stored dimensionless, and `get_value`
    compares dimensions strictly — so where the schema declares a dimension
    with a NON-EMPTY signature, a bare number is a runtime abort, not a
    default. Where the signature is empty it is not: `angle` resolves to no
    exponents at all, so `check_unit` accepts a bare number against "rad".
    Three-plus words store nothing and `get_value` aborts with "key is not a
    real". The old rule here ("a bare number is legal and dimensionless") was a
    false pass; testing the dimension's NAME rather than its signature was the
    over-correction.
    """
    declared = kspec.get("dimension")
    # Whether a bare number is legal is decided by the SIGNATURE, not by the
    # spelling. check_unit compares nothing but the dimension code, and `angle`
    # carries none — unit_to_si gives "rad" scale 1.0 and no exponents — so
    # `phase : 0` passes in the C++ exactly as `0 deg` does. Testing the name
    # against a two-word denylist reported a legal line in examples/2D_Tapestack
    # as a runtime abort.
    signature = units.expected(declared) if declared is not None else None
    dimensioned = (declared is not None
                   and declared not in ("dimensionless", "-")
                   and signature != {})
    words = value.split()

    if len(words) < 2:
        if kspec.get("unit_mandatory") or dimensioned:
            res.tick("unit dimensions")
            res.add(ERROR, f"`{stmt.key} : {value}` needs a unit", line,
                    sec.path,
                    f"a bare number is stored dimensionless and get_value "
                    f"expects {declared or 'a unit'} — runtime abort")
        return

    if len(words) > 2:
        res.tick("unit dimensions")
        res.add(ERROR,
                f"`{stmt.key} : {value}` must be `<number> <unit>`", line,
                sec.path,
                "with three or more words the C++ stores no value at all and "
                "get_value aborts")
        return

    token = words[1]
    res.tick("unit dimensions")
    actual = units.dimension_of(token, root)
    if actual is None:
        res.add(ERROR, f"unknown unit `{token}`", line, sec.path,
                "unit_to_si would abort with 'Unknown unit'")
        return
    _compare_dimension(stmt, kspec, res, sec, line, value, token, actual)


def _compare_dimension(stmt, kspec, res, sec, line, value, token, actual) -> None:
    declared = kspec.get("dimension")
    if declared is None:
        return
    want = units.expected(declared)
    if want is None:
        # a dimension spelling this checker cannot resolve must not pass in
        # silence — that is the "never looked" failure mode
        res.add(NOTE,
                f"`{stmt.key}`: schema dimension {declared!r} is not known "
                f"here, unit not checked", line, sec.path)
        return
    if actual != want:
        res.add(ERROR,
                f"`{stmt.key} : {value}` has the wrong dimension", line, sec.path,
                f"expected {declared} ({units.describe(want)}), "
                f"`{token}` is {units.describe(actual)}")


# ------------------------------------------------------------ cross-checks

def _check_cross_references(doc: Document, res: Result, root: str,
                            units_ok: bool) -> None:
    materials = doc.root.section("materials")
    known_materials = {s.type for s in materials.sections} if materials else set()
    if not known_materials:
        res.skipped.append(
            "layer materials resolving in `materials` — the deck has no "
            "materials section to resolve against")

    topo = doc.root.section("topology")
    if topo:
        for shell in [s for s in topo.sections
                      if s.type in ("thinshell", "tape", "shell")]:
            res.tick("thinshell/layers pairing")
            line = _line_of(doc.text, shell.header_offset)
            if doc.root.section("layers", shell.label) is None:
                res.add(ERROR,
                        f"thin shell `{shell.label}` has no `layers : "
                        f"{shell.label}` block", line, shell.path,
                        "a missing layer stack is fatal")

        periodic = topo.section("periodic")
        if periodic:
            for key in ("source", "target"):
                res.tick("periodic node triples")
                raw = periodic.get(key)
                if raw is None:
                    res.add(ERROR, f"`periodic` has no `{key}`",
                            _line_of(doc.text, periodic.header_offset),
                            periodic.path)
                    continue
                # counted EXPANDED, because the C++ reads these via get_ids
                # and `1:3` is three node ids there, not one token
                n = len(expand_ids(raw))
                if n != 3:
                    res.add(ERROR,
                            f"`periodic {key}` needs exactly 3 node ids, has {n}",
                            _line_of(doc.text, periodic.header_offset),
                            periodic.path, "three points define the plane")

    for layers in doc.root.find_all("layers"):
        for stmt in layers.statements:
            res.tick("layer materials")
            line = _line_of(doc.text, stmt.key_offset)
            if known_materials and stmt.name not in known_materials:
                res.add(ERROR,
                        f"layer material `{stmt.key}` is not defined in `materials`",
                        line, f"layers:{layers.label}")

            # thickness: read as tWords(2)=number, tWords(3)=unit with NO
            # bounds guard, then check_unit against "m" — a missing unit walks
            # off the word list, and a wrong dimension is a hard error
            words = stmt.value.split()
            res.tick("layer thicknesses")
            if len(words) < 2:
                res.add(ERROR,
                        f"layer `{stmt.key} : {stmt.value}` needs a thickness "
                        f"unit", line, f"layers:{layers.label}",
                        "the reader indexes the unit word unguarded — "
                        "a bare number is fatal")
            elif units_ok:
                actual = units.dimension_of(words[1], root)
                if actual is None:
                    res.add(ERROR, f"unknown unit `{words[1]}`", line,
                            f"layers:{layers.label}")
                elif actual != {"length": 1}:
                    res.add(ERROR,
                            f"layer thickness `{stmt.value}` is not a length",
                            line, f"layers:{layers.label}",
                            f"expected m or same dimension, `{words[1]}` is "
                            f"{units.describe(actual)}")


def _check_against_mesh(doc: Document, res: Result) -> str | None:
    """Check the deck's group ids against the mesh it names.

    Returns a reason string when the check could not run, which the caller adds
    to the skipped list. A missing mesh is NORMAL — the examples ship `.geo`
    and generate the `.msh` — so it never produces an error, only an honest
    statement that a class of check was unavailable.
    """
    from . import mesh as mesh_mod

    if doc.path is None:
        return "mesh checks — the deck was not read from a file"

    section = doc.root.section("mesh")
    filename = section.get("file") if section else None
    if not filename:
        return "mesh checks — no `mesh { file }` to resolve"

    info = mesh_mod.inspect(doc.path, filename)

    if info.kind == "missing":
        return (f"mesh checks — {info.path.name} is not present "
                f"(normal: the examples generate it from .geo)")
    if info.kind == "bfm":
        return f"mesh checks — {info.detail}"
    if info.kind == "unreadable":
        return f"mesh checks — {info.path.name}: {info.detail}"
    if not info.usable:
        return f"mesh checks — {info.detail or 'mesh not parsed'}"

    res.tick("mesh id existence")

    # A sibling .bfm is preferred by the factory when its stamps match, so the
    # ids actually used may come from there rather than from the .msh named.
    if info.sibling_bfm is not None:
        res.add(NOTE,
                f"a sibling `{info.sibling_bfm.name}` exists and is preferred "
                f"when its stamps match", 1, "mesh",
                "ids below are checked against the .msh, which may not be what runs")

    topo = doc.root.section("topology")
    if topo is not None:
        for sub in topo.sections:
            for key, pool, what in (("block", info.block_ids, "block"),
                                    ("blocks", info.block_ids, "block"),
                                    ("sideset", info.sideset_ids, "sideset"),
                                    ("sidesets", info.sideset_ids, "sideset")):
                raw = sub.get(key)
                if raw is None:
                    continue
                missing = [i for i in mesh_mod.expand_ids(raw) if i not in pool]
                if missing:
                    res.add(ERROR,
                            f"{what} id(s) {_brief(missing)} are not in "
                            f"{info.path.name}",
                            _line_of(doc.text, sub.header_offset), sub.path,
                            f"mesh has {what}s {_brief(sorted(pool))}"
                            if pool else f"mesh has no {what}s at all")

    _check_dimension_rules(doc, res, info)

    for sec in doc.walk():
        if not _under_boundary_conditions(sec):
            continue
        for key in ("sideset", "sidesets"):
            raw = sec.get(key)
            if raw is None:
                continue
            missing = [i for i in mesh_mod.expand_ids(raw)
                       if i not in info.sideset_ids]
            if missing:
                res.add(ERROR,
                        f"sideset id(s) {_brief(missing)} are not in "
                        f"{info.path.name}",
                        _line_of(doc.text, sec.header_offset), sec.path,
                        f"mesh has sidesets {_brief(sorted(info.sideset_ids))}")

    return None


def _check_dimension_rules(doc: Document, res: Result, info) -> None:
    """Rules whose predicate is the MESH dimension.

    `length` on a 2-D voltage BC is the only one today. It sat in T3 until
    building the mesh reader made clear that 2-D is a property of the mesh, not
    of the deck — no predicate over input can answer it.
    """
    if info.dimension is None:
        return
    res.tick("dimension rules")

    for sec in doc.walk():
        if not _under_boundary_conditions(sec) or sec.type != "voltage":
            continue
        if info.dimension == 2 and not sec.key_exists("length"):
            res.add(ERROR,
                    "a 2-D voltage boundary condition requires `length`",
                    _line_of(doc.text, sec.header_offset), sec.path,
                    f"{info.path.name} is 2-D; the amplitude is scaled by "
                    f"1/length")


def _brief(ids, limit: int = 8) -> str:
    """A readable id list. A mistyped range can name a hundred ids, and
    printing all of them buries the one line the reader needs."""
    ids = list(ids)
    if len(ids) <= limit:
        return "[" + ", ".join(map(str, ids)) + "]"
    shown = ", ".join(map(str, ids[:limit]))
    return f"[{shown}, … and {len(ids) - limit} more]"


def _check_landmines(doc: Document, res: Result) -> None:
    """Traps the C++ does not catch and the reference documents as pitfalls."""
    for sec in doc.walk():
        if sec.type != "background":
            continue
        raw = sec.get("direction")
        if raw is None:
            continue
        res.tick("background direction")
        parts = [p for p in re.split(r"[,\s]+", raw) if p]
        if any(("." in p or "e" in p.lower()) for p in parts):
            res.add(ERROR,
                    f"`direction : {raw}` has non-integer components",
                    _line_of(doc.text, sec.header_offset), sec.path,
                    "components are read through std::stoi, so 0.5 truncates "
                    "to 0 — use integer ratios")


def _as_list(value) -> list:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]
