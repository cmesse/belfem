# Cross-review: python/ input-file checker (belfem_conf)

**Date:** 2026-08-13
**Purpose:** Three-AI jury audit (Claude Fable pre-registered; Codex + Grok blind) of the
Opus-written `python/belfem_conf` package before it is relied on.
**Exchange:** `tmp/ai_exchange/review_python_conf_tool.md` (full record; this is the distillation)

## Outcome

29 findings, 28 CONFIRMED by source trace, 1 PLAUSIBLE. No P0 — the tool is
read-only and has no data-loss path. Ten P1s, all mechanical divergences
against the C++ or the schema; no physics questions. All three voices agree the
core design is sound: the lossless span parser + probe gate (roundtrip 6/6
byte-identical), the YAML-tree schema walker, elementary-tag mesh ids, and the
"say what you skipped" policy.

## P1 list (fix candidates, pending approval)

1. `required: true` below top level is unenforced — a deck with no
   `mesh { file }` or no `simulation time` validates clean (Codex).
2. `value_runtime_status` ignored — `crc` / `crank-nicolson` / `galerkin` pass
   although the schema says a validator must reject them (Codex).
3. Unknown keys are silently skipped for boundary-condition subsections, and
   the whole `circuit` subtree is unvalidated — README claims otherwise (Grok).
4. `unit_token` (mesh.unit) and layer thicknesses never get a dimension/unit
   check (Grok + Codex).
5. drift's code→schema scan is line-local and misses the multi-line
   `get_string(\n "unit")` in cl_MaxwellFactory.cpp — the one key parsed only
   that way (Claude class / Grok instance).
6. BC section types unvalidated — a typo'd type or the rejected
   `background dirichlet` passes (Codex).
7. HTS `one_of` sees only one of the three builtin-selection paths of
   cl_MaterialFactory.cpp (Codex).
8. `IndexError` escapes `mesh.inspect()` on a truncated `.msh` — crash instead
   of the designed "unreadable" verdict (Claude).
9. README / `__init__.py` / validate.py docstring contradict the code
   (validator "planned"; T3 described as open though CLOSED 2026-08-11) (Grok).
10. The tool's own first run found real repo drift: `metis nodendp`
    (cl_SolverParameters.cpp:75) is absent from the schema — an Input Contract
    violation to fix in input_schema.yaml AND input_file_reference.md — and the
    two `'"read_signed_sidesets"'` anchors are quoted as key literals though
    they name a function (schema authoring fix).

P2s (18): unit-string preprocessing (µ/²/³, first-`/`-only) not mirrored;
parser is cursor-based where the C++ is line-oriented; gmsh element table
omits PYRA5/TRI10/LINE4/TET20 and uses a z-tolerance where the C++ compares
exactly; id ranges (`1:3`) not expanded in counting checks; CONSUMERS list is
handwritten and misses fn_mesh_config_tag.hpp; duplicate-section lookup is
first-match vs C++ last-wins; notes rendered as warnings; and assorted
hygiene. Full table in the exchange file.

## Tooling note

`cross_review.sh` cannot take a directory, and inlining the 2,863-line bundle
blew the argv limit in `ask_codex.sh` (`codex exec: Argument list too long`,
exit 126) — the Codex leg was rerun blind with a pointer prompt (file list
only), which worked. Consider having `ask_codex.sh` pass long prompts via file
or stdin rather than argv.

## Fixes applied (same day, approved scope: python/ + contract docs)

Christian granted autonomy within the Python tool, C++/Fortran untouched. All
ten P1s and the actionable P2s are fixed; the full per-finding list is the
17:10 entry in the exchange file. Highlights:

- Required sections/keys now enforced at every nesting level; runtime-refused
  enum values (`crank-nicolson` et al.) rejected; BC type names and the whole
  circuit contract validated; unknown-key checks live on BC/topology
  subsections and material bodies.
- Unit checking now mirrors the C++ exactly: `unit_to_si` preprocessing
  (µ/²/³, first-slash-only), `unit_token` values, layer-thickness units, and
  the `create_key`/`get_value` value shape — the old "a bare number is legal
  and dimensionless" rule was a false pass; `get_value` aborts on it.
- drift's code→schema scan is whole-file, so the two-line `get_string(
  "unit" )` call is finally seen; the schema anchor arithmetic is per-token.
- Mesh reader hardened (IndexError → "unreadable", exact `min(z)==max(z)`,
  2.2/4.1 only, binary declared unread, full element-id table).
- **New bug found by the regression probes, missed by all three reviewers:**
  a comment between two statements leaked into the next statement's key —
  the parser's comment mask guarded structure but not the `key`/`value`
  accessors. Fixed, with diagnostics now pointing at the key's real line.
- Contract docs: `metis nodendp` added to schema + reference (§4.1); the two
  quoted `read_signed_sidesets` anchors fixed. `drift` is now **clean**
  (124/124 anchors, 111 literals all known).

Gate: roundtrip 6/6 byte-for-byte, drift clean, all six examples ok — under
both python3.14 and the SCLS python3.9. Broken-deck probes confirm every new
check fires. Still out of reach (reported by `-v`): binary .msh, id_groups
arity, mesh-state gates, material-shape prose.

## Status

Review complete, fixes applied and gated. The `ask_codex.sh` argv overflow was
fixed 2026-08-14: the prompt now goes through a temp file and stdin
(`codex exec -`), matching `ask_grok.sh`'s `--prompt-file` pattern. Gated with
a stub binary at 10 MB (the old path died at ~128 KiB) and a live end-to-end
Codex call; the stdin redirect also removes the known hang-on-open-stdin
pitfall. Remaining follow-ups: the T4 items above.
