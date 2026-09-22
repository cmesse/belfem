---
name: literature
description: Route a BELFEM question to the right reference, find the passage, and cite it in house format. Use when a task needs the theory, algorithm, or validation basis behind the code (h-phi, thin shell, cuts and cohomology, MPFA, element technology, quench, transient solver control, material property data); when writing or checking a citation in a source comment, in doc/, or in a module guide ("how do I cite this", "what is the DOI for X", author-year format, is this reference right); when a kernel comment needs a symbol annotation translated from a paper's notation into BELFEM's unified notation; or when CLAUDE.md's literature-first rule applies before judging a formulation, an algorithm, or a validation choice.
---

# Literature

Three jobs, one workflow: **route** a question to a source, **read** the passage, **cite** what
you read. The routing tables and per-source navigation files already exist — this skill is the
procedure that uses them, not a copy of them.

## Before anything else: is the library present?

```bash
ls literature/README.md 2>/dev/null && echo PRESENT || echo ABSENT
```

`literature/` is a separate, proprietary repository, `.gitignore`d here. It is often absent.

| | Present | Absent |
|---|---|---|
| Routing | `literature/README.md` → per-directory `index.md` | `doc/literature_references.md` routing tables |
| Passage | the `.txt` extraction | **not available** |
| Citation | the nav file's `## Citation Format` section | `doc/literature_references.md` citation key index |

**Absent means you cannot read the source.** Say so. State that the claim rests on code analysis,
name the work and the section that would settle it, and stop — do not reconstruct an equation from
memory and present it as the paper's. A remembered equation attributed to a numbered equation in a
paper you did not open is the worst failure mode this skill exists to prevent.

## Routing: four steps, in order

1. **Classify the question.** CLAUDE.md's "Literature Routing Strategy" table maps question type to
   work. It is the fastest hop and it is in context already.
2. **Open the directory index**, not the paper: `literature/papers/fem/index.md`,
   `papers/fvm/index.md`, `papers/topology/index.md`, `books/fem/index.md`, `books/math/index.md`,
   `books/coding/index.md`, `books/magnets/index.md`, `lectures/index.md`. Each opens with a
   task → source routing table.
3. **Open the nav file** for the source (`<id>.md` beside the `.txt`). Read its YAML header
   (`covers`, `does_not_cover`, `primary_for`, `pairs_with`) before the prose — it answers "is this
   the right source" in one screen, and `does_not_cover` prevents the common miss of hunting a
   topic in a book that never treats it.
4. **Search the `.txt`** with the nav file's `## Search Patterns` section, which gives that
   source's own section, equation, and figure conventions.

Details, search recipes, and the stop rule: `references/searching.md`.

## Citing

Author-year, always with the exact section or equation, per CLAUDE.md:
`Messe et al. 2023, §2.7` · `Alves et al. 2022b, Eq. 19` · `Bathe, §8.6` · `Thomas & Hunt 2020, Topic 25`

A number lifted from a source carries its DOI next to it. The citation key index in
`doc/literature_references.md` is the source of record for the full reference, the DOI, and the
ISBN. Material property data has its own bibliography in
`src/physics/materials/doc/material_property_sources.md`.

`scripts/check_citations.py` verifies that every author-year token in `src/` and `doc/` resolves to
an entry. Run it after adding a citation.

Per-medium formats, the traps (Alves a/b, Badel's deposit year, Lucchini's misspelled file,
retired `paperN` aliases), and what to do for a work that is not in the index:
`references/citing.md`.

## Notation and symbol annotations

A kernel written in a paper's notation earns a symbol annotation — the one comment form
`doc/commenting_guidelines.md` §2 admits at the end of a declaration. When you write one, translate
the source's symbol into BELFEM's unified notation and say which source's notation the block
follows, so the reader knows whether `B` is a strain-displacement matrix or a flux density.

`literature/doc/notation_guide.md` holds the translation tables (N, B, E, C; field variables;
per-book warnings). `references/notation.md` bridges them to the naming rules in CLAUDE.md and
`doc/commenting_guidelines.md`.

## Rules that do not bend

- **Never cite a work you did not open.** Cite the passage you read, at the granularity you read it.
- **Page numbers come from the `--- Page N ---` markers** in the extraction, and only for sources
  whose nav file does not warn that its pages are reflow artifacts (Thomas & Hunt is cited by Topic
  and Tip, never by page).
- **Reviewed, not verified.** Reading a paper and agreeing with the code is a review. "Verified"
  needs an executable gate that ran (protocol §11).
- **The library is licensed, not ours.** Extended verbatim text from `literature/` must not go into
  `tmp/ai_exchange/` — that channel is routed through external AI vendors. Quote a line where a
  claim turns on the exact wording; otherwise cite section and equation and paraphrase.
- **Documents state facts inline.** A `doc/` file cites literature, never a devlog or a task ID.
