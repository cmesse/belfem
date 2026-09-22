# Literature Skill and the Citation Key Index

**Date:** 2026-09-21 (continued 2026-09-22)
**Purpose:** Record the first BELFEM skill (`.claude/skills/literature/`), the machine-checked
citation key index added to `doc/literature_references.md`, `scripts/check_citations.py`, and the
citation defects the checker surfaced — a wrong publisher, a wrong book authorship, a year label
contradicting the DOI beside it, and 49 works cited in the tree that resolved to nothing.
**Module:** meta (`.claude/skills`, `doc`, `scripts`) + `src/physics/materials` (two comment lines)
**AIs involved:** Claude Opus 5 (1M context) (primary); Codex `gpt-5.6-luna`/medium (language sweep
over the new prose). **No jury round ran** — the session produced skill text, documentation and one
tooling script, which the exemption in `code-with-jury` covers; `scripts/check_citations.py` has
had no code audit and that is Christian's call to make.
**Status:** COMPLETE — skill, index and checker landed; the checker reaches **0 unresolved**.
O1 and O3 both resolved. Open: R6 (whether the checker becomes a gate) and O2 (marking which
indexed works are held locally).
**Literature References:** Messe et al. 2023 (its own bibliography settled the MUMPS and STRUMPACK
citations). The index holds **122 works** under 234 accepted spellings; three further property-data
sources went into `src/physics/materials/doc/material_property_sources.md`
**Verification:** reviewed, not verified — no build or test gate applies. Two checkers ran:
`scripts/check_citations.py` (**509 resolved, 0 unresolved, exit 0**) and
`scripts/check_doc_claims.py` (38/38 at the time of the first commit). DOIs were confirmed against the Crossref REST API or transcribed from an in-tree source
corroborated by a second in-tree entry; none was composed.
**Exchange:** `tmp/ai_exchange/citation_index_sweep.md` — the Codex language sweep; its outcome and
the two declined suggestions are distilled below

---

## What prompted it

Christian asked what a skill could do for the literature workflow: how to search for a topic, and
how to cite what is found — books by author, title, year, publisher, chapter and equation; papers
by DOI. Mid-session he added two things: the literature also carries variable-naming guidance that
belongs in the commenting style, and the references themselves look inconsistent.

Four rulings, all taken on the recommended option: one skill with three reference files; checked
into `.claude/skills/` in this repository; `doc/literature_references.md` extended and given a
checker rather than a second machine-readable twin; and the notation half cites the existing rule
in `doc/commenting_guidelines.md` §2 instead of adding one.

## What the survey found

The routing machinery already existed and is good: `literature/README.md` → per-directory
`index.md` → per-source nav `.md` (YAML header with `covers` / `does_not_cover` / `primary_for`,
a `## Search Patterns` section, a `## Citation Format` section) → the `.txt` extraction with
`--- Page N ---` markers. `literature/doc/notation_guide.md` even closes with a section titled
"Instructions for Claude Code". None of it loads unless a session remembers to walk the chain.

Three gaps were real:

1. **Paper nav headers carry `doi`, `journal`, `volume`, `pages`; book nav headers carry only
   `title` and `author`** — no year, edition, publisher or ISBN. For a book, the only citation data
   was the prose in `doc/literature_references.md`, and the four `books/magnets/` titles (Iwasa,
   Russenschuck, Barth, Glukhikh) were not in it at all.
2. **Two prose entries were wrong.** `Finite Element Procedures` 2nd ed. was attributed to Prentice
   Hall; the copyright page reads "© 2014 by Klaus-Jürgen Bathe, 2nd edition", ISBN
   978-0-9790049-5-7, published by the author, and states that Prentice Hall published the *first*
   edition. `Finite Element Exterior Calculus` was credited to Arnold, Falk and Winther; the book's
   own LOC record names Arnold alone (SIAM, CBMS-NSF 93, 2018) — Falk and Winther are co-authors of
   the Acta Numerica articles.
3. **`agelas2008.txt` is not the work the alias table cites.** The retired `F1` alias records
   Agélas, Di Pietro & Masson 2008 (FVCA V, DOI:10.1515/JNUM.2008.006); the extraction held locally
   is the two-author C. R. Acad. Sci. note, DOI:10.1016/j.crma.2008.07.015. Both are real. Recorded
   as a note under the index rather than silently resolved.

## What landed

**`.claude/skills/literature/`** — `SKILL.md` is a thin dispatcher (87 lines): is the library
present, the four-step routing chain, the citation format in one line, the notation hop, and five
rules that do not bend. The bulk is progressive-disclosure references:

| File | Covers |
|---|---|
| `references/searching.md` | the chain, which index, the page-aware `awk` search recipe, extraction reality (broken ligatures, search by stem), reading order, the stop rule, cross-checking a claim backwards |
| `references/citing.md` | per-medium formats, where the canonical data lives, twelve traps, what to do for a work not in the index |
| `references/notation.md` | the three naming systems, what the notation guide governs (equations in prose) and what it does not (kernel variable names), the symbol-annotation form, the collision table, the sweep checklist |

The rule that the licensed library must not be quoted at length into `tmp/ai_exchange/` is new. It
follows the existing nonfree posture — that channel goes to external vendors — and nothing stated
it before.

**Citation key index** in `doc/literature_references.md`, between
`<!-- BEGIN CITATION INDEX -->` / `<!-- END CITATION INDEX -->`: 78 rows, each `Key | Cite as |
Reference | DOI / ISBN`, the key cell carrying every accepted spelling separated by `;`. Paper rows
were generated from the nav YAML headers; book rows from the publication data in the extractions'
own front matter; the timestepping, Anderson, ghost-penalty and homology rows were harvested from
the module documents that already carried full references. No DOI was composed — entries without
one say `DOI not confirmed locally`.

**`scripts/check_citations.py`** parses that index plus the bibliography in
`src/physics/materials/doc/material_property_sources.md`, scans `src/` and `doc/` for author-year
tokens, and exits 1 on anything unresolved. Bibliographies are excluded from the scan: a table of
definitions is not a set of citations. `devlog/` and `todo/` are excluded too — records, kept as
written.

## What the checker found

First run: 137 distinct unresolved. After tightening the matcher (dates are not citations;
composite bibliography keys such as `Kamm-Alers 1964` register the `Kamm & Alers 1964` and
`Kamm 1964` spellings citations actually use; co-author surnames from the Authors column) and
adding the harvested rows, **402 citations resolve and 24 do not**. The residue is in
`todo/citation_index_completion.md`.

Spelling and year drift the run exposed, each now carried as an alias rather than a source edit:

| In the tree | In the bibliography | Handling |
|---|---|---|
| `Lüthi 1960` | `Luethi 1960` | alias |
| `Launey 1959` | `De Launay 1959` | alias |
| `O'Neal 1965` (code says 1964) | `O'Neil 1965` | alias; the author is O'Neal |
| `Freedman 2010`, `Hirani & Krishnamoorthy 2011`, `Geuzaine & Legros 1999` | indexed under their first authors | alias |
| `Gerritsen 1952` | `Gerritsen 1956` | **resolved 2026-09-22 — see below** |

`src/physics/materials/doc/resistivity_laws.md` ends with "Full citations:
`doc/literature_references.md`", and none of Rhyner 1993, Duron et al. 2004, Plummer & Evetts 1987
or Riva 2021 were there. All four are now indexed.

## O1 closed the next day: three Gerritsen papers, not a year conflict

Christian supplied the three DOIs, all confirmed against Crossref:

| DOI | Work | Year |
|---|---|---|
| 10.1016/S0031-8914(52)80222-8 | Gerritsen & Linde, "The electrical resistances of alloys of a noble metal and a transition metal", Physica 18:877-890 | 1952 |
| 10.1016/S0031-8914(53)80007-8 | Gerritsen, "The magnetoresistances of alloys of a noble metal and a transition metal at low temperatures", Physica 19:61-73 | 1953 |
| 10.1016/S0031-8914(56)90035-0 | Gerritsen & Linde, "Thermal conductivity of some dilute silver alloys", Physica 22:821-831 | 1956 |

There was never a conflict between two sources. `cl_Material_Silver.cpp:64` carried the **1956 DOI
on a line labeled 1952**, inside the thermal-conductivity block — and 1956 is the dilute-silver-alloy
conductivity paper, with `sample AG2T, RRR 407` one of its samples. The DOI was right; the year
label was wrong, at `:37` and `:64`. Both corrected on Christian's approval (comment-only, exempt
from the jury round). The bibliography row gained Linde as second author and the journal locus.

The 1952 and 1953 papers are deliberately **not** indexed: nothing in the tree cites them. Silver's
resistivity comes from Matula (Bloch-Grüneisen, 273.15 K anchor) and its magnetoresistance from
Lüthi, Strom-Olsen, Stout and Iwasa. They are recorded here so the three works can be told apart
the next time one is reached for.

The checker now reads 404 resolved, 23 unresolved.

## The completion pass (2026-09-22): 23 residual works down to 4

Worked in the planned order — materials and solvers, then Maxwell and element theory, then the gas
models. The method throughout: harvest whatever the citing document already states, then confirm
against Crossref, and accept a match only when title, authors, journal and year agree with what the
tree independently says.

**Nineteen resolved.** Most of the work was transcription, not research, because the tree already
held the answer somewhere other than where the citation sat:

- Every gas-model DOI was already in its own module header — `10.1063/1.3160306`,
  `10.1021/je301273j` (plus the 2022 erratum), `10.1063/1.3606499`,
  `10.1023/B:IJOT.0000022327.04529.f3`, `10.1063/1.555828`. `Leachman 2009` and
  `Penoncello & Lemmon 2009` are the **same paper**, matched twice out of one four-author line.
- `nedelec_derivation.md:487-489` carries full references for Zaglmayr's thesis and the
  EXODUS II report, neither of which has a DOI.
- `cl_Material_Lead.cpp:55` already carried `10.1179/mtlr.1961.6.1.1` — Köster & Franz, *Poisson's
  Ratio for Metals and Alloys*. The code spells the name `Koester`; both spellings are indexed.
- The MUMPS and STRUMPACK DOIs were cross-checked against references [24] and [26] of the BELFEM
  paper's own bibliography (`messe2023.txt:901-906`) rather than trusted from a search.

Crossref confirmed the rest, several of them satisfyingly exact: Giebe & Blechschmidt 1931 is
*"Über den Einfluss der Magnetisierung auf den Elastizitätsmodul bei Dehnungsschwingungen
ferromagnetischer Stäbe"*, which is precisely the Delta-E measurement the nickel header cites;
Hayasaka & Ito 2019 is the contact-probing current-transfer-length paper the contact-impedance
table describes.

One near-miss worth recording: a Crossref query for Grigorieff also returned **Calvo, Grande &
Grigorieff 1990**, *"On the zero stability of the variable order variable stepsize BDF-Formulas"*
(10.1007/BF01386395). That is arguably a closer source for the **order-dependent** ratio bounds
that `bdf_timestepping_theory.md` tabulates than the 1983 paper the docs cite. Not changed — the
docs cite 1983 and 1983 does cover variable-grid multistep stability. Flagged here only.

**Four could not be identified** and are logged as O3 in the todo: `Herrmann 1978`,
`Kuczmann & Iványi 2008`, `Meunier 2008`, `Robinson 1978`. Each is a bare surname and year with no
fuller reference anywhere in the tree. Two of them are described in their own document as coming
"from the notes", which is not a source this repository can reach.

**A checker blind spot, closed.** The matcher required `Author YYYY` and so never saw
`Author (YYYY)`. Widening it to accept a parenthesized year raised the resolved count to 475 and
exposed a second population of about 26 works, nearly all sitting in module reference lists that
already carry title, journal, volume and pages. Several of their DOIs are in the gasmodels headers
but in a different file from the list that names the work. That is R7, and it is mostly
transcription again.

The live count therefore reads 30, not 4: 4 from the original residue and ~26 newly visible. The
widening is kept because it is more correct, not because the number is nicer.

## R7 (2026-09-22): the parenthesized population, and three matcher defects

Harvesting the ~26 works the widened matcher exposed went the same way as the first pass: the tree
already knew the answers, filed apart from the citations that needed them. The module reference
lists (`gasmodels_usage_guide.md:1139-1162`, `gasmodels/doc/README.md:829-841`,
`math/graph/doc/README.md:202-209`, `interpolation_usage_guide.md:1877-1878`,
`iwg/doc/README.md:260-261`, `maxwell/doc/README.md:538-542`) carry title, journal, volume and
pages for nearly every one. Five DOIs were sitting in module headers in a *different file* from the
list naming the work — Soave 1972 in `cl_GM_EoS_AlphaFunctionFactory.cpp:33`, Span et al. 2000 in
`cl_GM_EoS_Nitrogen.hpp:25`, Schmidt & Wagner 1985 in `cl_GM_EoS_Oxygen.hpp:25`, Setzmann & Wagner
1991 in `cl_GM_EoS_Methane.hpp:24`, and both Zienkiewicz & Zhu papers in
`postprocessor_recovery_theory.md:542-552`. Each was accepted only because a second, independent
in-tree entry agreed on journal, volume and pages.

Six needed Crossref: Peng & Robinson 1976, Wilke 1950, Cuthill & McKee 1969, Liu & Sherman 1976,
Gibbs, Poole & Stockmeyer 1976, Dunavant 1985. Nine are books and carry edition and publisher
instead of a DOI.

**Three matcher defects, all found by reading the output rather than by testing:**

| Reported | Actually | Fix |
|---|---|---|
| `Nédélec's 1980` | a possessive of a work already indexed | strip a trailing `'s` |
| `J.-C 1980` | the initials in "Nédélec, J.-C. (1980)" | require three alphabetic characters in a surname |
| `SVEHLA 1994` | a line of `share/fluid/trans.inp` quoted inside a fenced block | skip fenced blocks in documents |

The last is the one worth keeping in mind: a fenced block in a document quotes a data file, a deck
or program output, and whatever looks like a citation inside it belongs to the quoted thing. That
change also dropped ten previously counted citations, which is the correct direction.

Two spellings became aliases rather than source edits: `Sousa 2022b` (the tree elsewhere writes
"de Sousa Alves") and `Zunino 2009` (indexed under its first author, Ern).

**The checker now reads 505 resolved, 4 unresolved**, and all four are the O3 works no document in
the tree identifies. Every author-year citation in `src/` and `doc/` that names a findable work now
resolves to a full reference.

## O3 closed the same day: the checker reaches zero

Christian supplied the two books and ruled on the fourth; the first settled on evidence.

**`Herrmann 1978` does not exist.** A Crossref author search over 1970-1982 returns no Herrmann
superconvergence paper — the only FEM Herrmann of 1978 is *Finite Element Analysis of Contact
Problems*, unrelated. The theorem at `postprocessor_recovery_theory.md:603` (derivatives
superconvergent at the Gauss points, `O(h^(p+1))` there against `O(h^p)` elsewhere) is **Barlow
1976**, *Optimal stress locations in finite element models*, 10.1002/nme.1620100202. Attribution
corrected, Barlow added to that document's reference list as item 5 with items 5-10 renumbered to
6-11. The near-miss that probably produced the wrong name is Herrmann **1972**, *Interpretation of
Finite Element Procedure as Stress Error Minimization Procedure* (10.1061/JMCEA3.0001684) — the
variational companion to the result, not the `O(h^(p+1))` statement itself. Not indexed, because
nothing cites it.

`Kuczmann & Iványi 2008` is *The Finite Element Method in Magnetics*, Akadémiai Kiadó, Budapest,
ISBN 978-963-05-8649-8. `Meunier 2008` is *The Finite Element Method for Electromagnetic Modeling*,
ISTE/Wiley, 10.1002/9780470611173 — confirmed against Crossref, which also returned the print ISBN.

`Robinson 1978` is indexed as **PR78** on Christian's ruling that PR78 is the standard name. It has
no DOI and will not get one: the revision is grey literature rather than a journal article. The row
describes the m(ω) branch for ω > 0.491 that `cl_GM_EoS_AlphaFunctionFactory.cpp:53-57` implements
and points at the base equation, Peng & Robinson 1976. No source edit — `Robinson 1978`,
`Peng Robinson 1978` and `PR78` are all accepted spellings on the row, per the minimal-rename rule.

**`scripts/check_citations.py` now reports 509 resolved, 0 unresolved, and exits 0.** Every
author-year citation in `src/` and `doc/` resolves to a full reference. That makes R6 — whether the
checker joins the session checklist or the post-commit hook — a real choice for the first time,
since it no longer fails by design.

## Status

Reviewed, not verified. The skill files, the index and the devlog are text; the checker is a new
Python script and has had no jury round — it is tooling, like `check_doc_claims.py`, and Christian
decides whether that class takes one. `scripts/check_citations.py` **exits 0**: 509 citations
resolve and none do not.

A Codex language sweep ran over the new prose (luna/medium, slug `citation_index_sweep`, exchange
entry stamped; `git status` identical before and after, so the round stayed read-only). It reported
no suspected errors and returned four rewrites. Three were applied whole or in part — "parses" →
"reads", a tighter apposition for `agelas2008.txt`, "Prentice Hall published only the first
edition", and the Arnold parenthetical dissolved into a clause. Two suggestions were declined:
"They live in" → "They are in" flattens without clarifying, and "published by the author in
Watertown MA **in 2014**" asserts a publication year where only the copyright line was verified —
the entry keeps "copyright 2014".

Not done: the four magnets books are in the index but still absent from the prose sections of
`doc/literature_references.md`, and no sweep has run over the prose added after the first commit.
