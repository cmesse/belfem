# Citing

**Purpose:** The house citation format, where the canonical data lives, and the traps that have
already produced wrong citations in this tree.

---

## The format

Author and year, then the **exact** section or equation. The exactness is the whole point — a
citation without a section is a gesture, not a trail.

| Medium | Form | Example |
|---|---|---|
| Paper | `Author et al. YYYY, §X` or `, Eq. N` | `Alves et al. 2022b, Eq. 19` |
| Paper, one author | `Author YYYY, §X` | `Aavatsmark 2002, §3` |
| Paper, two authors | `A & B YYYY, §X` | `Klausen & Winther 2006, §4` |
| Book | `Author, §X.Y` or `, Ch. N` | `Bathe, §8.6` · `Hughes, Ch. 4` |
| Book by Topic | `Author & Author YYYY, Topic N` | `Thomas & Hunt 2020, Topic 25` |
| Lecture notes | `Author YYYY, §X` | `Evans 2017, §3.2` · `Felippa IFEM, Ch. 11` |

Books are cited by author and section; the year lives in the full reference, and is added
(`Bathe 2016, Ch. 8`) only where the surrounding text also cites papers and the bare surname would
read as one. Either form is acceptable — what is not acceptable is a year that contradicts the
citation key index.

**A number carries its DOI.** A value lifted from a source — a coefficient, a measured property, a
benchmark result — gets the DOI beside it, so the next reader can reach the table it came from:

```cpp
// Kohler magnetoresistance fit: Iwasa 1993, 10.1016/0011-2275(93)90199-X
```

In a source comment the citation is a fragment, never a sentence, and never a derivation:

```cpp
// Static condensation (Messe et al. 2023 — avoids zeros on the diagonal)
// Interface conditions from Alves et al. 2022b, Eq. 19
// Newton-Raphson convergence strategy (Bathe, §8.6)
```

`doc/commenting_guidelines.md` §2 admits the citation as one of the six things a comment may say.
It also forbids the history that often travels with one: no dates, no reviewer names, no task IDs.

## Where the canonical data lives

| Need | Source of record |
|---|---|
| Full reference, DOI, ISBN, cite-as key | `doc/literature_references.md` → citation key index |
| Material property data provenance | `src/physics/materials/doc/material_property_sources.md` → Bibliography |
| Per-source citation examples | `literature/<dir>/<id>.md` → `## Citation Format` |
| Paper metadata (journal, volume, pages, DOI) | the nav file's YAML header |

The citation key index is machine-checked. `scripts/check_citations.py` reads it, scans `src/` and
`doc/` for author-year tokens, and reports any that do not resolve. Run it after adding a citation,
the same way `scripts/check_doc_claims.py` is run after editing a convention document.

Book nav files carry only `title` and `author` in their YAML headers — no year, publisher, edition,
or ISBN. For a book, the citation key index is the only source for those fields.

## Traps

These have each produced a wrong citation in this tree or in an AI review of it.

| Trap | What goes wrong | Rule |
|---|---|---|
| **Alves 2022a / 2022b** | Two papers, same authors, same year, different models (3D h-based vs 2D thin shell in h-φ) | Always carry the letter |
| **Badel** | Local file is `badel2021` — the HAL deposit year | Cite **Badel et al. 2019** |
| **Lucchini** | Local file is `luccini2025`, a misspelling | Cite **Lucchini 2025** |
| **Arsenault 2023** | The 2026 erratum corrects its air-domain coupling | Cite the erratum alongside it |
| **Messe** | `Messe et al. 2023` (five authors, SUST) and `Messe 2022` (single author, conference) are different works | Never drop `et al.` from the 2023 paper |
| **Iwasa** | `Iwasa 1993` is a Cryogenics paper on magnetoresistance; `Iwasa 2009` is the magnets textbook | Same author, unrelated works |
| **Riva** | `Riva 2021` is the EPFL thesis (8754); `Riva et al. 2023` is the Sparselizard paper | Year distinguishes them |
| **Zienkiewicz** | Vol. 1 (2013) and Vol. 2 (2014) are separate books | Name the volume |
| **Bathe printing** | The 2nd edition is self-published by the author, © 2014; the copy here is the 2016 4th printing. Prentice Hall published the **first** edition only | Never attribute the 2nd edition to Prentice Hall |
| **Hughes pagination** | The extraction is the 1987 Prentice-Hall printing; the Dover 2000 reprint is unabridged, so chapter and section numbers agree | Cite by chapter and section, not page |
| **Thomas & Hunt pages** | The extraction's page markers are reflow pages, not printed pages | Cite by Topic (1–53) and Tip (1–100) |
| **`paperN` aliases** | Retired 2026-08-11 | Never write a new one; the decoder in `doc/literature_references.md` exists only to read old devlogs |

## A work that is not in the index

Common, and not a reason to leave a number uncited. In order:

1. **Is it material property data?** Then it belongs in
   `src/physics/materials/doc/material_property_sources.md`, not in the general reference list —
   that document's bibliography already carries about sixty such sources with DOIs.
2. **Is it in the library but not the index?** Take the full reference from the nav file's YAML
   header and add a row to the citation key index.
3. **Is it neither?** The work is behind a paywall or simply not held. Ask before spending effort
   on it; purchased or retrieved PDFs land in `tmp/papers/`. Until then, cite what you can verify
   and mark the rest as uncited — an honest gap beats a plausible DOI.

**Never compose a DOI.** DOIs are opaque strings; a reconstructed one resolves to the wrong paper
or to nothing, and it looks exactly as authoritative as a real one. If you do not have the DOI in
front of you, write the citation without it and say the DOI is missing.
