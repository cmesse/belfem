# Clean Code, second edition, added to the literature library

**Date:** 2026-09-16
**Purpose:** Index Martin's *Clean Code*, 2nd ed. (2025) in `literature/books/coding/` and move the main repository's citations from the 2008 first edition to it
**Module:** `literature/` (separate, proprietary repository), `doc/`

## What was done

- `literature/books/coding/martin2025.md` written from the book navigation template: YAML header, a chapter table carrying printed page and extraction line for every chapter (sections are unnumbered, so the citation form is "Martin 2025, Ch. N 'Section', p. P"), question routing, reading paths, and a two-part BELFEM mapping — thirteen rules the book backs, five it argues against (prefix encodings, output arguments, status returns, kernel function size, test-first).
- `literature/books/coding/index.md` rewritten for three books; `literature/README.md` tree, routing and version history (v2.6); the one out-of-scope sentence in `literature/doc/book_fem_coverage_matrix.md`.
- **Extraction facts recorded:** printed page = PDF page − 33; effectively no ligatures (unlike the two Cambridge extractions); running headers carry a spurious space (`Cl ean Code`); chapter titles are glued to the preceding sentence; code listings collapsed onto one line.
- **Main repository, on Christian's approval ("update every to the new edition"):** `doc/commenting_guidelines.md` Sources entry moved to the 2nd edition — Ch. 4 → **Ch. 5**, publisher corrected, the good/bad taxonomy re-listed from the 2025 contents. Two substantive edition changes recorded: TODO comments moved from Martin's good list to his bad list ("TODO means Don't Do"; he no longer checks them in), so §7's owned-TODO rule is now stated as the narrow exception it is; the closing-brace category no longer exists in the book, so the §3 row is a house rule, not a citation. The opening quotation corrected to the book's wording ("to compensate for our failure to express ourselves in code"); "inobvious connection" → "unobvious connection" (the 2nd-edition section title). The dangling last sentence of the Sources block completed. `doc/coding_philosophy.md` External References gained a Martin 2nd-ed. bullet beside Ousterhout.

## Not done

- `doc/commenting_guidelines.md` is still untracked and still marked proposed; the jury round and Codex sweep it owes are unchanged by this.
- (Follow-up the same session, on Christian's call: because nothing was committed yet, the earlier devlog `dl20260916_commenting_guidelines.md` and its `devlog/README.md` line were updated too, so no document in either repository cites the first edition by chapter number. The `tmp/martin/` extract itself is first-edition text and is ephemeral.)

Reviewed, not verified — documentation only, no code touched.
