# The Pragmatic Programmer (20th anniversary edition) added to the literature library

**Date:** 2026-09-16
**Purpose:** Index Thomas & Hunt 2020 in `literature/books/coding/` and map it onto BELFEM's conventions
**Module:** `literature/` (separate, proprietary repository)

## What was done

- `literature/books/coding/thomas2020.md` written from the book navigation template: YAML header, a chapter-and-topic map carrying the extraction line of each of the 53 topics, question routing, reading paths, the full list of the 100 tips by topic, and a two-part BELFEM mapping.
- `literature/books/coding/index.md` extended to four books (fourth division column, flowchart branches, topic-table rows, tags, divergence items 7–8, per-file pitfalls row); `literature/README.md` tree, routing and version history (v2.7); the out-of-scope sentence in `literature/doc/book_fem_coverage_matrix.md`.
- **Finding for `CLAUDE.md` "Error Handling":** this is the first book in the library that states BELFEM's three-tier policy directly. Topic 25 draws the assertion / error-handling boundary ("Don't use assertions in place of real error handling"), Topic 24 argues crash-early for `BELFEM_ERROR`, and the "turn off only those assertions that really hit you" carve-out in Topic 25 is the published form of compiling `BELFEM_ASSERT` out on hot paths. No BELFEM document cites the book yet; a citation in `doc/coding_philosophy.md` would be the natural place (not added — read-only outside `literature/` this round).
- **Extraction facts recorded:** ebook build (P1.0, 2019-09-13), so `--- Page N ---` markers are reflow pages with no relation to the printed book — the guide forbids page citations and anchors on Topic and Tip numbers. Lines break mid-sentence, italic words sit on their own lines, inline mathematics was dropped (Topic 39's Big-O material is incomplete), chapter headings split across two lines. No ligatures.

## Not done

- No citation added to the main repository; `doc/coding_philosophy.md` External References and the `CLAUDE.md` error-handling section are the candidates.

Reviewed, not verified — documentation only, no code touched.

## Follow-up the same session: `oliveria2006` → `oliveira2006`

On Christian's request the misspelled stem was corrected in the literature repository: both files `git mv`ed, every path and `book_id` in the library retargeted, the two "filename misspells the author" caveats deleted, and the one pointer in this repository (`doc/commenting_guidelines.md`, Sources) updated. The library's version history keeps the old stem in its v2.2 record; v2.8 records the rename.

## Follow-up the same session: inline formulas restored

Christian supplied the original PDF (`tmp/thomas2020.pdf`). It is a calibre build of the publisher's ebook and renders every inline formula as an empty frame, so the PDF could not supply the text; `pdftotext` reproduces the same gaps. What it did supply was the *location* of every gap: `pdfimages -list` shows the conversion's missing images as 16×16 placeholders, and the pages whose placeholder count exceeds their `images/*.png` figure lines are exactly the pages with dropped formulas (95, 104, 174, 343–353, 397, 495; 275, 276 and 420 are lost figures). Each page was rendered and read, and the 83 gaps filled from the surrounding sentences — Topics 10, 11, 20, 39, 44 and Answers 29–30, the Figure 3 Big-O table and the displayed `O(n²/2 + 3n)` identity included. Line count and every line number in the guide unchanged; editorial note on line 2 of the `.txt`. Confidence high for the Big-O forms, medium for two spellings (`m × n`, `O(log₁₀ n)`), stated in the guide. Library README v2.9.

