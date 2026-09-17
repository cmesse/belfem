Language sweep of rewritten source comments (read-only; return edits, do not write files).

You are given a unified diff of a comment-only commit in BELFEM. The rules the new comments
must follow are `doc/commenting_guidelines.md`; read §2, §3, §8 and §10 first.

Audience of the comments: a developer who reads English as a second language and did not sit
in the session that wrote them. Only the lines that the diff ADDS are in scope. Do not comment
on deleted lines, and do not comment on code.

Return a list of edits, each as: file, the exact added comment line(s), the proposed
replacement, and one clause on why. Order by how much the edit helps the reader. If an added
comment is fine, do not list it. Do not return a rewritten file.

You MAY:
- shorten a sentence, split one that carries a decision and a reason and a consequence
- replace an idiom or a word a non-native reader may not know
- fix punctuation, articles, agreement, US spelling
- flag an added comment that still carries a date, a reviewer name, an audit round label, an
  incident or debt-register ID, a `tmp/` path, or an experiment log, under a separate heading
  "Guideline violations" — those are not language edits
- flag an added comment whose meaning you could not determine (say so; do not guess)

You may NOT:
- change a technical claim, a unit, a citation, a symbol, or the side a comment takes — if you
  believe a claim is wrong, list it under "Suspected errors" with your reason and leave it out
  of the edit list
- touch code, identifiers, macro names, or anything outside comment text
- add a comment where the diff removed one
- propose Doxygen markup, HTML, or a longer comment than the one given

Keep the house voice: direct, one idea per sentence, full sentences for a reason or a
contract, fragments for a citation or a unit, no exclamation, no hedging.
