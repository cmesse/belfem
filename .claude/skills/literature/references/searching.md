# Searching the Library

**Purpose:** How to get from a question to a passage you can cite, and when to stop looking.

---

## The chain

```
question type          CLAUDE.md "Literature Routing Strategy" table
   ↓
directory              literature/README.md  →  <dir>/index.md
   ↓
source                 <dir>/<id>.md   (YAML header: covers / does_not_cover / primary_for)
   ↓
passage                <dir>/<id>.txt  (searched with that nav file's patterns)
   ↓
citation               references/citing.md
```

Skipping to the `.txt` is the common mistake. The nav file costs one read and tells you whether
the source treats the topic at all — `does_not_cover` is there precisely so you do not spend twenty
minutes hunting edge elements in Bathe.

## Which index

| Directory | Holds | Index |
|---|---|---|
| `papers/fem/` | h-φ, thin shell, HTS material models, quench, transient control | `papers/fem/index.md` |
| `papers/fvm/` | MPFA, transmissibility, coercivity, mixed FVM | `papers/fvm/index.md` |
| `papers/topology/` | cohomology bases, cut generation, optimal-cut complexity | `papers/topology/index.md` |
| `books/fem/` | element technology, locking, inf-sup, solvers, Maxwell, FEEC | `books/fem/index.md` |
| `books/math/` | formulas and identities, graph algorithms, topological EM theory | `books/math/index.md` |
| `books/magnets/` | magnet physics, conductors and cables, accelerator fields, fusion magnets | `books/magnets/index.md` |
| `books/coding/` | the published rationale behind BELFEM's own conventions | `books/coding/index.md` |
| `lectures/` | functional-analysis foundations, step-by-step FEM | `lectures/index.md` |
| `doc/` | coverage matrices, `notation_guide.md`, nav templates | — |

The coverage matrices in `literature/doc/` answer the comparative question — *which of the eleven
FEM books treats this best* — faster than opening eleven nav files.

## Searching a `.txt`

Every extraction carries `--- Page N ---` markers. This recipe returns the page with the hit, which
is what a citation needs:

```bash
awk -v pat="checkerboarding" '
    /^--- Page [0-9]+ ---/ { page = $3 }
    $0 ~ pat { printf "p.%s\tline %d\t%.90s\n", page, NR, $0 }
' literature/papers/fem/messe2023.txt
```

For a case-insensitive or regex search, `rg -n` first and recover the page after:

```bash
rg -n 'inf-sup' literature/books/fem/bathe.txt | head
awk 'NR<=17422 && /^--- Page [0-9]+ ---/ { p=$3 } END { print "page", p }' literature/books/fem/bathe.txt
```

**Extraction reality.** These are ASCII extractions of typeset PDFs. Ligatures break
(`ﬁnite`, `satisﬁes`), math becomes soup, hyphens split across lines, and accented names mangle
(`Nicol \`o Riva`). Consequences:

- Search for a **stem**, not a word: `convergen`, `condens`, `checkerboard`.
- Search for the **equation number** rather than the equation: `(14)`, `Eq. 19`.
- A zero-hit search proves nothing until you have tried the stem and the nav file's own suggested
  term. The nav file's `## Search Patterns` section lists terms that are known to work in that file.
- Never paste an equation out of a `.txt` into code or documentation without reading it back
  against the surrounding prose — the extraction loses superscripts silently.

## Reading order for an unfamiliar topic

Two angles beat one, and the second angle is cheap once the first has given you the vocabulary:

| First | Then | Why |
|---|---|---|
| the BELFEM paper (Messe et al. 2023) | the formulation paper it builds on | what we chose, then why it works |
| a book chapter | the paper that specializes it | general mechanism, then the HTS-specific form |
| a symptom (Bathe Ch. 8) | the BELFEM paper's §2.7 | robustness theory, then our knob |

## The stop rule

Stop when you can name the source, the section, and the sentence that answers the question — or
when you can state precisely which source and section *would* answer it and why you cannot reach
it. Both are reportable results. What is not reportable is a confident answer assembled from
recollection plus code reading, presented as though the library backed it.

When the library is absent, that second form is the only honest output:

> The controller's integral step matches the form in Badel et al. 2019 §4.3 on code reading alone;
> `literature/` is not present in this checkout, so the coefficient could not be checked against
> the paper. DOI:10.1088/1361-6668/ab181f. Confidence: medium.

## Cross-checking a claim you already have

Going the other way — a devlog, a comment, or an AI reviewer asserts "X, per Y et al. Z" — is the
same chain run backwards, and it is worth running. Open the nav file, confirm the work covers the
topic at all, then find the section. Citations drift: a section number that was right for a
preprint is wrong for the published version, and a reviewer that cannot open the library will
produce a plausible section number rather than admit it. Treat an uncheckable cite as unverified,
per `doc/ai_collaboration_protocol.md` §11.
