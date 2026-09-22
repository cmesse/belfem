# Notation, Naming, and Symbol Annotations

**Purpose:** How a symbol from a paper becomes a variable name and a comment in BELFEM code, and
which document governs which half.

---

## Three naming systems, already in force

| Where | System | Governed by |
|---|---|---|
| Framework code | `a` argument, `t` temporary, `m` member, `g` global, then camelCase | CLAUDE.md, "Naming Conventions" |
| `src/math`, `src/physics` kernels | plain mathematical names — `a, b, c`, `p, q, r`, `D`, `lambda`; `m` and `g` still apply | CLAUDE.md, "Math/physics module exemption" |
| `physics/gastables`, `gasmodels` | the fixed thermophysical symbol table — `T`, `p`, `rho`, `v`, `h`, `s`, `u`, `A`, `alpha` | CLAUDE.md, "Canonical thermophysical symbols" |

The literature adds a fourth layer, and it is the one that is easy to miss: **which symbol the
source used**. That is what a symbol annotation records.

## What the notation guide is for

`literature/doc/notation_guide.md` defines BELFEM's unified operator notation — **N** node
interpolation, **E** edge interpolation, **B** gradient, **C** curl (B is to grad as C is to curl;
N is to nodes as E is to edges) — plus field-variable conventions, per-book translation tables, and
a warnings section for the books whose notation collides with it.

Its scope is **equations in prose**: module documentation, usage guides, devlog derivations, and
anything presented to a reader as "the weak form we implement". Translate into the unified notation
there, and say which source you translated from.

Its scope is **not** variable names in a kernel. Code that implements a published equation keeps
the published symbols, because that correspondence is the point — `doc/coding_philosophy.md`
protects it, and `doc/commenting_guidelines.md` §11 records why splitting a stencil into named
helpers destroys it. The bridge between the two is the annotation.

## The symbol annotation

`doc/commenting_guidelines.md` §2 lists **symbol annotation** among the six things a comment may
say: *the mathematical symbol a statement or a variable stands for, inside a kernel written in
domain notation*. §4 makes it the one comment form allowed at the end of a line, on a declaration
or an enumerator.

Two existing forms, both already in the tree:

```cpp
real tRho ;                          // [Ω·m]
aPoints.set_row( 1, tPoints.row( 0 ));   // eta
```

When the kernel follows a specific source, the annotation carries the source too, so the next
reader knows whose symbol table is in force:

```cpp
// h-φ interface, Alves et al. 2022b Eq. 19 notation
Vector< real > & mHt ;               // h_t, tangential field
Vector< real > & mB ;                // b, flux density  — NOT the gradient operator
```

The block citation above the kernel states the source once; the trailing annotations name the
symbols. One idea per comment, per §4.

## Collisions: annotate, never unify

The same letter carries different meanings in different modules, and every one of these is correct
in its own context:

| Symbol | One meaning | Another meaning |
|---|---|---|
| `B` | flux density **b** in the Maxwell calculator (`mB = mu * h`) | the gradient operator **B** in the notation guide |
| `u` | velocity in the gas flow routines | internal energy in the caloric accessors |
| `J` | current density | the Jacobian |
| `N` | node interpolation operator | the layer count in a thin-shell model |
| `n` | outward normal | the power-law exponent |

CLAUDE.md already rules on `u`: both meanings are standard, the symbol is context-bound, **do not
unify them**. The same holds for the rest. The fix for a collision is an annotation at the
declaration, never a rename that breaks the correspondence with the paper.

## Per-source notation warnings

Every book nav file has a `## Notation Warnings` section with a translation table — Bathe writes
**H** for the interpolation matrix where the unified system writes **N**, and puts the
configuration on a left superscript; Monk and Boffi use the edge-element conventions; Hughes and
Zienkiewicz each differ again. Read that section before transcribing an equation out of a book, and
before arguing that a code symbol contradicts a source.

## What the sweep will catch

`doc/commenting_guidelines.md` §3 lists **mumbling** — "a symbol above a statement it does not
describe" — and notes that a legitimate symbol annotation and mumbling look alike, so the sweep has
to tell them apart. Keep yours on the right side of that line:

- [ ] The annotation sits on the statement or declaration it names, not near it
- [ ] The symbol matches the source cited above the block, or the annotation says which source
- [ ] A unit annotation carries the unit in brackets: `// [Ω·m]`, `// [A/m²]`
- [ ] No date, no name, no task ID travels with the citation (§8)
- [ ] Prose and documentation use the unified notation and name the source they translated from
- [ ] The comment is edited in the same diff as the code it describes (§12)
