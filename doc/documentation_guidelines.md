# Documentation Guidelines {#doc_documentation_guidelines}

**Purpose:** Standards for organizing and naming markdown documentation in BELFEM.

---

## Directory Structure

```
belfem/
├── README.md                    # Project overview
├── CLAUDE.md                    # Claude Code instructions
│
├── todo/                        # Task planning and implementation plans
│   ├── README.md               # Index of tasks
│   └── *.md                    # Task files
│
├── devlog/                      # Session logs of changes made
│   ├── README.md               # Index of devlog entries
│   └── dlYYYYMMDD_topic.md     # Individual session logs
│
├── doc/                         # General project documentation
│   ├── README.md               # Index of documentation
│   └── *.md                    # Documentation files
│
└── src/[module]/doc/           # Module-specific documentation
    └── *.md                    # Co-located with source code
```

---

## File Categories

### **Task Planning** → `./todo/`

**Purpose:** Implementation plans, refactoring tasks, feature designs, work-in-progress

**When to use:**
- Planning a new feature
- Designing a refactoring approach
- Documenting implementation steps
- Tracking technical debt items

**Examples:**
- `thin_shell_volume_conductor_coupling.md` - Implementation plan
- `solver_performance_optimization.md` - Optimization task
- `mesh_format_migration.md` - Migration plan

### **Session Devlogs** → `./devlog/`

**Purpose:** Records of what was changed and why during AI-assisted development sessions

**When to use:**
- End of any debugging or investigation session
- After completing a refactoring task
- When resolving a significant exchange thread

**Naming:** `dlYYYYMMDD_topic.md`

**Examples:**
- `dl20260318_ghost_thinshell.md` - DG edge decoupling session
- `dl20260320_vertex_capacity_refactoring.md` - Refactoring session

**Distinction from `./todo/`:** Devlogs are backward-looking (what changed). Task files are forward-looking (what to do).

### **General Documentation** → `./doc/`

**Purpose:** Architecture explanations, analysis, comparisons, guides

**When to use:**
- Explaining how existing code works
- Performance analysis and benchmarks
- Design decision comparisons
- Cross-cutting architecture documentation
- Build system guides

**Examples:**
- `fem_dof_manager_and_hphi_formulation.md` - Architecture explanation
- `solver_performance_comparison.md` - Benchmark analysis
- `build_system_guide.md` - Developer guide

### **Module Documentation** → `./src/[module]/doc/`

**Purpose:** Module-specific technical documentation

**When to use:**
- Documenting algorithms specific to one module
- Mathematical theory for a component
- Module-specific API guides

**Examples:**
- `src/homology/doc/cohomology_algorithms.md`
- `src/fem/iwg/doc/iwg_usage_guide.md`
- `src/fem/maxwell/doc/maxwell_weak_forms.md`

### **Root Directory**

**Keep minimal - only essential files:**
- `README.md` - Project entry point
- `CLAUDE.md` - AI assistant instructions
- `LICENSE` - Legal information
- Configuration files (`.gitignore`, etc.)

---

## Naming Convention

### **Standard Format**

Use **lowercase_with_underscores.md** for all documentation files:

```
lowercase_descriptive_name.md
```

**Rules:**
- All lowercase letters
- Words separated by underscores `_`
- No spaces (filesystem compatibility)
- Descriptive but concise
- Acronyms in lowercase: `fem`, `iwg`, `hdf5`

### **Examples**

✅ **Good:**
- `fem_dof_manager_architecture.md`
- `maxwell_thin_shell_formulation.md`
- `build_system_guide.md`
- `solver_performance_analysis.md`

❌ **Bad:**
- `FEM_DOF_Manager.md` (mixed case)
- `MaxwellThinShell.md` (camel case)
- `Build System Guide.md` (spaces)
- `DOF-Manager.md` (hyphens, all caps)

### **Exceptions**

**Root-level convention files use CAPS:**
- `README.md` - Universal convention
- `CLAUDE.md` - Special configuration
- `LICENSE` - Legal convention

---

## Creating New Documentation

### Decision Tree

```
Is this about future work, a session record, or existing code?
│
├─► Future work (task planning)
│   └─► Create in: ./todo/task_name.md
│
├─► Session record (what changed)
│   └─► Create in: ./devlog/dlYYYYMMDD_topic.md
│
└─► Existing code (documentation)
    │
    ├─► Specific to one module?
    │   └─► Create in: ./src/[module]/doc/topic.md
    │
    └─► Cross-cutting or general?
        └─► Create in: ./doc/topic.md
```

### Workflow

1. **Choose directory** using decision tree above
2. **Name file** using lowercase_with_underscores
3. **Update README.md** in that directory
4. **Include metadata** at top:
   ```markdown
   # Title

   **Date:** YYYY-MM-DD
   **Purpose:** Brief description
   ```

---

## Publishing to the Doxygen Site

Everything under `doc/` and `src/<module>/doc/` is rendered into the generated
documentation by `make doc`, so these files have two audiences: people reading the
repository on the web, and people browsing the generated site.

Two pieces of Doxygen wiring are **generated, never hand-written**:

- the `{#anchor}` on each document's top heading, which fixes its page URL
- `doc/doxygen_nav.dox`, which holds the whole page hierarchy as `@subpage` directives

Run `scripts/update_doc_index.py` after adding, renaming or removing a document. It is
idempotent, and `make doc` runs it automatically. `--check` reports pending changes and
exits non-zero, for a hook or CI job.

Consequences for authors:

- **Give every module doc directory a `README.md`.** It becomes the module's index page;
  a directory without one cannot be nested and the script will say so.
- **Start every document with a single `# ` heading.** The anchor and the page title both
  come from it.
- **Link to files, not directories.** `](../../mesh/doc/)` renders as a dead link in the
  generated site; `](../../mesh/doc/README.md)` works on both surfaces.
- **Do not write `@subpage` or `@ref` into a module document.** Doxygen directives show
  up as literal text when the file is viewed on the web. `doc/mainpage.md` is the sole
  exception, because it exists only for Doxygen.
- **Avoid an inline code span inside bold text.** Doxygen 1.9.1 mis-parses it and leaves
  the emphasis open for the rest of the file — one such line was holding a `<strong>`
  open across 370 lines of `circuit_usage_guide.md`. Close the bold before the code span:

  ```
  good:  **the last node index** (`n - 1`) **is the ground node**
  bad:   **the last node index (`n - 1`) is the ground node**
  ```

Warnings go to `<build>/doc/doxygen_warnings.log`, and `make doc` reports the count
against a known baseline. Anything above the baseline is a regression worth fixing.

## Index Files (README.md)

Each documentation directory should have a `README.md` index file.

### **`./todo/README.md`**

```markdown
# Task Planning

## Active Tasks
- [task_name.md](task_name.md) - Brief description

## Completed Tasks
- [old_task.md](old_task.md) - Brief description

## Superseded
- [obsolete_design.md](obsolete_design.md) - Reason superseded
```

### **`./devlog/README.md`**

```markdown
# Devlog

## Entries
- [dlYYYYMMDD_topic.md](dlYYYYMMDD_topic.md) - Brief description
```

### **`./doc/README.md`**

```markdown
# BELFEM Documentation

## Architecture
- [architecture_doc.md](architecture_doc.md) - Description

## Analysis
- [performance_analysis.md](performance_analysis.md) - Description

## Guides
- [build_guide.md](build_guide.md) - Description

## Module-Specific
- [Module Name](../src/module/doc/)
```

### **`./src/[module]/doc/README.md`**

Required for any module whose docs should appear in the generated navigation: `update_doc_index.py` warns and drops a README-less directory from nesting rather than failing the run.

An index: what the module is, what documents it has, and a quick-reference table of its key
classes and pitfalls. It does **not** need a references section — see "Where citations live"
below. As of 2026-08-31 exactly one of the 27 module READMEs carries a DOI (`homology`), and
that is the expected state, not a gap: there the algorithm-to-paper mapping is the module's
subject matter, so the list earns its place — and it now points at the project citation list
as authoritative, so the two cannot drift apart silently.

---

## Best Practices

### **Naming**

- Be descriptive: `fem_dof_hanging_implementation.md` not `dof_stuff.md`
- Be concise: `solver_guide.md` not `complete_guide_to_solver_architecture.md`
- Use prefixes for related docs:
  - `fem_dof_manager.md`
  - `fem_dof_hanging.md`
  - `fem_element_types.md`

### **Organization**

- **One topic per file** - don't create mega-documents
- **Link between documents** - use relative links
- **Keep README.md updated** - it's the entry point
- **Date your documents** - include creation date in header
- **Mark obsolete docs** - add "SUPERSEDED" or "OBSOLETE" header

### **Content**

- **Code references** - include file paths and line numbers
  ```markdown
  See `cl_FEM_DofManager.cpp:123-145` for implementation
  ```
- **Diagrams** - ASCII art or saved images in `doc/images/`
- **Examples** - include code snippets
- **Links to literature** - see the rule below; link to the citation list, do not
  copy citations into module docs

### **Where citations live (one place, on purpose)**

`doc/literature_references.md` is the **single source of truth** for full citations and DOIs.
A module document that rests on published work should **name the work and link to that file**,
not carry its own DOI:

```markdown
Cuts follow the cohomology basis of Pellikka et al. — see
[literature_references.md](../../../doc/literature_references.md) for the full citation.
```

**Why not a reference section in every module README.** Duplicated citations rot
independently: the same paper ends up with three DOIs, two spellings and one broken link, and
nothing tells you which is current. One list, many links, is the arrangement that survives.

**This is not a requirement to invent references.** Most modules implement no published
method — `core`, `containers`, `comm`, `io` and friends have nothing to cite, and a
"References" heading with nothing under it is worse than its absence. Cite where a module
genuinely implements someone's method; stay silent elsewhere.

*Recorded 2026-08-11: an earlier version of this list read simply "reference papers in
`literature/`", which was read as an expectation on every module README. Measured against the
tree, 22 of the 23 module READMEs carried no reference at all — so the convention existed only
on paper. Rather than manufacture 22 reference sections, the convention is now stated as what
is actually useful and actually followed.*

### **Maintenance**

- **Update when code changes** - keep docs synchronized
- **Archive obsolete docs** - don't delete, mark as superseded
- **Review periodically** - remove truly obsolete content

---

## Module Documentation Structure

For modules with substantial documentation, use subdirectories:

```
src/module/doc/
├── README.md                    # Module documentation index
├── architecture.md             # Overall design
├── algorithm_name.md           # Specific algorithms
├── theory_and_implementation.md # Mathematical background
└── images/                     # Diagrams and figures
    └── *.png
```

---

## Examples

### Creating a New Task

```bash
# 1. Create file in todo/
$ cat > todo/maxwell_adaptive_timestepping.md
# Maxwell Adaptive Timestepping

**Date:** 2026-01-16
**Purpose:** Implement adaptive time-stepping for Maxwell solver

## Background
...

# 2. Update todo/README.md
$ vim todo/README.md
# Add: - [maxwell_adaptive_timestepping.md](...)
```

### Creating Module Documentation

```bash
# 1. Create module doc directory
$ mkdir -p src/solver/doc

# 2. Create documentation
$ cat > src/solver/doc/newton_raphson_implementation.md
# Newton-Raphson Implementation

**Date:** 2026-01-16
**Module:** solver

## Algorithm
...

# 3. Update doc/README.md to link to module
$ vim doc/README.md
# Add module reference
```

---

## Migration from Old Conventions

If you find files that don't follow these conventions:

1. **Check git history** - is it still relevant?
2. **Rename file** to lowercase_with_underscores
3. **Move to appropriate directory** (todo/ or doc/)
4. **Update any references** to the old filename
5. **Update README.md** in destination directory
6. **Commit with message**: `docs: organize [filename] per documentation guidelines`

---

## Rationale

### Why lowercase_with_underscores?

- **Readability:** Easier to read than ALLCAPS or CamelCase
- **Typing:** No shift key needed, faster in terminal
- **Compatibility:** Works on case-insensitive filesystems
- **Convention:** Matches Unix/Linux documentation standards
- **URLs:** Cleaner in GitHub Pages and web links

### Why Separate todo/ and doc/?

- **Clear intent:** "What we're building" vs "how it works"
- **Discoverability:** Easy to find active work items
- **Maintenance:** Mark tasks complete without deleting
- **Review:** Easy to see what documentation needs updating

### Why Module-Specific doc/ Folders?

- **Proximity:** Documentation near the code it describes
- **Modularity:** Self-contained modules
- **Ownership:** Module maintainers own their docs
- **Clarity:** No ambiguity about scope

---

## Quick Reference

```bash
# Task planning
./todo/feature_name.md

# Session devlogs
./devlog/dlYYYYMMDD_topic.md

# General documentation
./doc/topic_name.md

# Module documentation
./src/module/doc/algorithm_name.md

# Always lowercase_with_underscores.md
# Update corresponding README.md
# Include date and purpose in header
```

---

## See Also

- Project README: `../../README.md`
- Claude Instructions: `../../CLAUDE.md`
- Literature Organization: `../../literature/README.md`
