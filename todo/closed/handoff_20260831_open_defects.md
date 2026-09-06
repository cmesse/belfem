# Handoff: Defects Found During the Documentation Campaign, Still Open

**Date:** 2026-08-31
**Purpose:** The findings the 2026-08-31 documentation sweep and repair campaign turned up that are
**not** documentation problems and were therefore not fixed by it — each with what was verified, how,
and what decision it needs
**Module:** cross-cutting (`math/graph`, `fem/maxwell`, `homology/doc`)

---

## Why this file exists

The campaign's scope was documentation and comments only. Everything it found that was a *documentation*
defect is fixed (`devlog/dl20260831_doc_repair_campaign.md`). Four things it found are not, and one of
them is an error of mine that a ticked checkbox is currently hiding.

One code defect **was** fixed mid-campaign on Christian's ruling — `Vertex::mOwner` initialising to
`gNoID` instead of `gNoOwner` — and it is not repeated here except for the residue in F3, which is a
consequence of the fix rather than of the bug.

**Evidence posture:** every claim below was re-checked at session close with `command grep` (the shell
`grep` here wraps `ugrep --ignore-files` and silently skips `archive/`, `nonfree/`, `literature/` and
`tmp/` — see F5). Nothing here was compiled or run. **Reviewed, not verified.**

---

## F1 — `build_pargraph_adjacency` indexes an array by an unassigned owner sentinel

**Severity: LOW as it stands, because the function is unreachable (F2). Would be HIGH if wired.**

`src/math/graph/graphtools.hpp:136` sizes the counter to the rank count:

```cpp
Vector< T > tCount( tCommSize, 0 );
```

and `:142-146` guards exactly one unassigned marker before using the owner as the index:

```cpp
// bugfix for non-assigned dofs e.g. from circuit model
if ( tVertex->owner() == tCommSize )   // catches the comm_size() marker only
{
    tVertex->set_owner( 0 );
}
++tCount( tVertex->owner() );          // <-- unguarded for either sentinel
```

The guard catches the deliberate `owner() == comm_size()` marker set at `cl_FEM_Kernel.cpp:195,203`.
It does **not** catch `gNoOwner`. A vertex still holding the unassigned sentinel reaches `:146` and
indexes `tCount` at **2147483647** against a length of `comm_size()`.

`tCount` is a `Vector<T>`, so a debug build asserts. **A release build performs an out-of-bounds
increment — a write, not a read.**

Note the interaction with the `mOwner` fix that landed this session: before it, the stale value was
`-1` → `2^64-1` under `Vector`'s `size_t` subscript, which always faults. After it, the sentinel is
`2^31-1`, which is a *plausible* offset and therefore **less** certain to fault. The fix is right for
the reasons in F3, but it moves this particular site from loud to quiet.

**RESOLVED 2026-08-31, Christian: F2 is to be wired, so this must be fixed.** It becomes **D1** of
`todo/parmetis_ptscotch_wiring.md` and is a precondition for the first caller, not follow-up work —
wiring is exactly what makes the out-of-bounds write reachable.

## F2 — `parmetis_nd`, `ptscotch_nd` and `build_pargraph_adjacency` are dead code

**Severity: none today. It is a maintenance and honesty question, not a bug.**

Verified with `command grep` across `src/`, `tests/`, `archive/` and `nonfree/`: both entry points are
**defined and declared, and called from nowhere.**

| Symbol | Definition | Declaration | Call sites |
|---|---|---|---|
| `parmetis_nd` | `fn_Graph_ParMETIS.cpp:31` | `fn_Graph_ParMETIS.hpp:23` | **none** |
| `ptscotch_nd` | `fn_Graph_PTSCOTCH.cpp:28` | `fn_Graph_PTSCOTCH.hpp:23` | **none** |
| `build_pargraph_adjacency` | `graphtools.hpp:126` | — | only the two above |

(The `ptscotch_nd` hits under `tmp/STRUMPACK/` are a vendored third-party tree and a different symbol,
`sep_tree_from_ptscotch_nd_tree`.)

**Decision needed — three ways, and they are not equal:**

1. **Wire it.** Then F1 must be fixed first, and the path needs a test; nothing in `make check`
   currently exercises it.
2. **Move to `archive/`.** Matches what was done with `fn_Graph_tarjan`, and keeps the prototype
   findable — which mattered this session, see F5.
3. **Delete.** Cheapest, and loses the work.

Recommendation was **(2)**. **RULING 2026-08-31, Christian: (1) wire it** — parallel nested
dissection is wanted. Planned in `todo/parmetis_ptscotch_wiring.md`, D1 first, then a test that
actually runs multi-rank, then a benchmark. Scheduled 2026-09-01.

## F3 — Residue of the `mOwner` fix: failure got quieter, and that should be on the record

**Severity: informational. Not a defect. Do not "fix" it back.**

The fix (`cl_Graph_Vertex.hpp:42`, `gNoID` → `gNoOwner`) is correct, and the strongest evidence is not
the naming but the algorithms: the `std::min` ownership sweeps at `cl_Mesh_Partitioner.cpp:253-260` and
`cl_FEM_Kernel.cpp:427` need the sentinel to behave as **+∞**. Under `gNoID`-narrowed-to-`-1` it was an
absorbing element, and every facet would have ended unassigned. Several sites are **repaired** by it —
`cl_Mesh.cpp:1669,1684` and `cl_FEM_Kernel.cpp:379,918,1103` all did signed comparisons that `-1`
slipped through.

Two things stated during the campaign were **overstated and are corrected here**, both downward:

- **Blast radius.** `mesh::Basis::Basis()` sets owner to `0` at `cl_Mesh_Basis.cpp:25`, so **no mesh
  entity ever observed the default.** Only direct `graph::Vertex` users could.
- **Failure mode.** The old value converted to `2^64-1` under `Vector`'s `size_t` subscript and
  therefore *always* faulted. The new one gives `2^31-1`. This is a small regression in failure
  **loudness**, which is what makes F1 worth fixing rather than leaving.

**Still owed:** `make check` has not run since the fix and its test companion landed.

## F4 — ~~A doc residue I ticked as done, and it is not done~~ — **CLOSED 2026-08-31**

**Severity: MEDIUM. It teaches a pipeline that does not exist.**

> **CLOSED later the same day**, in a follow-on session under explicit authorization to edit
> homology module documentation (source untouched). All four sites in `homology_usage_guide.md`
> now carry the correction: the glossary rows, the Internal Workflow step, the "Manifold Cleanup"
> section (which gained the same design-note banner the sibling documents use), and the
> troubleshooting entry at `:2205-2208`. **That fourth site was not in the list below** — it told
> readers to tune the nonexistent function's thresholds, which is the most actionable form the
> error can take. Nothing was deleted; the design text stands, marked, for Gregory's ruling.
>
> Two corrections to what this section asserted:
> - the guide's residue was **four** sites, not two
> - "a three-phase Tarjan ... cleanup that no shipped code applies" is right about the pipeline,
>   but Tarjan itself is not fictional — `archive/graph/fn_Graph_tarjan.{hpp,cpp}` implements it,
>   outside the build. The glossary now says so rather than implying it ships.
>
> **Found in the same pass, and also fixed:** the advertised debug outputs `thick_cut_*.vtk`,
> `thin_cut_*.vtk` and `manifold_*.vtk` are written by no code; `save_debug_meshes()` emits
> `cut_<index>.vtk` and belongs to `CutProcessor`, not `CutFactory`. A sweep of all 44 code
> identifiers asserted across the six homology documents turned up no further missing symbols.

D6 in `todo/doc_currentness_fixes.md` named three files describing a manifold-cleanup workflow built on
symbols that do not exist. Two were fixed: `cohomology_theory_and_implementation.md:47-48` now carries
the correction banner, and `homology_usage_guide.md:839` now names the live path
(`clean_spfa()` → `remove_cut_pockets()`). **`src/homology/doc/homology_usage_guide.md:1329-1350` and
`:335` were not**, and I ticked the box anyway.

Those lines still present, with no banner, as the live pipeline:

- `manifold_filter_3d()` — **zero hits tree-wide**, `archive/` included
- `check_surface()` — **zero hits tree-wide**
- a three-phase Tarjan / region-growing / BFS cleanup with specific metric thresholds
  (cycle density > 0.3, compactness < 2.0, size < 20 faces) that no shipped code applies
- "**Implementation:** `manifold_filter_3d()` in CutProcessor" — `CutProcessor` is real (55 hits in
  `src/`), the method is not, which is the most misleading form this can take

The live path is `Cohomology::clean_spfa()` → `remove_cut_pockets()` (`cl_Cohomology.hpp:135,145`).

**This is documentation, so it is in the campaign's scope and I simply missed it.** It is listed here
rather than fixed because `src/homology/doc/` belongs with the module owner and is already routed
to Gregory Giard through the internal handoff. The
`§7.1` ban covers the cohomology **source**, not its docs, so whoever picks this up may edit these
files — but the content decision is Gregory's.

## F5 — The `grep` hazard: blast radius measured, and it is small

**Severity: informational, and the reason it is here is that the measurement is the useful part.**

Shell `grep` in this environment is a function wrapping `ugrep --ignore-files`, which honours
`.gitignore` — so every repository-root search silently excluded `archive/`, `nonfree/`, `literature/`
and `tmp/`. It caused one wrong claim this session: that no Tarjan implementation existed, when
`archive/graph/fn_Graph_tarjan.{hpp,cpp}` does. That claim reached a module document *and* the handoff
to the module owner before Codex caught it. Both are corrected, the archive file is now cited by path,
and the mechanism is saved as a memory.

**What was then measured, because one instance says nothing about the class:** every absolute-negative
claim in user-facing documentation was enumerated and re-checked with `command grep`. There are five.
The two load-bearing ones both **hold**:

| Claim | Location | Re-checked |
|---|---|---|
| `BELFEM_PHDF5` is defined nowhere | `src/io/doc/io_usage_guide.md:1006` | holds — no hit in `src/`, `config/`, `tests/`, `archive/`, `CMakeLists.txt` |
| `add_abstract_dof` — no such symbol | the internal cohomology handoff for Gregory Giard (not in the published tree) | holds — no hit in `src/`, `tests/`, `archive/`, `nonfree/`; live call is `extract_abstract_dofs_from_mesh()` (`cl_FEM_DofManager.hpp:205`) |

So the exposure in shipped documentation was **one claim, and it is fixed**. No audit sweep is needed.
Devlogs are a different matter and are deliberately **not** covered: they are dated records kept as
written, and a false negative in one stays as history.

---

## Open question, not a defect

**Q1 — Is the absence of AC-loss postprocessing deliberate?**

`src/fem/maxwell/doc/README.md:20` used to list **losses** among the postprocessed fields. It was
removed this session because a case-insensitive search for `loss` across every `.cpp`/`.hpp` under
`src/fem/maxwell/` returns **zero hits**, while B/H/J/JJC are all real
(`cl_MaxwellPostprocessor.cpp:113,117,146,150`).

For an HTS code this is the headline engineering quantity, so the correction removed an advertised
feature rather than a typo. **The question for Christian is which of these is true:** loss is intended
to be computed downstream from the exported `J` and `E` fields and the documentation was simply
aspirational — in which case `doc/` should say so, because a reader will look for it; or it is a real
gap worth a `todo/`.

**RULING 2026-08-31, Christian: a real gap, not aspirational documentation.** The heat imposed on the
structure is already computed — `rho * norm_j * norm_j` per integration point, already weighted by
`w(k)` and `dV(k)` (`mt_thermal_h.cpp:49,124`) — so summing it is roughly a half-hour task. Planned in
`todo/ac_loss_postprocessing.md`, which logs the three things that could make the sum quietly wrong
(resistivity source and isothermal runs, the `[0, 1e10]` clamp, thin shells). Scheduled 2026-09-01.

---

## Summary

| ID | What | Severity | Needs |
|---|---|---|---|
| F1 | Owner sentinel indexes `tCount` out of bounds in release | LOW now → **HIGH, F2 is being wired** | **D1** of `parmetis_ptscotch_wiring.md`, fix before first caller |
| F2 | ParMETIS/PT-Scotch nested dissection is dead code | none today | **RULED: wire it** → `parmetis_ptscotch_wiring.md` |
| F3 | `mOwner` fix made failure quieter; two claims corrected downward | informational | `make check` |
| F4 | ~~Homology guide still teaches a nonexistent cleanup pipeline~~ | ~~MEDIUM~~ **CLOSED** | fixed 2026-08-31 (4 sites, +debug-vtk defect); Gregory's content call still open |
| F5 | `grep` false negatives — blast radius measured at one, fixed | informational | nothing |
| Q1 | No AC-loss postprocessing | **RULED: real gap** | → `ac_loss_postprocessing.md` |

**Nothing in this file was compiled or run. Reviewed, not verified.**
