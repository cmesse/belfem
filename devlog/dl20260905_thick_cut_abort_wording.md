# Thick-cut abort: honest wording, and why "every thick cut is a thin cut" does not apply to a fixed mesh

**Date:** 2026-09-05
**Purpose:** Reword the `clean_spfa` abort in the cohomology core (message text only, authorized by
Christian as main developer and Gregory's advisor), and answer the question of which theorem the
thin-cut conversion seems to contradict.
**Module:** homology

## 1. The edit

`src/homology/cl_Cohomology.cpp:533-543`, the `BELFEM_ERROR` fired when the SPFA difference-constraint
system is infeasible. The old text said "no thin cut exists on this mesh ... the mesh is too coarse".
The new text states the actual limit first — BELFEM only supports thin cuts (a surface of element
faces with a single potential jump); thick cuts as basis functions are not implemented — then the
obstruction (the certificate loop winds around the conductor more often than it has edges) and the
remedy (refine along the loop). Format arguments (`%u`, `%lu`, `%s`) and their order unchanged; the
"betwen" typo fixed.

The cohomology core is closed to AI edits (protocol §7.1). This one edit was explicitly authorized
by Christian for the message text; no functionality was touched.

## 1.1 Codex language sweep and the `CutProcessor` guards

Codex gpt-5.6-luna medium (`tmp/ai_exchange/thick_cut_abort_wording.md`) swept the message plus the
two `collect_facets` guards of 2026-09-02 (`cl_CutProcessor.cpp`, empty and collapsed thin cut).
Taken: its rewrite of the core message ("No representative of this generator has only -1, 0, and 1
as coefficients"), and the up-front sentence "The solver supports only thin cuts; thick cuts are not
implemented yet" on both 3D guards. Rejected: its closing "Use a mesh or generator that produces a
tight one-sheet cut" — the user cannot choose the generator and Christian's 2026-09-02 ruling records
no remedy; replaced by "No user-side remedy exists yet". Its flag that the 2D branch has its own
empty-cut guard (`cl_CutProcessor.cpp:384`, conjugated edges) was right; that message got the same
treatment. Placeholders checked programmatically (same tokens, same order) on all four edits; both
files pass `-fsyntax-only` with the `build/` tree's compiler, defines and flags. Not built, not run.

## 2. What we missed (answer to Christian's question)

The theorem is Gross & Kotiuga, Mathematical Appendix MA-I: `H¹(Ω;ℤ) ≅ [Ω, S¹]`, every class is the
pullback of the generator under a map to the circle, and the preimage of a regular value is an
embedded orientable surface — the cut. With Poincaré–Lefschetz, the thick cut (a surface in the dual
complex) and the thin cut (a surface of primal faces) represent the same class. That is the
"every thick cut is a thin cut" statement, and it holds **in the manifold**.

Its mesh realization is Gross & Kotiuga Ch. 6, eq. 6-31 and the paragraph after it (p. 171): the
level set of the harmonic map is perturbed onto element faces, "unambiguous if the mesh is fine
enough to ensure that, over an element, θ does not go more than one third of the way around the
circle". That is the hypothesis a fixed mesh can violate. Its discrete form is already in
`src/homology/doc/thin_cut_nonunit_rectification.md` §2–3: for every closed edge loop `z`,
`|⟨c, z⟩| ≤ length(z)` is necessary for a unit representative, the pairing is a class invariant,
and the SPFA negative cycle is a loop that breaks it. Refinement always restores feasibility since
the pairing is invariant while the loop length grows. Gross & Kotiuga even recommend computing the
cut on the coarsest mesh and refining it with the mesh — the topology is mesh-independent, the
facet representative is not.

Second gap (confidence medium): Kotiuga's representative is tight for free — a level set of a
smooth function is locally the boundary of a sublevel set. BELFEM's representative comes out of the
Smith normal form and recombination as *some* cocycle in the class; unit coefficients are necessary
but not sufficient for the facet push. That is the CORC sheath of 2026-09-02. Pellikka 2013's thick
cut (the cocycle used directly as an edge basis function, eq. 4.5) sidesteps both; it is option 2
of `todo/cut_representative_options.md` and what "not implemented yet" refers to.

Reviewed, not verified: no build, no run.
