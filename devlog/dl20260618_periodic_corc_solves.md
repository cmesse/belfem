# corc periodic thin-cut reproducer SOLVES — Step 9 (Option 2) validated end to end

**Date:** 2026-06-18
**Purpose:** Record the milestone: the z-periodic CORC reproducer now runs the full pipeline
(cuts → cohomology → thin-cut duplication/relink → periodic rebuild → assembly → nonlinear
solve) and produces a physical HTS current distribution. This closes the functional φ-continuity
blocker that only appeared when cohomology thin cuts meet periodic boundary conditions.
**Module:** src/homology, src/mesh, src/fem/maxwell

---

## The moment

corc computes. The output (`cmake-build-debug/corc.png`) shows the CORC cable's helically-wound
HTS tapes colored by **J/Jc magnitude** (≈0.44–0.65), with the characteristic diamond/stripe
current-sharing pattern along the tapes — a physically sensible distribution. The run cleared
`create_facet_map` ("Node N not flagged"), `set_entity_dependencies`' edge wiring, assembly, and
reached the magnetic solve.

## What finally unlocked it (Step 9 — Option 2)

The last blocker was the Maxwell-pass periodic rebuild rejecting 178 seam cut-duplicates. The
chain of insight that resolved it:

1. **Diagnosis (measured, not guessed):** the 178 offending seam nodes are thin-cut duplicates
   `dup = original + cutDOF` — sourced from abstract cohomology cut DOF(s) + a *periodic*
   original. Abstract-aware classification: **178 induced / 0 genuine-jump / 4 abstract cut
   DOFs**. They were absent from the restored periodic node lists because `restore_node_pairs`
   replays only the *paired* backup and these one-bit dups carry no pair.
2. **The principle (Christian):** periodicity and hanging are orthogonal. A hanging dup's
   periodicity is *induced* through its sources — `slave = master + cutDOF` follows once the
   original's periodicity is applied. So the dup never needs to be a matched periodic node.
3. **Option 2 — key by original identity.** `create_facet_map`/`create_edge_map`/`match_edges`
   key periodic facets/edges by `node->original()` (which returns the geometric original for a
   cut dup and `this` for thin-shell dups / non-dups — automatically scoped). The value-
   periodicity rides the hanging/source chain; the DOF manager adopts it.
4. **Edge-dependency pairing fix.** `set_entity_dependencies` had been pairing master/slave
   edges *by list position*, but those lists are collected independently per side. Switched to
   the authoritative `A->periodic()` pointer (set by `match_edges`) + an original-aware
   orientation assert + an always-active `BELFEM_ERROR` null guard. Added release-active
   degenerate/collision guards to `create_edge_map`.

Steps 1–8 (periodic sideset routing, deterministic node-pair backup/restore, the Step-4
period-direction-jump diagnosis, Step 5a/5f) were the necessary groundwork; Step 9 was the
final key turn.

## Validated vs remaining

- **Validated:** the run reaches and completes the magnetic solve on a z-periodic CORC. The
  functional blocker is gone.
- **Remaining (Step 7 acceptance — correctness, now runnable):** φ continuity across the seam
  *up to the intended cohomology jump* (7g/3j); cut reaches the seam in debug VTK (7d/1g); full
  time-stepping. Then **Step 8 cleanup**: remove the `#DIAG 5f–5k` probes, `check_health`, the
  parked `match_nodes_and_edges` block, and fix the stale `periodicity.md:84` positional-pairing
  doc.

## Acknowledgements — a genuinely joint effort

This was solved collaboratively, with Christian (Prof. Messe) driving the physics and the key
conceptual calls, and three AI collaborators contributing distinct strengths:

- **Claude Opus** — broad diagnosis, probe instrumentation (#DIAG 4b–5k), synthesis, and the
  implementation/audit loop.
- **Codex** — precision source audits and the catch that "M always exists" was overstated, the
  edge-list positional-pairing flaw, and the release-safety push on the guards.
- **Grok** — independent third-voice verification (high-confidence confirmations + citation
  corrections), and corroboration of the inducible-dup classification.
- **Fable** — last week's Step-5 emission/rebuild drafts (`tmp/periodic/`), several of which
  were salvaged into the plan (the jump-side filter signature, the allocate-vs-source-fallback
  framing).

The repeated pattern — measure before fixing, cross-check non-trivial claims with at least one
independent auditor, and reverse conclusions when a probe contradicted them — is what kept the
diagnosis honest through several wrong turns (the cohomology "regression", the index-roundtrip
hypothesis, "M always exists"). Each was caught and corrected by the cross-check discipline.
