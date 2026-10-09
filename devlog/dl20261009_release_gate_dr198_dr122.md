# Devlog 2026-10-09 — Release gate: the shipped decks on f681f98 vs f9b9421 (DR-198) and 0e1e913 vs f681f98 (DR-122)

**Date:** 2026-10-09 (overnight, Christian's autopilot ruling of 2026-10-08)
**Topic:** Two deck comparisons requested by the Server session on Christian's word ("make sure this doesn't break any of our examples"): (A) the decks whose connectivity DR-198 changes, fixed tip f681f98 against the pre-DR-198 f9b9421; (B) the thin-shell decks with terminals on DrClean's 0e1e913 (DR-122) against f681f98
**AIs involved:** Claude Fable 5.1 (runs, analysis); Server session (request, release audit); DrClean session (0e1e913)
**Claude Confidence:** high on what ran; the numbers are executed, not reasoned
**Literature References:** N/A

## Method

Serial, debug, MKL flavor, two OpenMP threads, each run in a fresh directory with the same
regenerated mesh for both binaries and no `.bfm` (a matching cache skips the enrichment, Grok's
finding of 2026-10-08). Simulation times shortened to 80–100 ms (3-D) or kept (2-D, rlc); each run
capped at 25 min wall time, so "capped" below is the gate, not a failure. Compared: the setup log's
base-mesh counts (nodes, elements, faces, 1-cochains), the conductor edge count after cut and
interface processing, the free-unknown count, the number of accepted steps and the first rows of
the terminal history (`iv_results.csv` / `CircuitResults.txt`). Script and summaries in Claude's
scratchpad (`deck_gate.sh`, `gate_summary.py`, `gate_A_summary.txt`, `gate_B_summary.txt`).

## (A) f681f98 vs f9b9421 — what DR-198 changes

Base mesh identical on every deck. "edges" = conductor edges after enrichment, "free" = free unknowns.

| deck | edges ref → fix | free ref → fix | run | result |
|---|---|---|---|---|
| helix | 96703 → 96703 | 99060 → 99060 | 23 steps to 150 ms, exit 0 both | I/U agree to 1.3e-10. PASS (no shared cond\|air line) |
| tapestack2d_racetrack | identical | identical | 26 steps, capped at 1.36 ms both | 3.6e-11. PASS (2-D untouched) |
| tapestack3d | 90539 → 89339 | 97095 → 95883 | 11 steps to 100 ms, exit 0 | I_1 identical, U_1 to 1e-6; I_2/U_2 are 1e-10-level noise. PASS, dof set changed |
| tapestack3d_quench | as tapestack3d | as tapestack3d | 2 steps to 100 ms, exit 0 | U_1 to 1e-6, U_2 noise. PASS |
| pancake | 815073 → 802821 | 2235524 → 2235524 | 2 steps to 100 ms, capped both | first save agrees to 2.7e-7. PASS (merged copies were hanging edges) |
| rlc_circuit | 47743 → 47127 | 60362 → 59946 | capped at ~22 ms of 80 both | I_L1 to ~1.5e-5, U_L1 2.5e-3 at the first row then ~1e-6. PASS with a real small shift: this deck had free doubled edges |
| racetrack_usermat | 1490854 → 1447135 | 1286964 → 1286964 | capped before the first 50 s save (ref 19 s, fix 39 s simulated) | setup PASS; no result comparable in 25 min |
| pancake_tape | 906842 → 906712 | — → 2174470 | fails both | f9b9421: setup abort "couldn't figure out how to assign dependencies of dof … ( edge … )" (the DR-196 mechanism); f681f98: setup completes, first solve MUMPS −10 (singular) at every Δt cut, abort. Consistent with DR-196; no row filed |

Reading: wherever two conductor blocks met the air along a line, the merged per-block copies
remove conductor edges; where those edges were free dofs (rlc_circuit, tapestack3d) the answer
shifts in the fifth to sixth digit, where they were hanging (pancake, racetrack) nothing but the
edge count moves. The 1e-10-level secondary cut currents differ run to run on one binary by the
same relative amount as between binaries (two-thread reductions), so their 5–10 % "differences"
are noise. Caveat from the Server session: equal results on this toolchain do not show that the
shared-copy append order (an `unordered_map` iteration) is platform independent.

## (B) 0e1e913 vs f681f98 — what DR-122 changes

| deck | counts | run | result |
|---|---|---|---|
| tapestack2d_racetrack | identical | 26 steps, capped | 6.3e-11. PASS |
| tapestack2d_layered | identical (11 cochains, 38655 free) | 185 / 180 steps to 9.0 / 8.7 s, capped | 11 terminals, U to 6e-7 (the reversed return-pass arclength moves intermediate samples). PASS |
| tapestack3d | identical | 11 steps to 100 ms, exit 0 | I_1 identical, U_1 to 1e-7, I_2/U_2 noise. PASS |
| tapestack3d_quench | identical | 2 steps to 100 ms, exit 0 | U_1 to 1e-7, U_2 noise. PASS |
| racetrack_usermat | identical | capped before the first save (39 / 44 s simulated) | setup PASS, no result comparable |
| pancake | identical (802821 edges, 2235524 free) | first pair lost to the OOM killer (three 1.35 M-element runs overlapped, 31 GB resident); rerun alone: 2 steps to 100 ms, capped both | first save agrees to 3.7e-8. PASS |

No new `BELFEM_ERROR` fired on any deck in (B); the two always-active refusals DR-122 adds did not
trigger on the shipped thin-shell decks.

## Not done, by design

No register row filed (Server's rule for this gate; Christian decides on pancake_tape). No code
change. `main` untouched (pinned at v0.9.2 by the Server session). Residues: racetrack_usermat and
pancake need longer than 25 min for a first save under debug; the pancake_tape failure is DR-196.
