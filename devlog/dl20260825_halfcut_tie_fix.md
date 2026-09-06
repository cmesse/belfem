# Pure Half-Cut Periodic Ties: the Cap-Corner Fix (Option B), Full Round

**Date:** 2026-08-25 (third arc; mechanism identification is
dl20260825_cap_corner_defect.md)
**Purpose:** Record the plan+audit → code+audit round that landed the
conditional tying of pure half-cut edges across the periodic seam, fixing
the corc cap-corner antisymmetry's root cause.
**Modules:** src/mesh (PeriodicityFactory, match_edges), src/fem/maxwell
(PART 2 hanging guards)
**Exchange threads:** `tmp/ai_exchange/corc_corner_fix_plan.md`,
`tmp/ai_exchange/corc_corner_fix_code.md`
**Status:** LANDED AND GATED same evening. Gate 1 GREEN (405 tied,
66 pure half-cut, 0 untied; census 66x tied-pure-halfcut; no aborts).
Gate 2 GREEN (through PART 2, no allocate assert). G1' PASS beyond its
promise: corner-table deviations collapsed from +-45% antisymmetric
(outer) / +-15% (inner) to <= ~6% near-symmetric — the cap-mirrored
circulation is dead. G2' PASS: jump plateau uniform at cap and interior;
the per-cap pinch survives cap-symmetrically (Option A backlog, as
booked). First gate attempt was void (stale build/bin/hphirun from the
second build tree — worth remembering as a gate-hygiene trap). Gate 3 GREEN (helix:
256 pure half-cut tied, 0 untied/displaced — the reference conductor
deck carried the same defect class silently). Gate 4 GREEN (full make
check, after restoring the betterdoc-merge-regressed rel-tol default
1e-9 -> 1e-10 caught by the DefaultParametersValid pin). Remaining: G3'
on a longer run; probe removal + stale-comment sweep.

## Root cause (measured; see dl20260825_cap_corner_defect.md)

The seam policy left half-cut edges (one duplicated endpoint) untied across
the periodic map. At the tape terminals these are the corner-adjacent and
sheet-adjacent conductor/shell edges — free edge_h dofs, independent per
cap — so the discrete problem was not periodic exactly at the tape corners
and each cap relaxed the corner pinch independently (the antisymmetric,
cap-mirrored, sigma-selected corner circulation; 2.4x/0.23x at the
immediate cap row; identical in serial and 10-rank runs, identical in both
homology lineages).

## Plan round

Pre-registered plan (tie pure half-cuts) audited blind by Codex + Grok.
Kills and hardenings adopted: endpointwise four-node is_periodic()
predicate (Codex's (orig,dup)x(dup,orig) counterexample); backup
containers are NOT available at the edit site (use Node::periodic());
Grok's exactness split (condensed phi-region: tie redundant; free edge_h:
tie is the missing constraint; 1:2 cut-trace twins: must never tie —
INC-034); the PART 2 allocate_source_container landmine; G-gates rewritten
(the fix SELECTS the periodic member — corner solution expected to move;
Option B kills the ANTISYMMETRY, the per-cap pinch survives and belongs to
Option A, deliberately out of scope). Blocking census probe then measured:
all 66 half-cuts on corc_solder are both-sides-half-cut, all-four-
endpoints-paired free-H edges (12 tape corners + 54 sheet-adjacent solder
edges), zero 1:2, zero twin-of-tied. Christian's routing: B.

## Implementation (as landed)

- `match_edges` half-cut branch: candidates {E, F, swap-decision} collected
  WITHOUT mutation when both sides are half-cut and all four raw endpoints
  are periodic-paired, with the endpointwise raw-carrier test (each
  endpoint's periodic() must be the corresponding endpoint after
  alignment) as final arbiter. DEFERRED application after the main loop —
  candidates tie only if neither edge was claimed — making pure-twin-wins
  order-independent by construction (Codex C1 fix, in-round). Displaced
  candidates stay untied (censused "halfcut-displaced"). New counter
  tNumHalfCutTied; log line extended.
- PART 1-style EDGE-source skip guard added to all four
  hang_thinshell_edges_* loops (newly tied slave shell edges keep their
  periodic tie; master side hangs normally; cascade closes the chain).
- Temporary env-gated census probe (BELFEM_PROBE_FUSED_ROWS) tags
  tied-pure-halfcut / policy-skip / halfcut-displaced / twin-of-tied.

## Code round

Codex: no-land until C1 (alias inversion: a later pure twin would be
aliased to an earlier half-cut tie) — fixed via the deferral, Codex's own
structurally-safe option. Grok on the corrected artifact: LAND, C1-C6 all
pass, incl. verification that +1.0 is the correct dependency weight (the
alignment is baked into node order; a sign would double-count the swap),
compute_edge_directions rebuilds Whitney signs post-swap, and the PART 2
node-hang skip is REQUIRED (hanging a tied slave on slave-cap phi would
destroy the identification). Residuals deferred to probe removal: stale
"untied by policy" comments, tNumHalfCut folding three sub-classes,
vector-vs-Cell nit.

## Gate ladder (pending)

1. corc_solder debug, probe on, fresh bfm, setup-only: 66 tied-pure-halfcut
   / 0 policy-skip / 0 displaced; falsifiers = compaction abort,
   orientation assert.
2. Through PART 2: falsifier = allocate_source_container assert.
3. helix / Twisted_Conductor_Periodic setup-only (conductor-watch):
   displaced>0 is deferral SUCCESS; falsifier = half-cut tied to a pure.
4. check-fast (non-regression).
Then physics: G1' corner table cap-periodic (antisymmetry dead; pinch
survives -> Option A backlog), G2' jump plateau + still-pinched corners,
G3' interior/IV within tolerance (not bit-identity), new-to-new memdump
and BFM only (old caches will not match).

## Cleanup (done, same evening, after all four gates green)

All session probes stripped: orient_terminal_curves print+flip,
probe_edge_edge_row (DofData), probe_halfcut_census (PeriodicityFactory)
— Christian's pre-existing BELFEM_PROBE_FUSED_ROWS edge-on-node hook in
DofData is untouched. Stale "untied by policy" comments in match_edges
(header, counter, classification block) refreshed to the conditional-tie
policy. All three files syntax-checked.
