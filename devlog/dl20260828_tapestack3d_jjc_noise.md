# Devlog 2026-08-28 — tapestack3d j/jc Noise Onset (Jury Investigation)

**Date:** 2026-08-28
**Topic:** Root cause of the sudden j/jc noise at t = 2.125 s in the fine-cadence
tapestack3d run (piecewise law); blind jury round + post-freeze buffer-φ probes
**AIs involved:** Claude (primary, pre-registered), Codex + Grok (blind jury)
**Claude Confidence:** high on the measured facts; medium on the composite mechanism
**Codex Audit Confidence:** medium (~80 %) · **Grok:** high on the controller path
**Literature References:** Messe et al. 2023 §4, Eqs. 10–13 (checkerboarding, ε = 1e-11 —
note: §4, not the §2.7 CLAUDE.md's routing table says); Rhyner 1993; Riva 2021
**Verification:** numeric probes (exodus/iv/h5 table, scipy) + static source traces;
no executable gate ran — reviewed, not verified. Gates owed: clean one-key A/B,
Gauss-point/accepted-state instrumentation, bearing trace.

## Summary

The "noise" is a violent, ratcheting reorganization, not display jitter: between
t = 2.100 and 2.150 s max|B| triples (3.3→9.1 mT) under a smoothly rising 10 A
transport current, the bottom two tapes swap current (tape-1 transverse |Jx| ×40),
and the free axial generator I₁ jumps ×40 — landing on the state the clean riva
comparison reaches only at t ≈ 2.70 s. Verdict (3 voices, reconciled): the user's
options "bad material model" and "error in piecewise law" are REFUTED as the
trigger (table clamped below its 10 mT B-origin and smooth in θ; the two laws agree
to ~1e-28 in the entire visited regime, and only ρ-only Picard ran through onset).
The mechanism is the "undiscovered defect" family, three interacting parts:
(1) **~9 unpinned gauge constants** — each buffer patch carries one free φ constant
(measured: single per-patch value, redrawn each solve, envelope growing 4.9e-3 →
−3.9e7 over the run), and the *air* constant floats too (φ at the bearing corner
node wanders −1.9e6…+3.4e4 despite `fixed = 1`) → exactly singular matrix, MUMPS
condest 1e17–1e21; (2) resistively unpinned physical modes (ρ_HTS ~ 1e-30 at
max j/jc = 0.05); (3) **always-accept Picard semantics** — the reported residual is
the PRE-update entry state; the committed (Anderson-mixed) post-update state is
never residual-checked, min-2-iteration exit at 1e-7 with the absolute escape inert
(documented design, greg3 rationale — its contractivity assumption fails here).
The riva-vs-piecewise correlation is confounded (build, backend, executable,
algorithm keys, BDF order, save cadence — cadence is numerical via save-point Δt
clamping). Christian's mid-round buffer-unpinned hypothesis was decisive for (1);
his "h-connection to pinned air pins it" assumption is refuted by measurement —
h-coupling pins the gradient, never the constant.

## Key Findings

- Full report: `todo/tapestack3d_jjc_noise_2125ms.md` (mechanism, cleared items,
  ranked mitigation, provenance). Exchange record:
  `tmp/ai_exchange/review_tapestack3d_jjc_noise.md` (frozen pre-registration,
  both audits, verification of every cited line, reconciliation table).
- Register rows filed: DR-126 (unpinned φ constants / ineffective bearing, P1),
  DR-127 (accepted-state guard ruling, P1, Christian), DR-128 (sp-ap B-axis
  starts at 10 mT — no B-dependence over this deck's whole operating range, P2).
- Brief errors caught by the jury and withdrawn: κ·eps arithmetic, Anderson
  "window of 8" (depth is 3), κ-spike/onset alignment (intermittent), deck-diff
  staleness (coarse input.conf edited to Picard at 02:38 mid-session).
- By-catch: ghost + gauge chi are default-ON although the deck's penalty blocks
  are commented out (misleading deck comments); Messe 2023 routing mis-cite in
  CLAUDE.md (§2.7 → §4); DR-125 (last night) confirmed as the unrelated cause of
  the coarse 4.3 s band.

## Changes Made / Proposed

- No source edits (read-only diagnosis). New: this devlog, the todo report,
  three register rows, todo/README + devlog/README index lines.
- Proposed (ranked, ⚑ = Christian's ruling): clean one-key A/B; pin the gauge
  (repair air bearing + one anchor per buffer patch) ⚑; guard the accepted
  state (post-update residual or min-iterations/absolute-tolerance) ⚑; tighten
  tolerance toward 1e-11 only together with a floor escape; extend sp-ap below
  10 mT ⚑. I₁ stays free.

## Open Questions

- What the single prescribed magnetic dof actually pins; why written air φ floats.
- Gauss-point envelope vs nodal output (instrumentation owed).
- Whether backend (Blaze/Armadillo) or BDF order shifts the onset (A/B ladder).

## Files Updated

- todo/tapestack3d_jjc_noise_2125ms.md (new)
- todo/debt_register.md (DR-126..128)
- todo/README.md, devlog/README.md (index lines)

## Continuation (same night) — the bearing root cause

Christian asked why the air floats despite `bearing { nodes : 9 }`, then authorized
autonomous execution and left. Chase, in order: (1) static trace of every link in
the bearing chain — BC parse, vertex→node map, enrichment survival, create/link
order, impose site, fix-flag survival — ALL correct individually; (2) run-data
refutation of each candidate (the one prescribed dof is λ₀ from
`IWG_Maxwell::set_currents`, not the bearing; φ(node 9) unpinned from frame 1 in
BOTH the fine and coarse runs — the coarse φ-garbage reaches −2.1e9 at t = 4.35 s,
i.e. the 4.3 s wall band was fought on top of nine-orders-above-physical gauge
noise); (3) gdb on a serial boot (`cmake-build-debug/bearing_probe/`, isolated
deck copy): the bearing FIRES and FIXES the right dof — which turns out to be
HANGING (`mIndex = gNoIndex`, 1 source, weight 1): **node 9 lies on the periodic
target face z = L, and the periodic condensation eliminates the pinned dof from
the system.** Bearing × periodicity composition defect; silent because
`impose_dirichlet` accepts a hanging dof without complaint. Mesh check: hanging
low-id nodes are exactly the z = L corner twins {7,8,9,10,14,15,16,20}.

**Reproducer gate PASSED** (verified by execution, serial 30 ms probe): bearing
moved to node 4 (the z = 0 partner) → magnetic dof table flips
`prescribed (fixed) : 1 → 2`, free count −1. Step-1 conditioning stays ~1.0e18 —
as predicted, since the 8 buffer-patch constants still float (they, not the air,
dominate the condest; matches Christian's no-buffer diagnostic). Deck fix:
`bearing { nodes : 4 }`. Code fix proposed (NOT applied — needs approval + the
standing audit round): ERROR on hanging/empty bearing targets or chase the
weight-1 source chain; print which dofs are prescribed when the count is small.
DR-126 updated in place. Report: Addendum 2 of
`todo/tapestack3d_jjc_noise_2125ms.md`.

**Final gate (verified by execution, keyframe t = 25 ms):** φ(node 4) = 0.000000
exactly; φ(node 9) = 0 as well — its condensed dof follows the pinned source
through the weight-1 periodic pair, the condensation mechanism confirming
itself; far-air φ at physical scale (max 0.134 A), floating constant gone.
Conditioning stays 3.5e17–7.2e18 (buffer constants still dominate, as
predicted). The air-bearing half of DR-126 is closed at reproducer level.

## Continuation (afternoon) — the bearing code fix, landed and gate-verified

Christian confirmed the deck fix's effect in production (condest 1e17 → 1e8 with the
buffer off + node-4 bearing) and directed the code fix: pin through the periodic
condensation, or error loudly — then sharpened it mid-plan ("periodic partners share one
dof; only one dof pinned"). Full plan+audit → code+audit loop, both vendors, exchange
`review_bearing_pin.md`. The plan round's decisive finding (both voices): the parallel
Dirichlet channel is `synchronize_dirichlet_bcs()` (master's flags redistributed in
init_matrices) and bearings only ever link dofs on the master — so the fix collapsed to
ONE site, `Bearing::impose_dirichlet`: reroute a pin on a condensed dof onto its single
unit-weight source (value / weight, Codex catch, so the alias lands exactly), loud
rank-0 BELFEM_ERRORs (with point/node attribution, Grok catch) for empty/unknown targets
and any other condensation shape. Maxwell, thermal, and gauge inherit through the shared
site. Input contract updated both artifacts. Code-round fixes applied: value exactness,
error attribution, source-node id in the info line, rank-0 message gate, null-source
guard, comment corrections. Accepted divergence: always-active ERROR on the flattening
invariant (Grok endorsed, Codex preferred ASSERT). Residues filed in DR-126: gauge
first-pin on a condensed node lands after the graphs freeze (needs a factory-time
bearing on the same point); SideSet::impose_dirichlet is the unaudited sibling.

**Gates (verified by execution, Christian's 14:11 static rebuild):** nodes:9 → reroute
line once, fixed 1→2, φ(9)=φ(4)=0 exactly, fields numerically identical to the direct
pin; same at np=2; nodes:4 byte-identical behavior, no reroute line; nodes:999 → loud
abort naming the deck list. `make check-fast` owed (production run occupies the machine).

## Continuation — jury audit of Christian's find_autopins() draft (probe stage)

Blind jury on the new auto-pin generator (exchange `review_find_autopins.md`; the
call-site exit(0) is Christian's deliberate probe stop, excluded by instruction).
Convergent verdict, all citations re-verified: the probe machinery runs, but the body
needs the §8 rework before it can ship. THREE P0s: (1) both periodic-duplicate loops
dropped the `d <` from their condition (infinite loop + unchecked OOB the moment a
periodic partner has duplicates — correct in the fn_Mesh_symrcm_nodes.cpp original it
was adapted from; fix worth making even in the probe); (2) find_connected_partitions
runs on LIVE mesh nodes, so dfs overwrites node owner()/level() — the MPI ownership
map — and the cleanup restores only indices; poisoned owners would persist into the
.bfm save once the probe exit goes (Christian's own §8 A1 ruling — mirror graph, idiom
at cl_CutFactory.cpp:1066-1100 — already prescribes the cure); (3) the symrcm pass
steers each component's representative onto a pseudo-peripheral node — the periodic-
plane class whose pin evaporates, violating A3/A7 (smallest non-periodic id). P1s: air
double-pin vs the deck bearing (A2 skip-if-pinned), .bfm-reload path skips the block
entirely (A5), conductor-mediated adjacency can merge separate φ components on other
decks, bearing map is vertex-keyed (consumer must not pass node ids). Benign-but-
fragile: dfs shares flag bit 0 with the membership mark (sequencing saves it; the
mirror graph removes the dance). Expected count gate on this deck: 9 (air + 8 buffer
patches), aPins(0) = air. Reviewed, not verified — nothing executed.
