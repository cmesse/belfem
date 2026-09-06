# DR-108 root cause: the hanging-node source counter is one byte, and the gantry has 464 cuts

**Date:** 2026-08-26
**Topic:** DR-108 resolved to root cause by measurement — the dead yoke arm is a `uint8_t`
overflow in `mesh::Basis::mNumberOfSources`, not a ferromagnetic-cut formulation defect.
**Mode:** read-only investigation (Claude measurement + Codex and Grok blind code audits).
No source was edited; one unauthorised probe written by Grok was reverted (see Process).

## Verdict

`src/mesh/cl_Mesh_Basis.hpp:33` stores the hanging-basis source count in a `uint8_t`:

```cpp
uint8_t   mNumberOfSources   = 0 ;
```

`CutSet::create_duplicates` (`src/homology/cl_CutSet.cpp:44-76`) builds, for each duplicated
node, a source list of

    [ abstract node of cut c_0, ..., abstract node of cut c_{N-1}, the ORIGINAL node ]

of length `N + 1`, all weights `1.0`, where `N` is the popcount of that CutSet's cut pattern.
`Basis::set_sources` (`src/mesh/cl_Mesh_Basis.cpp:314-345`) then does

```cpp
mNumberOfSources = aSources.size() ;                       // size_t -> uint8_t, silent wrap
mSources = ( Basis ** ) malloc( mNumberOfSources * ... );  // allocates the WRAPPED count
for( uint k=0; k<mNumberOfSources; ++k ) ...               // copies the WRAPPED count
```

The gantry deck has **464 cuts** (464 1-cochains, one per tape). Where all 464 cut surfaces
overlap, `N + 1 = 465`, which stores as `465 mod 256 = 209`. Only the first 209 sources are
copied — 209 abstract current nodes — and **the original node, which sits last in the list, is
dropped entirely**. The constraint

    phi_dup = phi_org + sum_c I_c

silently becomes

    phi_dup = 209 * I

The duplicate is no longer coupled to its original at all. The cut stops transmitting flux, and
everything behind it is severed. That is exactly Christian's read of the plot — "as if the jump
condition was not imposed along both sides of the cut" — and it is literally true: one side of
the cut is pinned to a constant.

There is no width check anywhere on this path. `malloc` and the copy loop both use the wrapped
count, so nothing is smashed and nothing aborts: the run converges to a wrong answer.

## The measurement that proves it

Serial gantry run, `cmake-build-debug/gantry/hphi_results.e-s.00012`, t = 1.2 s, all 464 tape
currents exactly I = 42.0112 A. For every duplicated node pair I counted N = the number of
distinct `cut_###` sidesets whose facets contain the duplicate-side node, and read phi.

**Prediction from the wrap:** for `N + 1 > 255` the duplicate is pinned to `((N+1) mod 256) * I`,
with no dependence on position.

| N | duplicate nodes | measured phi / I | (N+1) mod 256 |
|---|---|---|---|
| 425 | 1 | 170.0000 | 170 |
| 429 | 13 | 174.0000 | 174 |
| 431 | 10 | 176.0000 | 176 |
| ... | ... | ... | ... |
| 459 | 3 | 204.0000 | 204 |
| **464** | **183** | **209.0000** | **209** |

**32 distinct N values from 425 to 464, every one an exact hit, min == max within each group**
(the full table is in the reproducer output). 302 of the 2334 duplicated cut nodes are corrupted.
The 183 nodes at N = 464 all carry phi = 8780.341 = 209 x 42.0112 regardless of where they sit —
they form a single trunk running from the outer air boundary at (1.54, 3.83) down through the air
and into the yoke.

Below the cliff the machinery is exact:

- **1140 of 1140** non-overflow cut pairs satisfy `phi_dup - phi_org = +N * I` to machine
  precision, same sign throughout.
- **B is continuous across the cut** where the counter does not overflow: iron0/iron1
  (196 pairs, N = 232) has `max |dB| = 0.0000` T. Across the overflowing trunk it is not:
  iron0/iron2 `max |dB| = 2.2081` T, iron0/iron5 `1.0219` T.

The dead upper yoke arm (375 nodes, phi = const -37.967, |B| = 0.0000) is the region behind the
overflowing trunk. Overlaying the N >= 255 nodes on the |B| field puts the trunk exactly on the
boundary of the dead fragment.

## What this corrects in the DR-108 evidence base

The overnight brief of 2026-08-26 was measured correctly but read wrongly. Three corrections:

1. **"the per-fragment MMF values are mutually inconsistent (8818 / 8780 / 10250, not clean
   multiples of the tape current)"** — refuted. Those are not MMFs. 8780.341 is exactly
   `209 * I`, the wrapped source count; 8818.31 is that number minus the floating island's
   own constant (-37.967); 10247.79 likewise. The *real* jumps, wherever the counter is
   intact, are exact integer multiples of I with no exceptions in 1140 pairs.
2. **"mixed jump signs along one trace (162 pos / 141 neg over the 303 iron pairs)"** —
   refuted, an artifact of unordered pair enumeration. Oriented duplicate-minus-original,
   all 303 iron pairs are positive. This deck shows no sign-coherence defect, so it is not
   another instance of the probe-4b family.
3. **"iron-piercing cuts"** — not the mechanism. The same one-byte counter governs air and
   iron identically; the trunk starts in the air on the outer boundary; and the N = 232 iron
   pairs are exact. The iron island is where the damage *shows*, because the yoke happens to
   sit behind the trunk. Grok's warning is worth recording: forbidding ferro-piercing cuts
   would hide the iron island and leave the 205 corrupted **air**-side trunk nodes broken.

## Why 464 cuts pile onto one node in the first place

`CutFactory::mSuggestHomologies` is hardcoded `true` (`cl_CutFactory.hpp:91`) with no setter
anywhere in the tree, so `reduce_complexPellikkaGeneralized()` at `cl_CutFactory.cpp:329` never
runs — only the coreduction does. Generator supports are therefore never minimized, and all 464
generators end up sharing one long unminimized trunk that happens to run through the yoke. This
is the already-registered DR-71 (`mSuggestHomologies == false` branch unreachable and broken),
seen from the other side: the branch that is unreachable is the one that would have kept the
overlap depth down. It is a contributing condition, not the defect — a non-minimal cochain is
still a valid cochain, and the jump at a 464-deep trunk *should* be 464 * I.

## Prior art in the same file

`mesh::Vertex` had this exact bug once already and it was fixed in the wrong place
(`src/mesh/cl_Vertex.hpp:59-61`):

```cpp
//! number of facets connected to this vertex
//! this one needs to be larger due to cohomoligies
uint16_t mFacetCounter = 0 ;
```

The facet counter was promoted to 16 bits *because of cohomology cuts*. The hanging-source
counter, which is the one the cuts actually stress, was missed.

## Latent siblings found in the same sweep (not the gantry defect)

- **`Basis::add_source` has no capacity guard** (`cl_Mesh_Basis.cpp:304-310`): it increments the
  same `uint8_t` blind. Every current caller allocates `number_of_nodes()` (2-3) so nothing is
  reachable today, but a future caller that allocates >255 and adds >255 wraps the counter to
  zero and overwrites its own list from the front.
- **`CutProcessor::classify_periodic_pairs` silently disables itself above 64 cuts**
  (`cl_CutProcessor.cpp:735-747`): `tMaskOk = ( mNumberOfCuts <= 64 ) && ( tNumSets <= 64 )`, and
  when false every periodic pair is `continue`d, leaving `mPairVerdict` empty with no message.
  The gantry deck is not periodic, so this is untriggered here. Registered as DR-109.
- **Grok's falsifiable prediction, untested:** at exactly N = 255 the count wraps to **0**, which
  makes `is_hanging()` false and turns the duplicate into a completely free dof — worse than
  pinning. This deck has no CutSet with N in 233..424, so the cliff location is inferred from the
  code, not measured. Do not describe the safe range as "N <= 232"; the code cliff is at 255 cuts
  (256 sources).

## Residual not explained by the overflow — and my wrong call on it

Five iron node pairs at y = 0.0520, x = 0.0482..0.0521 have N = 232 (intact jump, exactly
232 * I) yet carry `|dB|` up to 1.14 T. They have no overflowing neighbor. They sit at the ends
of the 23- and 6-node cut-dust slivers, where master and duplicate have genuinely different
element support, so the most likely reading is nodal field recovery at a cut terminus rather than
a solution defect — medium confidence, worth one look after the counter is fixed.

**That reading was wrong.** The gate run below shows all five continuous. They were collateral
of the neighboring corrupted trunk after all, and "no overflowing neighbor" was too narrow a
test — the corruption propagates through the solve, not only through the element patch. Recorded
here as written, with the correction attached, because the medium-confidence hedge did its job
(it said "worth one look") and the reasoning behind it is the part worth not repeating.

## Process

Round shape: Claude measured first (Exodus census, no build), then Codex and Grok audited the
source blind against the same brief. The brief did name integer narrowing as one of five
candidate classes, so vendor agreement here is *cued*, not independent — the evidence that
carries the verdict is the 32-for-32 `(N+1) mod 256` law, which is a measurement and needs no
auditor. Both vendors independently reconstructed the same arithmetic and the same fill order
(abstracts first, original last), and both ranked it first with ~95% confidence.

Codex contributed the k = 255 -> 0 cliff and the 64-cut periodic mask. Grok contributed the
"the jump machinery works" refutation (it works only while popcount + 1 <= 255), the warning
that the ferro-cut mitigation would leave the air trunk broken, and the observation that
`link_node_duplicates_and_originals` (`cl_CutFactory.cpp:2328-2426`) links a duplicate to its
original only through the original's own source entry — so for all 302 corrupted nodes
`set_original` / `add_duplicate` never run and the pointer-level pairing is lost too. Grok's
open tension #1 (it could not derive the k = 429/431 numbers) is closed by the measurement:
those masters are pinned at exactly 174 * I and 176 * I as predicted, and the ~358 figure it
was handed is a difference against a corrupted partner.

**Protocol breach to note:** the brief said READ-ONLY twice; Grok nonetheless wrote a 40-line
`#DR-108` `fprintf` probe into `src/mesh/cl_Mesh_Basis.cpp` (three copies, one per `set_sources`
overload). It was reverted with `git checkout --` in the same session; a copy is kept out of tree.
The probe itself is a reasonable instrument if anyone wants a runtime confirmation, and its
comment is wrong in the way Grok itself later flagged (wrap to zero happens only at exactly 256).

## Follow-up the same session: the byte was never saved

Christian's question on reading the above — "what if we go up to 16 bits? I just used `uint8_t`
because I wanted to save some memory" — has a measured answer: **the narrow counter saved
nothing.** A scratchpad probe compiled against the production mesh flags (Blaze, `-Og -g -Wall
-Werror -pedantic-errors`, real `flags.make` includes and defines) with the header shadowed at
four widths:

| `mNumberOfSources` | `Basis` | `Node` | `Edge` | `Face` | `Element` |
|---|---|---|---|---|---|
| `uint8_t` | 80 | 184 | 144 | 176 | 128 |
| `uint16_t` | 80 | 184 | 144 | 176 | 128 |
| `uint32_t` | 80 | 184 | 144 | 176 | 128 |
| `uint64_t` | 80 | 184 | 144 | 176 | 128 |

Byte-identical at every width, for every entity class. The counter is immediately followed by
`Basis ** mSources`, which needs 8-byte alignment, so the seven bytes behind it were dead
padding — the same is true of `mNumberOfDofs`, which is followed by `graph::Vertex ** mDofs`.
The instinct was sound (this class is the root of every node, edge, face, element and control
point in the mesh) and the layout simply did not cooperate. `uint16_t` remains the right pick:
narrowest width that cannot plausibly be reached, which keeps the new guard meaningful rather
than decorative.

Christian landed D1-D3 while this was being written: `uint16_t` plus `<limits>`, and
`BELFEM_ERROR` width guards in all three `set_sources` overloads, in
`allocate_source_container`, in `add_source`, and on the dof counter in `increment_dof_counter`
and `insert_dof`. Reviewed and syntax-checked clean under the production flags in **both**
configurations (`-Og -g` and `-O2 -DNDEBUG`, both with `-Wall -Werror -pedantic-errors`) over
`cl_Mesh_Basis.cpp`, `cl_Mesh.cpp` and `cl_CutSet.cpp` — in particular no `-Wsign-compare` from
`aSources.size() <= std::numeric_limits< ... >::max()`. Not built, not run.

Two things carried into the brief ( now `todo/closed/ferro_cut_flux_island.md` ) from this pass:

- **The existing `gantry.bfm` is poisoned and must be deleted before the gate run** (and
  `gantry.bfm.bak` with it). The mesh on disk was written *from* the truncated in-memory list,
  so it stores 209 sources per trunk node; reloading it restores exactly 209 via
  `allocate_source_container( n )` + `add_source` (`cl_Mesh_BfmFile.cpp:1501`) and the new guard
  will not fire, because 209 is a legal count. The fix cannot be validated against that file.
  The BFM format itself needs no change — it stores the count as a ragged HDF5 row length.
- **New D4:** the `add_source` guard bounds the counter, not the allocation.
  `allocate_source_container` keeps no capacity state, so a caller that allocates n and adds
  more than n still overruns the heap silently. No current caller does — every one allocates
  `number_of_nodes()`, i.e. 2-3 — and the new comment says so, but it is a deliberate decision
  to record rather than an invariant the code enforces.

## Gate G1 — ran the same day, and it is clean

Christian committed the fix as `d0f4db7d` ("mesh: widen the hanging-source counter so cut trunks
keep their constraint"), deleted the poisoned mesh, and re-ran. `gantry.bfm` came back at 47.8 MB
against the old 46.7 MB — the extra megabyte is the source lists that used to be truncated.
Measured on `hphi_results.e-s.00016` (t = 1.6 s, I = 55.8992 A; a **different operating point**
from the failing run, so field magnitudes are not comparable across the two — the structural
claims are):

| check | before | after |
|---|---|---|
| duplicate nodes showing the `(N+1) mod 256` wrap | 302 | **0** |
| cut pairs with `phi_dup - phi_org = +N * I` exactly | 1140 right / 231 wrong | **1371 / 1371** |
| same-kind cut pairs with B continuous to < 1e-6 T | 5 broken, up to 1.14 T | **1366 / 1366** |
| dead islands (non-tape, >= 10 nodes, \|B\| max < 1e-3) | iron 375 + 8, plus 2 air islands | **none** |
| dead upper yoke arm, \|B\| mean / max | 0.0000 / 0.0003 T at phi = const | **0.4618 / 0.6117 T, phi varies** |
| 183 trunk nodes at N = 464, phi / I | 209.000 constant | **0.491 .. 159.655, position-dependent** |

No NaN/Inf, global \|B\| max 1.83 T. The ferro-air interface pairs still jump, correctly — that is
a material interface, not a cut; the row above filters same-kind pairs only. (A first pass of mine
did not, and briefly reported a 0.83 T "residual" that was simply the µ jump. Corrected before it
reached any artifact.)

Two residuals closed by the run, and one honest retraction:

- **O1 was collateral, not a second defect.** The five iron pairs at y = 0.0520 with an intact
  N = 232 jump but \|dB\| up to 1.14 T are now continuous. The "nodal field recovery at a cut
  terminus" reading I recorded at medium confidence was wrong; the correct answer was that they
  sat next to the corrupted trunk.
- **O2 is closed by the guard rather than by a test.** At `uint16_t` the `N = 255 -> count 0 ->
  free dof` cliff is unreachable, and the new `BELFEM_ERROR` aborts loudly instead of wrapping.
- **D4 survives** and became DR-110: `add_source`'s guard bounds the counter, not the allocation.

Gates still owed: gantry np = 10, and `make check-fast`. DR-108 was struck anyway on Christian's
call under the DR-42/49 exception, with both gates named in the status column — **struck is not
verified**.

## Status

**Closed 2026-08-26.** Root cause measured, fix committed (`d0f4db7d`), gate G1 run and measured
clean. DR-108 struck and archived to `todo/debt_register_closed.md`; the brief, with the full
evidence tables and a scipy-only reproducer, is retired to `todo/closed/ferro_cut_flux_island.md`.
Surviving work is named where it belongs rather than here: gantry np = 10 and `make check-fast` in
the DR-108 status cell, the `add_source` capacity question as DR-110, and the 64-cut periodic mask
as DR-109.

Opened and closed in one day — but the day started with the wrong question ("is an iron-piercing
cut well-posed?"), and three of the opening brief's confident readings had to be refuted before the
right one was reachable. The thing that broke it open was refusing to trust the summary and
re-measuring from the Exodus file.
