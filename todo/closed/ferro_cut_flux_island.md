# Ferro-cut flux island (gantry dead yoke arm) — investigation brief

**Date:** 2026-08-26 (opened), 2026-08-26 (root cause found)
**Register:** DR-108
**Status:** CLOSED 2026-08-26. Root cause measured, fix landed and committed (`d0f4db7d`), gate G1 run
and measured clean. Gates G2 (np = 10) and G3 (`make check-fast`) survive in the register;
DR-108 struck under the DR-42/49 exception. Residual D4 spun out as DR-110. The defect is NOT
ferromagnetic and NOT a cut-formulation question — it is a one-byte counter overflow in
`mesh::Basis::mNumberOfSources`, and widening it costs zero bytes (measured).

## Root cause

`src/mesh/cl_Mesh_Basis.hpp:33` — `uint8_t mNumberOfSources = 0 ;`

`CutSet::create_duplicates` (`src/homology/cl_CutSet.cpp:44-76`) gives each duplicated cut node
a source list of `N + 1` entries, all weight `1.0`:

    [ abstract node of cut c_0, ... , abstract node of cut c_{N-1}, ORIGINAL node ]

where `N` = popcount of that CutSet's cut pattern. `Basis::set_sources`
(`src/mesh/cl_Mesh_Basis.cpp:314-345`, and the two sibling overloads at `:192` and `:232`)
assigns `mNumberOfSources = aSources.size()` — a silent `size_t -> uint8_t` narrowing — then
`malloc`s and copies **the wrapped count**.

The gantry has 464 cuts. Where all 464 cut surfaces overlap, `N + 1 = 465 -> 209`. Only the
first 209 sources survive, all of them abstract current nodes; **the original node, last in the
list, is dropped**. The intended constraint

    phi_dup = phi_org + sum_c I_c

degenerates to `phi_dup = 209 * I`. The duplicate is decoupled from its original, the cut stops
transmitting flux, and the region behind it is severed — which is precisely Christian's read of
the plot ("as if the jump condition was not imposed along both sides of the cut").

No width check exists on this path. `malloc` and the copy loop both use the wrapped count, so
nothing is smashed and nothing aborts: the run converges to a wrong answer.

## The proof

Measured on `cmake-build-debug/gantry/hphi_results.e-s.00012` (serial, t = 1.2 s, all 464 tape
currents exactly I = 42.0112 A). For every duplicated node pair, N = number of distinct
`cut_###` sidesets containing the duplicate-side node.

Prediction from the wrap: for `N + 1 > 255`, phi is pinned to `((N+1) mod 256) * I`, independent
of position. Measured, per N group, min == max:

| N | dup nodes | phi / I measured | (N+1) mod 256 |
|---|---|---|---|
| 425 | 1 | 170.0000 | 170 |
| 429 | 13 | 174.0000 | 174 |
| 431 | 10 | 176.0000 | 176 |
| 432..459 | 3-4 each | 177 .. 204 | 177 .. 204 |
| **464** | **183** | **209.0000** | **209** |

**32 distinct N values, 32 exact hits.** 302 of 2334 duplicated cut nodes are corrupted. The
183 trunk nodes at N = 464 all read phi = 8780.341 = 209 x 42.0112 regardless of position; they
form one trunk from the outer air boundary (1.54, 3.83) down through the air into the yoke.

Below the cliff the machinery is exact:

- **1140 / 1140** non-overflow cut pairs satisfy `phi_dup - phi_org = +N * I` to machine
  precision, one consistent sign.
- **B is continuous across intact cut pairs**: iron0/iron1 (196 pairs, N = 232)
  `max |dB| = 0.0000` T. Across the overflowing trunk it is not: iron0/iron2 `2.2081` T,
  iron0/iron5 `1.0219` T, iron0/iron4 `0.3302` T.

The dead upper yoke arm (375 nodes, phi = const -37.967, |B| = 0.0000) is the region behind
the overflowing trunk.

## Corrections to the original evidence base

The 2026-08-26 opening brief measured correctly and read wrongly:

- **"per-fragment MMF values are mutually inconsistent (8818 / 8780 / 10250)"** — refuted.
  Those are not MMFs. `8780.341 = 209 * I` is the wrapped source count; `8818.31` is that minus
  the floating island's own constant `-37.967`; `10247.79` likewise. Every intact jump in the
  mesh is an exact integer multiple of I.
- **"mixed jump signs (162 pos / 141 neg over 303 iron pairs)"** — refuted, an artifact of
  unordered pair enumeration. Oriented duplicate-minus-original, all 303 are positive. This deck
  shows no sign-coherence defect; it is not another probe-4b instance.
- **"iron-piercing cuts"** — not the mechanism. The counter governs air and iron identically,
  the trunk starts in air on the outer boundary, and the N = 232 iron pairs are exact. **205 of
  the 302 corrupted nodes are air nodes.** Forbidding ferro-piercing cuts would hide the iron
  island and leave the air trunk broken.

## Why 464 cuts pile onto one node

`CutFactory::mSuggestHomologies` is hardcoded `true` (`cl_CutFactory.hpp:91`) with no setter in
the tree, so `reduce_complexPellikkaGeneralized()` (`cl_CutFactory.cpp:329`) never runs — only
the coreduction does. Generator supports are never minimized, so all 464 generators share one
long trunk. This is DR-71 seen from the other side. It is a contributing condition, not the
defect: a non-minimal cochain is still valid, and a 464-deep trunk *should* carry `464 * I`.

## Prior art in the same file

`src/mesh/cl_Vertex.hpp:59-61` already carries this exact fix, applied to the wrong counter:

```cpp
//! number of facets connected to this vertex
//! this one needs to be larger due to cohomoligies
uint16_t mFacetCounter = 0 ;
```

The facet counter was widened *because of cohomology cuts*. The hanging-source counter, which is
what the cuts actually stress, was missed.

## Root questions — status

- [x] Is the per-segment MMF assignment wrong, or the fragmented routing itself?
      **Neither.** The assignment is exact wherever the source counter does not overflow
      (1140/1140 pairs). The routing is non-minimal but valid.
- [x] Are the mixed jump signs the probe-4b sign-coherence defect surfacing here?
      **No.** The mixed signs were a measurement artifact; all iron pairs are sign-coherent.
- [x] Why does the cut fragment into slivers inside the iron — is `reduce_complexPellikka`
      bypassed by hardcoded suggestions? **Yes**, `mSuggestHomologies` is hardcoded `true`
      (DR-71). It explains the sprawl and the dust, not the dead arm.
- [x] Literature check first (Pellikka 2013, Alves 2022b, Gross & Kotiuga): is an iron-piercing
      cut well-posed? **Moot for this defect** — the failure is arithmetic, not formulational,
      and the same code path runs in air. The question stays open on its own merits for the
      mitigation idea below, and should be answered before that idea is implemented.

## Fix shape — D1/D2 landed by Christian 2026-08-26, uncommitted

**Widening is free.** Measured with a scratchpad probe compiled against the production mesh
flags (`-Og -g -Wall -Werror -pedantic-errors`, Blaze, DEBUG): `sizeof` is byte-identical at
`uint8_t`, `uint16_t`, `uint32_t` and even `uint64_t` for the counter —

| counter width | `Basis` | `Node` | `Edge` | `Face` | `Element` |
|---|---|---|---|---|---|
| uint8_t (old) | 80 | 184 | 144 | 176 | 128 |
| uint16_t (new) | 80 | 184 | 144 | 176 | 128 |
| uint32_t | 80 | 184 | 144 | 176 | 128 |
| uint64_t | 80 | 184 | 144 | 176 | 128 |

The counter is immediately followed by `Basis ** mSources`, which needs 8-byte alignment, so the
seven bytes after it were dead padding already. The byte saved by `uint8_t` was never saved.
`uint16_t` is the right choice anyway — it is the narrowest width that cannot plausibly be hit,
and it keeps the guard meaningful rather than decorative.

- D1 [x] `Basis::mNumberOfSources` is now `uint16_t` with the provenance comment
        (`cl_Mesh_Basis.hpp:38`). `mNumberOfDofs` deliberately stays `uint8_t` and is guarded
        instead — one dof per node per type is what is in play, and the guard makes a future
        surprise loud rather than silent.
- D2 [x] `BELFEM_ERROR` width guards in all three `set_sources` overloads, in
        `allocate_source_container`, in `add_source`, and on the dof counter in
        `increment_dof_counter` / `insert_dof`. `BELFEM_ERROR` is correct here (setup code,
        runs once) so the release build cannot swallow it. `<limits>` added to the header.
        Syntax-checked clean under the production flags in **both** configurations
        (`-Og -g -Wall -Werror -pedantic-errors` and `-O2 -DNDEBUG` + the same warnings) over
        `cl_Mesh_Basis.cpp`, `cl_Mesh.cpp` and `cl_CutSet.cpp`; no `-Wsign-compare`, the signed
        `numeric_limits` operand is a known-non-negative constant.
- D3 [x] Rest of the hanging-source path re-read after the widening: `cl_Mesh_BfmFile.cpp`
        (`:1452`, `:1501` and the edge/face/facet/control-point siblings) and
        `cl_Mesh_Distributor.cpp` (`:608`, `:667-745`, `:780-909`) both carry the count in `uint`
        and the BFM format stores it as a ragged HDF5 row length, so no format change is needed.
- D4 [ ] Residual in `add_source`: the new guard bounds the **counter**, not the **allocation**.
        `allocate_source_container` keeps no capacity state, so a caller that allocates n and
        adds more than n still overruns the heap silently. No current caller does — every one
        allocates `number_of_nodes()` (2-3) — and Christian's comment says as much. Decide
        deliberately whether to add a capacity member or leave it documented.
- G1 [x] Gate: rebuild, **delete `gantry.bfm` AND `gantry.bfm.bak` first**, rerun serial gantry.
        This is not housekeeping — the mesh on disk was written *from* the truncated in-memory
        list, so it stores 209 sources per trunk node. Reloading it restores exactly 209 through
        `allocate_source_container( n )` + `add_source` (`cl_Mesh_BfmFile.cpp:1501`) and the new
        guard will **not** fire, because 209 is a legal count. The fix cannot be validated against
        the existing `.bfm`; it has to re-run the cut pipeline. Expected afterwards: the 183 trunk
        nodes carry `phi_org + 464 * I`, `max |dB|` across every cut pair returns to ~0, and the
        375-node yoke arm carries field. Re-run the reproducer below and require 0 corrupted nodes.
- G2 [x] Gate: gantry np = 10 (the distributed path ships the already-truncated list, so it
        must be re-measured, not assumed).
- G3 [x] Gate: `make check-fast`.
- O1 [ ] Residual: five iron pairs at y = 0.0520, x = 0.0482..0.0521 have an intact N = 232 jump
        but `|dB|` up to 1.14 T, with no overflowing neighbor. They sit at the ends of the
        23- and 6-node cut-dust slivers, so nodal field recovery at a cut terminus is the likely
        reading (medium confidence). Re-measure after G1 before opening anything.
- O2 [ ] Untested and falsifiable: at exactly `N = 255` the count wraps to **0**, `is_hanging()`
        returns false, and the duplicate becomes a free unconstrained dof — worse than pinning.
        This deck has no CutSet with N in 233..424, so the cliff position is read off the code,
        not measured. **Do not describe the safe range as "N <= 232".** The cliff is 255 cuts.

## Gate G1 — run and measured 2026-08-26 (PASS)

Serial gantry, rebuilt after `d0f4db7d`, with `gantry.bfm` regenerated (47.8 MB against the
poisoned 46.7 MB). Measured on `hphi_results.e-s.00016`, t = 1.6 s, I = 55.8992 A. Note this is
a **different operating point** from the failing run (t = 1.2 s, 42.0112 A), so field magnitudes
are not comparable across the two; every structural claim below is.

| check | before | after |
|---|---|---|
| duplicate nodes showing the `(N+1) mod 256` wrap | 302 | **0** |
| cut pairs with `phi_dup - phi_org = +N * I` exactly | 1140 right / 231 wrong | **1371 / 1371** |
| same-kind cut pairs with B continuous to < 1e-6 T | 5 broken, up to 1.14 T | **1366 / 1366** |
| dead islands in the mesh (|B| max < 1e-3, non-tape, >= 10 nodes) | iron arm 375 + fragment 8 + 2 air | **none** |
| dead upper yoke arm, |B| mean / max | 0.0000 / 0.0003 T, phi = const | **0.4618 / 0.6117 T, phi varies** |
| 183 trunk nodes at N = 464, phi / I | 209.000 constant | **0.491 .. 159.655, position-dependent** |

No NaN/Inf; global |B| max 1.83 T. The ferro-air interface pairs still jump, which is correct —
that is a material interface, not a cut, and the check above filters same-kind pairs only.

**O1 is closed by the same fix.** The five iron pairs at y = 0.0520 that carried an intact
N = 232 jump but `|dB|` up to 1.14 T are now continuous. They were collateral of the neighboring
corrupted trunk, not a second defect — the "cut-terminus field recovery" reading recorded at
medium confidence was wrong, and the honest note is that it was never needed.

**O2 is closed by the guard, not by a test.** At `uint16_t` the old `N = 255 -> count 0 ->
free dof` cliff is unreachable, and the new `BELFEM_ERROR` aborts loudly instead of wrapping if
anything ever approaches the new width.

**Still owed and recorded in the register:** G2 (gantry np = 10 — the distributed path ships
whatever the mesh holds) and G3 (`make check-fast`). DR-108 was struck anyway under the
DR-42/49 exception, on Christian's call; struck is not verified.

**D4 survived and became DR-110:** `add_source`'s guard bounds the counter, not the allocation.

## Christian's mitigation idea (recorded 2026-08-26 — still NOT to be implemented)

An `input.conf` flag telling the cohomology algorithm whether a cut **may** run through a
ferromagnetic domain; when forbidden, treat the yoke like a coil block during cut generation so
its cut is always condensed out and the physical cuts route around the iron. Christian solved
this geometry manually with no iron-piercing cuts two years ago, so an iron-free routing exists
here; a flag rather than a hard rule because some geometries may have no iron-free option.

**Now known:** this would not have fixed DR-108. 205 of the 302 corrupted nodes are in the air,
and the trunk begins on the outer air boundary. The idea stands on its own merits (support
minimization, cut quality in a nonlinear mu(H) region) and still triggers the two-artifact
input-contract rule when implemented — but it must not be sold as the fix for the dead arm.

## Reproducer

Run in the gantry output directory. Needs `scipy` only; no build. Prints one row per N group
above the cliff; any row where `phi/I` differs from `(N+1) mod 256` means the wrap is gone.

```python
from scipy.io import netcdf_file
import numpy as np, collections
f = netcdf_file('hphi_results.e-s.00012', 'r', mmap=False); v = f.variables
nm = lambda n: [b''.join(r).decode().strip().replace('\x00','') for r in v[n][:]]
nb = len(v['eb_prop1'][:]); conn = {b: v['connect%d'%b][:]-1 for b in range(1, nb+1)}
phi = v['vals_nod_var%d' % (nm('name_nod_var').index('phi')+1)][-1]
I = v['vals_glo_var'][-1][0]                      # all tape currents are equal in this deck
off = {}; o = 0
for b in range(1, nb+1): off[b] = o; o += conn[b].shape[0]
g2b = np.zeros(o, int); g2l = np.zeros(o, int)
for b in range(1, nb+1):
    n = conn[b].shape[0]; g2b[off[b]:off[b]+n] = b; g2l[off[b]:off[b]+n] = np.arange(n)
sides = {3: [(0,1),(1,2),(2,0)], 4: [(0,1),(1,2),(2,3),(3,0)]}
k = collections.Counter()
for i, s in enumerate(nm('ss_names')):
    if not s.startswith('cut_'): continue
    seen = set()
    for e, sd in zip(v['elem_ss%d'%(i+1)][:]-1, v['side_ss%d'%(i+1)][:]-1):
        c = conn[g2b[e]][g2l[e]]
        for q in sides[len(c)][sd]: seen.add(int(c[q]))
    for n in seen: k[n] += 1
by = collections.defaultdict(list)
for n, kk in k.items(): by[kk].append(n)
bad = 0
for N in sorted(by):
    if N + 1 <= 255: continue
    p = phi[by[N]] / I; pred = (N + 1) % 256
    ok = abs(p.mean() - pred) < 1e-3 and p.ptp() < 1e-3
    bad += len(by[N])
    print(f"N={N:4d} nodes={len(by[N]):4d} phi/I={p.min():9.4f}..{p.max():9.4f} "
          f"pred={pred:4d} {'WRAPPED' if ok else 'not the wrap'}")
print("corrupted duplicate nodes:", bad, "of", len(k))
```

The companion checks (B continuity across every cut pair, the iron component census, the
overflow-trunk overlay plot) are in the session devlog
`devlog/dl20260826_dr108_uint8_source_overflow.md`.

## Process

Investigation and any fix ride the standing plan+audit -> code+audit round (Codex + Grok).
The root-cause round is done and recorded in the devlog; the fix round has not started.

## Fix as landed (2026-08-26)

`src/mesh/cl_Mesh_Basis.{hpp,cpp}`, +58/-1:

- `mNumberOfSources` widened `uint8_t` -> `uint16_t`. Both auditors confirmed
  `sizeof(Basis)` is unchanged (the member sits in the 8-byte slot before `Basis** mSources`
  on LP64). `uint16_t` chosen over `uint` deliberately: it keeps the allocation guard live
  and matches the `mFacetCounter` precedent, which had already been widened *because of
  cohomology cuts*.
- `BELFEM_ERROR` width guards before every narrowing write: the **three** `set_sources`
  overloads (the plan's "four" was a miscount both auditors caught), plus an allocation-size
  guard in `allocate_source_container` and a no-new-member wrap guard in `add_source`.
- `mNumberOfDofs` left `uint8_t` but guarded in `increment_dof_counter` and `insert_dof`.
  This was upgraded from judgment call to REQUIRED in audit: a wrapped dof counter makes
  `allocate_dof_container` skip its `malloc`, and `insert_dof` then writes through `mDofs`,
  which has no `nullptr` initializer.
- All limits via `std::numeric_limits<decltype(member)>::max()`; no hardcoded cliffs.
- No `.bfm` load guard, by design: a wrapped file is *below* any width check, so it cannot
  be detected there. Deleting stale `.bfm` is a gate condition instead.

### Gate results

- **G1 (serial gantry, fresh mesh, stale bfm deleted):** the brief's census reports
  **0 corrupted nodes of 2334**. Every N-group above the old cliff now reads "not the wrap".
  The 183 trunk nodes at N = 464, formerly pinned to exactly 209*I regardless of position,
  now span 0.4621..159.6016 *I with real spatial structure.
- **G2 (np = 10, fresh build):** zero "Could not find master" aborts and the residual
  sequence is **identical to serial** (0.00 / -71.33, -4.77 / -79.07, -4.77 / -77.34,
  -4.77 / -75.62 dB), so the distributed path carries the full 465-entry lists.
- **G3:** `make check-fast` 9/9.

### Unplanned corroboration

With the fix, each timestep converges in **2 Picard iterations** to -71..-85 dB. Before it,
every step needed Picard 1-2 plus Newton 3-6. The wrapped constraint had been making the
solver fight an inconsistent system; the conditioning improvement is independent evidence
that the constraint is now right.

Round record: `tmp/ai_exchange/dr108_uint8_source_overflow_fix.md`.
