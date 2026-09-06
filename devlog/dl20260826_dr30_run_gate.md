# DR-30 run gate: the guard is live but not triggered, and the nominated deck could never have shown it

**Date:** 2026-08-26
**Module:** fem/maxwell (postprocessing), fem/kernel (postprocessor element selection — read only)
**Follows:** [dl20260826_dr30_stale_thermal_peer.md](dl20260826_dr30_stale_thermal_peer.md), which landed the code half

## What was run

`hphiTrun` rebuilt from HEAD `2320d938` plus the uncommitted DR-30 diff, carrying a temporary
counter probe in `compute_superconductor` / `compute_superconductor_ts`. Once per postprocessing
pass per rank it reported how many times the thermal peer was read, and how many of those were
declined because the peer's linked element id differed from the maxwell calculator's current one.

Two decks, each serial and on 2 ranks, cold started with no `memdump.hdf5`:

- **corc**, thermal-coupled (`cmake-build-debug/corc2/input.conf`, 200 A ramp, bdf5, fully coupled
  thermal), shortened to `simulation time : 0.06 s` / `save every : 0.02 s`. This is the deck the
  register row nominated.
- **dr30_bulk**, derived from `tapestack3d`: volume conductor blocks 2:8 repointed from `solder`
  to a new power-law `hts` material, `simulation time : 0.02 s` / `save every : 5 ms`. 43 276
  nodes, 247 755 elements, 8 blocks. See "the substitute deck" below for why this exists and why
  it is deliberately unphysical.

## The measurements

Thermal reads per postprocessing pass, and the stale-element count:

| deck | path | serial | rank 0 | rank 1 | rank sum | aura overhead | stale | null |
|---|---|---|---|---|---|---|---|---|
| corc | TS | 380 160 | 190 080 | 190 080 | 380 160 | **0** | 0 | 0 |
| dr30_bulk | TS | 468 512 | 234 256 | 234 256 | 468 512 | **0** | 0 | 0 |
| dr30_bulk | bulk | 1 334 102 | 1 147 784 | 331 034 | 1 478 818 | **+144 716 (10.8 %)** | 0 | 0 |

Counts were stable across passes — every pass added exactly the same increment (corc 2-rank ran
three passes, dr30_bulk 2-rank ran four), so these are structural set sizes, not transients. Both
2-rank runs exited 0.

## What the numbers establish

**The aura pass is real on the bulk path and absent on the thin-shell path — measured, not
inferred.** `Postprocessor::select_elements_and_owned_nodes` adds a second, aura pass for
non-thin-shell blocks and skips it for thin shells (`cl_FEM_Postprocessor.cpp:258`,
`if ( a == 1 && tIsThinShell ) continue ;`). The consequence is visible directly in the table: the
thin-shell rank counts sum to exactly the serial count in both decks (zero redundancy, so no
element is recovered twice), while the bulk counts exceed it by 144 716 reads. Those extra reads
are aura elements being recovered on both ranks. The audit's central reachability claim is
therefore **confirmed**.

**The stale condition itself did not occur.** Over 1.48 M reads per pass, across four passes, on a
partition split roughly 78/22, the thermal peer's linked element matched the maxwell calculator's
every single time — including for all 144 716 aura reads. `element_exists()` never returned false,
so the silent-keep in `link_element_maxwell_thermal` never produced a stale temperature and the
guard never had anything to catch. On this deck the thermal group's element set covers the maxwell
group's aura.

**No regression.** corc terminal voltages agree between 1 and 2 ranks to 1e-8 .. 1e-6 relative on
all six tapes at all three save points — Newton-path level, not a discrepancy. The first residuals
of the serial and 2-rank bulk runs are identical to the printed precision (−105.92 dB thermal).

## Why the nominated gate could never have worked

The register row asked for "2-rank thermal-coupled corc". That deck **structurally cannot reach the
guarded code**, on either path:

- `MaxwellFactory` constructs a `SuperConductor` postprocessor only for `s == 0` — *volume*
  conductor blocks — whose material has `jc` (`cl_MaxwellFactory.cpp:2318`). corc declares only
  `air` (blocks 1:3) and `thinshell : tape`, so it has no volume conductor block and
  `compute_superconductor` is never constructed. No `bulk SC` probe line was ever printed there.
- The thin-shell variant is constructed, and ran 190 080 times per rank per pass, but its blocks
  are excluded from the aura pass, so its thermal peer is always current.

This is not bad luck in partitioning: no deck configuration of corc could trip it. Recording that
gate as passed would have certified nothing.

Nor is any *shipped* deck able to do better. Every thermal-coupled deck in the tree is thin-shell
tape only (`corc*`, `sidecoating`, `sidecoatings`) except the `tapestack3d` family, whose bulk
blocks 2:8 are `solder` (Pb38Sn62, no `jc`) and so build a `Conductor` postprocessor, which never
touches the thermal peer at all. The one deck with a bulk ybco block, `bigtape`, has no thermal
section and last ran 2026-05-07.

## The substitute deck

`dr30_bulk` exists only to make the factory build a bulk `SuperConductor` postprocessor, and it is
**not physical**: it makes the inter-tape solder superconducting. Recipe, from `tapestack3d`:

1. add a material `hts` — `builtin : ybco`, `jc : 1e10 A/m^2`, `n : 20`, `ec : 1e-4 V/m`,
   `resistivity type : power-law`. Power law rather than the `sp-ap.hdf5` table so the deck carries
   no measured-data dependency.
2. in `topology`, change the `conductor` block's `material : solder` to `material : hts`.
3. shorten `simulation time` and `save every` so every step postprocesses.

Nothing else changes. The magnetic system is badly conditioned as a result (6.1e16 reported), which
is expected from the conductivity contrast and did not prevent convergence. The deck is described
here rather than checked in, so no fictional physics ships in `examples/`.

## Evidence class

Executed gate, not review: four runs, two decks, serial and 2-rank, cold start, exit 0. The probe
was a temporary diagnostic and has been removed; `cl_MaxwellPostprocessor.cpp` is back to exactly
HEAD plus the DR-30 diff (28 insertions, 16 deletions), re-checked with
`g++ -fsyntax-only -std=gnu++17` under the target's own flags.

One methodological note worth keeping: the first serial control was invalid because the run
directory was seeded with a `memdump.hdf5` copied from the 2-rank run, so it warm started at BDF2.
`examples/scripts/Allrun` prints a warning for exactly this case; calling `mpirun` directly
bypasses it. The control was rerun cold.

## Residuals, unchanged by this run

1. T is not in `mPostprocessorSourceFields` — postproc correctness still rests on a preceding
   thermal solve+distribute.
2. `gTbulk` on maxwell-only aura cells would bias partition-boundary nodal J/Jc if the groups ever
   did diverge; the physics-complete fix is expanding the thermal aura, not a postproc guard.
3. Pre-existing interpolation mismatch: postproc `norm(N*q)` (no clamp) vs assembly
   `dot(Nvec,q)` + clamp.

## The silent-keep is kept, and now says so

The run showed the miss never occurring, which invites the reading that the branch is dead and
could be tightened into an assert. **Christian's ruling is that it stays**: it is the sanctioned
fallback for a maxwell element with no thermal counterpart, and coupled models added later may
need exactly that behaviour. Latent is not dead.

What was missing was the intent, not the code. The comment at `link_element_maxwell_thermal`
explained the hazard but read as a description of an accident, so a later reader — human or model
— could reasonably have "fixed" it. It now states that keeping the peer is deliberate, that the
obligation sits on the consumer rather than the linker, and that the consumer discharging it is
the superconductor postproc's element-id comparison. The measured result is recorded there too, so
the next person to find the branch cold knows it is latent rather than unreachable.
