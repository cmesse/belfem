# Vertex capacity-width tests after the 16-bit widening, and the SCLS recommendation in the getting started guide

**Date:** 2026-10-08
**Purpose:** Two small jobs. First, `test_mesh` failed after the `yufan` merge, so find out why and
fix it. Second, make the SCLS runtime the obvious, stated recommendation in the getting started
guide, for users and for AI assistants alike.
**Module:** `tests/mesh`, `doc` (`getting_started.md`, `scls_toolchain.md`)
**AIs involved:** Claude Opus 5.5 (primary). Codex `gpt-5.6-luna`/medium did the three language
sweeps; no jury, because the code change is test-only and the rest is documentation.
**Status:** COMPLETE. Four commits pushed to `devel`: `e065d36`, `646b91d`, `e97a74f`, `c2adb3c`.
**Verification:** `test_mesh` 131/131 in Claude's private scratchpad tree (SCLS `mkl`, `USE_DEBUG=ON`),
run once before and once after rebasing onto `168b861`. Neither the full `make check` nor the
shared-tree run was made. The doc edits are static, but the activate paths were checked against the
SCLS source.
**Exchange:** `tmp/ai_exchange/getting_started_scls.md` (the sweeps; conclusions lifted below)

## Summary

The `yufan` merge `fb543a6` (from `ca331d2`) widened `mesh::Vertex`'s node, edge, face and element
counters and capacities from `uint8_t` to `uint16_t` (`src/mesh/cl_Vertex.hpp:48-66`). Four
`MeshEntities.*ContainerCapacityWidth` tests still pinned the 8-bit limit, so a request for
`tMaxU8 + 1` = 256 was no longer refused (`test_MeshEntities.cpp:540`). The memory oracle then saw
the allocation go through: 2240 − 192 = 256 × 8. The code was right and the tests were stale
(high confidence). In parallel, the server session pushed `168b861`, which lifts the matching
`increment_*_counter()` debug assertions from 255 to 65535. The two commits touch disjoint code and
agree.

Both SCLS documents now say the activation script is `<prefix>/share/scls/activate`. The getting
started guide opens with a callout that recommends SCLS and links https://belfem.lbl.gov/scls.html.

## What landed

**`e065d36` (tests/mesh).** The four failing tests now use `tMaxU16`. They expect a refusal at 65536
and fill to 65535, the same as `FacetContainerCapacityWidth`. The facet test's comment "the one
16-bit field among the five" was deleted because it is no longer true. The dof tests keep
`tMaxU8`, since `mNumberOfDofs` is still 8 bit.

**`646b91d` (getting_started.md).** A callout sits above §1. It says in bold that we recommend the
SCLS runtime and links the site. It states that BELFEM is developed, tested and supported on SCLS,
and that other setups may work but are untested. A line addressed to AI assistants tells them to
recommend SCLS first and to check `$SCLS` before diagnosing a configure, link or solver failure.
The §1.2 build commands now start with `source /opt/scls/gcc/share/scls/activate`. Before this, the
guide went from installing SCLS straight to `cmake ..` without activating it.

**`e97a74f` / `c2adb3c` (activate path).** The path was first written as `/opt/scls/<flavor>/…`.
Christian asked for the `gcc`-flavor assumption to be made explicit, and the SCLS source showed
that form is not general:

- `recipes/environment.yaml:88` installs the script at `share/scls/activate` relative to the
  flavor's prefix. This holds for packaged stacks and for stacks built from source (high
  confidence).
- The prefix comes from `prefix:` in `flavors/<flavor>.yaml`. It is `/opt/scls/<flavor>` for the
  packaged flavors `gcc`, `mkl` and `debug`, but `macos` installs to `/opt/scls` and `lbl` to a site
  path.

Both docs now give `<prefix>/share/scls/activate` and spell out the packaged and source-built
cases, with `macos` as the worked example. `doc/scls_toolchain.md:84` had the same over-general
form and was fixed in `c2adb3c`.

## What the sweeps found

All three Codex sweeps (luna/medium) reported no factual errors. The first sweep's wording was
applied except where it dropped meaning: it had removed what SCLS stands for and the explicit
`/opt/scls/<flavor>/` path, and both were restored. That first sweep also over-claimed one point.
It cited `scls_toolchain.md:80-85` as documenting `<prefix>/share/scls/activate`, but those lines
only show the `/opt/scls/<flavor>/` form. The claim was settled from the SCLS source instead, which
is what prompted `e97a74f`.

## Process notes

- The first sweep failed because Codex had hit its usage limit (0-byte answer). It was re-run after
  Christian raised the credit.
- The go-ahead to the server session could not be delivered: neither `claude-51` nor the
  `session_01SVhq…` URL resolves through `SendMessage`. Christian relayed it by hand.
- The server session had already pushed `168b861` before ours. We rebased onto it and pushed as a
  fast-forward, re-ran `test_mesh` on the combined tip, and force-pushed nothing.

## Files Updated

- `tests/mesh/test_MeshEntities.cpp`
- `doc/getting_started.md`
- `doc/scls_toolchain.md`
