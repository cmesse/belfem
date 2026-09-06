# Devlog 2026-08-05 — Cross-Review Tooling

**Date:** 2026-08-05
**Topic:** Freeze the manual three-AI review protocol into `/cross-review`, an opt-in
post-commit light review, and a public `METHODOLOGY.md` draft. Plus one real jury round on
the gesvd/posv/Spline WIP diff as the acceptance test.
**AIs involved:** Claude (implementation + primary review), Codex + Grok (jury on the WIP
diff; Codex as the quick-mode demo auditor)
**Claude Confidence:** high — every mode exercised against stubs and once for real
**Literature References:** none (tooling errand)

## Summary

The five-month manual protocol (pre-registration → headless read-only audits →
citation verification → reconciliation) is now tooling. No protocol logic was changed;
the scripts orchestrate what sessions already did by hand, reusing the existing
`ask_codex.sh` / `ask_grok.sh` wrappers unchanged. Plan: `todo/cross_review_tooling.md`.

## What landed (new files, uncommitted)

- `scripts/cross_review.sh` (~170 lines) — `--jury` (parallel, blind, per-auditor temp
  slugs merged in fixed order), `--relay` (sequential on the main slug), `--quick <sha>`
  (one auditor by commit-hash parity, flock-guarded append to
  `tmp/ai_exchange/autoreview.md`, P0 flag file). Diff capped at 400 lines
  (`DIFF_CAP`) with an explicit truncation marker + file list. Mode-split priming:
  the "find what the prior reviewers missed" clause is relay-only; `--quick` mandates
  P0/P1/P2 labels. Skips merge commits and mid-rebase states.
- `scripts/install_autoreview_hook.sh` — installs the post-commit hook into
  `git-common-dir/hooks` (one install covers all worktrees); refuses to overwrite a
  foreign hook; `--uninstall`. Hook is a no-op unless `BELFEM_AUTOREVIEW=1`, resolves
  the commit sha at commit time (no race with the next commit), and backgrounds the
  review via `nohup nice -n 19` — measured commit overhead 0.01–0.11 s.
- `scripts/review_status.sh` — "commits since last auto-review" + open P0 flags.
- `.claude/commands/cross-review.md` — the slash command; Claude does the judgment
  steps (pre-registration frozen BEFORE dispatch, per-citation verification labeled
  CONFIRMED/REFUTED/UNVERIFIABLE with quoted evidence, reconciliation table
  `finding | raised by | verdict | severity | agreement`, single-raiser findings routed
  to Christian, never auto-fix).
- `METHODOLOGY.md` — DRAFT at repo root, NOT committed; mapping table + deviations list
  await Christian's sign-off.
- `CLAUDE.md` — usage-lines section under AI Cooperation.

## Verification

- All modes smoke-tested with stub auditors in a scratch repo (jury merge order, relay
  sequencing, quick flock + flag, merge/rebase skip, file-target mode). `shellcheck` is
  not installed on this host; `bash -n` + stub tests + one real run stand in.
- **Acceptance 1 (jury, real):** full round on the WIP diff, record in
  `tmp/ai_exchange/review_gesvd_spline_wip.md`. Outcome: no P0; 3/3 clean on gesvd
  packing arithmetic, posv ldb, spline formulas; two P1 coverage gaps (gesvd
  provided-buffer/reuse branch untested — `test_lapack.cpp:526,570` always pass an empty
  `tWork`; new spline `aCol` overloads untested — `tests/math/test_Spline.cpp` is
  single-arg only), both single-raiser (Grok), Claude-verified, routed to Christian.
  P2: "no internal allocation" docstring overstatement (`fn_gesvd.hpp:322`), the
  accepted real-hosts-complex Work aliasing convention (geev/gees precedent), mixed WIP
  tree (tracked `.claude/ai_exchange_pos.txt` cursor churn riding with LAPACK code).
- **Acceptance 2 (hook, throwaway worktree):** default-OFF commit is a no-op; opt-in
  commit returned in 0.02 s; Codex findings landed in the worktree's `autoreview.md`
  under `## commit <sha> <ts> auditor=codex`; `review_status.sh` counted correctly.

## Defects found and fixed during acceptance

- **D1 (fixed 2026-08-05, Claude):** the P0-flag grep matched the token inside
  negations — Codex's clean verdict "**No** P0/P1 defect is introduced" raised the
  flag. Fix: negation filter on the grep + a priming line telling quick-mode auditors
  not to write the token P0 except for actual P0 findings. Both polarities unit-tested.
- **Operational caveat (not a defect):** the hook's `[ -x $TOP/scripts/cross_review.sh ]`
  guard means worktrees whose checkout predates the scripts silently skip the review
  (commit never blocked). Self-resolving once the scripts are committed.

## Round 2 (same day): METHODOLOGY sign-off, environment pin, commit slicing

- METHODOLOGY.md approved with edits: session-header block dropped (living root document);
  its `doc/literature_references.md` reference made real — that file now carries the
  authoritative **paper alias index** (paper0–paper9 + paperA → source txt → one-line
  citation, plus F0–F4 and the unaliased newer author-year papers). The paper9/paperA
  drift (both bind to schnaubelt2023: paper9 in Apr–May devlogs, paperA since July) is
  documented in place, not silently repaired.
- `scripts/scls_env.sh`: idempotent toolchain pin — prepends `/opt/scls/gcc/bin` and the
  generation-captured Python prefix (`/home/christian/Applications/python`) only when
  absent; guards (cmake/mpicxx under `/opt/scls/`, python3 under the pinned prefix)
  hard-fail when sourced non-interactively, warn-only interactively. Sourced first by
  `cross_review.sh` and the post-commit hook (which prints a failed pin but still exits 0
  — never blocks the commit). Layout facts recorded: SCLS root `/opt/scls` with `debug/`,
  `gcc/`, `mkl/` subtrees; only `gcc/bin` belongs on PATH; `debug/` and `mkl/` are
  CMake-discovered library stacks; Python is not part of SCLS.
- **Deviation needing Christian's confirmation:** the requested "g++ under /opt/scls/"
  guard is impossible on this host — SCLS ships no compiler binary; `mpicxx -show` proves
  it wraps the bare system `g++` with SCLS include/lib flags. Implemented as an
  existence-only check for g++, with the evidence quoted in the script comment.

- **Bare-environment acceptance (final):** in a fresh worktree of the tooling commit, an
  `env -i HOME=… PATH=/usr/bin:/bin` commit took 0.01 s; the background `--quick` job
  survived the launching shell's exit (reparented to `systemd --user`), `scls_env.sh`
  pinned the toolchain, and Grok caught a deliberately planted `strlen`-without-`+1`
  off-by-one, labeled it P0 — the flag file captured the genuine finding heading,
  proving the D1 negation filter passes real P0s through. `review_status.sh` counted
  correctly. Both hook demo paths (clean verdict, real P0) are now exercised.

## Open items

- Christian's sign-off on `METHODOLOGY.md` (mapping table + deviations list) before it
  is committed; commit of the tooling files themselves.
- Codex prose pass on `todo/cross_review_tooling.md` (deferred until after the §1
  summary was approved; approval arrived — pass still to run).
- The two P1 coverage gaps from the jury round belong to the LAPACK/spline WIP, not
  this task; recorded here and in the exchange file for Christian.

## Files Updated

- scripts/cross_review.sh, scripts/install_autoreview_hook.sh,
  scripts/review_status.sh, .claude/commands/cross-review.md (new)
- METHODOLOGY.md (new, draft), CLAUDE.md, todo/cross_review_tooling.md, todo/README.md
- .git/hooks/post-commit (installed, shared across worktrees)
