# Jury Wrappers: Concurrent Temp-File Collisions and bash 3.2 Portability

**Date:** 2026-11-08
**Purpose:** Fix four defects in the jury tooling. (1) Two `cross_review.sh` rounds running at once
collide on one literal temp filename. (2) Two rounds using the same slug delete each other's
per-seat receipts. (3) The `--quick` path uses bash-4 syntax. (4) `ask_copilot.sh`, and
`ask_grok.sh` when run directly under `/bin/bash`, abort under macOS's bash 3.2 when an argument
array is empty. After the fix, these three scripts run under any bash ≥ 3.2 and no round shares a
temp path with another.
**Module:** meta (tooling: `scripts/`, `.claude/scripts/`)
**AIs involved:** Claude (plan); Codex `gpt-5.6-terra`/high + Grok `grok-4.7`/high (plan jury round 1)
**Status:** ✅ COMPLETE (2026-11-08). R1–R5 landed, plus three amendments from the code juries, each approved by Christian: RUNTAG lowercased (Grok, round 1: on a case-sensitive filesystem the receipts were missed); the jury merge appends under a `mkdir` lock taken before `collect`, failing closed on a 120 s timeout (Codex + Grok, rounds 1–2); `SLUG` canonicalized as the wrappers do (all three voices, pre-existing). Gates executed: G1–G6 on the real tree, and on a case-sensitive APFS image same-slug concurrency, a non-canonical slug, a stale lock and SIGTERM during the lock wait, with positive controls for the mktemp collision and the RUNTAG case defect. Filed: DR-204 (`--quick` needs `flock(1)`, absent on stock macOS)

> **Scope guards:**
> - In scope: temp and receipt paths in `scripts/cross_review.sh`; bash-4-only syntax in that file;
>   empty-array expansion under `set -u` in `.claude/scripts/ask_copilot.sh` and `.claude/scripts/ask_grok.sh`.
> - Out of scope: prompts, depth handling, exchange-file format, the auditors' CLI flags, the shared
>   `_cross_review.log` (append-only diagnostics), and same-slug appends to the final exchange file
>   (one slug per task is the protocol's rule).

---

## 1. Current Behaviour and How It Fails

| Failure | Mechanism | Evidence |
|---|---|---|
| A second concurrent round fails at start: `mktemp: mkstemp failed on …/crossrev_subject_XXXXXX.md: File exists` | `cross_review.sh:95-96` pass `…_XXXXXX.md` to `mktemp`. BSD `mktemp` replaces only **trailing** X's, so every round on this Mac uses the same two paths. A round killed by SIGKILL also leaves them behind (the trap at `:97` is EXIT-only), and the next run fails even when it is alone. GNU `mktemp` probably rejects the template outright (Grok; unverified, no GNU mktemp here) | observed 2026-11-07 20:39 (DR-202 code round against the DR-203 round) |
| Two rounds with the same slug delete each other's seat receipts | jury receipts are `${SLUG}_jury_{codex,grok}.md` (`:184`), `rm -f`'d at start (`:185`) and again by `collect` (`:177`); `SLUG` defaults to `scratch` (`:40`). The `--quick` receipt `QT` is named by commit SHA alone and removed at start (`:205-207`) | source trace (Codex) |
| `--quick` under bash 3.2: `bad substitution` | `${AUDITOR^^}` (`:220`) is bash ≥ 4.0. `bash -n` passes the file, because the error is raised only when the line runs. `${SHA: -1}` (`:201`) is fine on 3.2 | executed 2026-11-08: `/bin/bash -c 'A=codex; echo ${A^^}'` → bad substitution; `${S: -1}` → `f`; `/bin/bash -n scripts/cross_review.sh` → rc 0 |
| `ask_copilot.sh` aborts with `EFFORT_ARGS[@]: unbound variable` when `COPILOT_MODEL=auto` | `#!/bin/bash` (`:1`) is 3.2.57 on macOS; `set -u` (`:36`); `EFFORT_ARGS=()` (`:88`); `"${EFFORT_ARGS[@]}"` (`:169`). Bash < 4.4 treats an empty array's `[@]` as unset under `set -u` | `/bin/bash -uc 'a=(); echo "${a[@]}"'` → unbound variable, rc 127 |
| Same latent defect in `ask_grok.sh`, masked by its shebang | `PERM_ARGS=()` (`:268`) stays empty unless `GROK_PERMISSION_MODE` is set (`:189` defaults it empty); `tDisallow=()` (`:415`) stays empty when the caller sets `GROK_DISALLOWED_TOOLS=` explicitly (`:139-142` defaults it only when unset). Both are expanded bare at `:424,:430`. `#!/usr/bin/env bash` resolves to bash 5.3 here. The comment at `:264-265` calls `"${arr[@]}"` "the set -u-safe expansion", which holds only for bash ≥ 4.4 | `/usr/bin/env bash --version` → 5.3.0 |

**Bottom line:** one temp template that does not randomize, receipt names that are not unique per
round, one bash-4 expansion, and an empty-array idiom that needs a newer bash than macOS ships.
(high confidence)

## 2. Mechanism

- **Per-round work directory.** `WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/cross_review.XXXXXX")"`,
  then the trap `trap 'rm -rf "$WORKDIR"' EXIT`, installed only after `mktemp -d` has succeeded
  (`set -e` at `:22` exits on a failed assignment before any trap exists). Never write
  `"$WORKDIR/"`. `SUBJECT_FILE="$WORKDIR/subject.md"`, `PROMPT_FILE="$WORKDIR/prompt.md"`. The
  template ends in X's, so BSD and GNU `mktemp` both randomize it; do not shorten it to a bare
  `mktemp -d`, which is GNU-only. Nothing depends on the `.md` names (the files are only written,
  `cat`ed and fed on stdin, `:111-160`).
- **Unique receipt names.** `RUNTAG="${WORKDIR##*.}"` (the random suffix). The jury receipts become
  `${SLUG}_jury_codex_${RUNTAG}.md` and `${SLUG}_jury_grok_${RUNTAG}.md`, and the quick receipt
  slug becomes `autoreview_${SHA:0:12}_${RUNTAG}`. The final exchange file and its `flock` are
  unchanged.
- **bash-4 expansion.** `${AUDITOR^^}` → `$(printf '%s' "$AUDITOR" | tr '[:lower:]' '[:upper:]')`.
- **Empty arrays.** Every array that can be empty is expanded as `${arr[@]+"${arr[@]}"}`, outer
  expansion unquoted. Under bash 3.2 an empty `[@]` counts as unset, so `+` substitutes nothing and
  the inner expansion never runs. Under ≥ 4.4 the inner, quoted expansion runs, and an empty array
  is legal there. Applied at `ask_copilot.sh:169` (`EFFORT_ARGS`) and `ask_grok.sh:424,430`
  (`tDisallow`, `PERM_ARGS`); the comment at `:264-265` is corrected. `READ_TOOLS`
  (`ask_copilot.sh:174`) and `ask_grok.sh`'s `EFFORT_ARGS` (`:266`, `:431`) are never empty and stay.
- `ask_copilot.sh`'s shebang moves to `#!/usr/bin/env bash`, the form the other three scripts use.
  That is for consistency; the idiom alone is the fix, and the script also runs under `/bin/bash`.

Rejected: a shebang that points at `/opt/bash/bin/bash`, or a version check that refuses to run.
The first hardcodes a path that exists only on this machine (Christian's standing preference keeps
machine-specific tool paths out of tracked files). The second still fails on stock macOS, while
the idiom costs nothing.

## 3. Steps

- [x] **R1** `scripts/cross_review.sh:95-97`: `WORKDIR` via `mktemp -d`, then `SUBJECT_FILE` and
  `PROMPT_FILE` inside it, then the trap (`rm -rf`). `RUNTAG` follows.
- [x] **R2** `scripts/cross_review.sh:184` and `:205`: receipt names carry `${RUNTAG}`.
  (after: R1)
- [x] **R3** `scripts/cross_review.sh:220`: `${AUDITOR^^}` → `tr`.
- [x] **R4** `.claude/scripts/ask_copilot.sh`: shebang `#!/usr/bin/env bash`; `:169` to the `+` idiom.
- [x] **R5** `.claude/scripts/ask_grok.sh:424,430`: the `+` idiom; comment `:264-265` corrected.
- [x] **R6** gate (§4), run by Claude in the scratchpad (no build, no solver). (after: R1–R5)
- [x] **R7** (amendment, code round 1, approved 2026-11-08) `RUNTAG` lowercased; jury merge under a `mkdir` lock;
  the GNU-mktemp sentence softened. Gate: G2 on a case-sensitive APFS sparse image, with a positive control.
- [x] **R8** (amendment, code round 2, approved 2026-11-08) lock taken before `collect`, fail closed on timeout
  (exit 3, receipts kept); `SLUG` canonicalized with the wrappers' `sanitize_slug`. Gates: same-slug
  concurrency, `Review-Jury` slug, stale lock, SIGTERM during the wait.

## 4. Gate

Stub CLIs are scratchpad scripts selected with `CODEX_BIN` (`ask_codex.sh:38-39`), `GROK_BIN`
(`ask_grok.sh:95`) and `COPILOT_BIN` (`ask_copilot.sh:39-40`). Each records its argv one argument
per line and returns a canned audit of more than 200 characters starting with `## Audit & Verdict`.
The Codex stub writes it to the `--output-last-message` file; the Grok stub emits the JSON
`parse_grok_json` expects. Every gate uses a throwaway slug (`AI_EXCHANGE_SLUG=wrapper_gate_<n>`),
and its exchange files are deleted afterwards.

| ID | Check | Pass |
|---|---|---|
| G1 | `bash -n` and `/bin/bash -n` on all three scripts | rc 0 under 5.3 and 3.2 |
| G2 | the **real** `cross_review.sh --jury` started twice concurrently with stub CLIs, the **same** slug and no `TMPDIR` | both runs exit 0; the exchange file holds four stub entries (two per run); no `cross_review.*` directory and no `*_jury_*` receipt remains afterwards |
| G3 | `/bin/bash scripts/cross_review.sh --quick HEAD` with stub CLIs | exit 0, no `bad substitution`; the auditor heading is upper case |
| G4a | `/bin/bash .claude/scripts/ask_copilot.sh` with `COPILOT_MODEL=auto` (empty array) | no `unbound variable`; recorded argv has no `--reasoning-effort` |
| G4b | the same with `COPILOT_MODEL=gpt-5.6-luna COPILOT_EFFORT=high` (non-empty array) | `--reasoning-effort` and `high` are two separate argv entries |
| G5a | `/bin/bash .claude/scripts/ask_grok.sh` with `GROK_PERMISSION_MODE` unset and `GROK_DISALLOWED_TOOLS=` (both arrays empty) | no `unbound variable`; argv has neither `--permission-mode` nor `--disallowed-tools` |
| G5b | the same with `GROK_PERMISSION_MODE='read only'` (a value containing a space; whether the real CLI accepts it does not matter to a stub) and the default disallowed tools | `--permission-mode` and `read only` are two separate argv entries, as are `--disallowed-tools` and `search_tool,use_tool` |
| G6 | `python3 scripts/check_practice_drift.py`, then separately `--self-test` | no new finding; self-test green. Housekeeping only: neither command looks at the array fix |

## 5. Open Questions

- **O1** Gate seams. RESOLVED 2026-11-08: `SCRIPTS` is hardcoded at `cross_review.sh:27`, but all
  three wrappers take a CLI path from the environment, so the real script runs end to end against
  stubs (Grok, plan round 1). The first version of this plan wrongly fell back to a probe copy of
  `:95-97`.

## 6. Audit Trail

- `tmp/ai_exchange/review_jury_wrapper_portability.md`: plan jury round 1 (Codex `gpt-5.6-terra`/high,
  Grok `grok-4.7`/high), with verification and reconciliation. v2 folds in:
  - `${AUDITOR^^}` (both auditors)
  - unique receipt names (Codex)
  - the JSON Grok stub (both)
  - empty and non-empty argv gates, including `tDisallow` (Grok)
  - the real-script gate (Grok)
  - the trap order (Grok)

  Refuted:
  - my "trap deletes the other round's file" inference (Grok)
  - "`${SHA: -1}` needs bash 4.2" (executed on 3.2)

  The round's parent process was stopped after dispatch, so its merge step did not run; both seat
  files were complete and were merged by hand.
