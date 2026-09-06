# Cross-Review Tooling — Freeze the Three-AI Protocol into Commands

**Date:** 2026-08-05
**Purpose:** Turn the existing manual Claude + Codex + Grok review protocol into a
`/cross-review` slash command, an opt-in commit-triggered light review, and a public
`METHODOLOGY.md`. Freeze rather than redesign: the scripts orchestrate exactly what the last
five months of sessions already did by hand.
**Module:** meta (`.claude/commands/`, `scripts/`, repo root)
**AIs involved:** Claude (recon + plan + implementation), Codex + Grok (dry-run subjects;
Codex prose pass on this file after the protocol summary is approved)
**Status:** ✅ COMPLETE 2026-08-05 — R1–R11 done, tooling commit landed on `sideconnectors`
(commit "cross review: freeze the three-AI protocol into tooling"). Both hook demo paths
exercised (clean verdict + real planted P0); D1 negation-filter defect found in acceptance and
fixed; METHODOLOGY.md approved (header dropped, alias index added to
doc/literature_references.md). One open confirmation for Christian: the scls_env g++ guard is
existence-only, not under-/opt/scls/ (SCLS ships no compiler; `mpicxx -show` wraps system g++
— see R10). The two P1 coverage gaps from the Acceptance-1 jury belong to the LAPACK/spline
WIP and stay with Christian. §1 was approved as written; design corrections applied: commit
sha resolved at commit time, flock on autoreview.md, mode-split priming; O1 jury=blind,
O2 P0/P1/P2 labels in --quick, O3 skip merges+rebases, O4 "commits since last auto-review".
**Closed 2026-08-09** (todo currentness sweep): every box is ticked and the whole
deliverable is committed in `bd73583b` — `scripts/cross_review.sh`,
`scripts/install_autoreview_hook.sh`, `scripts/review_status.sh`,
`.claude/commands/cross-review.md` and `METHODOLOGY.md` are all in the tree. The single
remaining line — Christian's confirmation that the `scls_env` g++ guard is
existence-only — is a one-question confirmation, not work in flight; it is carried in
`falsification_tooling.md`'s awaiting-Christian block so it does not get lost.

> **Scope guards (from the task brief):**
> - Freeze the protocol as practiced; change nothing about its logic.
> - Bash + existing CLIs only; no new dependencies, no config files/YAML — flags and env vars only.
> - Well under ~300 lines of new shell total. If a part wants to grow past that, stop and ask.
> - `tmp/ai_exchange/` stays gitignored; scripts and command file are committed.
> - Never auto-apply fixes. Auditors stay read-only.
> - `METHODOLOGY.md` mapping table and deviations list need explicit sign-off before commit.

---

## §1 — The de-facto protocol, as observed (CORRECT ME HERE)

Extracted from `doc/ai_collaboration_protocol.md`, `.claude/scripts/ask_codex.sh`,
`.claude/scripts/ask_grok.sh`, `tmp/ai_exchange/sideconnector_bfm_gaps.md` /
`greg2_top_tie_flip_codex.md`, and the audit sections of
`devlog/dl20260804_sideconnector_bfm_persistence.md` and
`dl20260730_lapack_unified_interface_review.md`.

1. **Channel.** One topic = one file, `tmp/ai_exchange/<slug>.md` (gitignored, ephemeral).
   Entries are `# AI_NAME YYYY-MM-DD HH:MM:SS TZ` + `## Section` headings +
   `Confidence: high|medium|low (~N%)`, separated by `---`. Wrappers append the voice
   header themselves; the auditor writes only the `##` body. Slug comes from
   `$AI_EXCHANGE_SLUG` (else a session tag, else `scratch`).
2. **Sequence per round.** Claude writes its own `# CLAUDE` entry first (findings/claims
   with per-claim confidence and file:line citations), then dispatches auditors, then
   writes a `# CLAUDE … ## Resolution` entry: per-claim verdicts ("confirmed / refuted /
   partial, high|medium"), what was applied, and what stays open. Devlogs show this as
   "3/3 confirmed", "two P0s caught", per-finding attribution.
3. **Auditor invocation.** Both headless, read-only, working root = repo:
   - Codex: `codex exec --sandbox read-only -C <root> --output-last-message <tmp>`
     (via `ask_codex.sh`). Preamble: read `AGENTS.md` + protocol + the exchange thread;
     independent verification over agreement; `## Audit & Verdict` + `Confidence:`;
     file:line for every claim.
   - Grok: `grok --sandbox read-only --tools read_file,grep,list_dir` + JSON output
     (via `ask_grok.sh`), with the shell-cancel footgun mitigations: no-shell tool
     allowlist, quality gate (min chars, `##` heading required, reject
     `stopReason=Cancelled`), retries with backoff, `--resume` salvage. Preamble adds the
     refutation mandate ("refutation and blind-spot detection, not agreeable synthesis").
   - Both wrappers append the result to the exchange file AND echo it to stdout.
4. **Severity & verdict vocabulary.** Findings carry P0/P1/P2 (or CRITICAL/HIGH/…);
   Claude's verification labels are per-claim confirmed/refuted/partial with evidence
   citations; false positives are retracted in place, never deleted (plan-template rule).
5. **Adjudication.** Reviewer agreement settles code-level facts; physics/design
   questions go to Christian by name (e.g. B7 deprioritized "decided, Christian").
   Literature outranks vote where routing applies.
6. **Distillation.** Exchange threads are ephemeral; conclusions are lifted into the
   devlog (and todo trackers) before sweep. Devlog cites the thread path.
7. **Known footguns encoded in the wrappers** (must survive the freeze): Grok headless
   shell-cancel; Codex background invocations need `< /dev/null` (stdin hang) and must
   not be nested in compound pipe commands.

Confidence: high on 1–3 and 6–7 (read directly from scripts/protocol/threads); medium
(~80%) on 4–5 being *invariant* practice rather than recent style — correct me.

---

## §2 — Design

### D1 — `/cross-review` (`.claude/commands/cross-review.md` + `scripts/cross_review.sh`)

Division of labor — **the script talks to auditors; Claude does the judgment steps**:

The slash command instructs Claude (the session) to:
1. Resolve target: no arg → `git diff HEAD` (tree clean → `git show` of last commit);
   arg → that file. Pick a slug `review_<topic>`.
2. **Pre-register:** write a complete `# CLAUDE` findings entry (per-claim confidence,
   file:line) to `tmp/ai_exchange/<slug>.md` BEFORE any auditor runs. Never edited afterward —
   synthesis is appended below.
3. Run `scripts/cross_review.sh --jury|--relay [target]` with `AI_EXCHANGE_SLUG=<slug>`.
4. **Verify:** re-check EVERY auditor file:line citation against source; label each
   finding CONFIRMED / REFUTED / UNVERIFIABLE with the evidence line quoted.
5. Append the reconciliation table:
   `finding | raised by | verdict | severity P0–P2 | agreement`.
   ≥2 independent raisers = high confidence; single-raiser → flagged for Christian;
   physics/design questions routed to Christian explicitly, never voted.
6. Print to the user: short summary + the P0/P1 list only. Never apply fixes.

`scripts/cross_review.sh` (target ~150 lines):
- Builds the grounded auditor prompt: target content (diff or file), minimal context
  (branch, target description), and the fixed priming "You are an independent reviewer.
  Find what the others missed. Cite file:line for every claim. State a confidence level
  per finding."
- **Diff cap:** ≤400 lines inline; beyond that, send the changed-file list + instruct the
  auditor to read the files directly (both CLIs read the repo), with an explicit
  `[diff truncated at 400 lines — file list follows]` marker.
- **Reuses `ask_codex.sh` / `ask_grok.sh` verbatim** (they already encode the footgun
  mitigations, quality gates, house-format append). Codex calls get `< /dev/null`.
- `--jury` (default): both auditors in parallel, blind — each wrapper writes to a private
  temp slug (`<slug>_jury_codex`, `<slug>_jury_grok`) so neither sees the other or
  Claude's pre-registration; on completion the driver appends both bodies, in fixed
  Codex-then-Grok order, to the main exchange file and removes the temp files. (Parallel
  appends to one file could interleave; temp slugs avoid it.)
- `--relay`: sequential on the main slug — auditor 1 sees Claude's pre-registration,
  auditor 2 additionally sees auditor 1's entry (the wrappers already instruct "read the
  thread first"). Never mixed with jury in one round.
- `--quick <commit>`: single pass for the hook — ONE auditor chosen by commit-hash parity
  (last hex digit even → Codex, odd → Grok), diff-only prompt, findings appended to
  `tmp/ai_exchange/autoreview.md` under a `## commit <hash> <timestamp>` marker. If the
  auditor body mentions a P0, touch `tmp/ai_exchange/AUTOREVIEW_P0.flag` with a one-line
  summary. No pre-registration, no verification pass (that is the light tier; escalation
  is a manual `/cross-review`).

### D2 — Commit-triggered light review

- `scripts/install_autoreview_hook.sh` (~30 lines): writes a post-commit hook into
  `$(git rev-parse --git-common-dir)/hooks/` so ONE install covers every worktree (hooks
  run with cwd = the committing worktree's root, and the hook re-resolves paths from
  `$PWD`, so it is worktree-correct by construction). Refuses to overwrite a pre-existing
  foreign post-commit hook. `--uninstall` flag removes it.
- Hook body (~10 lines, heredoc in the installer): exit 0 immediately unless
  `BELFEM_AUTOREVIEW=1`; else
  `nohup nice -n 19 scripts/cross_review.sh --quick HEAD >/dev/null 2>&1 &` — never
  blocks the commit, never touches the tree, disowned.
- `scripts/review_status.sh` (~30 lines): commits since the last
  `## commit <hash>` entry in `autoreview.md` (`git rev-list <hash>..HEAD --count`),
  plus contents of `AUTOREVIEW_P0.flag` if present.

### D3 — `METHODOLOGY.md` (repo root, committed)

One page, factual: intro paragraph (Claude primary, Codex + Grok cross-vendor auditors,
ephemeral exchange → durable devlog); mapping table (cross-model adversarial review /
N-version programming–design diversity (Avizienis), jury/relay, pre-registration =
blinded analysis, blackboard architecture, hypothesis-driven debugging, engineering
daybook + ADRs); deviations list (pre-registered findings; source-grounded verification
of every file:line claim; physics/literature outranks reviewer agreement; durable
devlog). Internal names kept verbatim. **Shown to Christian before commit.**

### Line budget

cross_review.sh ~150 + installer ~30 + hook ~10 + review_status.sh ~30 ≈ 220 < 300. The
command file and METHODOLOGY.md are markdown, not shell. `shellcheck` is not installed on
this machine — will run it if available, otherwise note that in the devlog.

---

## §3 — Steps

- [x] R1 — `scripts/cross_review.sh` (`--jury` / `--relay` / `--quick`, diff cap, parity pick) —
      2026-08-05, 231 lines total incl. R3 scripts; all modes stub-tested (jury merge order,
      relay sequencing, quick flock+P0 flag, merge/rebase skip, file target); `shellcheck`
      unavailable on this host, `bash -n` + stub smoke used instead
- [x] R2 — `.claude/commands/cross-review.md` (pre-registration → dispatch → verification →
      reconciliation → stdout summary; never auto-fix) — 2026-08-05
- [x] R3 — `scripts/install_autoreview_hook.sh` + hook body + `scripts/review_status.sh` — 2026-08-05;
      hook installed into git-common-dir, foreign-hook guard + `--uninstall` verified by inspection
- [x] R4 — Acceptance 1 complete 2026-08-05: full jury round in
      `tmp/ai_exchange/review_gesvd_spline_wip.md` (frozen pre-registration, blind Codex+Grok,
      per-citation verification, reconciliation table). No P0; two single-raiser P1 coverage
      gaps routed to Christian
- [x] R5 — Acceptance 2 complete 2026-08-05: throwaway worktree, default-OFF no-op proven,
      opt-in commits 0.01–0.02 s (non-blocking), Codex findings landed in the worktree's
      `autoreview.md`, `review_status.sh` counts correct; worktree + branch removed, hook kept.
      D1 issue found and fixed: P0-flag grep matched the negation "No P0/P1 defect" — negation filter +
      quick-priming token rule, both polarities unit-tested. Caveat: hook is a silent no-op in
      worktrees whose checkout predates `scripts/` (self-resolving once committed)
- [x] R6 — `CLAUDE.md`: usage-lines-only section added under AI Cooperation — 2026-08-05
- [x] R7 — `METHODOLOGY.md` drafted 2026-08-05 (repo root); mapping table + deviations list
      approved by Christian in R9 and **committed** in `bd73583b` together with the scripts
      (box closed 2026-08-09 — the "NOT committed" caveat is obsolete)
- [x] R8 — Devlog `dl20260805_cross_review_tooling.md` written + indexed 2026-08-05; the
      demonstration thread `review_gesvd_spline_wip.md` is deliberately NOT swept — it is the
      Acceptance-1 exhibit; sweep after Christian has inspected it
- [x] R9 — METHODOLOGY.md approved with edits (2026-08-05, Christian): header block dropped;
      `doc/literature_references.md` extended with the authoritative alias index
      (paper0–paper9 + paperA, F0–F4, one line each: alias → txt → citation; paperA/paper9
      drift to schnaubelt2023 documented rather than hidden; newer author-year papers listed
      unaliased)
- [x] R10 — `scripts/scls_env.sh` (2026-08-05): idempotent SCLS pin (prepend
      /opt/scls/gcc/bin + generation-captured Python prefix only if absent; no reordering),
      hard-fail guards non-interactive / warn interactive; sourced first by cross_review.sh
      and the hook (hook prints the failure but never blocks). DEVIATION flagged for
      Christian: the g++ guard is existence-only, not under-/opt/scls/ — SCLS ships no
      compiler binary; its mpicxx wraps the bare system g++ (`mpicxx -show`), so the
      as-specified guard can never pass. Verified: bare-env source PASS, double-source
      idempotent, guards pass on this host
- [x] R11 — Sliced tooling commit landed 2026-08-05 (scripts/, command file, CLAUDE.md, devlog,
      METHODOLOGY.md, literature index, this plan + todo/README.md entry, + `.gitignore`
      whitelist for `.claude/commands/cross-review.md`; WIP excluded). Bare-env re-demo PASSED:
      `env -i HOME=... PATH=/usr/bin:/bin` commit in a fresh worktree took 0.01 s; the --quick
      job survived the parent shell's exit (reparented to `systemd --user`), scls_env pinned the
      toolchain, Grok found the planted off-by-one and labeled it P0, the P0 flag captured the
      real finding heading (D1 negation filter holding), `review_status.sh` counted correctly.
      Worktree + branch removed; hook remains installed

## §4 — Open questions

- O1 — Jury blinding: brief says auditors get "target + minimal context"; the manual
  habit sometimes allowed auditors to read Claude's exchange entry. Plan follows the brief
  (jury = fully blind, relay = sees prior findings). Confirm.
- O2 — `--quick` P0 detection is a grep for `P0`/`CRITICAL` in the auditor body — cheap
  and slightly over-sensitive. Acceptable for a glance-flag?
- O3 — Hook target `--quick HEAD` reviews the commit that just landed; merge commits and
  rebases produce noisy diffs — skip merges (`git rev-parse HEAD^2` test)? Proposed: yes.
- O4 — `review_status.sh` counts against `autoreview.md` only (not manual devlog audits).
  Good enough for "unreviewed commits since last audit entry"?
