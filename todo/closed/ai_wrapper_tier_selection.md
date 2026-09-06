# Explicit Model + Effort Selection for the AI Auditor Wrappers

> **CLOSED 2026-08-30 — SOLVED, moved out of the active set by the register/todo currentness sweep.**
> All thirteen steps landed and the R10 gate ran green; the depth table lives in
> `doc/ai_collaboration_protocol.md` §9.1 and both wrappers stamp `(model=…, effort=…)` into the exchange
> header. Zero open checkboxes.

**Date:** 2026-08-30
**Purpose:** Make the backing model and the reasoning effort an explicit, validated, and
**recorded** parameter of every `ask_codex.sh` / `ask_grok.sh` invocation, instead of a silent
inheritance from each vendor's user config. Mechanism: two env vars per wrapper, an allowlist
validator, a tier table in the protocol, and the chosen tier stamped into the exchange-file
header line so any finding is attributable to the depth that produced it.
**Module:** `.claude/scripts`, `scripts/cross_review.sh`, `.claude/commands/cross-review.md`,
`doc/ai_collaboration_protocol.md`, `CLAUDE.md`
**AIs involved:** Claude (plan), Codex (audit), Grok (audit)
**Status:** ✅ COMPLETE 2026-08-30. All thirteen steps landed and the R10 gate ran green (a)-(j),
including the two checks that prove the design rather than the code: an 821 KB prompt through the
newly-flagged argv path (6.4× the historical E2BIG breaking point), and `--quick` stamping
`luna`/`medium` while the invoking shell exported `sol`/`xhigh`. Two jury rounds ran: a plan round
(Codex REVISE, Grok APPROVE-WITH-CORRECTIONS) and a code round (ten findings, all confirmed and
fixed — see §4.2). O1-O6 decided by Christian, all six as proposed, on the principle that the change
which starts recording depth must not also change it. **Residual:** the `sol` rung is allowlisted
but named by no row and stays that way until somebody A/Bs it against `terra/xhigh` on a real diff;
the wrapper defaults do not track vendor-config drift and will silently diverge if
`~/.codex/config.toml` or `~/.grok/config.toml` moves.

> **Scope guards:**
> - OUT of scope: the Grok quality gate, the retry/salvage control flow, the sandbox posture, and
>   the tool allowlists. This plan touches argv construction, one `printf` per wrapper, the
>   dispatch driver, and documentation.
> - OUT of scope: teaching the wrappers to *choose* a tier. The tier is chosen by the caller per
>   the §5 table and merely enforced and recorded by the script.
> - **Partial exception to "no prompt-text changes"** (flagged by Grok): `ask_grok.sh:166` tells
>   the model "the wrapper adds the `# GROK <timestamp>` header itself". R5 changes that header,
>   so that one preamble line is corrected. Nothing else in either preamble is touched.
> - **Compatibility means "does not fail", not "preserves depth".** R1 and R3 both change the
>   effective tier of existing unset call sites. That is the subject of O4 and O5, not a
>   side effect to be waved through.

**Evidence labelling.** Claims about vendor model lists, default efforts, and the local vendor
configs were read from `~/.codex/models_cache.json`, `~/.codex/config.toml`,
`~/.grok/models_cache.json`, `~/.grok/config.toml`, `codex exec --help` and `grok --help` on
2026-08-30. Those files are outside the repository, so both auditors correctly returned
CANNOT VERIFY on them. They are marked **[locally observed]** below and are not repository-verifiable;
R10 is the only step that can settle them.

---

## 1. Current Behaviour and How It Fails

`ask_codex.sh:100-105` builds the whole invocation as:

```
"$CODEX" exec --sandbox read-only -C "$BELFEM" --output-last-message "$TMPOUT" - < "$TMPPROMPT"
```

No `-m`, no `-c model_reasoning_effort`. Every Codex audit dispatched through this wrapper
therefore inherits `~/.codex/config.toml`, currently `model = "gpt-5.6-terra"` and
`model_reasoning_effort = "medium"` **[locally observed]**.

`ask_grok.sh:197-200` does plumb effort (`EFFORT_ARGS=(--effort "$GROK_EFFORT")`), but the
variable defaults empty (`:115`), so the flag is omitted and Grok audits inherit
`~/.grok/config.toml` `default_reasoning_effort = "xhigh"` **[locally observed]**. There is no
model knob at all.

| Failure | Mechanism | Evidence |
|---|---|---|
| Silent tier drift | the effective depth of every audit is a vendor config file outside the repo; editing it retunes the whole jury with no trace | `ask_codex.sh:100-105`; `ask_grok.sh:115,197-200` (in-repo, verified). The *current values* are **[locally observed]** |
| Unattributable finding | the exchange records vendor + timestamp only, so a thin round-1 finding is indistinguishable from a thin deep one | `ask_codex.sh:124`, `ask_grok.sh:479` |
| Escalation that downgrades (Codex) | `gpt-5.6-sol` carries `default_reasoning_level: "low"`, so `-m gpt-5.6-sol` alone lands *below* terra/medium | **[locally observed]** `~/.codex/models_cache.json`. This is a hypothesis from the repo's point of view; it is the reason effort must always be passed *with* the model |
| Escalation that downgrades (Grok) | the config default is already `xhigh`, so an explicit `high` would be a reduction | **[locally observed]** `~/.grok/config.toml` |
| Wrong effort value list in the docs | `ask_grok.sh:28` advertises `none\|minimal\|low\|medium\|high\|xhigh\|max`; `grok-4.6` offers `low\|medium\|high\|xhigh` **[locally observed]** | `ask_grok.sh:28` (in-repo, verified) |
| Asymmetric wrappers | Codex has a model choice and no plumbing; Grok has effort plumbing and no model knob | as above. Description of today, not itself a failure |

~~Stale documented quirk: `ask_grok.sh:28` claims `GROK_EFFORT` is unsupported.~~
**Struck 2026-08-30, refuted by both auditors.** Line 28 makes no such claim, and the durable
in-repo record `doc/lessons_learned_evidence.md` INC-487 already states the HTTP 400 came from the
retired backing model `grok-build` and that "the constraint was model-specific and expired". The
stale belief survives only in the extra-repo memory file. Replaced by the value-list row above.

**Bottom line:** the two wrappers spend an unknown and externally-mutable amount of reasoning per
audit, and the record cannot tell you which. Every downstream claim about "the jury agreed" is
missing the one variable that most directly controls how much the jury actually looked.

## 2. Architecture: Enforce and Record, Do Not Decide

The wrappers stay dumb. They validate what they are handed, pass it to the CLI, and stamp it into
the record. The *choice* lives in the §5 table, keyed on two things observable without
self-assessment: the **scope** of the subject and the **round number** on the same problem.
Both auditors endorsed this spine and told me not to redesign it.

Rejected alternative: keying the tier on Claude's stated confidence in the claim under audit. The
claim's confidence is the thing being audited; when it is confidently wrong is exactly when a cheap
round would be ordered.

Rejected alternative: a shared `belfem_ai_tier.sh` sourced by both wrappers. Two env-var blocks of
about fifteen lines each do not justify a third file in the chain, and the wrappers are
deliberately standalone.

**Known limit of the round-number key (Grok B1, accepted):** "round >= 2 on the same problem" is
not machine-observable — identity is the slug, and a new slug on the same bug never escalates. The
table is therefore advisory for the common path and enforced only in that the *values* must be
supplied. A narrow-looking claim that is really an MPI or lifetime bug can sit at the low rung
while both auditors agree. Mitigation: escalation triggers on a safety-boundary subject
(ownership, lifetime, MPI collectives, ABI) as well as on round number, and stays cheap to invoke.

## 3. Gap Table

| # | State | Needed for | Handled today? | Class | Citation / rationale |
|---|---|---|---|---|---|
| 1 | Codex model | tier ladder | no | (c) explicit | `ask_codex.sh:100-105` has no `-m` |
| 2 | Codex effort | tier ladder | no | (c) explicit | `codex exec` has no `--effort` flag; effort is `-c model_reasoning_effort=<x>`. Key name confirmed against `~/.codex/config.toml:2`, which sets it **[locally observed]** |
| 3 | Grok effort | tier ladder | partly | (c) explicit | plumbed at `ask_grok.sh:197-200`, never set (`:115`) |
| 4 | Grok model | symmetry + a name to stamp | no | (b) see O3 | `-m, --model` confirmed in `grok --help` |
| 5 | Tier in the record | attribution | no | (c) explicit | `ask_codex.sh:124`, `ask_grok.sh:479` |
| 6 | Tier for `--quick` hook | cost control | no | (c) explicit | `cross_review.sh:144`, reached from `install_autoreview_hook.sh:40` |
| 7 | Tier for jury/relay | jury discipline | no | (c) explicit | `cross_review.sh:121-132` |
| 8 | Invalid-value handling | typo safety | n/a today | (c) explicit | `CODEX_EFFORT` does not exist yet and `GROK_EFFORT` is unset, so the unchecked-typo path is **prospective, not observed** (Grok correction) |
| 9 | `max` / `ultra` reachability | cost + sandbox posture | n/a | (b) see O2 | `ultra` is described as "maximum reasoning with automatic task delegation" **[locally observed]**, against the `--no-subagents` already passed at `ask_grok.sh:323` |
| 10 | `none` / `minimal` / `max` in the Grok doc list | doc accuracy | no | (b) see O2 | `ask_grok.sh:28` lists them; R3's allowlist would newly reject them, which must be a logged decision rather than an accident |
| 11 | `/cross-review` dispatch line | the command keeps working | no | (c) explicit | `.claude/commands/cross-review.md:17` |
| 12 | Permission allowlist match | audits run without a prompt | unknown | (b) see O6 | `.claude/settings.local.json:64` is start-anchored on `ask_codex.sh` and has no `ask_grok.sh` entry |

### 3.1 Cross-cutting findings

- **The record stamp (row 5) is load-bearing.** Rows 1-4 without it produce a policy that cannot be
  checked after the fact.
- **Environment does not reach the child by assignment alone.** Corrected from the first draft:
  `cross_review.sh:103` uses a *prefix assignment*, which preserves the inherited environment, so
  "bare environment" was wrong. The real constraint is the opposite one, and it was verified by
  experiment: a plain (unexported) assignment inside `cross_review.sh` is **invisible** to
  `ask_*.sh`. R6 must prefix-assign onto the child, or `--quick` silently keeps running at the
  wrapper default — the exact failure this plan opens by describing.
- **Validation must fail loudly and in the right order.** `ask_codex.sh` has no `set -e`
  (deliberately — it is why the `TMPERR` dump at `:108-111` is reachable), so a validator must
  `exit 1` itself; a `return 1` would be ignored. `cross_review.sh` *does* use `set -euo pipefail`
  (`:18`), so its unset checks must read `${VAR:-}` before any bare expansion or bash reports
  "unbound variable" instead of the intended message.

## 4. Ordered Steps

- [x] **R1** `ask_codex.sh`: add `CODEX_MODEL` (default `gpt-5.6-terra`) and `CODEX_EFFORT`
      (default `medium`, per O5). Treat an empty
      string as unset (`${VAR:-}`). Order: detect supplied-vs-unset, default, validate, exec.
      Allowlist `gpt-5.6-sol|gpt-5.6-terra|gpt-5.6-luna|gpt-5.5|gpt-5.4|gpt-5.4-mini` and
      `low|medium|high|xhigh`; on a miss print the allowlist and `exit 1` (not `return 1`).
      Models outside the §5 table are an escape hatch, not a tier, and say so in the usage header.
- [x] **R2** `ask_codex.sh`: insert `-m "$CODEX_MODEL" -c model_reasoning_effort="$CODEX_EFFORT"`
      as two option pairs placed **before the positional `-`** (after: R1). *(Amended 2026-08-30
      after the code round: the draft over-specified the position as "after
      `--output-last-message`"; the code places them after `-C`, which satisfies the only
      constraint that matters. Grok flagged the mismatch and was right that the box had been
      ticked against prose the code did not follow.)* `-c` and its `key=value` stay two separate words. The prompt
      must remain on stdin via `- < "$TMPPROMPT"`; nothing moves to argv (the E2BIG history is at
      `:93-97` and INC-492). Keep `CODEX_EXIT=$?` immediately after the invocation.
- [x] **R3** `ask_grok.sh`: add `GROK_MODEL`, allowlist `grok-4.6|grok-4.5` (**not** the Codex
      list), default `grok-4.6` per O3. Validate `GROK_EFFORT` against `low|medium|high|xhigh`
      and default it to `xhigh` per O4. Pass `-m "$GROK_MODEL"` **inside `run_grok`**, not at the call site — `run_grok` is
      reused for the `--resume` salvage pass (`:412`), and a model flag added only at the first
      call site would be silently dropped on salvage. Update the `:28` value list in the same edit.
- [x] **R4** Both wrappers: track supplied-vs-defaulted **per knob** and emit one stderr warning on
      a default. Annotation is per knob, e.g. `model=gpt-5.6-terra, effort=high [defaulted]`, never
      a trailing comma clause that reads as part of the effort value (Grok).
- [x] **R5** Both wrappers: change the exchange header to
      `# CODEX <ts>  (model=<m>, effort=<e>)` with per-knob `[defaulted]` markers, and the Grok
      equivalent (after: R4). The `# CODEX` / `# GROK` token stays first on the line. Correct the
      one preamble line at `ask_grok.sh:166` that describes the old header, and the comment block
      at `ask_grok.sh:54-57`.
- [x] **R6** `cross_review.sh` (after: R1-R3 — it needs wrappers that *honour* the vars, not the
      stamp): resolve the four names per mode, then **prefix-assign all four onto the child** in
      `run_auditor` alongside the existing `AI_EXCHANGE_SLUG`. `--quick` sets them itself so the
      unattended hook is pinned and a user-exported value cannot leak in. `--jury`/`--relay`
      require all four from the caller — including `GROK_MODEL`, missing from the first draft
      (Codex) — checking `${VAR:-}` first, and exit 1 naming the §5 table.
- [x] **R7** `doc/ai_collaboration_protocol.md`: add §9.1 "Depth Selection" with the §5 table;
      update the §2 sample (`:74-97`) **and the rule bullet at `:101`** ("Always start with
      `# AI_NAME YYYY-MM-DD HH:MM:SS TZ`", which the stamp would otherwise violate); update §9
      step 3 (the argv description) and step 5 at `:316`.
- [x] **R8** `CLAUDE.md`: update the prose-sweep example at `:545` and the cross-review cheat-sheet
      at `:23-25`. Then run `scripts/check_doc_claims.py` (mandatory after any `CLAUDE.md` edit;
      note it probes flags/targets/executables and will not itself check env-var examples).
- [x] **R9** Correct the stale `GROK_EFFORT`-is-broken belief in the extra-repo memory file
      `feedback_grok_third_voice_rules.md`. **Do not** edit the dated devlogs
      (`dl20260612…`, `dl20260723…`) — they are immutable records, and INC-487 in
      `doc/lessons_learned_evidence.md` already records the expiry correctly. The `:28` value-list
      fix rides with R3, not here.
- [x] **R10 (gate)** See §4.1. The only executable gate in the plan.
- [x] **R11** `.claude/commands/cross-review.md:17`: add the tier export to the dispatch line and
      teach which §5 row applies. **Must land in the same change as R6** or `/cross-review` breaks
      on arrival (Grok).
- [x] **R12** Document the new variables in both wrapper usage headers (`ask_codex.sh:5-9`,
      `ask_grok.sh:22-36`).
- [x] **R13** `.claude/settings.local.json`: verify whether a leading `VAR=value` prefix still
      matches the start-anchored `Bash(/…/ask_codex.sh *)` rule at `:64`, and add the missing
      `ask_grok.sh` entry. If prefixes do not match, prefer `export` on a preceding line. See O6.

### 4.1 R10 — the executable gate

**Run 2026-08-30.** Nine of ten green; (i) rides on the code-audit jury round. Evidence:

| Check | Result |
|---|---|
| (a) | `bash -n` clean on all three scripts |
| (b) | `(model=gpt-5.6-luna, effort=low)` and `(model=grok-4.5, effort=low)` stamped in `tier_gate.md` — both vendors accept the flags |
| (c) | all four invalid values exit 1 with the allowlist printed, before stdin is read |
| (d) | model unset, effort `xhigh` → `(model=gpt-5.6-terra [defaulted], effort=xhigh)`; only the unchosen knob is marked |
| (e) | 821 036-byte prompt accepted, 6.4× the ~128 KiB argv breaking point — the flags did not disturb the stdin path |
| (f) | `codex exec` accepts `xhigh`; this settles one of the **[locally observed]** claims |
| (g) | shell exported `sol`/`xhigh`/`grok-4.5`/`xhigh`; `--quick` still ran and stamped `luna`/`medium`. **This is the proof for the defect Grok found** — a plain assignment would have leaked the exported value through |
| (h) | `--jury` and `--relay` exit 1 naming each unset variable and the table, with no vendor call |
| (j) | Grok `--model grok-4.5` accepted and stamped |


- [x] (a) `bash -n` on every changed script.
- [x] (b) Each wrapper at a non-default tier on a throwaway slug: exit 0, no vendor 400 on the
      effort parameter, and the header line carries the requested model and effort.
- [x] (c) An invalid value exits 1 **with no network call made** (not merely exit 1).
- [x] (d) Unset variables produce the `[defaulted]` stamp and the stderr warning.
- [x] (e) A large prompt still runs, confirming the flags did not disturb the stdin path.
- [x] (f) The round-2 pair specifically (`CODEX_EFFORT=xhigh`), since `xhigh` acceptance by
      `codex exec` is **[locally observed]** only.
- [x] (g) `--quick` produces a header stamped with the pinned cheap tier and **not** `[defaulted]`,
      including when the invoking shell exports a conflicting value.
- [x] (h) `--jury` with any one of the four unset exits 1 **before** contacting a vendor, naming
      the table.
- [x] (i) `/cross-review` still dispatches after R11 — the code-audit jury round was itself dispatched through `scripts/cross_review.sh --jury` at `terra/xhigh` + `grok-4.6/xhigh` and both entries carry the stamp.
- [x] (j) Grok `-m` accepted (flag name confirmed in `grok --help`, but the wrapper path is untested).

## 5. The Tier Table (the policy this plan enforces)

| Situation | Codex | Grok |
|---|---|---|
| Prose sweep of an ordinary guide or README; citation and doc-claim mechanics | `gpt-5.6-luna`, `medium` | not used |
| Prose sweep of a dense technical document (input reference, coding philosophy) | `gpt-5.6-terra`, `medium` | not used |
| Single narrow claim, round 1 | `gpt-5.6-terra`, `medium` | `high` |
| Plan audit or code-diff audit, round 1 | `gpt-5.6-terra`, `high` | `high` |
| Round >= 2, a round-1 split verdict, or a safety-boundary subject (ownership, lifetime, MPI collectives, ABI) | `gpt-5.6-terra`, `xhigh` | `xhigh` |
| Unattended post-commit `--quick` | `gpt-5.6-luna`, `medium` | `grok-4.6`, `medium` |

**The second rung raises effort, not model.** Changed 2026-08-30 on both auditors' advice. The
first draft escalated `terra -> sol` at round 2; Grok's objection decides it — moving both knobs
confounds which one helped, and if sol's `low` default is real then `-c` becomes the only thing
standing between "escalate" and a downgrade. `gpt-5.6-sol` is demoted to an allowlisted escape
hatch pending a measured head-to-head, which R10 does not provide.

Effort is not a substitute for a sharp prompt. The Grok narration-stub history (INC-489) was a
prompt-shape defect, and no effort setting would have touched it.

## 6. Open Design Questions — all RESOLVED 2026-08-30 by Christian

All six were decided as proposed. The governing principle he selected across O4 and O5 is
**record current depth first**: the change that starts recording depth must not also change it,
or no new finding is comparable to any old one.

- **O1 — RESOLVED: the split, as proposed** (decided 2026-08-30, Christian). Proposal: annotated default in the
  wrappers, hard fail in `cross_review.sh --jury|--relay`. Both auditors accept the split but warn
  it is toothless unless (i) every documented invocation carries the vars, (ii) the stamp says
  *which knob* was defaulted, and (iii) R11 lands with R6. Grok adds a sharp point: `run_auditor`
  sends wrapper stderr to `_cross_review.log` (`cross_review.sh:121-122,144`), so the warning is
  invisible on the jury path and the stamp is the only real signal.
- **O2 — RESOLVED: cap at `xhigh`** (decided 2026-08-30, Christian). `max` and `ultra` are
  rejected, as are `none`/`minimal`. Proposal: cap at `xhigh`; reject `max` and `ultra`, and also
  `none`/`minimal` which `ask_grok.sh:28` currently advertises. `ultra` delegating to subagents is
  **[locally observed]** and neither auditor could confirm it, so this is a cost-and-conservatism
  decision, not a proven sandbox requirement.
- **O3 — RESOLVED: add it, defaulted `grok-4.6`** (decided 2026-08-30, Christian). Only `grok-4.6` and `grok-4.5` exist **[locally observed]**.
  Proposal: add it, default `grok-4.6`, never vary in practice — the stamp needs a model name to
  print regardless.
- **O4 — RESOLVED: default `xhigh`, i.e. record current depth** (decided 2026-08-30, Christian).
  **Proposal reversed on both auditors' advice.** The first draft
  defaulted to `high`, a reduction from today's effective `xhigh`. Grok's argument decides it:
  reducing depth in the same change that starts *recording* depth destroys comparability with
  every past round while calling the work attribution. **New proposal: default `xhigh`**, i.e.
  record current behaviour first. The ladder still works, because `--quick` supplies the cheap
  rung and Codex retains luna/terra and medium/high/xhigh.
- **O5 — RESOLVED: default `medium`, option (a), symmetric with O4** (decided 2026-08-30,
  Christian). Found by Grok. The first draft defaulted `CODEX_EFFORT=high`
  against an inherited `medium` **[locally observed]** — a silent *upgrade* of every existing call,
  the same unlogged behaviour change I objected to in O4 and had not noticed in my own plan.
  Options: (a) default `medium`, recording current behaviour, and let the table raise it per
  invocation; (b) default `high`, accepting a deliberate across-the-board deepening (and a cost
  increase) as policy. Proposal: **(a)**, symmetric with O4.
- **O6 — RESOLVED: add the `ask_grok.sh` entry; prefer `export` if prefixes do not match**
  (decided 2026-08-30, Christian). `.claude/settings.local.json:64` allows
  `Bash(/…/ask_codex.sh *)`, start-anchored, with no `ask_grok.sh` entry. If an env-var prefix
  stops it matching, every audit starts prompting. Proposal: add the `ask_grok.sh` entry and, if
  prefixes do not match, `export` on a preceding line instead.

### 4.2 Code round, 2026-08-30

Jury round dispatched through `cross_review.sh --jury` at the round-2 depth. Codex REVISE (one P1,
two P2), Grok APPROVE WITH CORRECTIONS (two P1, five P2). **Ten findings, ten confirmed against the
tree, ten fixed** — full detail in `tmp/ai_exchange/review_ai_wrapper_tier_code.md`.

- [x] **D1 (P1, Codex)** `cross_review.sh` checked depth *presence* but not *validity*, so a typo
      passed the dispatcher, failed one wrapper, let the other leg bill, and exited 0 with a
      synthetic failure entry indistinguishable from a vendor outage. Fixed by validating values
      before dispatch; retested (`CODEX_EFFORT=hgih` and `GROK_MODEL=grok-5` both exit 1 with no
      vendor call). The allowlist now appears in three places on purpose: the wrappers must reject
      a bad value when called directly, the driver must reject it before paying for the other leg.
- [x] **D2 (P1, Grok)** protocol §9 "How to invoke" and workflow step 2 still taught the bare,
      defaulting invocation while §9.1 said not to default. Both recipes now carry the depth.
- [x] **D3 (P1, Grok)** the permission allowlist covered none of: the documented
      `scripts/cross_review.sh` dispatch line, effort-first prefixes, or prefixed *absolute*
      wrapper paths — the last a regression against the pre-existing absolute rule. Nine rules added.
- [x] **D4 (P2, both)** `collect()`'s synthetic failure header carried no depth, in exactly the
      case that most needs it. Fixed.
- [x] **D5 (P2, both)** the scripts and §9.1 restated **[locally observed]** vendor claims as fact.
      Requalified as one machine's config on one day, with the note that wrapper defaults do not
      track config drift, and the `xhigh` cap restated as policy rather than proven constraint.
- [x] **D6 (P2, Grok)** R2's specified flag placement was not what the code does — see the amended
      R2 above; code kept, prose corrected.
- [x] **D7 (P2, Grok)** the `-m`/`-c` comment said "two short argv words"; it is four.
- [x] **D8 (P2, Grok)** the usage header's "everything below `gpt-5.6-luna` is an escape hatch"
      excluded `gpt-5.6-sol`, the one model the plan actually demoted.
- [x] **D9 (P2, Grok)** `run_grok`'s comment named `-m` where the flag used is `--model`.
- [x] **D10 (P2, Grok)** the slash-command depth table carried Codex-only prose rows into a driver
      that always runs both legs. Rows removed with a pointer to `ask_codex.sh`; Grok cells now
      name the model.

**Split worth noting for future rounds:** Codex found the finding that costs money; Grok found the
three that make the policy toothless in practice; only D4 and D5 were found by both. The jury's
value here was vendor diversity, not voice count — consistent with why the round is run at all.

## 7. Definition-of-Done Checklist

- [x] Every gap-table row maps to a step or an open question.
- [x] Each claimed gap carries an in-repo citation, or is explicitly marked **[locally observed]**.
- [x] Ordered steps with dependencies noted (R6 after R1-R3, not after R5).
- [x] O1-O6 decided by Christian before R1/R3/R6 land, not silently chosen in code.
- [x] R11 lands in the same change as R6.
- [x] R10 (a)-(j) executed and pasted into the devlog. The plan is not "verified" without it —
      two auditors agreeing is the weakest rung of the evidence ladder (protocol §11).
- [x] `scripts/check_doc_claims.py` clean after R8.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/ai_wrapper_tier_selection.md` — pre-registration, both audits,
  and the verification/reconciliation table (18 rows).
- Plan audit 2026-08-30, Codex + Grok blind. Codex: REVISE BEFORE IMPLEMENTATION. Grok: APPROVE
  WITH CORRECTIONS. Both endorsed the architecture and told me not to redesign it.
- **Codex-only catch:** R6 omitted `GROK_MODEL` from the required set.
- **Grok-only catches:** the R6 child-environment defect (an unexported assignment never reaches
  `ask_*.sh`, so `--quick` would have stayed at the wrapper default forever — confirmed by
  experiment, not by agreement); the `/cross-review` breakage; the model flag needing to live
  inside `run_grok` for the salvage path; O5; the per-knob annotation ambiguity.
- **Both, against my draft:** the second rung should raise effort rather than switch model, and
  Grok's depth must not be reduced in the same change that starts recording it.
- **Refuted and recorded:** my `ask_grok.sh:28` claim, my "bare environment" wording, three wrong
  line numbers, and Grok's own assertion that reassigning an exported variable does not reach the
  child (false on bash 5; tested).
- Every auditor citation was re-checked against the tree before adoption; the vendor-config claims
  neither auditor could reach are relabelled **[locally observed]** rather than dropped or asserted.
