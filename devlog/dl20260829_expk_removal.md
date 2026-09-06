# Removing `expk`: a deck key that was the same knob as `fuzzyness`, wired to nothing

**Date:** 2026-08-29
**Purpose:** record the removal of the `expk` source-function key ( the second half of DR-146 ),
the ABI break it required, and the migration advice that was nearly wrong
**Module:** `src/numerics/sources`, the three BC/circuit factories, `doc/` input contract

## What it was

`function_sigmoid` builds its logistic rate from `fuzzyness` and `period`. `set_sigmoid( aExpK )`
stored `log( expk )` in a slot whose accessor had zero callers. Codex found it during a language
sweep of `circuit_usage_guide.md`; the question it left open was whether `expk` was a *broken*
knob or a *superseded* one.

The algebra answers it. With `fuzzyness = 1/(1+expk)`,

    log( (1-f)/f ) == log( expk )

and the two defaults are the same real number written twice: `log(99) = -log(0.01/0.99) = 4.59512`.
Grok took it one step further — `fuzzyness` **is** the normalised sigmoid evaluated at `t = offset`,
and `expk` is the odds at that same edge. They are one parameter in two spellings. `fuzzyness`
won; the conversion was never wired up; the old slot kept being written and never read.

## Why removal, and why now

Christian ruled removal over an alias, then set the timing: a design freeze the following evening,
so an ABI-breaking change either lands now or waits past 1.0.

## The break had to be total

I had offered a softer fallback — keep the five-argument `set_sigmoid` as a deprecated forwarder,
no ABI break. Codex showed that is the one genuinely unsafe option. `sigmoid()` was an **inline**
accessor, so a plugin compiled against the old header carries slot index 7 inside itself. Removing
the slot while keeping the forwarder gives that plugin an out-of-bounds read on a seven-element
`Cell` — a memory error dressed up as compatibility. The coherent choices are a clean break or a
full compatibility period, and nothing in between.

So: the fifth parameter, the private one-argument setter, the `sigmoid()` accessor, the
`BELFEM_BCVAL_SIGMOID` macro and slot 7 are all gone, and `mValues` is seven long. Plugins get
rebuilt. A stale five-argument call now fails to link, which is the C++ analogue of refusing the
key.

My stated reason for the break being safe was also wrong, and Codex corrected it:
`cl_SourceFunction.hpp` was already installed by the whole-tree header glob and is listed in
`todo/user_api_header_install.md`. It has been plugin API for a long time; the new umbrella header
only formalised it. The break is right for Christian's reasons, not the one I gave.

### What the break actually cost

Measured, not estimated. Christian rebuilt the user library and restarted the running
`tape_quench_usermat` quench case: "that's an easy one." `usermat.so` and `userdefect.so` were
recompiled against the 4-argument `set_sigmoid` at 22:41 and the run went straight on, past
`dlopen` and past `MaterialFactory`. The cost of the clean break, for the only person in the
project holding a real plugin, was one recompile. That is the outcome the design assumed and it is
worth recording as fact rather than as the prediction it was when the decision was taken.

The same run is also the first production exercise of DR-146's plugin branch -- `usermat { }`,
the `have(jc)` gating, `critical temperature` -- which no shipped example covers, and it confirms
that removing the never-opened `file : sst-1.hdf5` from that deck was behaviour-neutral.

## The correction that mattered most

My first draft of the error message told the user to migrate with `fuzzyness = 1/(1+expk)`.

Grok refuted it, and the refutation is the whole point of the row. That formula is the
*intended-equivalence* map, not the *preserve-your-run* map:

| deck today | keep current results | apply the formula |
|---|---|---|
| `expk : 1000`, no `fuzzyness` | delete `expk` — still runs at `0.01` | `fuzzyness : 1/1001` — **changes the waveform** |
| `fuzzyness : 1e-4` and `expk : 99` | delete `expk`, keep `1e-4` | **overwrites** `1e-4` with `0.01` |

Because `expk` never did anything, every deck carrying it has always run at whatever `fuzzyness`
said. "Removal changes no numerical result" is true only along the delete path — and my message
pointed the other way. It now says DELETE to preserve results, and offers the conversion only as
what `expk` was *meant* to express, explicitly flagged as changing the waveform.

That is the second time in one session an audit caught a defect that would have altered numbers
quietly rather than failing loudly.

## Placement

Both vendors required the check **before** the function-type dispatch, not inside the sigmoid
branch. `expk` was only ever read on the sigmoid path, so a sigmoid-only refusal would have left it
silently ignored on `ramp`, `sine`, `constant`, `userdefined`, and on any section with no `type` at
all — the same defect, relocated. `bearing` is correctly excluded: it reads no source function.

## Three YAML mistakes, all of which parsed

The schema edit went wrong three times and each version loaded without complaint:

1. an unclosed quoted scalar ( caught by Grok on the materials round );
2. `refused_keys` inserted mid-list, which silently swallowed `file`, `label` and `units` as
   properties of the removed key, deleting three valid `userdefined` keys from the accepted
   contract ( caught by Codex, blocking );
3. the corrected block indented one level too shallow, making it a sibling of
   `source_function_keys` instead of a child.

"It parses" is not verification for a machine-readable artifact. The check is now a structural
assertion — accepted keys must contain `file`/`label`/`units`/`fuzzyness` and must not contain
`expk`; `refused_keys` must hold exactly `expk` with exactly its five fields; no stray sibling may
exist — plus a cross-check that every source-function key the factory actually reads appears in the
accepted set. Same shape as the stale-binary trap: a green signal measuring the wrong thing.

## Documentation

The migration path is documented rather than erased, on Codex's advice: deleting every mention of
`expk` would hide the migration from exactly the person who needs it. `input_file_reference.md`
carries a struck row with the delete-don't-translate warning, the schema carries a `refused_keys`
block, and `circuit_usage_guide.md` says the same in prose.

## Status

Both audits satisfied — Codex REQUEST-CHANGES on the YAML defect, fixed and re-verified; Grok
APPROVE on the working tree, having noticed the diff it was sent predated the fix and re-read the
tree instead.

**`make check` has NOT run on this change.** The suite that passed at 22:14 predates it: the test
binary contains none of the three new cases and the built library still exports the five-argument
`set_sigmoid`. A rebuild is required before this half is verified.

Known coverage gap, left visible: the three tests exercise `ElectricalCircuitFactory` only. The
Maxwell and Thermal guards are confirmed statically by both vendors but not locked by a test, since
neither suite has a deck-driven fixture to hang one on.
