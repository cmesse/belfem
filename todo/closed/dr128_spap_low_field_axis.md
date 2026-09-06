# DR-128 — extend the sp-ap B axis below 10 mT from the measured self-field data

**Date:** 2026-08-29
**Purpose:** Remove the sp-ap table's 10 mT B-axis floor, which clamps jc/n to
B-independent values across the entire operating range of mT-scale decks and puts a
∂/∂B tangent discontinuity exactly where the field crosses the edge.
**Module:** physics/materials (table data only — no C++ change)

**Status:** **CLOSED 2026-09-03 by the tree-first `todo/` sweep — the transition function is
IMPLEMENTED, contrary to this plan's own "nothing implemented" framing.** The low-field
transition lives at `src/physics/materials/cl_JcFunction_Database.hpp:440-480`, carrying the
measured rationale in its own comment (the historical hard clamp was wrong by up to 23 % near
Tc, and its tangent jumped from 0 to the spline slope the instant |B| crossed `mBmin` — a
Newton tangent); `eval` and `deval_dB` apply it at `:580-671` behind `mHaveSelfField && normB <
mBmin`. Christian's §7 ruling was followed: no clamp, a Kohler-pattern transition.
**Residue: validation gates only** — G1–G4 are not evidenced by any test or run. Struck is not
verified.

*Superseded status text:* "PLAN v2 — drafted 2026-08-29. **§3 (axis extension) is SUPERSEDED by
§7**, Christian's ruling: follow the Kohler pattern the metals already use — do not clamp, build
a transition function." §3 is kept as the fallback if §7 is rejected. Round-1 jury was dispatched against §3 and its verdict is read in that light.

---

## 1. The defect, restated with measurements

The sp-ap table's B axis is `log10 B = -2 .. 1` (`build_spap.py:9,29`), i.e. 10 mT … 10 T.
`JcFunctionDatabase::eval` clamps `normB` into `[mBmin, mBmax]`
(`cl_JcFunction_Database.hpp:180`) and `deval_dB` returns exactly `0.0` outside that window
(`:203`). Both are internally consistent — a clamped function really does have zero
derivative — but the consequences are:

1. **No field dependence at all below 10 mT.** The tapestack3d deck runs at
   `max|B| = 2.7 … 5.6 mT` through the 2.0–2.4 s window, entirely under the floor.
2. **A C1 kink at the floor.** `∂jc/∂B` jumps from 0 to the spline's one-sided slope the
   moment `|B|` crosses 10 mT, and that is a Newton tangent. The tapestack3d field grows
   ×1.047/frame and crosses 10 mT at roughly t ≈ 2.7 s, so the kink is exercised in
   production, not hypothetically.

**How wrong is the clamped value?** Measured from `profiles.npz` — the *reduced anchors the
table is actually built from* — comparing the `B = 0` slice against the `B = 10 mT` slice at
the exact setpoints:

| T [K] | θ | jc(0)/jc(10 mT) | n(0) | n(10 mT) |
|---|---|---|---|---|
| 65 | 0° | 1.0063 | 27.46 | 27.01 |
| 75 | 0° | 1.0232 | 27.06 | 25.62 |
| **77.5** | **0°** | **1.0326** | **26.35** | **25.05** |
| 77.5 | 45° | 1.0130 | 26.41 | 25.87 |
| 77.5 | 90° | 1.0077 | 26.36 | 26.47 |
| 80 | 0° | 1.0598 | 26.19 | 23.03 |
| 85 | 0° | 1.2346 | 23.08 | 16.90 |

So at the operating point the clamp **under**estimates jc by ~3.3 % and n by ~5 %; at 85 K it
is 23 % and 37 %.

> **Correction (round-1 jury, Codex finding 4 — accepted).** An earlier draft of this table
> quoted raw-`sp-ap_htsdb.json` rows averaged over ±1.5 K and ±2° matching windows, which
> smears across setpoints and *damps* every ratio (it gave 1.001 / 1.013 / 1.022 for the
> 65 / 75 / 77 K rows above, and named a 77 K setpoint that does not exist — the axis carries
> 77.5 K). The reduced anchors are the correct source because they are what the builder
> consumes. **The defect is therefore slightly worse than first reported, not better.** The
> §7 validation was already computed entirely against `profiles.npz` and is unaffected.

## 2. The unlock: the primary data already contains B = 0

The DR-128 row records "htsdb source has no sub-10 mT scans", and that is true as written:
there is **nothing strictly between 0 and 10 mT**. But the source is not silent below the
floor — it contains **946 rows at exactly B = 0**, spanning 35 temperatures with full angle
coverage, and the reduction already carries them:

- `profiles.npz` has `BSET[0] = 0.0` and `LNJC[:, 0, :]` / `LNN[:, 0, :]` are **100 % finite**
  over all 16 setpoint temperatures × 181 angles.
- `build_spap.py:43` deliberately drops that slice from the B-axis interpolation
  (`lbs = np.log10(BSET[1:])`) — because `log10(0) = -∞` cannot sit on a log axis.
- The same slice is already trusted elsewhere in the build as physics: the high-T backbone
  uses the measured self-field sweep as the T-collapse at the low-B edge (`:149-155`).

**Therefore the axis can be extended from measured data at both ends.** This is not the
"inventing a decade for numerical convenience" that the earlier jury rejected: the bottom
anchor is a measurement, and only the *shape* of the connection between two measured
endpoints is modelled — which is what the existing table already does between every pair of
its 20 field setpoints.

## 3. Design

Extend the target grid from `LBg = -2 .. 1` to `LBg = -4 .. 1` (step 0.1, `NB` 31 → 51),
i.e. down to 0.1 mT, and fill the new nodes with a monotone C1 connection:

- **Bottom (`log10 B = -4`):** the measured B = 0 value, `LNJC[:, 0, :]` / `LNN[:, 0, :]`,
  with **slope zero**. Physically, jc(|B|) flattens into the self-field plateau; numerically,
  a zero slope at the new floor is what makes the surviving clamp below it C1-consistent, so
  the kink is removed rather than relocated.
- **Top (`log10 B = -2`):** the existing value and the existing one-sided slope of the
  in-range interpolant, so nothing above 10 mT moves by construction.
- **Between:** a cubic Hermite in `log10 B` matching both endpoint values and both slopes,
  checked for monotonicity in B and falling back to a monotone (PCHIP-style) limiter where
  the Hermite would overshoot. Most of the variation lands in the top decade, which is the
  Kim-like behaviour the measured setpoints already show above 10 mT.
- The same construction is applied to `n`, using its own B = 0 slice and its own validity
  mask.

**No C++ change is required.** `mBmin` and `mBmax` are derived from the file's own axis
(`cl_JcFunction_Database.hpp:99-100`), so a rebuilt table with a longer axis is picked up
with no code edit and no format change. Old tables keep working; they simply keep their
old floor.

## 4. Steps

- [ ] **D1** — extend `LBg`/`NB` in `build_spap.py`; thread the widened grid through the
      anchor arrays `AJ`/`AN` and every consumer of `NB` (the Hilton correction field, the
      high-T backbone `g`, `enforce_monotone_tb`).
- [ ] **D2** — implement the two-sided Hermite fill with the monotone fallback, as its own
      helper, applied to jc and n independently.
- [ ] **D3** — re-run the existing monotonicity guarantees (`enforce_monotone_tb`, the PAVA
      pass) over the widened grid; they must hold on the extension, not just in range.
- [ ] **D4** — C1 gate **on the B-spline control net**, not the nodal values and not the
      element Bernstein net (this is the standing rule for these tables); the gate must pass
      across the new `-2` join specifically.
- [ ] **D5** — regenerate `sp-ap.hdf5`, re-embed provenance (`embed_provenance.py`), and
      diff the in-range surface against the current file: the region `log10 B >= -2` must be
      **bit-comparable or explainably tiny**, since nothing there is meant to change.
- [ ] **R1** — jury round on this plan before any rebuild.
- [ ] **R2** — jury round on the rebuilt table's diff + gate output.
- [ ] **O1** — decide whether sst-1 gets the same treatment (it shares the pipeline and
      **shares the floor — confirmed by the round-1 jury, not assumed**). Under §7 this is fixed for free, since the transition is in the reader and every table gets it at once.

## 5. Gates (from the DR-128 row, plus one)

- [ ] **G1** — table probe: jc and n continuous through 10 mT with matching one-sided
      `∂/∂B` (the row's own gate).
- [ ] **G2** — deck rerun: tapestack3d no longer parks `max|B|` at the edge, and the
      Newton behaviour through the t ≈ 2.7 s crossing is unchanged or better.
- [ ] **G3** — in-range regression: a deck whose field never goes below 10 mT produces the
      same answer as before the rebuild (this is the one that protects every existing HTS
      result).

## 6. Known gaps

| Gap | Why it is acceptable / what closes it |
|---|---|
| The shape between 0 and 10 mT is modelled, not measured | Both endpoints are measured; the interpolant is the same class already used between the table's other setpoints. Closed only by a sub-10 mT scan, which the source does not have |
| `∂jc/∂B` at exactly B = 0 is set to zero | jc(\|B\|) generally has a cusp at the origin; zero slope is the choice that makes the residual clamp C1-consistent, and the deck never evaluates below 0.1 mT anyway |
| 85 K error is large (21 %) and this fix inherits the same two-point limitation there | The fix strictly improves it (from a 21 % clamp error to a measured endpoint); it does not make near-Tc decks trustworthy on its own |
| sst-1 not addressed | O1 |


---

## 7. REDESIGN (supersedes §3) — a Kohler-style transition, no clamp, no rebuild

**Christian's ruling, 2026-08-29:** "Do we even have to clamp? Check out my Kohler functions
for the metals. There I simply create a transition function."

`Copper::create_kohler` (`cl_Material_Copper.cpp:260-340`) is the pattern. The measured
magnetoresistance fit is valid only above a switch field `x = e^2.3`; below it the code does
not clamp, it builds a **cubic in linear x** matched to the fit's value, first and second
derivative at the switch point, with the constant term forced to zero (`p0(3) = 0.0`) so that
the magnetoresistance vanishes at B = 0 exactly as physics requires.

The same move applies here, and it dissolves the whole `log10(0) = -inf` problem that forced
§3's axis extension: **the transition lives in linear B, where B = 0 is an ordinary point.**

### 7.1 The form

Below `Bmin`, replace the clamp with

```
    b      = min( (djc/dB)|Bmin / ( 2 Bmin ), 0 )       // monotone guard
    jc(B)  = b B^2 + ( jc(Bmin) - b Bmin^2 )
    djc/dB = 2 b B
```

- `djc/dB = 0` at `B = 0` **by construction** — jc has its maximum at zero field, and the
  quadratic has no linear term, exactly as Kohler's cubic has no constant term.
- Value and slope match the spline at `Bmin`, so the join is C1 — which is precisely the
  DR-128 gate ("continuous through 10 mT with matching one-sided ∂/∂B").
- Monotone decreasing on `[0, Bmin]` whenever `b < 0`; the `min(..., 0)` guard covers the
  low-T samples where the spline slope at `Bmin` is flat or noise-positive, and there it
  degenerates to today's clamp (whose slope is ~0 anyway, so C1 survives to within that noise).
- Same construction for `n`.

**Why quadratic and not Kohler's cubic:** the table is an **order-2** tensor mesh
(`build_spap.py:11`), so its second derivative is piecewise constant and discontinuous at
element boundaries. Matching `d²/dB²` at `Bmin` would import that element noise into the
extrapolation. C1 is the well-posed contract for this spline; Kohler's C2 is right there
because its fit is an analytic Bezier/polynomial.

### 7.2 Validated against data we do not have to ship

The form makes `jc(0)` a **prediction**, so it can be tested against the 946 measured
self-field rows. Over all 208 sampled (T, θ) pairs, predicted-over-measured `jc(0)`:

| form | mean | worst error | monotone |
|---|---|---|---|
| C2 cubic (Kohler's exact form) | 0.9941 | 6.43 % | 87 / 208 |
| **C1 quadratic + guard (chosen)** | **0.9956** | **8.88 %** | **208 / 208** |
| clamp (today) | 0.9855 | 20.63 % | 208 / 208 |

The chosen form has the best mean, is monotone everywhere, and roughly halves the worst-case
error against the clamp. At the 77.5 K operating point it recovers `jc(0)` to 1.1 % (0° ),
0.7 % (45°) and 0.02 % (90°); the clamp is 3.2 % / 1.3 % / 0.8 % low there and 19 % low at
85 K / 0°.

### 7.3 Why this is better than §3

| | §3 axis extension | §7 transition |
|---|---|---|
| table rebuild | required | **none** |
| file format | new axis, `NB` 31 → 51 | **unchanged** |
| existing tables (sst-1, bscco, every shipped deck) | keep their old floor until each is rebuilt | **all fixed at once** |
| risk to the in-range surface | real — the high-T backbone normalizes by `slp[0]`, the low-B-edge slope, whose index moves | **zero, nothing in range is touched** |
| C++ change | none | ~15 lines in one header |
| B = 0 measured data | must be shipped in the file | used only to *validate*, never shipped |

### 7.4 Revised steps

- [ ] **D6** — implement the transition in `cl_JcFunction_Database.hpp`: `eval`, `deval_dB`,
      and the θ/T derivative paths must all agree below `Bmin` (today they each clamp
      independently). The `b` coefficient depends on (T, θ), so it is computed per call from
      `evaluate_derivy` at `Bmin` — no state, no cached table.
- [ ] **D7** — the same treatment for `n`, via whichever JcFunction subclass serves it.
- [ ] **D8** — confirm every other consumer of `mBmin` is consistent (`:180, :203, :226, :246`).
- [ ] **R3** — jury round on the §7 design + diff.
- [ ] ~~D1–D5~~ — superseded; kept in §4 as the fallback path only.

### 7.5 Revised gates

- [ ] **G1** (unchanged) — probe: jc and n continuous through 10 mT with matching one-sided
      ∂/∂B. Now expected to pass **by construction**, so the probe is a check on the
      implementation rather than on the design.
- [ ] **G2** — deck rerun: tapestack3d through the t ≈ 2.7 s crossing of 10 mT, Newton
      behaviour unchanged or better, no new iterate cost.
- [ ] **G3** — in-range regression: any deck with `|B| > 10 mT` throughout must be
      **bit-identical**, since nothing above `Bmin` is touched. This is now a much stronger
      claim than under §3 and should be the first gate run.
- [ ] **G4** (new) — the `jc(0)` prediction check above, rerun as a unit-style probe against
      the htsdb rows, so the 8.88 % worst case is a recorded number rather than a session note.

---

## 8. STATE as of 2026-08-29 (read this first on resuming)

**Nothing is implemented. `src/physics/` is untouched by this work.** DR-128 is plan-only.

### 8.1 Round-2 jury on §7 — BOTH REQUEST CHANGES ("right shape, not implementation-ready")

Both voices accept the diagnosis and accept that a linear-|B| C1 quadratic below `Bmin` is
the correct shape (no rebuild, B = 0 is an ordinary point). Both refuse the current §7 as
written. Findings to fix before any code:

- **F1 (both, high) — the monotone guard breaks C1 exactly where it fires.** With `s > 0` the
  guard sets `b = 0`, so the transition-side slope is 0 while the spline side keeps `s`. My
  "C1 within noise" wording is wrong: it is an exact derivative discontinuity, just a small
  one. It fires in **58 of 208** sampled (T, θ). Choose openly: keep every boundary slope and
  lose monotonicity on those slices, or keep monotonicity and report the discontinuity as a
  measured number. Do not claim both.
- **F2 (both, high) — the θ and T derivative paths need the chain rule I never wrote.**
  `b = b(T, θ)`, so below `Bmin` the θ/T tangents pick up `∂b/∂θ · B²` and `∂b/∂T · B²`
  terms. `deval_dT` is a live thermal Newton tangent. `Database` exposes no mixed partials,
  so this needs analytic mixed spline derivatives, a deliberately different transition form,
  or an explicitly documented deferral with its numerical impact. **D6 is not a ~15-line
  change** — that estimate was wrong.
- **F3 (Codex, high) — `n` must not inherit the jc model.** Applying the same
  monotone-decreasing transition to `n` contradicts this file's own measurements. Treat `n`
  independently, validate over all anchors, add an `n(0)` error gate; G4 currently validates
  only jc.
- **F4 (Codex) — a reader-wide change is not justified by sp-ap-only evidence.** Either make
  the transition opt-in via table metadata / constructor policy, or validate jc and n
  independently for every affected table first.
- **F5 (Codex) — my error percentages use inconsistent denominators.** Report either
  "measured is X % above the clamp" or `(measured − clamped)/measured`, consistently.
- **F6 (Grok) — the Python rosetta-stone reader was omitted** from the change scope.
- **F7 (Codex) — the Kohler analogy is structurally useful but overstated** in §7's prose.

### 8.2 Peer session (lift-factor storage) — coordination, BINDING

A parallel session is converting all eight `share/material` tables from absolute jc to a
dimensionless lift factor `L = log10(jc/jc_ref)`, and edits **the same file**
(`cl_JcFunction_Database.hpp`: a `const real mLogRef` added inside `eval()` before the `pow`).

- **They go first. I hold `src/physics/` — no edits, no commits — until they ping.** Then I
  rebase onto what they land. Their edit does not touch the below-`Bmin` branch, which is
  mine; the math commutes (an additive constant in log space passes through a quadratic built
  from `jc(Bmin)` and its slope, both shifted by the same factor).
- Their caution for my rebase: `mLogRef` must **not** reach `min_value()`
  (`cl_JcFunction_Database.hpp:132-135`), which `cl_Material.cpp:678-681` calls on the *n*
  function for the n-floor warning — n is never shifted.

### 8.3 Peer data input — the cusp-vs-maximum question now has a DATA answer

They measured the shipping tables' embedded source rows at exactly B = 0, within 1.5 K of
77 K: 58 rows each for fesc / fysc / sp-ap / sst-1 / superox, spread across 55 stage angles
of only **0.25 – 0.73 %**; one row each for the two sch04 tables; bscco carries documented
digitised meta instead.

**Therefore measured jc at B = 0 is angle-independent to well under 1 %.** Any below-`Bmin`
form whose B → 0 limit keeps a visible angle dependence is contradicted by data, whatever the
literature says. This is a stronger constraint than the literature route CODEX-6 asked for,
and it should be checked against the chosen form before implementation — my current quadratic
inherits `jc(Bmin)`'s angle dependence at B = 0 and **may well fail it**. Open question, and
it is now the single most likely thing to change the form.

### 8.4 Free validation gate they hand us

`meta/Icw_77p5K_sf_A_per_m` already exists on all seven REBCO tables. Once `L` is referenced
to measured true self-field, "L → 1 as B → 0" becomes a physical statement about measurement
rather than a normalisation artifact. Under today's clamp L at the `Bmin` node is ≈ 0.97
(1/1.0326 = 0.968 at 0°, 1/1.0077 = 0.992 at 90°). Under the bridge it should reach 1.00.

- [ ] **G5 (new)** — after both changes land: evaluate `L` at 77.5 K, B → 0, on five tables
      independently, at every angle. Should be 1 to within the bridge's own accuracy. An
      absolute check against measurement, not self-consistency; overshoot shows up directly
      as `L(0) > 1`.

### 8.5 Next action on resuming

Answer 8.3 first (does the chosen form respect angle-independence at B = 0?), because it can
invalidate the form. Then F1–F7. Then round 3. Do not write code before the peer pings.

### 8.6 §7 FORM REFUTED 2026-08-29 — verified against the shipped table

The peer session (lift-factor storage) measured the shipped tables' B = 0 rows; I re-verified
against `share/material/sp-ap.hdf5` directly. Results:

- **The angle-spread failure is real and the form is refuted.** At T = 77 over 181 angles:
  angle spread of `jc(Bmin)` = 4.82 %, of the predicted `jc(0)` = 2.61 %, against a **measured
  0.33 %**. The form moves the right way and falls several times short, because
  `jc(0) = jc(Bmin)·(1 − s/2)` inherits `jc(Bmin)`'s angular shape by construction. Peer's
  independent numbers at 77.5 K: sp-ap 5.67 → 4.01 vs 0.33 measured; superox 8.28 → 5.88 vs
  0.73; sst-1 0.84 → 0.89 vs 0.25; fesc 4.12 → 2.14 vs 0.65; fysc 2.44 → 1.20 vs 0.68.
- **Repair direction:** make angle-independence at B = 0 *structural* — pin the B → 0 limit to
  an angle-averaged self-field value and let the angular dependence switch on with field —
  rather than letting it fall out of the `Bmin` slice.
- **Development vs validation target:** develop against sst-1 (its table reproduces its own
  measurements to −0.02 %), but **never validate on it alone** — its `jc(Bmin)` angle spread
  is only 0.84 %, so the structural weakness is nearly invisible there. Gate on sp-ap and
  superox.

### 8.7 The monotone guard does NOT fire at the operating point — F1 is smaller than feared

`d(log10 jc)/d(log10 B)` at the first B node, sp-ap, T = 77, all 181 angles: min −0.0722,
max −0.0104, mean −0.0378. **Negative at 181 of 181.** So the shipped table is internally
decreasing in B at its low end, `s` has the right sign, and `b = min(…, 0)` never fires there.
The peer's worry that the guard might be "silently doing the work" for sp-ap is refuted, and
Codex's F1 — while correct as mathematics — does not bite at the production operating point.
Still must be resolved honestly for the slices where it does fire.

### 8.8 A separate, larger table defect — NOT DR-128, routed to Christian

The peer found the shipped sp-ap table sits **~6.3 % high** against the raw rows it was built
from (77.5 K, 10 mT: mean +6.33 %, range +4.26 … +7.38 %, n = 57), roughly uniformly, which
reads as scale rather than shape. sst-1 is clean (−0.02 %); superox has +2.38 % mean but a
−2.68 … +6.49 % range, i.e. angular shape error — a *different* defect.

Their suspected origin is confirmed at least as intent, and it is written in the file itself —
`share/material/sp-ap.hdf5` `meta/jc_type`:

> "layer critical current density, jc = (Ic/w)/t_eff; t_eff from PCHIP-in-T match of
> jc(76 K, 0.01 T, 90 deg) to the previous table"

with `meta/t_eff_m = 9.8575e-07` back-derived from it. So sp-ap's absolute scale is
deliberately inherited from the *previous* table for deck compatibility; if that legacy anchor
was high, every absolute number in sp-ap inherits it. **Christian's call, not ours.**

**This does not propagate into the DR-128 design.** The transition is built from the table's
own value and slope at `Bmin`, so a uniform scale error multiplies both alike: `b` scales with
it and the ratio `jc(0)/jc(Bmin) = 1 − s/2` is exactly invariant (verified — rescaling by
1.063 leaves the predicted angle spread at 2.61 %, unchanged to the digit). Consequently the
peer's absolute "+5.9 % predicted jc(0)" for sp-ap is **mostly the table's scale offset, not
the form's error**; the scale-free spread number is the one that refutes the form. The two
must be attacked separately — fixing either does not fix the other.

### 8.9 CORRECTION to §8.8 — the sp-ap offset is the PROJECTOR, not the legacy anchor

**§8.8's attribution is withdrawn, and the error was partly mine.** The peer session proposed
the legacy `JC_OLD` deck-compat anchor as the origin; I endorsed it and added apparent weight
by quoting `meta/jc_type`, and I reported it onward to Christian in that form. Grok refuted it
and the refutation is correct: **the `JC_OLD` scale cancels in any table-vs-measurement
ratio**, because table jc and "measured jc" are both divided by the same `t_eff`. That ratio
therefore tests the interpolant against the data and is structurally blind to the anchor. The
`meta/jc_type` quote is real; the inference I drew from it was not. (It is the same
cancellation argument I had used on the peer one message earlier for my own form's
scale-invariance, which I should have recognised.)

**Actual cause: the peer's B-spline projection stage.** Table/measured at 77.5 K over
0.01–8 T, pre- vs post-projection: sp-ap 0.9920 → **1.0748**; sst-1 0.9937 → 0.9969; superox
0.9988 → 1.0255. Pre-projection all three sit on the measurements within ~1 %. The projection
moved sp-ap by +8.3 %, and it is not low-field-specific — it is the 70–85 K band, which is
merely where the 77.5 K measurements live. Independently corroborated by their earlier
per-band measurement: sp-ap 70–85 K median projection cost 0.0312 in log10 = 7.4 %.

Their reported remedy (Christian's call, not ours): refining **T** fixes it — sp-ap
table/measured at 77.5 K is 1.0748 at 2 K × 0.1, 1.0713 at 2 K × 0.05 (finer B, no help),
1.0023 at 1 K × 0.1. Caveat they state honestly: even at 1 K the p95 error is 9.7 % against
5.5 % pre-projection, so C1 still costs angular scatter, just far less at 1 K.

**Consequence for the development/validation split in §8.6:** the *reason* changes even though
the choice stands. sst-1 is not intrinsically clean — pre-projection all three tapes are within
1 %. sst-1 is the tape that SURVIVES the projection. If Christian rebuilds at 1 K that ranking
can move, so the split must be re-derived after any rebuild rather than inherited.

### 8.10 The refutation SURVIVES on clean data — verified, not inherited

The obvious hope was that §8.6's refutation was itself a projector artifact, since it was
measured on the shipped file in the band where the projector is worst. **It is not.** Fed the
pre-projection anchors directly (`profiles.npz`, sp-ap, T = 77.5 K, 181 angles):

| quantity | pre-projection | shipped table |
|---|---|---|
| spread of jc(10 mT) | 4.07 % | 4.82 % |
| **spread of predicted jc(0) (this form)** | **3.88 %** | 2.61 % |
| **spread of MEASURED jc(0)** | **0.36 %** | 0.33 % |
| slope `s` at Bmin | −0.1623 … −0.0084, positive at 0/181 | −0.0722 … −0.0104, 0/181 |
| mean predicted/measured jc(0) | **1.0150** | ~1.059 |

So on clean anchors the form still predicts **~11× too much angular spread**. The refutation is
intrinsic, not inherited. (Note the projection accidentally *reduced* the predicted spread,
3.88 → 2.61, while adding the +7.5 % level bias — so the shipped file flattered the form
slightly on the one axis where it fails.)

**The precise statement of the defect, now that scale and shape are separated:** the form is
right in MAGNITUDE and wrong in SHAPE. Its mean is only 1.5 % high on clean data — the peer's
+5.9 % on the shipped file was mostly the projection bias, exactly as the scale-invariance
argument in §8.8 predicted — but its angular shape is wrong by an order of magnitude.

**The physics this reveals, and the design requirement it sets.** Measured jc is ~4 %
angle-dependent at 10 mT and only ~0.36 % angle-dependent at B = 0. The anisotropy therefore
**collapses to nothing between 10 mT and zero field**, which is physically what one expects —
anisotropy is a statement about field *direction*, and at zero field there is no direction.
Any below-`Bmin` form must therefore drive the angular dependence to zero along with the
field, not carry `jc(Bmin)`'s angular shape down to the axis. Concretely the transition has to
act on the anisotropy and the isotropic level separately, e.g. an angle-averaged self-field
level plus an anisotropic part that switches on with B. That is the round-3 design, and the
guard question (F1) is moot for sp-ap either way: `s < 0` at 181/181 angles pre-projection too.

## 9. Round-3 design work (2026-08-29) — the anisotropy-collapsing form

Requirement from §8.10: the transition must drive the ANGULAR dependence to zero along with the
field, not carry `jc(Bmin)`'s angular shape down to the axis.

**Construction.** Introduce an angle-free zero-field level computed *from the table itself*, so
nothing new has to be shipped:

    iso(T) = < jc(Bmin, theta) * ( 1 - s(theta)/2 ) >_theta        ( s = dln jc / dln B at Bmin )

then interpolate from `iso` at B = 0 to the table's own value and slope at `Bmin`.

**Candidate A — cubic, no linear term:** `jc(B) = iso + b B^2 + a B^3`, with `b, a` fixed by
matching value and slope at `Bmin`. Satisfies: angle-free at B = 0 (by construction), C1 at
`Bmin`, and `djc/dB(0) = 0`.

| T [K] | spread of jc(0) | measured spread | mean predicted/measured | monotone |
|---|---|---|---|---|
| 70 | **0.000 %** | 0.43 % | 1.0058 | **no** |
| 75 | **0.000 %** | 0.39 % | 1.0079 | **no** |
| 77.5 | **0.000 %** | 0.36 % | 1.0150 | yes |
| 80 | **0.000 %** | 0.44 % | 1.0213 | yes |

The angular defect is **solved** — zero spread by construction, tighter than the measurement's
own 0.36–0.44 % scatter — and absolute accuracy is 0.6–2.1 % at the operating band. Monotonicity
is not guaranteed and fails at 70 and 75 K.

**Candidate B — power law:** `jc(B) = iso + (V − iso)·(B/Bmin)^p` with `p = S·Bmin/(V − iso)`.
Monotone by construction (`V < iso` holds at **all** temperatures tested) and mean accuracy is
better (1.0021 at 65 K … 1.0565 at 85 K). **Rejected:** `p` is determined pointwise and
degenerates where `V → iso`, ranging 0…12 across angles. `p < 1` gives an *infinite* `djc/dB` at
B = 0 — a diverging Newton tangent, which is worse than the defect being fixed.

### 9.1 The open question that now blocks the choice — and it is the one Codex asked for

Both candidates impose `djc/dB(0) = 0`. **That may be physically wrong.** A Kim form,
`jc = jc0/(1 + |B|/B0)`, has `djc/d|B|(0) = −jc0/B0`, i.e. finite and NON-zero. If the true
behaviour is a finite-slope cusp rather than a smooth maximum, then requiring zero slope
over-constrains the fit — and dropping it frees a degree of freedom that would likely buy
monotonicity back (a quadratic `iso + cB + bB²` matches `iso`, value and slope with no
zero-slope condition).

So the literature route CODEX-6 demanded is no longer optional or cosmetic: **it now selects the
functional form.** It must be answered before the form is fixed:

- [ ] **R4** — route the low-field jc(|B|) question: for a REBCO coated conductor, is
      `djc/d|B|` at B = 0 zero (smooth maximum) or finite (Kim-like cusp)? Kim 1962 and the
      Bi-2223/REBCO angular-dependence literature are the starting points; note `doc/literature_references.md`
      records that a Kim `B_eff` form structurally cannot fit the Bi-2223 angular data, so do not
      assume the REBCO answer transfers.
- [ ] **R5** — with that settled, re-choose between Candidate A, a zero-slope-free quadratic, and
      a limited variant; re-run the four-column table above plus a monotonicity gate.

**Status: the angular refutation of §8.6 is SOLVED in principle; the form is not yet fixed.**
DR-128 cannot close until R4 and R5 land, the peer session releases
`cl_JcFunction_Database.hpp`, and the F2 chain-rule work (§8.1) is written.

## 10. LITERATURE ROUTE (R4) ANSWERED, and it changes the design

**Question:** for a REBCO coated conductor, is `djc/d|B|` at B = 0 zero (smooth maximum) or
finite (Kim-like cusp)?

**Answer: FINITE and NON-ZERO.** Routed via `literature/papers/fem/index.md:302` ("Kim model
for J_c(B)"), which points at Messe et al. 2023 and Riva et al. 2023. Three sources, one
structure:

- **Denis et al. 2026, Eq. 11:** `jc(b) = jc0 / ( 1 + sqrt(kc² bx² + by²)/b0 )^α`
- **Riva et al. 2023, Eq. 2:** `Jc(B∥,B⊥) = Jc,0 / ( 1 + sqrt(k|B∥|² + |B⊥|²)/Bc )^b`
- **Messe et al. 2023 §2.6** names Kim 1962 as the model BELFEM intends to implement,
  "frequently used in the HTS community"

For the Kim family, with `u = |B_eff|`, `d jc/du` at `u = 0` is `−jc0·α/B0` — finite and
non-zero. **So the zero-slope condition imposed by both §9 candidates is physically wrong**, and
Christian's instinct for a quadratic that keeps its linear term is correct.

**The same form independently explains the angular constraint of §8.10.** The anisotropy enters
only through the angular factor multiplying B *inside* the norm, so as B → 0 it vanishes and
`jc(0) = jc0` is angle-free — exactly the measured 0.36 % spread. The literature form has
precisely the structure the data demands, which is a strong corroboration of both.

### 10.1 What was tested, and the negative result that matters

- [x] **Candidate C — Christian's quadratic** `jc(B) = iso + cB + bB²` (linear term kept),
      `iso` fixed and `c, b` from value+slope at `Bmin`. Angle-free at B = 0 by construction.
      **Fails on monotonicity, but for a reason that is not the form's fault:** with `iso` taken
      as an angle *average*, `iso < jc(Bmin,θ)` at the angles above the mean, so the curve must
      RISE from `iso` to `V(θ)` there — `c` comes out positive at those angles (measured range
      −9.2e10 … +9.8e10, angle-spread 190–400 %). The estimator is wrong, not the polynomial.
- [x] **Candidate D — fit the Kim exponent** so that `jc0 = V·(1+x)^α`, `x = −s/(α+s)`, is
      angle-free; one unknown, 181 equations. **Rejected.** The fit runs to the bound at every
      temperature and the residual angular spread of `jc0` is 1.7 % (65 K) → 3.6 % (77 K) →
      9.9 % (85 K), an order of magnitude worse than the 0.36 % target.

**The negative result is the useful one.** Any construction that *extrapolates* the zero-field
level from `(V, s)` at `Bmin` must cancel a ~4 % angular variation in `V` against a ~200 %
angular variation in `s` to land inside 0.36 %. That is ill-conditioned, and no one-parameter
family fixes it. **The zero-field level must be MEASURED, not extrapolated.**

### 10.2 The design this settles on

The measurement already exists and is already reduced — `profiles.npz` carries `SWEEP_T` /
`SWEEP_ICW`, a 34-point self-field Ic/w curve from 13 K to 95 K (77.5 K: 38074 A/m, matching the
shipped `meta/Icw_77p5K_sf_A_per_m` = 37995.8 to 0.2 %). So:

    jc(0,T)   = measured self-field level          ( from file metadata, angle-free BY DATA )
    jc(B,θ)   = jc(0,T) + c(θ) B + b(θ) B²         ( c, b from the table's own V and s at Bmin )

- angle-free at B = 0 **by measurement**, not by construction or by an estimator;
- finite non-zero `djc/dB(0) = c(θ)`, Kim-consistent, angle-dependent as the literature requires;
- C1 at `Bmin` by construction;
- monotone iff `jc(0) ≥ jc(Bmin,θ)` for all θ, which the measurements satisfy (77.5 K ratios
  1.0077 … 1.0326 across angles).

- [ ] **D9** — the self-field curve must be readable from the table file. The scalar
      `meta/Icw_77p5K_sf_A_per_m` exists on all seven REBCO tables but is ONE temperature; the
      full `SWEEP_T`/`SWEEP_ICW` curve is needed. **This is a file-format addition and therefore
      couples to the lift-factor session**, which is already touching metadata and referencing
      `L` to measured self-field — the curve should land once, in their format change, not twice.
- [ ] **D10** — implement C1 quadratic below `Bmin` in `cl_JcFunction_Database.hpp` with the
      F2 chain-rule terms (§8.1); `n` treated independently (F3).
- [ ] **G6** — monotonicity assert: refuse (or clamp with a loud warning) if the stored
      self-field level is below `jc(Bmin,θ)` at any angle, which would indicate an inconsistent
      table rather than a physical case.

**DR-128 is now DESIGNED and literature-grounded, and blocked only on D9** — the self-field
curve reaching the file. It is no longer blocked on a physics question.
