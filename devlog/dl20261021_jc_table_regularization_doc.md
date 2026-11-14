# Jc table regularization: why the tables are monotone and C1, and how to edit one

**Date:** 2026-10-21
**Purpose:** An email thread (Yufan, Giorgio, Christian) asked whether the rules behind the
shipped Jc tables are documented. The thread came from a user's old `sst-1.hdf5`, which the
lowest-field slope check refused; Giorgio flattened the rising part to get it to run; Yufan noted
that Jc can genuinely rise with B. Christian: the rules regularize noisy Robinson data, and
well-behaved material functions plus the B-spline projection are a significant part of solver
stability; an edited table must be re-projected or Newton dies during quench.
**Module:** `share/material`, `src/physics/materials/doc`
**AIs involved:** Claude Opus 5.5 (primary); Codex `gpt-5.6-luna`/medium (two language sweeps);
Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high (blind documentation jury)
**Status:** COMPLETE — documentation only, no code change
**Literature References:** Wee et al. 2013 (DOI 10.1038/srep02310, Crossref-verified, full text
read via Europe PMC PMC3741626)
**Verification:** probe — the HDF5 kit-extraction snippet, `project_table.py sst-1.hdf5 -o …` and
`verify_tables.py` on a scratch copy ("all tables pass every gate"); the prose is reviewed, not
verified. `check_citations.py`: 514 resolved, 0 unresolved, after adding `Linux` to
`NOT_A_SURNAME`
**Exchange:** `tmp/ai_exchange/review_jc_table_regularization_doc.md`,
`tmp/ai_exchange/jc_table_regularization_doc.md`

## Summary

What was already documented: `share/material/README.md` "What a deck may assume" stated the
contracts (jc non-increasing in T and |B|; n non-increasing in |B|, bounded; C1 from the
B-spline projection). What was missing: the reason, the fact that real conductors need not
obey the jc rule, the load-time refusal, and how to re-project an edited table. Two new
subsections were added to `share/material/README.md` ("Why the tables are monotone and smooth",
"Using an older or modified table"), together with a pointer in `materials_usage_guide.md` §3.4 and a `Wee 2013`
row in the `material_property_sources.md` bibliography.

## What the jury found (all confirmed against source; full record in the exchange)

- **P1, Codex + Grok.** The first guide pointer claimed every table is non-increasing in T. REBCO
  `n` has no dn/dT contract, by decision (`share/material/README.md:55-57`;
  `dl20260904_rebco_tables_n_monotone_b16.md:121`: in-plane n(T) rises 0.2–0.26 dex between
  55 and 70 K on all five tapes). Fixed.
- **P1, Codex.** The re-projection recipe omitted `--mirror`, which the bscco-2223 table needs
  (embedded `project_table.py` docstring). Fixed: separate REBCO and Bi-2223 commands.
- **P1, Claude + Grok.** The noise example was about `n` but sat in a `jc` argument, and the
  two operators differ (isotonic L2 fit for jc, running minimum along log10 B for n).
- **P1, Grok.** "Acceptable when the change is smaller than the noise" holds for jc in the
  measured band, but not for the n gate, which removes one-sided outliers much larger than the
  scatter (sp-ap n 13 → 58 from one bad fit; `dl20260904:16-29`).
- **P1, Grok.** The draft conflated an interior rise (passes `open_checked`, which bounds only the
  lowest log10 B face, `cl_JcFunction_Database.hpp:190-199`) with a rise at Bmin (refused).
- **P2.** C1 holds for the interpolant only: the jc continuation at beta = 2 has discontinuous T and
  angle tangents (`cl_JcFunction_Database.hpp:552`), and the n floor is C0
  (`resistivity_laws.md:188`). The check also fires on the clamp path, `continued()` at
  `cl_JcFunction_Database.hpp:693` (bscco). The Bernstein bound is a conservative certificate
  (`dl20260923_jc_c1_continuation.md:130`). `howtoread.py` clamps below Bmin. `project_table.py`
  has no switch for the jc constraint. "Usual cause" was unsupported. The bscco file has no
  "## Construction" heading (h5py probe).
- **Newton argument, Claude + Grok.** The power-law tangent contains djc/dB and djc/dT
  (`powerlaws.hpp:1986`, `cl_FEM_Calculator.hpp:3153ff`). The Bézier-blend field derivative is
  staged 0 (`powerlaws.hpp:2077-2079`). The only stall the code documents is the omitted
  knot motion (`powerlaws.hpp:2142`). The doc now states the element-jump stall as experience,
  not isolated in a reproducer.

## The physics question: can Jc rise with B?

`./literature` has nothing on it (only the Kim model and general magnet texts). Wee et al. 2013
(IBAD-MgO templates, H∥c, 77 K) show that BZO nanocolumns flatten Jc(H) (alpha 0.48 → ~0.33) and
that the optimal doping follows the matching field. They show **no rise**: alpha > 0 in every
sample. A first wording, "need not decrease", overclaimed and was cut. The doc cites Wee only
for the flattening and the matching field. It states that the tables cannot represent a rise and
that none of the five REBCO tapes shows one. Yufan's "a rise occurs in some conductors" stays
uncited and unasserted. Confidence medium that such rises exist in research samples
(search summaries only, a SmBCO+BZO plateau upturn, not checked in full text).

## Observation from the probe

Re-projecting the shipped `sst-1.hdf5` moves the n interpolant maximum from 35.31 to 36.05, so
the projection is not idempotent on n. An edited and re-projected table will therefore not
reproduce the shipped values exactly, even where it was not edited. Harmless for jc; not
investigated further.

## Open

- Whether to add a `--no-jc-monotone` style switch to `project_table.py` for conductors with
  a genuine rise; Christian's call, and the embedded copies cannot change without new files.

## Files Updated

- `share/material/README.md` — two new subsections
- `src/physics/materials/doc/materials_usage_guide.md` — §3.4 pointer
- `src/physics/materials/doc/material_property_sources.md` — `Wee 2013` bibliography row
- `scripts/check_citations.py` — `Linux` added to `NOT_A_SURNAME`. "Amazon Linux 2023"
  (`doc/getting_started.md:44`, `doc/scls_toolchain.md:63`, from a03526b) is a product name; the
  checker had read it as an author-year citation and failed on it
