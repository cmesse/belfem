# Homology: three items for Gregory from the test-gap run

**Date:** 2026-09-22
**Purpose:** Three findings in `src/homology` that the overnight test-gap run could not act on
because the module is Gregory's and its cohomology core is closed to AI edits. Christian asked
for them in one note so he can forward it. Everything here is reviewed by reading and by the
tests named; nothing in the core was changed.
**Module:** `src/homology`

## 1. The 2-D suggested generator's sense versus the core's own comment (O10)

`tests/homology/test_Cohomology.cpp` (`5834e8b`) builds the annulus with a fan conductor, runs the
suggest path the way `CutFactory::compute_cohomologies` does, and reads the sign of the suggested
1-chain against the counterclockwise inner ring. Two auditors and the plan agent derived `+1`
from the code (suggest adds the fan boundary with `-1`, `reorient_generators` multiplies by `-1`
again), and the first run in both builds gave `+1`. The test pins that as a named calibration
constant, `gInnerRingSense = 1`, with the comment that it is what the tree produces, not a
physical convention.

The comment in `cl_Homology.cpp:667-668` says the generator must run **clockwise** for a
positive current along `+z`, and the paragraph after it says the measurements, not that sentence,
justify the unconditional `-1`. The 2026-08-28 devlog on the current sign records the same split:
the Ampère ratio went from about `-0.99` to about `+1` after the extra minus while the prose still
called the desired loop clockwise. The question for you: which sentence is the convention? If the
comment is the stale one, the test's constant stays `+1` and the comment changes; if the comment
is right, `reorient_generators` is off by a sign in 2-D and the test goes to `-1`. The test cannot
decide that by itself: it never builds a cut or `h = -grad φ`, so it sees no physical current.

## 2. `clean()`'s unit guarantee, and two small things beside it (O18)

The T-44 row added a refusal in the open file `cl_CutData.cpp`: `collect_coefficients` now
`BELFEM_ERROR`s on any cochain coefficient outside `{-1, 0, +1}`, naming the cut and the edge,
where before a `±2` silently became weight 0 with the edge still flagged (in a 2-D release run
that was a cut of the wrong sign with no message). The consumer sweep found that
`Cohomology::clean()` already carries two release-active unit gates after the one
non-unimodular producer (`cl_Cohomology.cpp:605-613` and `:883-890`, after Step 7 of
`updatekGeneratorsFromHomology`), so the new refusal is a third line of defence and no shipped
deck can reach it. What is missing is the contract: `clean()`'s header (`cl_Cohomology.hpp:98`)
does not say it guarantees units, and `CutData` documents no precondition. A sentence next to
`clean()` and a `@pre` on `CutData::collect_coefficients` would make the backstop a documented
one. Comment-only.

Two small things seen on the way: the new error says "Cut N" with `mID = index + 1` while
`clean_spfa` says "generator N" 0-based, so the same object is numbered two ways in two
messages; and `CutProcessor::compute_edge_bitsets()` fills `mCohomologyEdgesPlus/Minus`
(`cl_CutProcessor.cpp:219`, `:223`) that nothing reads — two bitset allocations per cohomology
edge for no consumer (grep-based, medium-high confidence).

## 3. `Chain::getSimplicesMap` is declared with the wrong key type (DR-177)

`cl_Chain.hpp:107-108` declares `OrderedMap< id_t, int > & getSimplicesMap()` and the inline body
at `:157-160` returns the member `OrderedMap< index_t, int > mSimplicesMap` (`:44`). The two
typedefs coincide on the default build (`uint32_t` and `unsigned int`), so it compiles; under
`BELFEM_INT64` (`USE_MKL_64BIT_API`, `CMakeLists.txt:155`) `index_t` widens to `uint64_t` and the
return is ill-formed. `Cochain::operator()` (`cl_Cochain.cpp:76`) and the new tests call it. One
word in the declaration, in a closed file, so it is yours.

## What was not touched

`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`, `cl_Chain`, `cl_Cochain`, `fn_Smith`:
read only. The T-44 change is in `cl_CutData.cpp` and its module note
`src/homology/doc/thin_cut_nonunit_rectification.md`, both open files, both juried.
