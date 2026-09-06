# Devlog: B4a Quaternion and Tensor Whitepaper

## Context

Read-only documentation pass for Tier B4a math helpers, split from the broader B4 "other features" bucket. Scope was limited to quaternion helpers and fourth-order tensor helpers.

## Changes

- Added `tmp/whitepaper/B4a_quaternions_tensors.md`.
- Confirmed quaternion code lives under `src/math/quaternion`, not `src/numerics` or `src/linalg`.
- Confirmed fourth-order tensor code lives under `src/math/tensor`, not `src/numerics` or `src/linalg`.
- Documented quaternion storage and rotation convention from source: scalar-first `(w,x,y,z)`, Hamilton multiplication, and implemented `q * v * q*` vector rotation.
- Documented tensor representation from source: full fourth-order storage internally, 6x6 Voigt conversion at the API boundary, engineering shear convention rather than Mandel scaling, and symmetry assumptions at conversion/inversion boundaries rather than raw storage enforcement.
- Recorded implemented-vs-tested coverage and open questions, including stale quaternion borrowed-buffer tests, the apparent stale FVM three-argument `Tensor` construction, and the `rotate42` comment/index wording.

## Verification

- Rechecked path citations for cited source line ranges.
- Checked the whitepaper is ASCII-only.
- Did not run C++ tests; this was a source-reading documentation task.

## Follow-up (same task, user-approved source edits)

Independent re-read by Claude (the B4a draft had been mis-routed to Codex). Claude
confirmed every substantive convention fact against the same files and firmed up
two items Codex left open. Then, with user approval:

### Comment fix
- Corrected the `rotate42` doc comment in both backends from
  `A_mnop = B_ijkl * R_im * R_jn * R_ko * R_np` to `... * R_lp`. The summation
  index `l` of `B_ijkl` was missing and `n` was duplicated; the kernel itself is
  correct (tested against a reference loop), only the comment was wrong.
  Files: `src/math/tensor/armadillo/fn_TR_rotate42_arma.hpp:24`,
  `src/math/tensor/blaze/fn_TR_rotate42_blaze.hpp:24`.

### New usage guides (module docs, following documentation_guidelines.md)
- `src/math/quaternion/doc/README.md` + `quaternion_usage_guide.md`
- `src/math/tensor/doc/README.md` + `tensor_usage_guide.md`
- Both guides lead with the convention/notation table (the silent-bug sources),
  then API quick reference, common patterns, numerical aspects, and status.
- Recorded that quaternions have no production user yet (test-only), and that the
  tensor 6x6 boundary is engineering-Voigt stiffness (no Mandel factors), with
  `invert_symmetric` doing the factor-of-2 shear bookkeeping internally.

### Not changed (flagged, awaiting decision)
- Stale quaternion borrowed-buffer tests in `tests/math/test_Quaternion.cpp`
  (describe an obsolete malloc/free storage model) -- left for a test refresh.
- `src/fvm/cl_FVM_Factory.cpp` three-arg `Tensor` construction -- stale/unbuilt
  (`fvm` is not added in `src/CMakeLists.txt`); author to confirm intent.

## Verification (follow-up)

- Confirmed both `rotate42` comments now read `R_lp`.
- Confirmed the four new doc files are ASCII-only.
- Whitepaper `tmp/whitepaper/B4a_quaternions_tensors.md` updated to reflect the
  two-AI process (Codex drafted, Claude independently verified) and the resolved
  items.
- No C++ compiled or tests run; comment-only source edit plus new Markdown docs.
