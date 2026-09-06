# ShiftRegister: memory-safe support for owning element types (Vector/Matrix)

**Date:** 2026-07-14
**Purpose:** Make `ShiftRegister<T>` genuinely safe for owning element types
(`Vector<T>`, `Matrix<T>`) instead of relying on a false "malloc gives valid
empty objects" premise. Audited by Codex + Grok.
**Module:** containers

## Trigger

Christian asked whether `ShiftRegister` (`src/containers/cl_ShiftRegister.hpp`)
can hold `belfem::Vector`. It had an `is_shift_register_safe` specialization
whitelisting `Vector<T>`/`Matrix<T>`, but the storage model did not back that up.

## Problem (why the original was wrong, not just misleading)

`ShiftRegister` keeps a `malloc`'d buffer of `capacity+1` slots (+1 = revert
backup). It was written for trivially copyable `T`: `malloc` leaves slots
uninitialized, elements are populated via `operator=`, and `free` runs no
destructors. For `Vector<T>` — which wraps a non-trivial, heap-owning
`arma::Mat<T>` (`src/linalg/armadillo/cl_AR_Vector.hpp:31,38`) — this was UB:

- `malloc` slots hold garbage; `mData[0] = aValue` invoked `arma::Mat::operator=`
  which reads the destination's existing state → free of a wild pointer.
- `std::move_backward`/`std::move` move-assigned into garbage slots → UB.
- destructor/`free`/`reserve` freed without `~T()` → every element's heap buffer
  leaked.

The trait comment claimed the backends "default-construct to a valid empty
state that is equivalent to zeroed memory" — true of a *default-constructed*
Vector, but the code never default-constructed the slots, so the property was
real yet irrelevant. Real consumer: `src/numerics/ode/cl_BDF.cpp` uses
`ShiftRegister<Vector<real>>` for BDF time-stepping history, so the UB was live.

## Fix

Model: for owning `T`, the buffer's `capacity+1` slots are **always live,
constructed objects** for the whole lifetime of the buffer.

- `construct_slots()` — placement-`new T()` into every slot right after `malloc`
  in `reserve()`.
- `destroy_slots()` — `~T()` every slot right before `std::free` (reserve
  teardown, `free()`, destructor, move-assignment).
- Both guarded by `if constexpr ( !std::is_trivially_copyable<T>::value )`, so
  the scalar fast path compiles to byte-identical code (zero overhead).

With every slot always valid, the existing element ops (`push`'s
`move_backward`, `mData[0]=aValue`, `std::copy`, `revert`'s `std::move`, `fill`)
act on constructed objects — exactly the contract `arma::Mat` assign/move
expect; a moved-from `arma::Mat` is a valid empty object.

Also:
- `static_assert` that an owning (non-trivially-copyable) `T` is
  default-constructible.
- Corrected the trait comment to describe the two storage regimes honestly.
- Removed the dead `if (mData != nullptr) std::free` branch in the move ctor
  (mData is always null in a ctor; it would have leaked owning slots if reached).
- Added a `BELFEM_ERROR` `malloc`-failure check in `reserve()` (construct_slots
  now dereferences immediately, so OOM must be caught, per
  `doc/coding_philosophy.md` allocation-failure rule).
- Added `BELFEM_ERROR( aCapacity > 0, ... )` at the top of `reserve()` — the
  single choke point every constructor routes through — closing the pre-existing
  zero-capacity `push`-through-`nullptr` hole (Christian's call, this session).

Tests (`tests/containers/test_ShiftRegister.cpp`): fixed the stale header
comment that declared owning types unsafe; added a `ShiftRegisterOwning`
section over `Vector<real>` (push/shift, push-beyond-capacity + backup shift,
revert, deep-copy independence for copy-ctor/assign, move transfer, reserve
teardown, fill, free) — all `[valgrind]` candidates.

## Audit (Codex + Grok, `tmp/ai_exchange/shiftregister_owning_types*.md`)

Both independently approved the lifetime model for the positive-capacity,
successful-`malloc` path (the BDF use case): no use-after-free, double-construct,
double-destruct, leak, or non-live-slot assignment. Both verified the three
tricky spots — `reserve()` destroy-old/construct-new ordering (incl.
same-capacity early-return not double-constructing), copy touching only
`[0, mSize)` while leaving `[mSize, mCapacity]` default-constructed, and the
full-buffer `push` shift through the backup slot. Confidence: high.

## Deferred (pre-existing, not owning-type-specific)

- **Aliasing `push(reg(i))`:** logical bug, not lifetime; BDF doesn't do it.
- **Use-after-`free()` `push`:** calling `push()` on a manually `free()`d
  register still derefs `nullptr`; treated as caller misuse (like using any
  freed resource), not guarded. The zero-capacity *construction* path is now
  closed by the `reserve()` guard (see fix list).
- Over-aligned trivially-copyable `T` (`alignas(64)`) vs plain `malloc`;
  constructor exception-unwind under exceptions-enabled builds — both acceptable
  under the project's `-fno-exceptions` + current `Vector`/`Matrix` use.

## Files touched

- `src/containers/cl_ShiftRegister.hpp` — lifetime management, static_assert,
  malloc check, move-ctor cleanup, comment fix.
- `tests/containers/test_ShiftRegister.cpp` — header comment, owning-type tests.

Build/run handed to Christian (not run here).
