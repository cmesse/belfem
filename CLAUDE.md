# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) and OpenAI Codex when working with code in this repository. It is the session bootstrap: conventions, coding standards, and literature routing in the form a session needs them up front. It is not the top of the authority stack — `doc/ai_collaboration_protocol.md` governs collaboration process and `doc/coding_philosophy.md` carries the full HPC rationale; where either is more specific, it wins (see Instruction Precedence).

## AI Cooperation

This repository uses a two-AI model: **Claude Code** (primary, broad exploration) and **Codex** (secondary, precision audit). Cross-review rounds add Grok as a third voice. The detailed collaboration protocol — exchange format, audit checklist, confidence calibration, and devlog conventions — lives in:

**`doc/ai_collaboration_protocol.md`** — read this at the start of every session.

Key points:
- **Calibrated uncertainty:** attach confidence (high / medium / low) to non-trivial claims
- **Exchange channel:** `./tmp/ai_exchange/<slug>.md` (AI-only, ephemeral, per-task AI-to-AI communication; distilled into the devlog before it is swept)
- **Nonfree exception:** for work scoped to `./nonfree/` — or porting from `./tmp/manta/` into it — do **not** create, read, or append the exchange at all. It is routed through external AI vendors and is a leak surface for proprietary work; keep coordination in chat and record outcomes in `./nonfree/devlog/`.
- **Session records:** `./devlog/dlYYYYMMDD_topic.md` (persistent devlogs)
- **Edit safety:** read-only by default; source edits only with explicit user approval
- **Literature-first:** mandatory for algorithms/formulations/validation; not for trivial fixes

For role definitions and edit safety rules, see `AGENTS.md`.

### Cross-Review Tooling

```
/cross-review [--jury|--relay] [path]      # three-AI review round (default: working-tree diff; --jury = parallel blind, --relay = sequential)
                                           # needs an explicit depth: CODEX_MODEL/CODEX_EFFORT/GROK_MODEL/GROK_EFFORT
scripts/cross_review.sh --quick <sha>      # single-auditor light review (used by the post-commit hook; pins its own cheap depth)
scripts/install_autoreview_hook.sh [--uninstall]   # opt-in post-commit hook, all worktrees; enable with BELFEM_AUTOREVIEW=1
scripts/review_status.sh                   # commits since last auto-review + open P0 flags
scripts/check_doc_claims.py [--verbose]    # verify this file's checkable claims against the tree
scripts/check_wrapper_policy.py            # L-21: direct third-party calls outside their wrapper module
```

`check_doc_claims.py` guards the class of error a reader cannot catch: a convention document
asserting a build flag, make target, executable name, module build status, CMake default, or the
existence of a named symbol that the tree contradicts. **Run it after editing this file or
`doc/coding_philosophy.md`, and after changing `CMakeLists.txt` or the compiler config.** Prose and
rationale are out of its scope — those need the cross-review round.

Findings land in `tmp/ai_exchange/` (reviews per-slug; auto-reviews in `autoreview.md`, P0s flagged in `AUTOREVIEW_P0.flag`). See `METHODOLOGY.md` for the method in published vocabulary.

### Instruction Precedence

If guidance conflicts across collaboration docs, apply this order:

1. `doc/ai_collaboration_protocol.md`
2. `AGENTS.md`
3. `CLAUDE.md`

### Minimum Session Compliance Checklist

Before substantial work, confirm all of the following:

- Read `AGENTS.md` and `doc/ai_collaboration_protocol.md`
- Read **Layer 1** of `doc/lessons_learned.md` (the tripwire layer, ~220 lines): match its
  trigger index against what you are about to do, and its failure-signature table against any
  symptom you are chasing. Layer 2 holds the evidence — read a card only when Layer 1 sends you
  to it
- Decide whether literature-first is required for the current task
- Use `./tmp/ai_exchange/<slug>.md` (AI-only, ephemeral) for AI-to-AI findings and audit threads — **except** for `./nonfree/` work, where the exchange tier is suppressed entirely
- Writing to `./todo/` and `./devlog/` is always allowed, even under the read-only default
- State confidence (high / medium / low) for non-trivial claims
- Say **"reviewed"**, not "verified", for static-only work — "verified" requires an executable gate that actually ran (evidence ladder, protocol §11); AI reviewer agreement is the weakest rung, not proof
- Keep investigations read-only unless the user explicitly approves edits
- Never edit the cohomology core (`cl_Cohomology`, `cl_Homology`, `cl_SimplicialComplex`,
  `cl_Chain`, `cl_Cochain`, `fn_Smith`) — closed to AI regardless of session edit approval;
  read and report to Gregory instead (`doc/ai_collaboration_protocol.md` §7.1)
- If the session touches any `input.conf` key, update **both** `doc/input_file_reference.md` and
  `doc/input_schema.yaml` (see §"The Input Contract: Two Artifacts, One Rule")
- If the session adds or revises **user-facing** documentation, run a **Codex language sweep**
  over it before treating it as finished. Devlogs, `METHODOLOGY.md` and `doc/lessons_learned.md`
  are excluded — the sweep dilutes, and their density is the point
  (see §"Prose Gets a Language Sweep")
- If session is meaningful, create `./devlog/dlYYYYMMDD_topic.md` and update `./devlog/README.md`

## Build System

This project uses CMake. Build trees are created by the user (`cmake-build-debug/` and `build/` both exist and are shared across worktrees).

**The user runs builds.** Do not invoke `make` or launch solver binaries unless explicitly asked — a partial or missing object file usually means a build is already in progress elsewhere in the shared tree.

### Common Build Commands

- **Build project**: `make` (or `make -j<cores>`) from the build directory
- **Clean build**: `make reset` (removes `bin lib test` and recreates them)
- **Run tests**: `make check` (all tests) or `make check-fast` (fast-labeled subset, whole set <~3 min).
  Tests are **on by default** (`USE_TEST=ON`, release policy since 2026-08-14). In a tree
  configured with `-DUSE_TEST=OFF`, `check-fast` does not exist, `check` prints a
  "re-run cmake with -DUSE_TEST=ON" error, and `test` is overridden to point at `make check`.
- **Generate documentation**: `make doc` (requires Doxygen)

The generator is pinned: CMake aborts on anything but Unix Makefiles unless `-DALLOW_NINJA=ON`.

### Build Configuration

The generator is not the only pinned part of the environment: **Open MPI is the only
supported MPI.** MPICH and Intel MPI (MPICH-derived) are untested — PETSc has been observed
to crash on that path — so `config/system/find_mpi.cmake` aborts the configure when it
detects them, overridable with `-DALLOW_UNTESTED_MPI=ON`. Using the Intel compilers is
fine; build Open MPI with them rather than substituting Intel MPI. The Open MPI library
names hardcoded in the MUMPS and MKL configs are correct, not a portability gap — see
`doc/mpi_support.md`.

The project has extensive configuration options in `CMakeLists.txt`:
- MPI support (USE_MPI=ON by default); OpenMP flags are enabled (USE_OPENMP=ON) because the third-party solvers need them. BELFEM's **own** parallelism is MPI: the C++ sources carry no `omp` pragmas, and the `!$omp` directives in `src/sparse/splinalg.f90`, `arpacktools.f90` and `parpacktools.f90` are gated on `BELFEM_OMP` (`USE_BELFEM_OPENMP`, **OFF** by default — DR-155). The plain `OMP` define is a different thing and still guards the third-party thread queries (`hatch_turtle()`, the banner's thread count, PARDISO's `gParameters(3)`); do not unify the two. See `doc/parallel_execution.md`
- Matrix libraries: Armadillo (default everywhere except Apple) or Blaze (default on Apple)
- Linear algebra: SuperLU, MUMPS, STRUMPACK, PETSc, ARPACK on by default; MKL, SuiteSparse, PARDISO opt-in
- I/O: HDF5, Exodus formats
- Optional: Maxwell modules (ON), gas models and combustion (OFF), tests (ON), coding examples (OFF)

## Architecture

BELFEM is a finite element framework with modular architecture. All paths below are relative to `src/`:

### Core Libraries
- **core/**: Fundamental utilities (Logger, Timer, Arguments, etc.)
- **comm/**: MPI communication layer
- **containers/**: Data structures (Bitset, Map, Queue, etc.)
- **linalg/**: Linear algebra abstraction (supports Blaze/Armadillo backends)
- **sparse/**: Sparse matrix solvers and interfaces
- **mesh/**: Mesh data structures and I/O (supports Gmsh, HDF5, Exodus, VTK)
- **io/**: File I/O utilities (HDF5, XML, ASCII)

### FEM Components
- **fem/kernel/**: Core FEM classes (DofManager, Element, Domain, etc.)
- **fem/iwg/**: Integral weak forms for different physics
- **fem/maxwell/**: Electromagnetic field solvers
- **fem/thermal/**: Thermal physics and boundary conditions
- **fem/interpolation/**: Shape functions and integration
- **fem/postproc/**: Post-processing utilities

### Specialized Modules
- **circuit/**: Lumped electrical circuit coupling
- **homology/**: Topological analysis and cut algorithms — its **cohomology core is closed to AI edits** (protocol §7.1)
- **math/**: Mathematical utilities (graph algorithms, tensors, quaternions)
- **numerics/**: Numerical methods (Bezier curves, splines, integration, ODE, source functions, and an NLOPT-backed `opt` optimizer)
- **physics/**: Material, gas, and database property definitions
- **visualizer/**: Rendering utilities — built only under `USE_VTK`, which defaults OFF

### Applications
- **executables/**: Solver application (`belfem`); `hphirun` and `hphiTrun` are retired and no longer built

The `nonfree/` directory is a separate, non-open-source repository (boundary layer, combustion, and other proprietary extensions) that plugs into the build when present. Its devlogs belong in `./nonfree/devlog`, never in the open-source `./devlog`.

## Key Design Patterns

- Template-based linear algebra with backend abstraction
- Factory pattern for element creation and solver selection
- Modular physics through IWG (Integral Weak Form) classes
- Parallel computing via MPI with distributed matrix operations
- Mesh-agnostic design supporting multiple file formats

## Coding Standards and Philosophy

**CRITICAL:** BELFEM is an HPC framework optimized for performance, not modern C++ safety conventions. Read `doc/coding_philosophy.md` to understand the deliberate design choices.

### Core Principles

**Performance First, Safety Through Testing**
- Release builds (`USE_DEBUG=OFF`): `-O2 -DNDEBUG` (GCC Fortran additionally `-march=native` on Linux x86, `-mtune=native` on Apple x86; ICC uses `-O1 -xHost` instead)
- Debug builds (`USE_DEBUG=ON`): `-Og -g` plus the full assertion set (GCC Fortran `-O0 -fcheck=bounds -fbacktrace`, ICC `-O0 -g`)
- **Which one you get by default depends on the SCLS flavor** (`config/system/find_scls_flavor.cmake`,
  read before the `option()` calls). With `$SCLS` set, `$SCLS_FLAVOR` supplies the defaults:
  a flavor whose name contains `mkl` → `USE_MKL=ON`, `USE_PARDISO=ON`; `debug` → `USE_DEBUG=ON`;
  every other flavor, and no `$SCLS` at all → all three OFF. `USE_PARDISO` always defaults in
  parity with `USE_MKL`. Nothing is checked or refused: an explicit `-D` or a ccmake edit always
  wins, and every vendor other than MKL is assumed to behave like reference BLAS/LAPACK.
  With no `$SCLS`, or a prefix whose flavor is not one of those three, both default OFF and the
  build behaves like any ordinary Linux build. These are *defaults only* — an explicit
  `-DUSE_DEBUG=…` / `-DUSE_MKL=…` always wins, and a value already in a build tree's cache is never
  overwritten (so repointing `$SCLS` under an existing tree does not change its knobs; the
  configure summary prints the flavor so the mismatch is visible)
- `CMAKE_BUILD_TYPE` is derived from `USE_DEBUG` (Debug/Release; a conflicting
  command-line value is overridden with a notice), and the per-config flag sets are
  blanked — the hand-set flags above are the only optimization source
- Zero abstraction penalty: code compiles to same machine instructions as hand-written C

**Honest testing posture:** since 2026-08-30 a GitLab CI server runs `make check` nightly in a
debug/release × Armadillo/Blaze matrix, so `make check` on a workstation is no longer the only thing
that ever executes the suite. `make check` includes the Tier 2 MPI suites (`commmpi`, `sparsempi`),
each registered at 2 and 4 ranks; the runner has two cores, so the 4-rank cases run oversubscribed,
which is a correctness run and not a timing one. Four limits still hold. The pipeline definition,
`.gitlab-ci.yml` at the repository root, runs only on a scheduled or manually started pipeline and on
a runner tagged `belfem-local`, so a clone carries the definition but cannot reproduce the run. The
job configures without `-DUSE_GASMODELS=ON`, so the gas-table and gas-model tests, most of the
physics test code by line count, have never run there (DR-157). There is still no sanitizer
configuration (`-fsanitize` appears in no CMake file), and Valgrind is run by hand. And a nightly
that passes says only that the suite passed: of 539 cataloged incidents (the locked catalog; addenda run to INC-564), 38 were closed by adding a test, so the
suite is not yet what finds our defects. `USE_TEST` defaults ON, so a fresh tree builds it.
Do not describe a change as validated on the strength of the build alone.

**Explicit Over Implicit**
- Manual memory management (`malloc`/`free`) on critical paths
- No hidden allocations, no reference counting, no automatic conversions
- Memory ownership is explicit via naming and documentation

**Mathematical Readability**
- Code should resemble equations from papers
- Relaxed naming rules in numerical kernels

### Naming Conventions (Strict in Framework Code)

| Prefix | Meaning | Example | Use For |
|--------|---------|---------|---------|
| `a` | Argument | `aTarget`, `aNumRows` | All function parameters |
| `t` | Temporary | `tCount`, `tBuffer` | Local variables |
| `m` | Member | `mMatrix`, `mRank` | Class data members |
| `g` | Global | `gLog`, `gComm`, `gTbulk` | Framework globals |

After the prefix, variables are camelCase (`aNumRows`, `tNodeList`); free functions and methods are
snake_case (`create_facet_map`, `comm_rank`). Capitalized mathematical symbols keep their case
(`aT`, `mB`).

**Mathematical Exceptions:** In numerical kernels implementing published algorithms, use standard notation:
- Matrices: `A`, `B`, `K` (stiffness), `M` (mass), `J` (Jacobian)
- Vectors: `x`, `y`, `f` (force), `u` (displacement)
- Indices: `i`, `j`, `k`
- This makes code match papers (e.g., Bathe, Hughes)
- **Math/physics module exemption (2026-07-23):** in `src/math` and
  `src/physics` kernels, the `a`/`t` prefixes may be dropped entirely —
  arguments and locals use plain mathematical names (`a, b, c, d`, `p, q, r`,
  `D`, `lambda`, …) so the code reads like the derivation. The `m` (member)
  and `g` (global) prefixes still apply. Example: `fn_cardano.hpp` /
  `fn_ferrari.hpp`.
- **Canonical thermophysical symbols (2026-08-06):** `src/physics/gastables`
  and `src/physics/gasmodels` apply the exemption above with a fixed symbol
  table, in both the API and the implementation:

  | quantity | symbol | | quantity | symbol |
  |---|---|---|---|---|
  | temperature | `T` | | mass specific enthalpy | `h` |
  | pressure | `p` | | mass specific entropy | `s` |
  | density | `rho` | | velocity | `u` |
  | mass specific volume | `v` | | cross section | `A` |
  | angles | `alpha`, `beta` | | | |

  Numbered states carry the index: a shock reads
  `shock( T1, p1, u1, alpha, T2, p2, u2, beta )`.

  Two rules make this safe, and both must be kept:
  1. **Always call member functions through `this->`.** The accessors are named
     `T`, `p`, `v`, `h`, `s`, `u`, `cp`, `mu`, … , so a parameter of the same
     name shadows them. Every self-call in these modules is already qualified;
     an unqualified one would silently become a compile error.
  2. **`u` is velocity in the flow routines and internal energy in the caloric
     accessors.** Both meanings are standard, so the symbol is context-bound —
     do not "unify" them.

  Not renamed, because the letters are taken by non-physical objects:
  `aA`/`aB` (comparator pointees, species labels in `RefGasFactory`) and `aC`
  (coefficient vector in `fn_GT_create_glue_poly.hpp`).

### Container Selection (Critical for Performance)

| Need | Use | Don't Use | Why |
|------|-----|-----------|-----|
| Dynamic arrays | `Cell<T>` | `std::vector<T>` | Debug bounds checking, FEM utilities, MPI-safe |
| Linear algebra | `Vector<T>`, `Matrix<T>` | `Cell<Cell<T>>` | BLAS/LAPACK, alignment, contiguous storage |
| Large boolean flags | `DynamicBitset` | `Cell<bool>` | 1 bit per flag in 64-bit words, bitwise ops |
| Time-stepping history | `ShiftRegister<T>` | `Cell<T>` | Fixed-depth, no reallocation |

**Critical:** `Vector`/`Matrix` are for linear algebra ONLY, not general arrays. Use `Cell` for node lists, IDs, etc.

**Size accessor split:** `size()` on `Cell` and `DynamicBitset`, `length()` on `Vector` and `Matrix`.
Historical, not yet unified — use whichever the container actually has.

### Matrix Memory Layout (Critical)

BELFEM matrices are **column-major** (Fortran convention) for BLAS/LAPACK interop. Column `j` is contiguous; row `i` is strided.

Inner loop must vary the **row** index:

```cpp
Matrix<real> A(nrows, ncols);

// FAST — inner over rows (contiguous):
for (uint j = 0; j < ncols; ++j) {
    for (uint i = 0; i < nrows; ++i) {
        A(i, j) = ...;
    }
}
```

Inverting the loops is correct but slow (cache thrash on every inner step).

**Backend caveat:** Under Armadillo, `A(i, j) == data[j*nrows + i]` exactly. Under Blaze, columns are padded for SIMD alignment, so the inter-column stride may exceed `nrows`. **Always use the `A(i, j)` accessor for element access** — never compute offsets from `data()`. `data()` is for bulk ops (BLAS calls, MPI transfers); MPI works because all ranks share the same backend. For MPI sends of a whole matrix, the payload length is `spacing() * n_cols`, never `capacity()`.

**Watch out:** many C++ libraries default to row-major (Eigen, NumPy, the default C-array mental model). BELFEM does not. See `doc/coding_philosophy.md` for the full rationale.

### MPI Distribution: share/receive vs broadcast (Critical)

| Payload | Use |
|---------|-----|
| Scalar | `broadcast(T&)` |
| Small fixed-size payload (~≤ 8 entries, e.g., a properties tuple or a 3×3 matrix) | `broadcast(...)` |
| Large or variable-size `Vector<T>` or `Cell<T>` (ephemeris columns, mesh quantities) | `share` / `receive` (chunked) |
| `Matrix<T>` | `broadcast(...)` — there is no `share( Matrix )`; for payloads too large for one unchunked broadcast use a manual `send`/`receive` loop |

`broadcast` is collective. It packs the payload into a single `MPI_Ibcast` without chunking — fine for small messages, but observed to fail on large datasets. For large variable-size data, use the **asymmetric** `share`/`receive` pair instead:

```cpp
if ( comm_rank() == 0 )
{
    belfem::share( vec );        // root only — chunked via comm_split
}
else
{
    belfem::receive( vec );      // non-root only — receives chunks
}
```

The `if/else` rank guard is required — `share` and `receive` are **NOT collective**. Calling either on the wrong rank is a deadlock. See `doc/coding_philosophy.md` for the full rationale.

### Memory Management Rules

**Use manual memory when:**
- Inside loops or on critical path (assembly, solve, communication)
- Interfacing with C libraries (MPI, BLAS, Fortran)

**Smart pointers OK when:**
- High-level utilities (I/O, post-processing)
- One-off allocations (setup, initialization)
- Outside hot paths

**Ownership patterns:**
- `Cell<T>` owns its data
- `Cell<T*>` non-owning (does not delete)
- Mesh owns Nodes/Elements (documented in header)
- Solver borrows Mesh (takes `Mesh*`, does not delete)
- Factories **produce and hand over**: a factory moves its data into the consumer and dies. No `mOwn` flags, no factory that outlives its product. (Holding a `shared_ptr` to the object under construction until handover is fine — several factories do.)

### Error Handling (Two-Tier System)

```cpp
// Debug-only checks (compiled out in release)
BELFEM_ASSERT(index < size, "Index %lu exceeds size %lu", index, size);

// Always-active runtime checks
BELFEM_ERROR(file.is_open(), "Failed to open: %s", filename.c_str());
```

**Use `BELFEM_ASSERT` for:** Logic bugs, bounds checks, preconditions, invariants in hot paths
**Use `BELFEM_ERROR` for:** File I/O, MPI failures, allocation failures — things that can fail in a correct program

Failed checks also land in the system log (identity `belfem`, with the MPI rank in the
payload) — `journalctl -t belfem` on Linux, `log show --predicate ... --last 1h` on macOS,
recovers the message of a crashed run even when the terminal output is gone. The hook is
POSIX `syslog(3)`, deliberately unguarded by platform. See
`src/core/doc/core_usage_guide.md`, "System Log Integration".

The tier follows the call rate: setup and initialization code, which runs once, may use `BELFEM_ERROR`
freely so the message survives a release build. Per-element and per-iteration code keeps `BELFEM_ASSERT`.

**Third category — expected algorithmic failure, never an abort.** Solver non-convergence, a failed
line search, a rejected timestep: these return a status and trigger the retry path above them
(timestep reduction, relaxation, controller fallback). `BELFEM_ERROR("did not converge")` inside an
algorithm that has a retry policy above it is a design error — it turns a recoverable state into a
run abort. See `doc/coding_philosophy.md`.

### Thread Safety

**BELFEM is deliberately NOT thread-safe internally:**
- MPI for parallelism (distributed memory), not threading
- No internal mutexes (adds latency to every operation)
- Users needing OpenMP: protect calls externally with `#pragma omp critical`

### Quick Reference Patterns

**Avoid in loops:**
```cpp
// BAD - hidden allocation
for (int i = 0; i < n; ++i) {
    std::vector<real> temp = compute();  // Allocation every iteration!
}

// GOOD - preallocate
Cell<real> temp;
for (int i = 0; i < n; ++i) {
    compute(temp);  // Reuse buffer
}
```

**No temporary `Vector`/`Matrix` in frequently-called member functions:**
The loop rule above extends to any class method expected to run repeatedly at
runtime (i.e., not tied to construction/initialization): the loop may live in
the *caller* — a driver or optimizer can invoke the method millions of times.
Allocate scratch objects as member variables, size them once in the
constructor, and mark them as scratch:

```cpp
// BAD - heap allocation on every call
void MyClass::update() {
    Vector<real> tAxis(3);          // Allocates each call!
    ...
}

// GOOD - member scratch, sized in the constructor
//! scratch for update() (header)
Vector<real> mAxis;                  // ctor: mAxis.set_size(3, BELFEM_QUIET_NAN);
void MyClass::update() {
    mAxis(0) = ...;                  // Reuse, no allocation
}
```

Expression-template assignments into an already-sized member
(e.g. `mInv = trans(mFwd)`) are fine — they evaluate in place. Temporaries
remain acceptable in constructors, initialization, and one-off setup paths.

**MPI-safe container access:**
```cpp
Cell<real> data;  // Even if empty
MPI_Send(data.data(), data.size(), ...);  // data() is well-defined when empty
```
`data()` on an empty container is well-defined but **possibly `nullptr`** (Armadillo returns null
deterministically, Blaze may not) — that composes with MPI, which accepts any pointer at count 0.
Never null-test `data()` as an emptiness proxy; test the size accessor instead.

**Alignment is backend-provided, not hand-rolled.** `src/` contains no `aligned_alloc` or
`posix_memalign` call, and nothing in BELFEM guarantees 64 bytes: Blaze aligns to the build's SIMD
width (32 under AVX2), Armadillo to 32 bytes for allocations ≥ 1 KiB and never 64. Do not write a
hand-aligned buffer — allocate through `Vector`/`Matrix` and let the backend align it.

**Marking and dedup:** entities carry an 8-slot flag bitset. Use it instead of building a
`Map`/`std::set` for "have I seen this node/edge/facet" work — and clear the flags both before
and after the pass.

**See `doc/coding_philosophy.md` for complete rationale and examples.**

## Development Notes

- C++17 standard with Fortran support for numerical libraries
- Extensive use of CMake for cross-platform compatibility
- Support for major linear algebra libraries (MKL, BLAS, etc.)
- Integration with scientific computing ecosystems (PETSc, MUMPS, etc.)
- The project builds as **one** library, `libbelfem`, assembled from per-module OBJECT libraries (`config/scripts/Add_Library.cmake`, `Add_BelfemLibrary.cmake`). It is STATIC by default and SHARED under `-DUSE_SHARED_LIBS=ON`; the per-module compile gate is `make belfem_<module>`. The plugin build templates (`src/physics/materials/User{Material,Library}Template.cmake`) and the example decks' plugins all produce MODULE libraries; they are dlopen()ed, never linked. `make install` deploys the library with the executables named in `BELFEM_INSTALL_EXECUTABLES` (`config/globals.cmake`), headers, `share/`, `examples/`, the plugin templates plus the annotated `example_user_*.cpp` under `share/belfem/templates/`, and `LICENSE` at the prefix root (see `todo/closed/shared_library_and_install_plan.md`)

## The Input Contract: Two Artifacts, One Rule

**CRITICAL — this rule is triggered by editing the C++, not by editing the docs.**

The `input.conf` contract is described in two places that must stay in lockstep:

| Artifact | Audience | Role |
|---|---|---|
| `doc/input_file_reference.md` | humans | prose reference: what a key means, why, and its pitfalls |
| `doc/input_schema.yaml` | programs | machine-readable contract: types, units, defaults, enums, aliases |

**If you add, remove, rename, or change the behavior of any `input.conf` key, you MUST update
BOTH files in the same session that lands the code change.** This applies to:

- a new key, or a key no longer read
- a changed default, unit dimension, or enum value set
- a new or removed alias, or a change in alias precedence
- a change in required-vs-optional, or in conditional requiredness
- a change in case-sensitivity of a value

Then run a Codex prose pass over the touched sections of the `.md` (read-only; apply what it returns).

### Why this rule keeps failing, and the one thing that fixes it

The rule previously lived only inside `doc/input_file_reference.md` itself — invisible to anyone
editing `cl_FEM_Controller.cpp`. Measured on 2026-08-10: of the 87 `file:line` citations in that
document, **26 were stale, 24 shifted by exactly +342** because the controller grew. Worse, the
document and the controller were edited **in the same commits** — the policy was being followed —
yet in commit `bc578b5e` the doc cites `tolerance` at `cl_FEM_Controller.cpp:2569` while that same
commit's controller has it at `:2911`. A line number read from a working copy is already wrong by
the time it is committed.

**Therefore: `doc/input_schema.yaml` anchors every key by a searchable token, never a line number.**

```yaml
anchor: '"tolerance"'                 # greppable, survives any edit
# NOT  site: cl_FEM_Controller.cpp:2569
```

The `.md` may keep `file:line` citations for human navigation, but they are advisory, not truth.

### The check that catches what discipline cannot

Human discipline reliably keeps the *prose* current — it has. What it cannot catch is a key added
to a factory and to no document, because whoever adds the key is precisely the person unaware the
document exists. That needs a mechanical pass in the **code → schema** direction: inventory every
`key_exists( "…" )` / `get_real( "…" )` / `get_value( "…" )` / `section_exists( "…" )` string
literal in the input consumers and diff it against `doc/input_schema.yaml`.

Note when writing that check: `Section::get_string` and `key_exists` both call `string_to_lower`,
so key lookup is case-insensitive from the caller side. Normalize case before comparing, or
`get_reals( "Direction" )` will look like an unknown key.

## Documentation Organization

**IMPORTANT:** BELFEM has a comprehensive three-tier documentation system. Follow `doc/documentation_guidelines.md` when creating or organizing documentation.

### Three-Tier Documentation Structure

```
belfem/
├── todo/                        # Task planning and implementation plans
│   ├── README.md               # Index of active/completed tasks
│   └── *.md                    # Future work, refactoring plans
│
├── devlog/                      # Session logs of changes made
│   ├── README.md               # Index of devlog entries
│   └── dlYYYYMMDD_topic.md     # Individual session logs
│
├── doc/                         # General project documentation
│   ├── README.md               # Master documentation index
│   ├── coding_philosophy.md    # CRITICAL: HPC design patterns
│   ├── documentation_guidelines.md  # How to organize docs
│   ├── literature_references.md     # Citations and DOIs for all referenced works
│   ├── input_file_reference.md / input_schema.yaml  # The input contract
│   └── *.md                    # Architecture, analysis, guides
│
└── src/[module]/doc/           # Module-specific technical documentation
    ├── README.md               # Module documentation index
    ├── [module]_usage_guide.md # Comprehensive user guide
    └── *.md                    # Algorithms, theory, implementation
```

**Documents state facts inline.** A `doc/` or `src/*/doc/` file must not cite a devlog, a todo file,
or a task ID as its source — devlogs may cite documentation, not the reverse.

### Navigation Strategy

**When you need to understand a module:**
1. Start with `src/[module]/doc/README.md` for overview and quick reference
2. Read the `*_usage_guide.md` it names — the stem follows the subject, not the directory
   (`container_usage_guide.md`, `dof_manager_usage_guide.md`, `lapack_usage_guide.md`)
3. Check specialized docs for algorithms, theory, implementation details
4. Cross-reference with `doc/README.md` for related modules

**Modules with their own `doc/` directory:**
`circuit`, `comm`, `containers`, `core`, `executables`, `fem`, `fem/interpolation`, `fem/iwg`,
`fem/kernel`, `fem/maxwell`, `fem/postproc`, `fem/thermal`, `homology`, `io`, `linalg`,
`math/graph`, `math/quaternion`, `math/tensor`, `mesh`, `numerics/opt`, `numerics/spline`,
`physics/database`, `physics/gasmodels`, `physics/gastables`, `physics/materials`, `sparse`,
`visualizer`.

**A module README typically provides:**
- Quick reference tables (key classes, files, operations)
- Common code patterns
- Development notes and common pitfalls
- Sometimes external references — but only one of the 27 currently carries a DOI, so do not expect
  the literature trail to live there

### Documentation Decision Tree

```
Need to create documentation?
│
├─► About future work? → ./todo/task_name.md
│
├─► Recording what changed? → ./devlog/dlYYYYMMDD_topic.md
│
├─► About existing code?
│   ├─► Module-specific? → ./src/[module]/doc/topic.md
│   └─► Cross-cutting? → ./doc/topic.md
│
└─► Root-level convention file? → Use CAPS (README.md, CLAUDE.md)
```

Work on the `nonfree/` tree logs into `./nonfree/devlog`, never into the open-source `./devlog`.

### Prose Gets a Language Sweep

**Whenever you add or revise prose written for a *user* to read — a module README, a `doc/*.md`
guide, a usage guide — hand it to Codex for a readability pass before treating it as finished.**
Codex is markedly better at plain-language editing than at the audit work it is otherwise given
here, and the drafts it improves are usually accurate but denser than they need to be.

```bash
CODEX_MODEL=gpt-5.6-luna CODEX_EFFORT=medium AI_EXCHANGE_SLUG=<slug> .claude/scripts/ask_codex.sh - < prompt.md
```

Pick the depth from the table in `doc/ai_collaboration_protocol.md` §9.1 — luna for an ordinary
guide, terra for a dense technical document. The wrappers stamp the chosen model and effort into
the exchange entry, so a finding is always attributable to how hard the auditor was asked to look.

Tell the sweep what it may **not** change, or it will smooth away meaning:

- technical claims — a better sentence that alters a fact is a regression, not an improvement; ask
  for suspected errors to be flagged separately rather than silently corrected
- file paths, code identifiers, class names, CMake option names, and `file:line` citations
- the document header block (`**Date:** / **Purpose:** / **Module:**`)
- any closing status paragraph that records what a document does *not* yet cover — that honesty is
  the point of it

Apply what comes back, but read each edit against the code before keeping it.

**The sweep is for user-facing documentation only.** These are never swept:

- `devlog/` entries and `tmp/ai_exchange/` files — records rather than documentation
- `METHODOLOGY.md` and `doc/lessons_learned.md` — descriptions of the method and its rules,
  written for us rather than for a user of the code
- `./todo/` files — working artifacts for the team. The earlier standing rule that swept every
  todo file was withdrawn on 2026-08-20 as not worth its cost
  (`doc/ai_collaboration_protocol.md` §6). Technical audits of a todo file's content are
  unaffected

The last two exclusions are there because **the sweep dilutes.** Smoothing costs information per
line. That is a good trade for a guide someone reads once to learn a module, and a bad trade for
a document whose whole value is how much it packs into a line that a session loads every time.
Where readability and density conflict in those files, keep the density.

This is the general form of the narrower rule in §"The Input Contract", which asks for the same
pass over `doc/input_file_reference.md`.

### Naming Convention (Strict)

**Standard format:** `lowercase_with_underscores.md`

✅ Good:
- `fem_dof_manager_architecture.md`
- `maxwell_thin_shell_formulation.md`
- `solver_performance_analysis.md`

❌ Bad:
- `FEM_DOF_Manager.md` (mixed case)
- `MaxwellThinShell.md` (camel case)
- `Build System Guide.md` (spaces)

**Always include in document header:**
```markdown
# Title

**Date:** YYYY-MM-DD
**Purpose:** Brief description
**Module:** [if module-specific]
```

### When Creating Documentation

1. **Choose directory** using decision tree above
2. **Name file** using lowercase_with_underscores
3. **Update README.md** in that directory
4. **Include code references** with file:line format: `cl_DofManager.cpp:123-145`
5. **Link to literature** using author-year citations: "Messe et al. 2023, §2.7" or "Bathe, §X.Y"
6. **Add quick reference tables** for classes, functions, common operations
7. **Document common pitfalls** and debugging tips

Task plans in `./todo/` follow `todo/plan_template.md` (Dn/Rn/On step IDs, gap table, living Status
line) and are registered in `todo/README.md`. Every step carries a live checkbox; obsolete steps are
struck through, not deleted, and boxes are ticked in the same session the step is finished.

**See `doc/documentation_guidelines.md` for complete workflow and examples.**

## Literature and References

BELFEM development is literature-first: formulations, algorithms, and validation choices are traced
to published work. Two artifacts support this:

| Artifact | What it is | Availability |
|---|---|---|
| `doc/literature_references.md` | Full citation list with DOIs — FEM papers, FVM papers, topology papers, textbooks and lecture notes — plus routing tables ("which reference answers which question") and the retired `paperN` alias decoder | Part of the open-source repository |
| `./literature/` | Curated local library of ASCII text extractions plus per-source navigation files and index/routing guides | **Proprietary, separate git repository**, `.gitignore`d here |

**Copyright:** the `literature/` directory is a completely separate repository, held under
institutional licenses for internal development use only. It **MUST NOT** be committed to the BELFEM
git, shared, or distributed. BELFEM itself is open source and fully usable without it. Its layout and
routing tables are documented in `literature/README.md` — start there when the directory is present,
and do not mirror its file inventory into this repository.

The `./nonfree/` repository, when present, carries an additional reference library; see
`./nonfree/doc/literature_routing.md`.

### CRITICAL: When to Consult Literature

**You MUST proactively consult the literature before making assumptions or relying solely on code analysis for:**

1. **Understanding algorithms** - How thin-shell, h-φ formulation, cohomology, MPFA work
2. **Resolving ambiguities** - Code behavior unclear or inconsistent
3. **Implementing new features** - Need theoretical foundation
4. **Debugging** - Results don't match expectations, need validation benchmarks
5. **Documenting code** - Need to explain *why* something is implemented a certain way
6. **Reviewing implementations** - Verify correctness against published methods

### Literature Routing Strategy

**Step 1: identify the question type**

| Question Type | Start With |
|---------------|------------|
| "How does BELFEM implement X?" | Messe et al. 2023 — BELFEM's own choices are documented there |
| "What's the theory behind h-φ?" | Arsenault et al. 2023 and 2021, **with** the 2026 erratum (it corrects the air-domain coupling) |
| "How do thin-shell models work?" | Messe et al. 2023, Alves et al. 2022a/2022b; Denis et al. 2026 for homogenized stacks |
| "What are cuts for, and how are they generated?" | Alves et al. 2022b, Schnaubelt et al. 2023; Pellikka et al. 2013 and Gross & Kotiuga for the cohomology |
| "How is a circuit coupled to the field problem?" | Dular et al. 1999 (natural H-formulation/circuit coupling); Dular et al. 1997 for source fields |
| "Why isn't this converging?" | Messe et al. 2023 §2.7, then Bathe 2016 Ch. 8 |
| "How do I control the timestep in a nonlinear transient?" | Badel et al. 2019 §4.3 (step controller, NR warm-start) |
| "What does a quench simulation look like end to end?" | Wozniak et al. 2025, then Badel et al. 2019 |
| "How does MPFA work?" | Aavatsmark 2002 (O-method fundamentals) |
| "Will my FVM mesh converge?" | Agélas et al. 2008 (coercivity), Klausen & Winther 2006 (rough grids) |
| "3D hexahedral FVM?" | Ingram et al. 2010 (enhanced BDDF₁), Wheeler et al. 2011 |
| "How do I implement element X?" | Hughes 2000 Ch. 5, Bathe 2016 Ch. 5 |
| "Why is my element locking?" | Hughes 2000 Ch. 4, Bathe 2016 §4.3–4.5 |
| "What's the formula for Z?" | Bronshtein |

**Step 2: use the index files** (when `./literature/` is present). `literature/README.md` holds the
master decision flowcharts; `papers/fem/index.md`, `papers/fvm/index.md`, `papers/topology/index.md`,
`books/index.md`, and `lectures/index.md` hold the detailed per-topic routing tables.

**Step 3: search the source text.** Each index file names the `.txt` extraction for the work it
routes you to. Search it with `rg -n` (or the equivalent search tool, keeping file+line output).

**Step 4: cross-reference.** Most topics benefit from two angles — theory plus implementation
(Brenner & Scott with Hughes), general plus specific (Zienkiewicz & Taylor with the HTS papers),
or robustness plus symptom (Bathe with the BELFEM papers).

### Key BELFEM Implementation Choices (from the literature)

| Decision | Choice | Source | Rationale |
|----------|--------|--------|-----------|
| Interface coupling | Static condensation (preferred) over Lagrange multipliers | Messe et al. 2023 | Avoids zeros on diagonal, better solver performance |
| h-φ variant | Magnetodynamic (H-φ/D) over quasi-static | Arsenault et al. 2023 (+ 2026 erratum) | Better convergence, natural coupling via Faraday |
| Thin-shell resolution | N > 1 for stacked tapes | Alves et al. 2022b | N=1 misses top/bottom losses (similar limitation as T-A) |
| Nonlinear tolerance | deck choice; the class default is `1e-6` (`cl_FEM_Controller.hpp`). Messe et al. 2023 recommend driving ε < 10⁻¹¹ after the Picard stage | Messe et al. 2023 | Prevent checkerboarding in HTS simulations |
| Solver choice | STRUMPACK first, MUMPS for robustness | Messe et al. 2023 | STRUMPACK ~2× faster but MUMPS more robust |
| Cut type | Thick cuts (cohomology) in BELFEM | Alves et al. 2022b, Schnaubelt et al. 2023 | Automatic generation vs manual thin cuts |

### Citation Format

Cite by **author and year**, with the exact section or equation:

```cpp
// Static condensation (Messe et al. 2023 — preferred method to avoid zeros on diagonal)
// Interface conditions from Alves et al. 2022b, Eq. 19
// Magnetodynamic coupling (Arsenault et al. 2023, Section II, Eqs. 5-7)
// Cohomology basis following Pellikka et al. 2013, Eq. 10

// General FEM: Mixed formulation stability (Bathe, §4.4.3)
// Newton-Raphson convergence strategy (Bathe, §8.6)
// Locking remedies (Hughes, Ch. 4)
```

- Papers: "Author et al. YYYY, Section X" or "Author et al. YYYY, Eq. Y"
- Books: "Author, §X.Y" or "Author, Ch. N"
- Always give the exact section/equation number for traceability

**Retired `paperN` aliases:** papers used to be cited as `paper1`, `paper6`, `F0`, … . Those aliases
are **retired** (2026-08-11) and were converted out of the source comments and the module
documentation; they survive only in dated devlog entries, which are kept as written. Do not write new
ones. `doc/literature_references.md` keeps the alias → work mapping so those records remain readable.

### Common Pitfalls (from the literature)

| Pitfall | Symptom | Solution | Source |
|---------|---------|----------|--------|
| Loose convergence tolerance | Checkerboarding | Drive ε < 10⁻¹¹ | Messe et al. 2023 |
| N=1 for stacked tapes | Missing top/bottom losses | Use N > 1 | Alves et al. 2022b, Table 1 |
| Identical FE orders at interface | Spurious oscillations | Hierarchical enrichment | Dular et al. 2021 |
| Missing cuts | Wrong current flow | Add cuts (thin or thick) | Alves et al. 2022b, Schnaubelt et al. 2023 |
| Lagrange multipliers | Zero diagonal blocks | Use static condensation | Messe et al. 2023 |

### If the Literature Directory Is Absent

1. Use `doc/literature_references.md` — full citation list, DOIs, and routing tables
2. State clearly when a claim rests on code analysis rather than literature
3. Flag the uncertainty and name the specific work that would settle it

---

## Quick Start Guide

**When starting a new session, remember:**

### 1. This is HPC Code — Performance First

- Manual memory management is intentional (not legacy code)
- Prefix-based naming (`a`, `t`, `m`, `g`) is strict in framework code
- `Cell<T>` for arrays, `Vector`/`Matrix` for linear algebra ONLY
- Matrices are column-major — the inner loop varies the row index
- `BELFEM_ASSERT` for debug checks, `BELFEM_ERROR` for runtime errors
- Not thread-safe by design (MPI-based parallelism)

### 2. Navigation Strategy

**Understanding a module:**
```
1. Read src/[module]/doc/README.md (quick reference)
2. Read [module]_usage_guide.md (comprehensive)
3. Check code with understanding from docs
```

**Answering a question:**
```
BELFEM-specific (HTS, h-φ, thin-shell, cuts)? → literature/papers/fem/index.md
FVM / MPFA?                                   → literature/papers/fvm/index.md
Cuts, cohomology, topology?                   → literature/papers/topology/index.md
General FEM theory?                           → literature/books/index.md
Functional-analysis foundations?              → literature/lectures/index.md
Formula lookup?                               → literature/books/index.md (Bronshtein)
No literature/ available?                     → doc/literature_references.md
```

**Creating documentation:**
```
Future work? → todo/task_name.md
Session log? → devlog/dlYYYYMMDD_topic.md
Module-specific? → src/module/doc/topic.md
Cross-cutting? → doc/topic.md
(Always lowercase_with_underscores.md)
```

### 3. Literature-First Approach

**Key FEM papers (HTS electromagnetics):**
- **Messe et al. 2023** — BELFEM core: architecture, h-φ, thin-shell, solver strategy — PRIMARY
- **Messe 2022** — BELFEM framework overview
- **Arsenault et al. 2023** — magnetodynamic h-φ coupling — RECOMMENDED for h-φ; read with the **Arsenault et al. 2026 erratum**, which corrects the air-domain coupling
- **Arsenault et al. 2021** — h-φ implementation in COMSOL
- **Alves et al. 2022a / 2022b** — 3D and h-φ thin-shell theory, cohomology cuts
- **Alves et al. 2024** — 2D thin-shell in COMSOL
- **Schnaubelt et al. 2023** — no-insulation coils, thin-shell approximation
- **Dular et al. 2021** — mixed formulation stability (inf-sup)
- **Dular et al. 1997 / 1999** — source-field computation and H-formulation/circuit coupling
- **Riva et al. 2023** — h-φ with domain decomposition in Sparselizard
- **Denis et al. 2026** — homogenized multi-scale h-φ thin shell for stacked coils
- **Lucchini 2025** — magnetization losses in 3D CORC tapes
- **Wozniak et al. 2025** — Ic defects, coupled magneto-thermal h-φ quench simulation, quench-detection margins
- **Badel et al. 2019** — dissipative zones in REBCO coils; adaptive time stepping (integral controller) and Newton-Raphson warm-starting (local file: `badel2021`)

**Key FVM papers (MPFA diffusion):**
- **Aavatsmark 2002** — MPFA O-method fundamentals, transmissibility — PRIMARY
- **Ingram et al. 2010** — enhanced BDDF₁ for 3D hexahedra — CRITICAL
- **Klausen & Winther 2006** — physical vs reference space, robust convergence
- **Agélas et al. 2008** — coercivity conditions for convergence
- **Wheeler et al. 2011** — MFMFE overview

**Key topology references (cuts and cohomology):**
- **Pellikka et al. 2013** — computation of cohomology bases for cuts
- **Gross & Kotiuga** — topological approach to electromagnetics
- **Dey et al. 2011, Chen & Freedman 2010, Dunfield & Hirani 2011** — optimal-cut theory and complexity
- **Cormen et al. (CLRS)** — the underlying graph algorithms

**Key books for general FEM:**
- **Bathe 2016** — robustness and reliability (Ch. 8 for solvers)
- **Hughes 2000** — linear FEM algorithms and implementation
- **Brenner & Scott 2008** — mathematical proofs and convergence theory
- **Zienkiewicz & Taylor 2013/2014** (Vol. 1/2) — comprehensive encyclopedia
- **Belytschko et al. 2014** — nonlinear FEM, XFEM
- **Arnold 2018** — topological foundations, why FEM spaces must exist (FEEC)
- **Boffi et al. 2013** — mixed FEM theory, saddle-point problems, inf-sup stability
- **Monk 2003** — Maxwell equations, edge elements (critical for electromagnetics!)
- **Evans 2017, Felippa** (lecture notes) — functional-analysis foundations and step-by-step FEM

Full citations and DOIs for all of the above: `doc/literature_references.md`.

---

**Ready to start? First action in any session:**
1. Identify if the question is BELFEM-specific or general FEM
2. Check the appropriate documentation (module READMEs, literature routing)
3. For algorithmic or validation-heavy tasks, understand the "why" from the literature before judging the code
4. Follow BELFEM coding patterns when writing or reviewing code
