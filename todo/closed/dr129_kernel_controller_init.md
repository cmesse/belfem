# DR-129: `Kernel::mController` is never initialized

**Date:** 2026-08-30
**Purpose:** `fem::Kernel` declares `Controller * mController` with no initializer and no
constructor-initializer-list entry, so every read of the member on a kernel that was never handed a
Controller — whether through `Kernel::controller()` or through the implicit copy constructor — is a
read of an indeterminate value. Give the member a `nullptr` default, add a `has_controller()`
predicate, guard the one accessor read that legitimately runs without a Controller
(`IWG_Maxwell::link_to_group`), and remove the production copy of an indeterminate `Kernel` in
`ThermalFactory`.
**Module:** `src/fem/kernel` (primary), `src/fem/maxwell` + `src/fem/thermal` (consumers),
`tests/fem` (gate)
**AIs involved:** Claude (exploration + plan + verification), Codex (audit), Grok (third voice)
**Status:** **CLOSED 2026-08-30 — all three run gates verified, register row STRUCK and archived
by Christian's ruling.** Debug `make check` 17/17 with `GhostElementContract [ OK ]`; release
(`USE_DEBUG=OFF`, `-O2 -DNDEBUG`) `make check` 17/17 with the same test, release-ness proven from
the linked binary rather than the cache; and a release coupled h-ɸ/T smoke (`tape_quench_usermat`,
8 ranks) that carried the R8 construction site at startup and then ran 483 converged coupled
timesteps with zero DR-129-class failures. Only R7 bookkeeping remains.
`make check` RAN GREEN 2026-08-30 in both modes with `GhostElementContract` passing in the rebuilt
binary — the re-enabled controller-less `link_to_group` executed and the `have_thermal()`
contract held. The release run (17/17 suites, `USE_DEBUG=OFF`, `-O2 -DNDEBUG`, ASSERT strings
confirmed absent from the linked binary) is the decisive one: it is the configuration in which the
pre-fix code dereferenced the indeterminate member with no assert to catch it. Previously: CODE-AUDITED — R1/R2/R3/R4/R5/R8 LANDED 2026-08-30 (Christian's
same-day rulings: O1 = (c), R1 = all three members, implementation approved). Code-audit round ran
same day at `terra`/`high` + `grok-4.6`/`high`: **unanimous APPROVE, zero logic defects**; both
of Grok's comment nits applied in the same session (the "both paths" comment overclaimed vs the
empty-group early return; the accessor comment said "garbage" post-R1). All four touched TUs pass
`g++ -fsyntax-only` with the build tree's defines — syntax gate only, **reviewed not verified.**
**Remaining: the R6 run gates (Christian's runs) and R7 closure.** Plan-audit round 1 also 2026-08-30,
at the §9.1 safety-boundary depth (`terra`/`xhigh` + `grok-4.6`/`xhigh`), reconciled below.

> **Scope guards:**
> - **OUT of scope:** any change to when a Controller is constructed, to `Controller`'s own
>   lifetime, or to the `Kernel`/`Controller` back-pointer topology. This row is about
>   initialization and the read contract, not about ownership.
> - **OUT of scope:** the other six `controller()` reads (§3) are audited here and left
>   *unguarded* by design unless O2 rules otherwise — they are ordered-safe on the production path
>   and guarding them would paper over a real ordering defect if one ever appeared. Both auditors
>   concur (round 1).
> - **Compatibility:** `has_controller()` is additive. No existing signature changes, no default
>   changes, no input-key changes — so this touches neither `doc/input_file_reference.md` nor
>   `doc/input_schema.yaml`.
> - The dangling-`mController`-after-`~Controller` teardown pattern (`belfem.cpp` declares
>   `tKernel` before `tControl`, so the Controller dies first) is **out of scope and pre-existing**:
>   `~Kernel` never calls `controller()`, `~Controller` does not clear the back-pointer, and none
>   of R1–R8 claims to police lifetime. `has_controller()` means "non-null association", not
>   "live Controller". (Codex, verified.)

---

## 1. Current Behaviour and How It Fails

`Kernel` declares the member with no default member initializer:

```
src/fem/kernel/cl_FEM_Kernel.hpp:70      Controller       * mController ;
```

Contrast the siblings four lines up and four lines down, which *do* carry one —
`Mesh * mSubMesh = nullptr ;` (`:67`) and `bool mOwnParameters = false ;` (`:73`). The omission is
an oversight in a class that otherwise defaults its optional members.

The constructor's initializer list (`cl_FEM_Kernel.cpp:46-52`) sets `mCommRank`, `mCommSize`,
`mParams`, `mMesh` and `mFieldOffset`. `mController` is not among them, and the constructor body
never assigns it. The only *named* writer in the tree is `Kernel::set_controller`
(`cl_FEM_Kernel.cpp:1127-1130`), called from exactly three places, all inside `Controller`
(`cl_FEM_Controller.cpp:87`, `:90`, `:4601`) — plus one *unnamed* writer the original filing
missed: the compiler-generated copy constructor (D3 below).

The accessor guards with the wrong tier:

```
src/fem/kernel/cl_FEM_Kernel.cpp:1135-1139
        Controller *
        Kernel::controller()
        {
            BELFEM_ASSERT( mController != nullptr, "Controller is not set for this kernel." );
            return mController ;
        }
```

`BELFEM_ASSERT` compiles out under `USE_DEBUG=OFF` entirely, and even in a debug build it cannot
see the failure mode that matters: an uninitialized pointer is overwhelmingly likely to be
**garbage-non-null**, which passes the assert and is then dereferenced.

| Failure | Mechanism | Evidence |
|---|---|---|
| SEGV, battery-order dependent | `link_to_group` dereferences the uninitialized member; whether it faults depends on what the allocator left in those 8 bytes | proven 2026-08-28: the DR-46 ghost fixture's `link_to_group` call passed solo and SEGV'd deterministically after 21 sibling tests |
| Guard cannot fire | assert is debug-only and null-only, the defect is release-live and garbage-valued | `cl_FEM_Kernel.cpp:1137` |
| Indeterminate-value copy on **every coupled production run** (D3) | `cl_ThermalFactory.cpp:138` copy-constructs the thermal kernel from a temporary whose `mController` and `mMyCommIndex` are indeterminate; formal UB, no dereference | `belfem.cpp:157`, `hphiTrun.cpp:80` both reach it; Codex round 1, re-verified |

> ~~Original filing also claimed a "silent wrong answer" mode: garbage `->thermal_kernel()`
> setting `mHaveThermal` wrong and mis-dispatching physics.~~ **REFUTED in round 1 by both
> auditors independently, re-verified:** `have_thermal()` has **zero consumers** in the tree —
> only its declaration and inline getter exist (`cl_IWG_Maxwell.hpp:90-115`) plus the write at
> `:266`. The write is a dead store today; the live failure is the dereference itself. R3
> survives (the read is still UB without a Controller), but "wrong physics" is not this row's
> severity story. (D5a.)

The fixture already documents the trap in a comment rather than fixing it, and skips the call
(`tests/fem/support/cl_TS_TestStack.hpp:892-899`); the consuming test is
`tests/fem/test_InterfaceOrientation.cpp:673` in the gtest binary `test_fem`.

**Reachability, re-verified 2026-08-30 and independently re-derived by both auditors.**
Every `IWG_Maxwell::link_to_group` call comes from a `DofManager` assembly routine
(`cl_FEM_DofManager.cpp:471`, `:518`, `:669`, `:709`, `:772`, `:865`, `:944`, `:991`, `:1224`),
and on the production path every executable builds the kernel before the controller whose
constructor writes the back-pointer (`belfem.cpp:132`→`:144`, `hphirun.cpp:69`→`:82`,
`hphiTrun.cpp:65`→`:77`); `load_memdump` comes after in all three, controller creation is not
rank-guarded, and `MaxwellFactory::create_controller` calls `link_maxwell` only *after* the
Controller ctor has run (`cl_MaxwellFactory.cpp:979-997`), which also covers the `MaxwellData`
member-init read (row 2). So: **no shipped path reaches a controller-less `Kernel::controller()`
read.** Confidence: high — three independent enumerations agree.

The claim must be worded that narrowly (D5b): controller-less kernels that *assemble* are common —
`poisson.cpp:40-65`, `staticheat.cpp`, `hangingnodes.cpp`,
`src/homology/cl_CutProcessorManual.cpp:83-114` all construct a `Kernel`, never a Controller, and
call `compute_jacobian` — they are safe only because their IWGs' `link_to_group` never reads
`controller()`. (Grok C5, spot-verified on `poisson.cpp`.)

**Bottom line:** two distinct reads of indeterminate state are live — a heap-lottery dereference
in every controller-less unit context, and a formal-UB member copy on every coupled production
run. Production ordering hides the first; the second executes on every thermal deck. Both are
cured by initialization plus one construction-site fix, at a cost of a handful of lines.

## 2. Architecture: Why a `nullptr` Default Plus One Guarded Read

Three things are being fixed and they should not be confused:

1. **The member must be `nullptr` by default** (with `mMyCommIndex` handled in the same stroke —
   D2 revised). Unconditional, not a design question. It makes the existing `BELFEM_ASSERT`
   able to fire, converts a garbage-valued dereference into a deterministic null dereference,
   and makes the D3 copy well-defined. Every other optional member of this class already does it.

2. **The `link_to_group` read needs a defined answer when there is no Controller.** Genuine
   choice; the recommended shape after round 1 is a both-path assignment (R3), not a skip.

3. **The production copy of an indeterminate `Kernel` should not exist at all** (R8). The
   temporary-plus-copy at `cl_ThermalFactory.cpp:138` buys nothing; `make_shared` can construct
   in place.

The reason the guarded read is cheap is that the target member already has the right default:
`bool mHaveThermal = false ;` (`cl_IWG_Maxwell.hpp:53`). A kernel with no Controller has no
thermal kernel by construction — `thermal_kernel()` is `Controller::mKernel2`, written only by the
Controller ctor and `set_thermal_kernel` (verified, Grok Q3). Rejected alternatives (unchanged
from the draft): a `BELFEM_ERROR`-only fix (kills the fixture), a dummy Controller for tests
(disproportionate), reading the thermal kernel from elsewhere (topology surgery). Precedent:
`cl_IWG_MaxwellThermal.hpp:23-24` documents the same guarded-at-the-consumer pattern.

## 3. Gap Table — every read of the indeterminate members in the tree

Rows 1–7: the complete set of `Kernel::controller()` call expressions, confirmed complete by
three independent enumerations (Claude grep; Grok including `parent()->parent()` chains and
aliased locals; Codex including headers, tests and `nonfree/`). Row 8 is the non-accessor read
Codex found.

| # | Site | Reads | Ordered-safe in production? | Class | Rationale |
|---|---|---|---|---|---|
| 1 | `cl_IWG_Maxwell.cpp:266` | `controller()->thermal_kernel()` | yes — assembly-time only (§1) | **(c) handle explicitly** | the one site with a legitimate controller-less caller (unit fixtures, future library consumers). Guard it → R3 |
| 2 | `cl_FEM_Calculator.cpp:78` | `controller()->time()` in a **member-init list**, bound to `const real & mTime` (`cl_FEM_Calculator.hpp:143`) | yes — `link_maxwell` runs from `create_controller()` after the Controller ctor (`cl_MaxwellFactory.cpp:986-994`) | (a) no change | **cannot** grow a local guard without turning the reference into a copy/pointer — this is the deciding argument for O1(c), both auditors |
| 3 | `cl_FEM_Calculator.cpp:1985` | `controller()->kernel()->dofmgr()->sideset(...)` | yes — field evaluation after the solve (used from `cl_FEM_Calculator.hpp:2564`, `cl_MaxwellPostprocessor.cpp:541`, `:852`) | (a) no change | |
| 4–7 | `cl_MaxwellPostprocessor.cpp:723`, `:784`, `:813`, `:903` | `thermal_kernel()` / `time()` | yes — postproc after the solve | (a) no change | |
| 8 | `cl_ThermalFactory.cpp:138` | implicit copy ctor reads indeterminate `mController` **and** `mMyCommIndex` | **no — executes on every coupled run** | **(c) handle explicitly** | `Kernel` declares a destructor and no copy/move members, so C++17 suppresses the implicit move ctor and `make_shared<Kernel>(Kernel(mKernelParameters))` copy-constructs from the temporary → R8 |

### 3.1 Cross-cutting findings

- **The register row understated the blast radius twice.** One consumer → seven accessor reads
  (D1), and an eighth read that is not an accessor call at all (D3). Severity direction differs:
  rows 2–7 are ordered-safe, row 8 executes in production.
- **Rows 2–7 are safe by ordering, not by construction.** After R1, a future ordering mistake at
  any of them produces a deterministic null dereference — and, if O1(c) lands, a named
  `BELFEM_ERROR` in release. That is the actual safety win.
- **The D3 copy is benign on the live path beyond the two indeterminate scalars** (Claude trace,
  round-1 calibration of Codex's finding): both production callers use the
  `ThermalFactory(file, Kernel*)` ctor, so `KernelParameters::mKernel != nullptr`
  (`cl_FEM_KernelParameters.cpp:49`) and both the Kernel ctor and `distribute_mesh` take their
  non-owning else branches — at the moment of copy the temporary has `mOwnCommTables`,
  `mOwnMesh`, `mOwnSubmesh`, `mOwnParameters` all false and empty
  `mDofManagers`/`mIWGs`/`mMaterials`/`mBoundaryConditions`, so its destructor deletes nothing:
  **no double-free, no use-after-free.** Confidence: high.
- **Latent hazard beside row 8:** the mesh-based `ThermalFactory` ctor
  (`cl_ThermalFactory.cpp:17-25`, `mKernel = nullptr`) has **zero callers** in the tree — but if
  it were ever used in parallel, the temporary at `:138` would *own* its comm tables
  (`cl_FEM_Kernel.cpp:696`), the copy would inherit both pointer set and ownership flag, and the
  temporary's destructor would free the tables under the live kernel: use-after-free plus double
  delete. R8 removes the whole class. Same implicit-copy family as DR-141
  (`ElectricalCircuit`/`Solver`, closed 2026-08-29).
- **Recorded future hazard (Grok, verified):** `create_thermal_kernel` runs
  `mThermalField->initialize()` (`cl_ThermalFactory.cpp:242`) *before* `set_thermal_kernel` is
  called by the executables — a real wall-clock window in which the thermal kernel has no
  Controller. Harmless today: thermal calculators are never `link_maxwell`'d in that window, so
  `allocate()` skips `MaxwellData`, and `IWG_MaxwellThermal::link_to_group` does not read
  `controller()`. **Do not "fix" this by constructing `MaxwellData` earlier** — that would create
  exactly this row's defect on the thermal kernel in every coupled run. Listed under O1's legal
  states.

## 4. Ordered Steps

- [x] **R1 — default the members.** ( LANDED 2026-08-30: `mController = nullptr` + `mMyCommIndex = 0` at `cl_FEM_Kernel.hpp:70`/`:76`, and `mFunMKF = nullptr` per Christian's O4 ruling — all three members ride the diff. ) `cl_FEM_Kernel.hpp:70`: `Controller * mController = nullptr ;`
      and (D2 revised — the member is *read* by the D3 copy on every coupled run, so this is no
      longer optional) `cl_FEM_Kernel.hpp:76`: `proc_t mMyCommIndex = 0 ;`. Removing
      `mMyCommIndex` outright stays a separate ruling (O3).
- [x] **R2 — add the predicate**  ( LANDED 2026-08-30: declared beside `controller()`, defined inline beside `is_master()` per its precedent. ) *(after: R1)*. `bool has_controller() const ;` on `Kernel`,
      returning `mController != nullptr`. `has_controller` does not exist anywhere in `src/`
      today. Naming precedent: `has_edge_dofs` (`cl_IWG_Maxwell.hpp:79`). Inline in the header
      (like `is_master`, `cl_FEM_Kernel.hpp:384-388`) or out-of-line beside `controller()` —
      implementer's choice, both conform.
- [x] **R3 — guard the one legitimate read**  ( LANDED 2026-08-30: both-path assignment through a `tKernel` local; `Group::parent()` → `DofManagerBase*` → `parent()` → `Kernel*` chain type-verified. ) *(after: R2 — does **not** wait on O1; D5c)*.
      `cl_IWG_Maxwell.cpp:266` becomes a **both-path assignment** (Codex round 1 — assign, don't
      skip, so a relink can never leave stale state):
      `mHaveThermal = aGroup->parent()->parent()->has_controller() && aGroup->parent()->parent()->controller()->thermal_kernel() != nullptr ;`
      (short-circuit safe; hoist the parent chain into a local per taste). The existing
      empty-group early return at `:260-264` stays the first precondition.
- [x] **R4 — apply the O1 ruling**  ( LANDED 2026-08-30: O1 ruled (c) by Christian; `BELFEM_ASSERT` → `BELFEM_ERROR` at the accessor. ) *(after: O1)*. If O1 lands as (b)/(c), promote the
      `BELFEM_ASSERT` at `cl_FEM_Kernel.cpp:1137` to `BELFEM_ERROR`. Per-group call rate, so the
      tier rule allows it; both auditors recommend it. **Note (Codex): O2's "fail loudly in
      release" property only exists if R4 lands — with the assert retained, release builds return
      null and dereference later. R4 is therefore a prerequisite of the O2 no-change ruling, not
      an optional extra.**
- [x] **R5 — re-enable the fixture call**  ( LANDED 2026-08-30: skip comment replaced by the call in `cl_TS_TestStack.hpp`; `ASSERT_FALSE( tStack.iwg()->have_thermal() )` added in `test_InterfaceOrientation.cpp` — `have_thermal()`'s first consumer. ) *(after: R3)*. Replace the skip comment at
      `tests/fem/support/cl_TS_TestStack.hpp:892-899` with the actual
      `mIWG->link_to_group( mGhostGroup )` call; rewrite the comment to record what the call now
      proves. Domain-safe: `DomainType::Ghost` is a handled dispatch case
      (`cl_IWG_Maxwell.cpp:417-420`, sets `mFunMKF = &maxwell::h_ghost`), and the call also sizes
      `mTimeStepMatrices` (`:271-273`). Add an assertion that `have_thermal()` is false after the
      link (Codex — this is what pins R3's semantics, and it gives `have_thermal()` its first
      consumer).
- [x] **R6 — GATE** ( 3 of 3 RUN: all sub-gates VERIFIED 2026-08-30 — see below ) *(rewritten after round 1; the draft's battery order-shuffle gate was
      refuted, D5d)*. After R1 the pointer is deterministically `nullptr` — shuffle cannot make a
      defaulted member garbage, and no in-tree shuffle target exists anyway (DR-46's shuffled
      runs were ad-hoc private-binary runs). The gate that can observe this row:
      1. **Debug:** `make check` green, with the re-enabled `link_to_group` in `test_fem`
         (`test_InterfaceOrientation.cpp:673` consumes the fixture) **passing** — the acceptance
         criterion is *pass*, not "named error or pass": R3 prevents the accessor call entirely
         (Codex). **RAN AND PASSED 2026-08-30 (Christian's build): `test_fem` rebuilt 15:55,
         `InterfaceOrientation.GhostElementContract [ OK ] (3 ms)` in
         `cmake-build-debug/Testing/Temporary/LastTest.log:17041-17092`, suite green; the new
         assertion string confirmed present in the binary, so the pass is the new test, not a
         stale link. VERIFIED. Confirmed a second time 2026-08-30 in the independent `build/`
         tree: `test_fem` relinked 16:38, full `make check` 17/17 suites green
         (`build/Testing/Temporary/LastTest.log:17041-17092`), `have_thermal` present in the
         binary. Both trees are `USE_DEBUG=ON` — neither run touches sub-gate 2.**
      2. **Release (`USE_DEBUG=OFF`):** the same test passes. This is the run that would still
         SEGV if R3 were omitted and O1 stayed at (a). **RAN AND PASSED 2026-08-30 (Christian's
         build): `build/` reconfigured to `USE_DEBUG=OFF` at 17:06 (`CMAKE_BUILD_TYPE=Release`,
         `-O2 -DNDEBUG`), `test_fem` relinked 17:13, full `make check` 17/17 suites green,
         `InterfaceOrientation.GhostElementContract [ OK ] (3 ms)` in
         `build/Testing/Temporary/LastTest.log:13062-13113`. Release-ness confirmed against the
         linked binary, not the cache: an ASSERT-only string is ABSENT (0 hits), so
         `BELFEM_ASSERT` really is compiled out and the pass cannot be coming from a debug-only
         guard — this is precisely the configuration in which the pre-fix code dereferenced the
         indeterminate member. The R4 `BELFEM_ERROR` string ("Controller is not set for this
         kernel") and `have_thermal` are both PRESENT in the same binary, so O1(c)'s always-active
         guard and R5's contract assertion survive into release. VERIFIED.**
         *(First attempt at this sub-gate, the 16:38 `build/` run, was a debug binary — the
         `-DUSE_DEBUG=OFF` had not reached the configure. `option( USE_DEBUG ... ON )`
         (`CMakeLists.txt:90`) takes the cached value unless `-D` is on the `cmake` line itself,
         so a bare `make` reconfigure silently keeps the old mode. Check
         `CXX_FLAGS`/`CXX_DEFINES` in `flags.make`, or an ASSERT-only string in the binary —
         not the label on the tree.)*
      3. **Coupled thermal smoke** (any small coupled deck) to exercise the R8 construction site.
         **RUNNING 2026-08-30 17:19:57 (Christian): `tape_quench_usermat`, `prterun -np 8`, from a
         RELEASE build (`cmake-build-debug/` reconfigured `USE_DEBUG=OFF` at 17:17, `-O2 -DNDEBUG`
         in `flags.make`, `bin/belfem` relinked 17:18 — mode read from the compiler artifacts, not
         the tree's misleading name). Deck is genuinely coupled: "Thermal solver section found in
         the input file: solving the coupled h-ɸ/T problem", and every step prints both Magnetic
         and Thermal Picard blocks. Warm restart from `memdump.hdf5`, resuming t = 11.0000 ms at
         step 203. The R8 site is traversed regardless of restart, because `ThermalFactory`
         constructs the thermal kernel BEFORE `load_memdump` in all three executables (the same
         ordering established in §1's reachability analysis). At 17:21 the run had reached step 208
         with magnetic 8.0e-11 / thermal 2.0e-13 and no errors — so the construction site is past
         and the constructed kernel is live. **PASSED — closed 2026-08-30 21:35 on 4 hours of
         evidence, without waiting for the deck to finish.** That run completed **483
         coupled timesteps** (resumed at step 203, reached step 686; 17:20:32-21:07:41), converging steadily at ~8e-11 magnetic
         / ~3e-13 thermal, with **zero DR-129-class failures**: no SEGV, no
         `"Controller is not set for this kernel"`, no symptom of an indeterminate read. The
         remaining hours cannot add evidence, because the R8 construction site executes **once, at
         startup** — it was behind the run by 17:20:32, and everything after is the deck's physics.
         *By-catch, NOT this row's business and left unfiled:* 8 `MUMPS soft fail: error code -20`
         and 2 `PARPACK failed to compute the eigenvalue conditioning estimate` in the
         eigenvalue/conditioning instrumentation path; every one was recovered (26 further
         timesteps succeeded after the last), so they are handled soft-fail paths, not aborts.**
      4. Optional belt-and-suspenders: `./test_fem --gtest_shuffle` a few recorded seeds — the
         binary is gtest (`tests/fem/test_Integration.cpp:12`), so the flag exists even though no
         make target wraps it.
- [x] **R7 — register + devlog** ( DONE 2026-08-30 ). Boxes ticked; the DR-129 row's status cell
      now names all three gates and what each actually ran; `todo/README.md` and
      `devlog/README.md` entries updated; the session devlog
      (`devlog/dl20260830_dr129_kernel_controller_init.md`) records both audit rounds, the release
      gate, the debug-binary false start, and the early-close reasoning for the smoke.
      `scripts/check_doc_claims.py` reports 37/37. **The register strike itself is left for
      Christian's ruling** — everything else in this step is complete.
- [x] **R8 — remove the production copy**  ( LANDED 2026-08-30: `std::make_shared<Kernel>( mKernelParameters )`, constructs in place. ) *(independent of R1–R5; before R6.3)*.
      `cl_ThermalFactory.cpp:138`: `std::make_shared<Kernel>( mKernelParameters )` — construct in
      place, no temporary, no implicit copy. Kills the live indeterminate-value read *and* the
      latent owning-copy double-free (§3.1). Numbered R8 to keep R1–R7 stable across the audit
      thread; lands logically with R1.

## 5. Open Design Questions (not silently decided)

- **O1 — RESOLVED 2026-08-30 → (c), decided by Christian** ( `BELFEM_ERROR` in `controller()`
  plus `has_controller()` ), landed as R4. Original question and options kept below for the record.
  After R1/R3, is a controller-less `Kernel` legal, or is it a contract violation?
  - *(a) Legal:* keep `BELFEM_ASSERT`. Cost (both auditors): rows 2–7 stay debug-only — release
    null-derefs without a message.
  - *(b) Contract violation:* promote to `BELFEM_ERROR`.
  - *(c) Both — `BELFEM_ERROR` in `controller()`, `has_controller()` for callers that ask first.*
    **Recommended by Claude, Codex and Grok independently.** The deciding reason (both auditors,
    verified): row 2 binds `const real & mTime` in a member-initializer list and *cannot* grow a
    local guard, so the accessor is the only place a release-loud check can live. The semantic
    rule ("assert for internal bugs") leans (a); overridden because after R1 the release
    alternative is a messageless null deref, and `IWG_MaxwellThermal::link_to_group` already uses
    `BELFEM_ERROR` at the same tier and call rate (`cl_IWG_MaxwellThermal.cpp:105-108`).
    **Christian's call; not decided here.**

  Legal states the `BELFEM_ERROR` must *accept* (must not fire) — enumerated per the
  new-guard-new-failure-mode rule: (1) magnetic kernel after the Controller ctor returns, all
  ranks, all three Maxwell executables; (2) thermal kernel after `set_thermal_kernel` (`:4601`,
  before the `link_maxwell` relink at `:4640-4653`); (3) the two-kernel ctor with `aKernel2`
  passed (both back-pointers written at `:87-90`). States that stay legal *as objects* but on
  which `controller()` must die: controller-less kernels in tests (`cl_TS_TestStack.hpp`,
  `cl_EF_TestVolume.hpp`, `test_DofSeeding.cpp`), `poisson.cpp`, `staticheat.cpp`,
  `hangingnodes.cpp`, `CutProcessorManual.cpp`; and the thermal kernel inside the
  `create_thermal_kernel` → `set_thermal_kernel` window (§3.1) — no current caller reads there, a
  future read *should* abort. What the error will never catch: the dangling back-pointer after
  `~Controller` (scope guard, above).

- **O2 — Should rows 2–7 be guarded too?** No — both auditors agree, argument unchanged: a
  fallback would bind `MaxwellData::mTime` to a dummy clock or make the postprocessor silently
  report `gTbulk`/t=0 with no thermal field. Loud failure is the better mode. **Conditioned on
  R4** (see R4 note). Also recorded (Codex): R3's predicate must not be read as implying MPI
  synchronization — in a future partial-rank library setup `mHaveThermal` could diverge across
  ranks; shipped drivers construct the Controller on every rank, so not a current concern.
- **O3 — remove `mMyCommIndex`, or keep it initialized?** R1 initializes it (that part is no
  longer optional — D2 revised). Whether to *delete* the member is a separate L-17-checked
  ruling: no named use, no accessor, no `friend` — but deletion is a distinct review. Christian's
  call, any session. ( Init-only landed 2026-08-30; the deletion question stays open. )
- **O4 — RESOLVED 2026-08-30 → ride R1, decided by Christian.** `IWG_Maxwell::mFunMKF` lacked `= nullptr` too (`cl_IWG_Maxwell.hpp:46-51`; Grok round 1,
  verified). Same defect class, same function; the thermal sibling has the default
  (`cl_IWG_MaxwellThermal.hpp:26-28`). Production always assigns it in `link_to_group` before
  use; the fixture currently leaves it garbage and calls `h_ghost` directly, and R5 will start
  assigning it. A one-line `= nullptr` would ride R1 naturally — but it is a *third* member and
  scope is Christian's call. Not this row's 2026-08-28 crash (verified: the crash was the
  `controller()` read).

## 6. Defect Tracker

**Round 1 — 2026-08-30, plan audit, Codex (`gpt-5.6-terra`/`xhigh`) + Grok (`grok-4.6`/`xhigh`),
reconciled by Claude with tree re-verification of every kept finding.**

- **D1 — the register row names one consumer; there are seven accessor reads.** LOW
  (documentation accuracy). Found by Claude 2026-08-30, confirmed complete by both auditors
  independently. Resolution: register row widened 2026-08-30. Reviewed.
- **D2 — `Kernel::mMyCommIndex` is uninitialized — REVISED in round 1: it is NOT dead.** Was
  filed as "dead member, optional cleanup". Codex: the D3 implicit copy *reads* it (indeterminate
  `int` copy) on every coupled thermal run, so it is live UB. The "no named use, no accessor, no
  friend" half of the original claim stands (all three enumerations agree). Consequence: R1 must
  initialize it; only the *deletion* question remains optional (O3). Found Claude / escalated
  Codex, verified Claude 2026-08-30. Reviewed.
- **D3 — production copy of an indeterminate `Kernel`.** HIGH (formal UB executing on every
  coupled run; benign in practice on the live path — §3.1 calibration). `cl_ThermalFactory.cpp:138`
  `std::make_shared<Kernel>(Kernel( mKernelParameters ))`: user-declared destructor
  (`cl_FEM_Kernel.hpp:102`) + no copy/move members → C++17 selects the implicit copy ctor for the
  rvalue; memberwise copy of indeterminate `mController`/`mMyCommIndex`; reached from
  `belfem.cpp:157` and `hphiTrun.cpp:80`. Found by Codex round 1; reachability, ownership-flag
  state and the no-double-free calibration verified by Claude against
  `cl_FEM_KernelParameters.cpp:46-49`, `cl_FEM_Kernel.cpp:649-733`, and both executables.
  Fix → R8 (+ R1 for the members). Reviewed.
- **D4 — `IWG_Maxwell::mFunMKF` uninitialized function pointer.** LOW (production always assigns
  before use). Found by Grok round 1, verified by Claude. Disposition → O4, Christian's call.
  Reviewed.
- **D5 — defects in this plan's own draft, found by the round and fixed in this revision:**
  - **(a)** the "silent wrong physics dispatch" severity claim — refuted, `have_thermal()` has
    zero consumers (both auditors independently; verified). §1 rewritten.
  - **(b)** the reachability sentence overclaimed — controller-less kernels *assemble* routinely
    (`poisson.cpp` et al.); only the `controller()` read is unreached (Grok C5; spot-verified).
    Narrowed in §1.
  - **(c)** circular dependency: the draft had R3 "(after O1)" while O1 said "after R1/R3"
    (Codex; Grok C3 concurring). R3 now depends only on R2.
  - **(d)** the battery order-shuffle gate — refuted: post-R1 the member is deterministically
    null, and no in-tree shuffle target exists (both auditors). R6 rewritten; gtest-level
    `--gtest_shuffle` survives as optional only.
  - **(e)** (Claude, self-found during reconciliation) the plan-drafting caller grep used a
    pattern that missed stack declarations (`ThermalFactory tTFactory(...)`) and wrongly
    concluded `belfem.cpp` was the only `ThermalFactory` caller; both vendors had it right.
    Corrected in §3.1/D3. The enumeration lesson: match the *type name*, not `TypeName(`.

**Round 2 — 2026-08-30, code audit on the landed diff, Codex (`gpt-5.6-terra`/`high`) + Grok
(`grok-4.6`/`high`).** Verdicts: **APPROVE / APPROVE**, no blocking findings, no new defects
introduced. Both independently confirmed the two highest-risk items: the re-enabled fixture call
survives a static line-by-line walk (`mTimeStepMatrices` is allocated in the `IWG_Timestep` base
ctor before the `IWG_Maxwell` body; the Ghost `FieldList`/dof-table arms are wired; nothing after
the switch dereferences the controller — Grok at ~90 % with the residue explicitly named as the
unrun gate), and R3 is bit-identical on every attached-controller path (`parent()` chain
non-virtual and unshadowed). `BELFEM_ERROR` confirmed unconditional in release
(`src/core/assert.hpp:263-276` vs the `BELFEM_ASSERTIONS_ACTIVE`-gated assert at `:242-259`).
`mFunMKF` has no null-test consumer in `IWG_Maxwell` (the null tests are on
`IWG_MaxwellThermal`'s *own* member), so the default is behavior-neutral. Nits applied same
session: two comment corrections (Grok Q7). Recorded residuals, all deliberate/pre-existing: the
implicit `Kernel` copy ctor still exists (Rule-of-5 hygiene, not this row — after R1 a copy is
defined but shallow on owning Cells; the one production copy site is gone), no `EXPECT_THROW`
coverage of the bare `controller()` abort, empty-group return still skips the `mHaveThermal`
write (pre-existing skip, unobservable — no production consumer).

## 7. Definition-of-Done Checklist

- [x] Every gap-table row mapped to a step (rows 1 → R3, 8 → R8) or an explicit no-change ruling
      (rows 2–7 → O2).
- [x] Each claimed gap backed by a citation re-opened at its named file, not by inference.
- [x] Plan audited by Codex **and** Grok (round 1, 2026-08-30; `git status` after the Grok round
      ran clean — modified-file set matched the session baseline).
- [x] O1 ruled by Christian before R4 lands (R3 does not wait) — ruled **(c)** 2026-08-30, applied
      in R4 (`BELFEM_ASSERT` → `BELFEM_ERROR` at the accessor).
- [x] O3/O4 scope ruled by Christian (ride R1, or stay out) — 2026-08-30: **all three members ride
      R1** (`mController`, `mMyCommIndex`, `mFunMKF`).
- [x] Code audited by Codex **and** Grok after R1–R5/R8 land (2026-08-30, APPROVE/APPROVE;
      `git status` after the Grok round matched the expected set exactly).
- [x] R6 gate RUN, not assumed: name the build(s), the suite result, the release-tree result and
      the thermal smoke in the DR-129 status cell. **Struck is not verified**; a clean compile is
      not this row's gate. ( All three RUN 2026-08-30: debug `make check` 17/17 with
      `GhostElementContract [ OK ]`; release `make check` 17/17 with the same test, release-ness
      proven by an ASSERT-only string being ABSENT from the linked binary; coupled thermal smoke
      483 steps clean in a release build. )

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/dr129_kernel_controller_init.md` (both round-1 entries
  appended by the wrappers with model/effort stamps; ephemeral — distilled here and into the
  devlog before sweep).
- **Codex round 2, code** (`gpt-5.6-terra`, `high`): APPROVE, all eight brief questions clean,
  `git diff --check` clean. **Grok round 2, code** (`grok-4.6`, `high`): APPROVE with two comment
  nits (applied). Depth note: round 2 ran at the §9.1 code-diff-round-1 tier (`high`), one tier
  below the plan round's safety-boundary `xhigh` — chosen by Claude on the grounds that the
  lifetime analysis was already bought at `xhigh`, and flagged to Christian at dispatch.
- **Codex round 1, plan** (`gpt-5.6-terra`, `xhigh`): verdict "do not approve unchanged". Unique
  findings: D3 (the blocker), D2 escalation, the R3 both-path assignment form, the O2-requires-R4
  condition, the R6 "must pass, not named-error-or-pass" criterion, the O1/R3 circularity. Also
  independently confirmed: seven-reader completeness (searched `nonfree/` too), reachability, O1(c).
- **Grok round 1, plan** (`grok-4.6`, `xhigh`): verdict approve-with-corrections. Unique findings: D4
  (`mFunMKF`), the thermal-initialize window (§3.1), the C5 reachability wording, the
  no-in-tree-shuffle-target half of D5d, the O1 legal-state enumeration, the R5 domain-safety and
  `mTimeStepMatrices` side-effect check. Also independently confirmed: seven-reader completeness,
  D2's no-use half, O1(c), O2, dead-store `have_thermal()`.
- Every finding above was re-verified by Claude against the tree before inclusion; the D3
  ownership-flag calibration (no double-free on the live path) is Claude's own trace, not in
  either audit.
- Source state: shared checkout with active peer sessions. `git status` taken before each edit of
  this file; none of the target source files (`cl_FEM_Kernel.hpp/.cpp`, `cl_IWG_Maxwell.cpp`,
  `cl_ThermalFactory.cpp`, `cl_TS_TestStack.hpp`) is modified by another session as of the last
  check. Re-check before every edit.
