# Ferro Undulator: Non-Constant mu Crash and the Coil-Interface Backstop

**Date:** 2026-08-12
**Purpose:** Two defects reported by Gregory from a ferro undulator benchmark, with his two
workarounds assessed. The first is a real debug-only crash whose proposed fix introduces
undefined behaviour that silently turns the iron yoke into vacuum; the correct fix is a
short-circuit one-liner. The second is a hard error he commented out — the error is a backstop
and commenting it out is *quieter* than the abort, so the fix belongs wherever the sideset is
admitted, not at the error.
**Module:** `src/fem/kernel` (D1), `src/fem/maxwell` (D2)
**AIs involved:** Claude (exploration + plan), Codex + Grok (blind jury,
`tmp/ai_exchange/review_ferro_undulator_bugs.md`)
**Status:** OPEN — **D1 fixed in this tree 2026-08-12 (R1 applied); R6 fall-through closed the
same session**; round 2 (full jury, both auditors) confirmed both edits, no P0/P1 (§7a). R2
(telling Gregory about his D1 variant) still owed. **R3 ANSWERED by Gregory 2026-08-12 evening,
and the answer refutes the round-1 D2 reconstruction (§8): his error was never the coil-interface
aborts — it is `fn_check_facet_orientation.hpp:78`, and his model has no coils at all** (2D,
HTS thin-shell lines in the poles touching only air; blocks 2 and 3 ferro, rest air). The real
defect is **D4**: `check_facet_orientation` implements only the 3D shared-edge test and can
never succeed for two distinct 2D facets, so `fix_facet_masters`' BFS aborts on any 2D model
where a same-type facet (his tapes) is node-adjacent to a cross-type seed (his ferro-air
interfaces). **Fix applied 2026-08-12 on Christian's go-ahead (R10, §8.3), jury-audited in
round 3 (both auditors accept, no P0).** The R8 run then got past that abort and hit the next
one — **D5, §9**: the 2D cut pipeline silently requires a *uniform master side per tape*, and
R10's selective flips broke it ("No surface found on tape/boundary"). Fixed by **R13/option B**
(no same-type BFS propagation in 2D; 3D bit-identical, §9.4). §9.2 **corrects §8.2**: Gregory's
comment-out was *not* harmful in 2D — flip-all is a uniform reversal and preserves the
invariant, which is why his run passed; its real cost is 3D-scoped. The
round-1 code findings stand *as code analysis* (D2' staleness, D3 `InterfaceCondFerro` no-case,
§3.4 unreachability — corroborated by Gregory never hitting those errors), but none of them was
his bug; O3 is moot. Audit trail: `tmp/ai_exchange/review_ferro_undulator_bugs.md`.

> **Scope guards:**
> - Gregory's changes exist only in **his** working copy. Nothing here reproduces them in this
>   tree, and nothing should be committed to `BELFEM_Template` from this side.
> - D2's fix is explicitly **not** "delete or weaken the `BELFEM_ERROR`". The error is doing its
>   job; what is unknown is how a sideset got past the admission whitelist.
> - Not in scope: activating coil interfaces in assembly. `interface_node_duplication_coil_ferro.md`
>   records that coil interfaces are geometry/postprocessing-only by design.

---

## 1. Current Behaviour and How It Fails

| # | Failure | Mechanism | Evidence |
|---|---|---|---|
| D1 | **Debug-only abort** on any material with non-constant mu (i.e. every ferro with a B-H curve) | `constant_property()` is an assert-then-return, and the call site invokes it unconditionally | `cl_FEM_Calculator.cpp:96`; `cl_Material.hpp:1445-1449` |
| D1' | **Release is accidentally correct** — and that is why this survived | the assert compiles out, the getter returns the stored `NaN`, and `NaN == constant::mu0` is false, which is the right answer | `cl_Material.hpp:1180-1182` (`is_constant` *is* the not-NaN test) |
| D1'' | **Gregory's fix: silent wrong physics** — iron may behave as vacuum | `tIsConstantMu0` is left uninitialized exactly when mu is not constant, and the ferro path reads it | reads at `cl_FEM_Calculator.cpp:186,195,218,297` |
| D2 | **Hard error** reaching assembly for a coil/conductor-air interface | one of three `BELFEM_ERROR( false, "… must be disabled!" )` | `cl_IWG_Maxwell.cpp:349,354,359` |
| D2' | **Commenting it out is quieter than the abort** | an empty case body leaves `mFunMKF` holding the previously linked group's kernel; it is dereferenced unconditionally | `cl_IWG_Maxwell.hpp:41-46`; `cl_IWG_Maxwell.cpp:247` |

**Bottom line:** D1 is a one-line fix whose *proposed* remedy is more dangerous than the bug,
because the failure it introduces is silent rather than loud. D2 is not diagnosed — only the
workaround is — and the error that was removed is the only thing currently preventing a sideset
from being assembled with an unrelated weak form.

---

## 2. D1 — Non-constant mu at `cl_FEM_Calculator.cpp:96`

### 2.1 What the line does today

```cpp
bool tIsConstantMu  = mMaterial->is_constant( MaterialProperty::mu );          // :94
bool tIsConstantMu0 = mMaterial->constant_property( MaterialProperty::mu ) == constant::mu0 ;  // :96
```

`Material::constant_property` asserts and returns (`cl_Material.hpp:1445-1449`):

```cpp
BELFEM_ASSERT( this->is_constant( aProperty ), "Property is not constant" ) ;
return mConstantProperties( static_cast< size_t >( aProperty ) ) ;
```

`is_constant` is `! std::isnan( mConstantProperties( … ) )` (`cl_Material.hpp:1180-1182`), so a
non-constant property stores `NaN`.

**Why a ferro's mu is NaN in the first place** (traced by Codex; this repository's original
analysis asserted it without proof): `Material::load_bh_curve` calls `reset_dependencies( mu )`
and then `set_dependency( mu, normH )` (`cl_Material.cpp:799-806`), and **both** write
`BELFEM_QUIET_NAN` into `mConstantProperties` (`cl_Material.cpp:206-220`). Loading a B-H curve is
precisely what clears the constant slot.

Therefore:

| build | behaviour | verdict |
|---|---|---|
| debug (`USE_DEBUG=ON`, the default) | assert fires, run dies | **the crash Gregory saw** |
| release (`NDEBUG`) | assert compiled out, returns `NaN`, `NaN == mu0` is `false` | **accidentally correct** |

Confidence: **high** — mechanism read at the cited lines.

### 2.2 Why Gregory's replacement is unsafe

```cpp
bool tIsConstantMu0 ;
if (tIsConstantMu) tIsConstantMu0 = mMaterial->constant_property( MaterialProperty::mu ) == constant::mu0 ;
```

`tIsConstantMu0` is uninitialized whenever mu is not constant — the ferro case. The
`MaxwellData` constructor branches three ways:

| branch | citation | reads the flag? |
|---|---|---|
| `if ( tIsSideConnector )` | `:170` | yes — `:186`, `:195` |
| `else if ( tIsThinShell )` | `:214` | yes — `:218` |
| `else // bulk` | `:282` | **yes — `:297`** |

A ferro block is neither a side connector nor a thin shell, so it takes the bulk branch and
**reads the uninitialized flag at `:297`**. Confidence: **high** — branch structure traced, not
assumed.

What the two outcomes select (`:297-311`):

| garbage byte | `mFunMu` / `mFundMudH` | physical meaning |
|---|---|---|
| truthy | `compute_mu_0` / `compute_dmu_zero` | **the iron yoke behaves as vacuum**, Newton tangent dµ/dH zeroed |
| falsy | `compute_mu_h` / `compute_dmu_material` | correct |

"It works" on his machine means the stack byte happened to be zero. The failure is silent, can
differ between builds, machines and optimisation levels, and for an undulator it does not crash —
it yields a plausible field computed with no iron in it.

**The build will not stop it**, and the answer is compiler-dependent — the original claim here
was understated for Intel:

| toolchain | flag | catches Gregory's UB? |
|---|---|---|
| GCC (the tree's default) | `-Wall -Werror … -Wno-error=maybe-uninitialized` (`config_gcc.cmake:69`) | warns, does **not** fail |
| Clang | `-Wall -Werror=uninitialized` (`config_gcc.cmake:67`) | would fail — but only on *definite* uninitialized use, and this read is branch-dependent |
| Intel | `-Wall -Werror -Wno-uninitialized` (`config_icc.cmake:36`) | **disabled outright — worse than GCC** |

### 2.3 The fix

```cpp
bool tIsConstantMu0 = tIsConstantMu
        && mMaterial->constant_property( MaterialProperty::mu ) == constant::mu0 ;
```

Short-circuit, so `constant_property` is called only when its assert would pass; always
initialized; and the value matches today's release behaviour on every path, so nothing that works
now changes. One consequence worth noting rather than acting on: `BELFEM_ERROR( tIsConstantMu0 ||
tIsConstantMu, … )` at `:186` reduces to `tIsConstantMu`, since the new flag implies it.

---

## 3. D2 — The commented-out interface error

### 3.1 What removing the error does

`mFunMKF` is a raw function-pointer member (`cl_IWG_Maxwell.hpp:41-46`), assigned per case in
`IWG_Maxwell::link_to_group` and dereferenced unconditionally in `compute_mkf`
(`cl_IWG_Maxwell.cpp:247`), reached from assembly via `cl_IWG_Timestep.cpp:676-693`. No reset
between groups exists — confidence **high** after Codex looked independently and found none
either. An empty case body therefore leaves the pointer holding **the previously linked group's
kernel**, and that sideset is assembled with the wrong weak form, silently. Strictly worse than
the abort: it produces numbers.

**And it can be worse than stale.** The constructor's initializer list does not initialize
`mFunMKF` at all (`cl_IWG_Maxwell.cpp:36-47`), so a first-linked group hitting an empty case would
dereference an **indeterminate** pointer. In practice blocks are linked before sidesets
(`cl_FEM_DofManager.cpp:655`, then the sideset loop), so a sideset case normally inherits a real
block kernel — which is exactly what makes the failure quiet rather than a crash.

Independently of the mechanism, removing an admission-tier `BELFEM_ERROR` contradicts
`doc/coding_philosophy.md:612-627`, which puts setup and admission checks in the always-active
tier precisely so that release does not continue into a garbage state.

### 3.2 Why the error should have been unreachable — and what that implies

`MaxwellFactory::select_sidesets` (`cl_MaxwellFactory.cpp:1940-1988`) whitelists domain types
before calling `set_sidesets`. It admits `InterfaceCondFerro` and `InterfaceFerroAir`; the three
erroring types — `InterfaceCondAir`, `InterfaceAirCoil`, `InterfaceFerroCoil` — fall through
`default: // pass`. `set_sidesets` has exactly one caller in `src/` (`:1988`).

So on the normal path those types never become FEM groups, and the three errors are a **backstop**.
If one fired, something admitted that sideset anyway — and **the fix belongs at the admission,
not at the error.**

Leading hypothesis, confidence **low** and explicitly unverified: the domain type is re-derived
after selection. `Topology` derives it from the two adjacent blocks
(`cl_Topology.cpp:300-340`), so a coil adjacent to iron yields `InterfaceFerroCoil`
automatically, and a sideset admitted under a whitelisted type could later read as an
unwhitelisted one. **No re-derivation has been shown to occur.**

### 3.3 A latent hole of the same shape, found while tracing (not Gregory's doing)

`case DomainType::InterfaceFerroAir` (`cl_IWG_Maxwell.cpp:312-325`) assigns `mFunMKF` **only**
inside `if ( mFormulation == maxwell::Formulation::HPhi )`. There is no `else`, so any other
formulation falls through and inherits the previous group's kernel — the same silent-wrong-kernel
mechanism as D2'. Confidence **high** that the branch has no else; **low** on whether a non-HPhi
formulation reaches that case in practice.

A coil *block* (as opposed to interface) has no case at all and would hit the switch `default`,
which prints the domain type and raises `"Not implemented"` (`:428-432`) — distinguishable from
the three interface messages, which is useful when Gregory reports which one he saw. Note the
factory's **block** selection is whitelisted too (`cl_MaxwellFactory.cpp:648-665`: Air, Buffer,
Conductor, Ferro, ThinShell only), so a coil block reaching `link_to_group` would likewise require
a bypass of the factory.

Practical severity of the `InterfaceFerroAir` hole is **lower than it first appears**:
`maxwell_usage_guide.md:61-90` records HPhi as the solving formulation and the L2 variants as
postprocessing-only, so the unguarded path is unlikely to be reached today. It still merits
closing (R6), because nothing enforces that. *(Closed 2026-08-12; round 2 additionally
established the stronger structural argument — `IWG_MaxwellPostproc` derives from `IWG`, not
`IWG_Maxwell`, so the L2 formulations cannot reach this switch by type hierarchy.)*

### 3.4 R4 traced in-house (2026-08-12, prompted by the 2D hint) — and a new defect D3

The admission machinery, end to end, all source-traced (full citations in the exchange thread):

1. **Admission is by ID, filtered by type-at-snapshot** —
   `set_block_types_in_magnetic_equation` (`cl_MaxwellFactory.cpp:1924-1988`; round 1 misnamed
   it `select_sidesets`, which is a `Topology` method).
2. **The type `link_to_group` consumes is not the snapshot**: the fem group re-reads the mesh
   type at construction (`cl_FEM_SideSet.cpp:49`) and it is re-stamped unconditionally at
   `cl_MaxwellFactory.cpp:733-738`. The round-1 re-derivation hypothesis has real plumbing —
   two snapshots, no re-check.
3. **But no writer stamps the coil/cond-air trio after admission** (post-whitelist writers:
   `Periodic`, `Inactive` only), deck strings cannot produce interface types
   (`en_DomainType.cpp`), and the trio was **never whitelisted in any commit** since
   `a75e59a2` introduced it. **Conclusion: in this tree the three errors are unreachable —
   Gregory's abort implies version skew, or a different message than assumed.**

**D3 (new, found while tracing):** `InterfaceCondFerro` **is** whitelisted
(`cl_MaxwellFactory.cpp:1960`), gets a dof table (`cl_Maxwell_FieldList.cpp:317-322`), defaults
to full activation (`cl_IWG.cpp:2149-2161`; only Inactive/GeometryOnly/ThinShell are dormant),
is linked (`cl_FEM_DofManager.cpp:690-696`) — and has **no case** in
`IWG_Maxwell::link_to_group`, so it lands in `default:` and aborts `"Not implemented"` after
printing the type. The hanging-edge pass condenses its facets under HPhi
(`cl_MaxwellFactory.cpp:1366-1376`) but never deactivates the sideset, unlike ThinShell and
FerroAir. **A 2D undulator whose coils are declared `conductor` and sit on the iron poles hits
exactly this** — the only in-tree abort a coil-adjacent sideset can reach today. Whether this
*is* Gregory's error turns on whether his message text was "Not implemented" (D3) or verbatim
"must be disabled" (version skew); that is now the first R3 question. Mechanism confidence
**high**; reachability-in-practice **medium** (no conductor-touching-ferro deck has been run
from here).

---

## 4. Ordered Steps

- [x] **R1** — Apply the D1 short-circuit fix at `cl_FEM_Calculator.cpp:96`. One line; no
      behaviour change on any currently-working path. *Applied 2026-08-12.*
- [ ] **R2** *(after R1)* — Tell Gregory why his version is unsafe, so the uninitialized variant
      does not persist in his working copy or reach `BELFEM_Template`. **This matters more than
      R1**: our tree being right does not help if his build silently drops the iron.
- [x] **R3** — ~~Get from Gregory: which of the three error messages fired~~ **ANSWERED
      2026-08-12 evening, and the premise was wrong (§8): the error is
      `fn_check_facet_orientation.hpp:78`, not a coil-interface abort, and his model has no
      coils.** The sharpened sub-questions (message text, coil declaration, adjacency) are moot;
      his BELFEM commit hash is no longer needed for D2 — the defect reproduces from this tree's
      source by inspection.
- [x] **R4** — ~~Identify the admission route~~ **Closed via §3.4 + §8: there was no admission —
      the erroring call is `fix_facet_masters`' BFS on the cuts path, reached with perfectly
      ordinary sidesets. The admission machinery trace (§3.4) stands as documentation and
      corroborates that the coil errors never fired.**
- [x] **R9** *(independent of Gregory)* — Fix D3: give `InterfaceCondFerro` either a
      kernel case in `link_to_group` or a deactivation after its hanging-edge condensation
      (mirroring the ThinShell / FerroAir treatment) — whichever matches the intended physics.
      **Christian's ruling 2026-08-12: no second group is needed** — the h-side is coupled to
      the phi-side by the condensation, not by a weak form. *Applied 2026-08-12: **both**
      halves, since they are complementary rather than alternatives.* The sideset is now set
      `Inactive` once its condensation is done (mirroring `InterfaceFerroAir`), **and** it
      gains a `"Conductor-Ferro Interfaces must be disabled!"` case beside `InterfaceCondAir`,
      so a bypassed deactivation fails with the rule instead of `"Not implemented"`.
      Deactivation was chosen over dropping the whitelist entry because the whitelist runs
      *before* the condensation pass while the fem groups take the **live** mesh type — the
      existing idiom, with no dependence on admission ordering. The same session established
      that **there is no ferro-air interface weak form either, in 2d or 3d**, and retired the
      `phi_phi_2d`/`phi_phi_3d` stub (which called `exit( 0 )`) along with it.
      Syntax gate green; run gate is R8's successor with a conductor-touching-ferro deck.
      Session record: `devlog/dl20260812_ferroair_interface_removal.md`.
- [x] **R10** — Make `check_facet_orientation` dimension-aware: keep the shared-edge test for
      ≥3-node facets, add the 2-node chain rule on `original()->id()` (§8.3). Keep the error for
      the no-shared-node case — it stays loud and correctly flags caller bugs. *Applied
      2026-08-12 on Christian's go-ahead; `cl_MaxwellFactory.cpp` (the header's only consumer)
      passes `-fsyntax-only` under the tree's flags; **round 3 jury: both auditors accept, no
      P0** (§7c); post-round naming/comment nits applied and re-checked.*
- [ ] **R11** *(after R10; rationale CORRECTED — see §9.2)* — Tell Gregory to revert his
      comment-out. **The original rationale here was wrong and must not be sent as written:**
      in 2D his workaround flips *every* reached facet, which is a uniform reversal per
      connected component and therefore *preserves* the per-tape uniform-master invariant —
      that is why his run passed the cut pipeline. Its real defect is 3D-scoped (a
      no-shared-edge 3D pair silently flips instead of erroring) plus the loss of a genuine
      diagnostic.
- [x] **R13** *(new — D5, the regression the R8 run exposed)* — Disable same-type BFS
      propagation in 2D in `fix_facet_masters` (option B). *Applied 2026-08-12 on Christian's
      choice; four-site `tPropagate` guard; 3D bit-identical; syntax gate green (§9).*
- [ ] **R12** *(low, found while diagnosing)* — The 3D loop in `check_facet_orientation`
      compares raw `node(i)->id()`, not `original()->id()`; on meshes with duplicate nodes
      (thin shells, periodicity) twin edges may not match. Unverified; check when touching R10.
- ~~**R5** — Fix at the admission point.~~ **Struck 2026-08-12: predicated on the refuted
      admission hypothesis; superseded by R10.**
- [x] **R6** — Close the `InterfaceFerroAir` fall-through (§3.3) with either an `else` that
      errors or an assignment, whichever matches the intended formulation coverage. Independent
      of R3–R5. *Closed 2026-08-12 with an `else` + `BELFEM_ERROR`, matching the HPhi-only
      solving coverage recorded in `maxwell_usage_guide.md`.*
- [ ] **R7** — Consider a `mFunMKF = nullptr` reset at the top of `link_to_group`, converting
      every present and future fall-through from a silent wrong kernel into a null deref at the
      call site. See O2.
- [ ] **R8** *(gate)* — Run the undulator benchmark in a **debug** build, which is where D1
      aborts and where any new assert would fire.

---

## 5. Open Design Questions

- **O1** — Should `Material::constant_property` be hardened so this class of misuse is not
  per-caller? Today a release build silently returns `NaN` to any caller that forgets the
  `is_constant` guard, and D1 shows the guard is easy to forget. Options: leave as is (hot-path
  accessor, assert is the documented contract); return `NaN` deliberately with the contract
  documented at the declaration; or split into `constant_property()` (asserting) and
  `constant_property_or_nan()`. **Not decided — Christian's call**, and it touches every material
  consumer.
- **O2** — Is a `mFunMKF = nullptr` reset at the top of `link_to_group` worth the null deref it
  would produce? It converts a silent wrong-kernel assembly into an immediate crash at
  `cl_IWG_Maxwell.cpp:247`. Argues for: the current failure mode is undetectable. Argues against:
  a null function-pointer call is a segfault, not a diagnosable error box — a `BELFEM_ERROR` on
  `mFunMKF == nullptr` in `compute_mkf` would be better but adds a per-element branch, so the
  tier rule (`doc/coding_philosophy.md`) points at `BELFEM_ASSERT`, which vanishes in release
  exactly where the silence hurts. **Not decided.**
- ~~**O3** — Are coil interfaces expected to appear at all in an undulator deck?~~ **Moot
  (2026-08-12, §8): Gregory's model has no coils. The question was an artifact of the refuted
  D2 reconstruction.**
- **O4** *(reworded 2026-08-12 after round 3 — the original wording was REFUTED by both
  auditors independently)* — The claim "consistency within each chain holds regardless of which
  seed reaches it first" is true **only for endpoint contacts** (Gregory's geometry: tape tip
  on a pole node). For a seed touching an **interior** node Q of a same-type chain, the BFS
  orients both arms from itself — one arm flips, both end up starting at Q, masters on opposite
  sides — and the arms never compare against each other (both unflagged on first visit): **the
  chain is torn at Q**. Hand-traced against `cl_MaxwellFactory.cpp:1281-1300`. The 3D
  non-manifold analogy does not carry (3D adjacency is by shared edge; a 2D seed can bisect a
  chain in a way no 3D edge-neighbor can). This is a **latent P1 in the BFS + node-based
  adjacency**, pre-existing and merely *exposed* by R10 making 2D propagation reachable; not
  fixable inside `check_facet_orientation`. Candidate fixes: restrict 2D BFS adjacency to
  same-sideset neighbors, or post-BFS chain-consistency sweep per sideset. **Christian's call;
  no known deck hits it today** (needs a cross-type sideset ending on the interior of a
  same-type chain).

---

## 6. Definition-of-Done Checklist

- [ ] D1 fixed in this tree and communicated to Gregory (R1 ✓, R2 open).
- [x] D2's real identity established with citations (R3, R4): **D4,
      `check_facet_orientation` is 3D-only**, reached from `fix_facet_masters`' BFS in 2D (§8).
- [ ] D4 fixed (R10) and Gregory's comment-out reverted (R11).
- [x] D3 (`InterfaceCondFerro` admitted-but-no-case) fixed per Christian's direction (R9) —
      deactivated after condensation **and** given a "must be disabled!" case.
- [x] The three coil-interface `BELFEM_ERROR`s still present and still unreachable (confirmed
      §3.4; re-confirmed after R9/R10 — the whitelist still admits none of the trio, and R9
      added a fourth error of the same tier for `InterfaceCondFerro`). Note the count is now
      **five**, not three: `InterfaceFerroAir` joined them when R9 retired its dead kernel.
- [x] `InterfaceFerroAir` fall-through closed (R6).
- [ ] D5 fixed (R13 ✓) and the undulator gets past `orient_terminal_curves_2D` (R8).
- [ ] O1, O2, O4, O5 decided or explicitly deferred, not silently resolved (O3 moot).
- [ ] Undulator benchmark runs in a **debug** build (R8) — now gates **both** D1 and D4: it
      aborts at `fix_facet_masters` today, and after R10 it is the reproducer for the fix.

---

## 7. Audit Trail

### 7c. Round 3 (2026-08-12 evening) — audit of the applied D4 fix (R10)

- Pre-registration frozen (6 findings + 4 refutation targets, author bias declared); scoped
  subject = the one hunk; blind `--jury`, both legs returned. **Both auditors accept the hunk;
  no P0.** All three load-bearing properties confirmed independently three ways: the 2D truth
  table (incl. coincident twins reproducing the old rule), `node(0)→node(1)` = master-side
  traversal with `flip()` reversing it on a CCW pair, and `original()` as the right comparison
  space (usually identity at this call site — duplicates don't exist yet at `:937`).
- **The round's catch, raised independently by BOTH auditors: my O4 wording was wrong.** A seed
  at an *interior* node of a same-type chain tears the chain (arms oriented from the seed, never
  compared with each other). Endpoint contacts — Gregory's case — are fine. O4 reworded; latent
  P1 in the BFS, pre-existing, not part of R10.
- Post-round polish, labelled: `tA0/tA1/tB0/tB1` rename (both auditors: `a`-prefix collision)
  and comment wording ("two-node facets", Grok C2). Syntax gate re-run green.
- Verification pass + reconciliation table in the exchange thread. Highest evidence: source
  trace — **reviewed, not verified**; R8 is the gate.

### 7a. Round 2 (2026-08-12 PM) — audit of the applied R1/R6 fix

- New pre-registration frozen (8 findings + 4 refutation targets, author bias declared), scoped
  subject file excluding the unrelated `hcur`/`hhist` hunk; blind `--jury`, **both legs
  returned** (the round-1 Grok failure was brief-specific — the narrowed subject ran clean).
- **Both auditors confirm both edits; no P0, no P1.** Codex: value table over all four mu
  cases, no sibling unguarded `constant_property(mu)`, tier correct. Grok: same conclusions
  independently, plus the sharpest contribution of the round — **a refutation of the round-1/
  pre-reg framing**: the ferro *assembly* kernels (`phi_ferro_newton/picard`) call
  `Material::dmudH` directly and never consumed the two flags; R1 unblocks the *constructor*,
  and independently sets the right *helper* pointer (`compute_mu_h`) for postproc/thermal
  consumers. Accepted; recorded here rather than papered over.
- Grok also flagged the `:189` guard as now-tautological (optional cleanup, not wrong) and
  corrected round 1's function naming (`set_block_types_in_magnetic_equation`, not
  `select_sidesets`).
- Verification pass re-checked every load-bearing citation by content; none refuted.
  Reconciliation table in the exchange thread. Highest evidence rung: source trace — **R8 (the
  debug undulator run) is still the only executable gate.**

### 7b. Round 1 (2026-08-12 AM) — audit of the diagnosis

- Exchange thread: `tmp/ai_exchange/review_ferro_undulator_bugs.md` — pre-registration frozen
  before dispatch; blind `--jury` round dispatched 2026-08-12; verification pass and
  reconciliation table appended on return.
- **Codex:** confirmed all seven findings, refuted none, and added three things the original
  analysis lacked — the `load_bh_curve` trace explaining *why* mu is NaN (§2.1), the three-way
  compiler-flag picture including Intel's outright `-Wno-uninitialized` (§2.2), and the
  observation that `mFunMKF` is never initialized in the constructor, which upgrades the D2'
  confidence from medium to high (§3.1). It also raised the coding-philosophy argument against
  removing an admission-tier `BELFEM_ERROR`. One of its citations had drifted the same day
  (`assert.hpp` is now gated on `BELFEM_ASSERTIONS_ACTIVE`, not `#ifndef NDEBUG`, after
  `5cc9d6dd`); the conclusion is unaffected.
- **Grok: failed twice** — the jury leg and a solo retry, both exit 1 with an empty body, both
  stopping at `attempt 1/5` before the quality gate. A smoke test through the same wrapper
  returned cleanly, so this is specific to this brief or round, not an outage; the plausible cause
  is exceeding `GROK_MAX_TURNS` (30) across the ~10 files the brief points at, but the log does
  not say so and it is not asserted here. **This round therefore rests on one auditor plus the
  in-house trace.** That is weaker than a jury and is recorded rather than glossed — though every
  verdict is source-traced, which ranks above auditor agreement on the evidence ladder. If a third
  voice is wanted, a re-run with a narrowed brief (D1 only) is the cheap version.
- The finding most wanted from the round was a **refutation**: a demonstration that the bulk
  branch is not reached for ferro, or that the flag is never read uninitialized, which would make
  Gregory's patch safe. None was produced. Declared bias: the diagnosis was already reported to
  Christian before the round was dispatched.
- All auditor citations were re-checked against the source by content before inclusion here.
- Source of the report: Gregory's text message, not his diff. The code analysis stands on the
  citations above; the reconstruction of *his* intent does not, and R3 exists because of that.

---

## 8. 2026-08-12 evening — R3 answered: D2 was misidentified; the real defect is D4

Gregory's answers: the error is **`fn_check_facet_orientation.hpp:78`** ("Could not determine
orientation of facets."), apparently called for facet pairs whose blocks are not neighbors; he
commented the line out and the run proceeded. And: **his model has no coils** — 2D, the pole
lines are HTS thin shells touching only air, blocks 2 and 3 ferro, everything else air. The
round-1 D2 reconstruction (coil-interface aborts) is therefore refuted by the reporter; §1-§3
stay as written because D2'/D3/§3.4 are real code findings — they were just never his bug.

### 8.1 D4 — `check_facet_orientation` is 3D-only, and 2D always aborts [high, source-traced]

The function looks for a shared **directed edge** — a shared node *pair* — between the two
facets' corner-node cycles (`fn_check_facet_orientation.hpp:52-77`). A 2D facet is a line
segment whose "cycle" is the segment traversed both ways, so two **distinct** 2D facets would
need to share *both* nodes to match; adjacent segments share exactly **one** → fall-through →
the `BELFEM_ERROR` at `:78`, unconditionally.

Its only caller is the BFS in `MaxwellFactory::fix_facet_masters`
(`cl_MaxwellFactory.cpp:1298`; rank 0, cohomology path, `:937`, before cuts exist):

| stage | what | Gregory's model |
|---|---|---|
| seeds | facets whose master/slave **block** types differ, orientation from type priority (`:1262-1279`) | ferro-air interfaces of blocks 2/3 |
| flags | facets with equal types both sides (`:1275-1278`) | every tape line (air\|air) |
| BFS | propagate over facet-to-facet adjacency, which **in 2D is node-based** — one shared node, originals *and duplicates* (`cl_Mesh_ConnectivityCalculator.cpp:375-406`) | tape endpoint touches a pole surface node |
| step | `check_facet_orientation( seed, tape )` → distinct segments, one shared node | **abort** |

The node-based adjacency crossing chains at a point contact is also exactly his "blocks that
were not even neighbors" observation. Structural consequence: **any 2D cohomology model where
a same-type sideset facet is node-adjacent to a cross-type facet dies here** — nothing specific
to his deck.

### 8.2 Why the comment-out is silently wrong [high]

With `:78` removed the function falls through to `return false`, so the BFS **flips every
flagged facet it reaches** (each exactly once — unflagged at first visit). `Facet::flip()`
swaps master/slave (`cl_Facet.cpp:89-98`), and stage 2 of the thin-shell pipeline hangs the
duplicate nodes on the **slave** side (`thin_shell_facet_orientation.md`, stage 2). A uniform
flip of all reached facets *preserves* the relative orientation of the element-ID-based initial
assignment — consistency along each tape is **not established, merely inherited from gmsh
element numbering**. "It runs fine" means his mesher happened to number elements consistently:
the same it-works-by-accident shape as his D1 patch, mesh-dependent instead of stack-dependent.

### 8.3 Proposed fix (R10) — not applied, Christian's call

Dimension-aware `check_facet_orientation`: keep the shared-edge loop for ≥3-node facets; for
`n == m == 2` compare on `original()->id()` (the 2D adjacency pass matches
originals+duplicates, so originals are the right comparison space):

```cpp
// A = (a0 → a1), B = (b0 → b1), ids via original()
if ( a1 == b0 || a0 == b1 ) return true ;   // B continues A's traversal; also reversed twin
if ( a0 == b0 || a1 == b1 ) return false ;  // both start or both end at the shared node; identical twin
BELFEM_ERROR( ... )                          // no shared node: caller bug, keep loud
```

The two twin cases (coincident segments) reproduce the 3D semantics exactly; the chain cases
are the correct 2D traversal-consistency rule; junction nodes (degree ≥ 3) stay pairwise-defined,
as at non-manifold edges in 3D. Found alongside, filed as R12: the 3D loop compares raw
`node(i)->id()`, not originals — potential duplicate-blindness on enriched/periodic meshes,
unverified.

---

## 9. D5 — the R8 run's next failure, and what it revealed (2026-08-12, late)

The debug undulator run got past `fix_facet_masters` (R10 worked) and died further down:
`BELFEM_ASSERT` "No surface found on tape/boundary" (`cl_CutFactory.cpp:2210`,
`orient_terminal_curves_2D`). Progress, then a regression — and **both round-3 auditors had
flagged this exact class** (selective 2D propagation at multi-contact points). "No known deck
hits it today" lasted a few hours.

### 9.1 Mechanism [high, traced end to end]

**The hidden invariant:** the 2D thin-shell cut pipeline requires a **uniform master side per
tape sideset**. Nothing establishes it — it is inherited from element numbering
(`connect_facets_to_elements` takes the lower element ID as master), and gmsh writes elements
grouped per physical surface, so the flanking air blocks land in disjoint ID ranges.

**The consumer:** `relink_slave_elements_with_duplicate_nodes` (`cl_CutFactory.cpp:1763`)
collects the slave-side **blocks** over all facets of a tape. Uniform masters → one slave block
→ only that side is relinked to duplicates. **A single flipped facet puts both blocks in the
set → both sides get duplicates → the original nodes are stranded.** `duplicate_and_relink_facets`
(`:1813`) then rebuilds the live sideset from both sides, and `orient_terminal_curves_2D`
(`:2156`) searches it for a facet carrying the terminal segment's *original* nodes (flag-based
identity, `:2196-2207`). None exists → the assert.

**Why R10 triggered it:** the chain rule made 2D propagation *reachable*, and it flips
selectively — mixed sides on one tape are guaranteed once a tape or its feeding chains touch
seeds at more than one point. In an undulator (many air blocks, ferro poles, tapes inside the
poles) that is the normal case, not a corner case.

### 9.2 Correction to §8.2 — Gregory's workaround was not harmful in 2D

§8.2 said his comment-out left consistency "inherited from gmsh numbering". The observation was
right; the verdict was wrong. With the error removed the function returns `false` for **every**
distinct 2D pair, so the BFS flips **every** reached facet exactly once — a uniform reversal per
connected component, which **preserves** the per-tape invariant. That is precisely why his run
cleared the cut pipeline while a pairwise-*correct* rule broke it: flip-all and flip-none both
preserve uniformity; only selective flips can destroy it. His workaround's real cost is 3D
(a no-shared-edge pair silently flips instead of erroring) and the lost diagnostic. R11 carries
the corrected rationale.

### 9.3 The design fact underneath

In 2D **no consumer benefits from same-type propagation**: tapes need uniformity (which
propagation can only break), nothing reads air-air block-boundary masters, and cut sidesets are
built later by CutFactory under its own conventions. And no 2D model ever ran propagation
successfully — pre-R10 it aborted at the first propagation step, structurally. R10 opened a door
into a room that was never built.

### 9.4 Fix applied (R13, option B)

`fix_facet_masters` gains `tPropagate = ( number_of_dimensions() == 3 )`, guarding four sites:
the facet-to-facet connectivity build, the same-type flagging, the BFS loop, and the
connectivity teardown. Cross-type seed normalization by domain-type priority still runs in 2D —
the only part 2D consumers ever depended on. **3D is bit-identical**; 2D returns to exactly its
historically-working behaviour and additionally skips a connectivity build it never needed.
The R10 chain rule stays correct and in place, simply unreached in 2D, and would be required
if 2D propagation is ever deliberately designed.

Not sent to a jury round: protocol §11 points at the executable gate, not another review — this
is a four-site guard restoring prior behaviour, the mechanism is traced, and the auditors had
already named the failure class. **Reviewed, not verified; R8 is the gate.**

### 9.5 Left open deliberately

**O5 (new, robustness):** the per-tape uniform-master invariant is still implicit and still
rests on gmsh numbering. Option (C) — enforce or assert it in CutFactory before duplication —
would have turned today's cryptic downstream assert into a one-line diagnostic at the true
cause. Recommended as a follow-up, **not applied**: it is a real behavioural guard on the cut
pipeline and deserves its own round rather than being folded into a hotfix. **Christian's call.**
