# Reorder the Field-Dependent `rho` / `lambda` API to B-First

> **CANCELLED 2026-09-03** (todo/ currentness sweep, round 3): Christian's ruling 2026-09-03: keep T-first. The tree is uniformly `(T, B, beta)`, the document inconsistency that motivated B-first is fixed, and the reorder would have broken the overload family's prefix consistency at 17 one-argument call sites. Status lines and checkboxes below are as they stood at closure and are not maintained.

**Date:** 2026-08-28
**Purpose:** Make the three-argument field-dependent material accessors uniformly
**field-first** — `rho( B, angle, T )` and `lambda( B, beta, T )` — replacing today's
`( T, B, beta )`. The goal is a convention a reader cannot misremember per function; the
mechanism is a two-pass reorder in which the intermediate state **fails to compile** at
every positional call site, so the change is compiler-enforced rather than trust-based.
**Module:** `src/physics/materials` (primary), `src/fem/kernel` (four call sites),
`src/physics/materials/doc` + `src/fem/*/doc` (documentation)
**AIs involved:** Claude (inventory + plan), Codex (audit), Grok (third voice)
**Status:** PLAN — drafted 2026-08-28, **audit dispatched, pending approval; no source
modified.** Christian's decision (2026-08-27 evening): B-first for both, "I think this is
the safest". Explicitly deferred out of that session — "this is too hot for tonight".
**Claude's plan-stage recommendation is AGAINST the reorder** (Appendix B, written at
Christian's request): the tree is already uniform, the inconsistency is in two documents
rather than in the code, and B-first would break the overload family's prefix consistency
at 17 one-argument call sites. Recorded for the auditors to attack; the decision above
stands until Christian revisits it, and O4 gates any implementation regardless.
**Audit round 1 (Codex) returned 2026-08-28: NOT implementation-ready** — inventory
undercounted (8 public call sites not 4, four Kohler pointers not one, `UserDefined` a third
overrider), and the change is invisible to the type system, so `override`, the pointer type
and the plugin ABI give none of the protection §2 claimed. Findings D1–D4 in §3.2; R2/R5/R6
need rework before any go. Claude has **withdrawn** the against-reorder recommendation on
Codex's evidence; the surviving recommendation is the distinct-name variant (Appendix B).
**RE-PLANNED 2026-08-28 against both audits** (Christian: "go ahead"): §1 now carries the
true eight-function family and ten call sites, §2's method deletes the whole family at once,
§3's gap table gains the private Kohler/table layer and the `Database::evaluate` trap, §3.1
A and C are rewritten, and §4 is re-ordered around a blocking **R0** (freeze ruling) and a
new **R5'** — the hand-enumerated private-layer pass, which is the step INC-115 failed in
July. Remaining before implementation is approvable: **O4/R0 is Christian's**, O1 and O2
have auditor recommendations awaiting a ruling, and the re-planned §4 has not itself been
audited.
**Grok round returned 2026-08-28, same verdict independently** (§3.3): not
implementation-ready, inventory incomplete in a load-bearing way. Combined result — the
public family is **eight functions, not four**, with **ten** positional call sites, not four;
four Kohler pointer forwards plus a parallel table family that R2 structurally cannot reach;
and one Claude citation error (D7 — INC-286 for INC-115, the register's own trap #2). Both
auditors independently reject Appendix B's empirical half; both decline to treat it as
grounds to reverse the decision. **Next action is a re-plan, not an implementation:** R2/R5/R6
need rework against the true family, and O2 (powerlaw) and O4 (freeze) are still open.

> **Scope guards:**
> - **No code tonight.** This file is a plan. Implementation needs Christian's separate go.
> - **In scope:** the 3-arg `rho` / `drhodT` / `lambda` / `dlambdadT` family on
>   `Material` and its two overriders, their four positional call sites, the Kohler
>   function-pointer forward, and the documentation that teaches the order.
> - **OUT of scope, and must not be touched — two traps for any mechanical sweep:**
>   (i) `Gas::rho( T, p )` and `Gas::lambda( T, p )` in `src/physics/gasmodels` and
>   `src/physics/gastables` — a different hierarchy under the frozen thermophysical symbol
>   table (`CLAUDE.md`, "Canonical thermophysical symbols"); a textual sweep on `rho(` /
>   `lambda(` hits ~15 of them. (ii) `Database::evaluate( theta, log10B, angle )` inside
>   `Alloy::rho` and `Metal::rho_table` — three `real`s in **lookup-table axis order**, which
>   has nothing to do with the public API and must not be permuted to match it (D8).
> - **OUT of scope:** all 1-argument forms (`rho(T)`, `lambda(T)`), which have no
>   ordering to get wrong; and the `powerlaw` / `piecewise` derivative family, pending O2.
> - This **reverses** the 2026-07-11 option (c) ruling that unified the family on T-first.
>   `rho_lambda_argument_convention.md` and its audit trail stop being current the day
>   this lands.

---

## 1. Current Behaviour and How It Fails

**The family is eight public functions, not four** (corrected 2026-08-28 after D1/D5; the
original four-function table was the error that made the rest of the plan wrong). All are
`( T, B, beta )` today, all three parameters `real`:

| symbol | base decl | base def | Metal | Alloy | live caller |
|---|---|---|---|---|---|
| `rho( T, B, beta )` | `:670` | `:2015` | `hpp:183` / `:576` | `hpp:80` / `:145` | `Calculator.hpp:2984`, `main.cpp:88,385` |
| `drhodT` | `:673` | `:2027` | `hpp:186` / `:594` | `hpp:83` / `:157` | `Calculator.hpp:3431` |
| `drhodB` | `:676` | `:2042` | `hpp:189` / `:600-602` | `hpp:86` / `:173` | `Calculator.hpp:3480` |
| `drhodbeta` | `:679` | `:2051` | `hpp:192` / `:606-608` | `hpp:89` / `:190` | `Calculator.hpp:3499` |
| `lambda( T, B, beta )` | `:616` | `:1881` | `hpp:204` / `:662` | `hpp:102` / `:207` | `Calculator.hpp:2910`, `main.cpp:94,402` |
| `dlambdadT` | `:619` | `:1894` | `hpp:207` / `:668` | `hpp:105` / `:213` | `Calculator.hpp:2929` |
| `dlambdadB` | `:622` | `:1907` | `hpp:210` / `:684-690` | not overridden | **none** |
| `dlambdadbeta` | `:625` | `:1920` | `hpp:213` / `:694-700` | not overridden | **none** |

Plus the 4-arg HTS pair `lambda`/`dlambdadT( T, B_par, B_perp, J )` (`:636,:639` / `:1936,
:1948`), which has no caller anywhere — see O1.

**Ten public positional call sites**, all verified by direct read 2026-08-28: six in
`cl_FEM_Calculator.hpp` (`:2910, :2929, :2984, :3431, :3480, :3499`) and four in
`src/physics/materials/main.cpp` (`:88, :94, :385, :402`).

**Behind them, a private layer R2 cannot reach** (D6): four Kohler member pointers
(`cl_Material.hpp:387-390`) forwarded at `cl_Material_Metal.hpp:578, 596, 602, 608`, their
pointees `rho_kohler`/`drhodT_kohler`/`drhodB_kohler`/`drhodbeta_kohler`
(`cl_Material.hpp:1496-1505`, defs `cl_Material.cpp:560-584`), a parallel `*_table` family
(`:1508-1517`), and `UserDefined`'s `rho_kohler`/`lambda_custom` overrides
(`cl_Material_UserDefined.hpp:310, 324`). Every one of them is `( T, normB, angle )` and
every one is three `real`s.

The failure this plan addresses is not a live wrong answer — it is a **standing hazard
plus a documentation split**:

| Failure | Mechanism | Evidence |
|---|---|---|
| Silent argument swap | every parameter is `real`, so a wrong order compiles clean and computes garbage | INC-114: `Metal::lambda(B,beta,T)` silently received T in the \|B\| slot from four call sites for months (`doc/lessons_learned_evidence.md`, card INC-114) |
| Two conventions in the tree's prose | the materials docs teach `(T,B,beta)`, two fem guides teach `(B,angle,T)` | `dl20260825_materials_doc_overhaul.md` §1 vs `dof_manager_usage_guide.md:3118` |
| Per-function rule is unmemorable | an asymmetric convention (rho field-first, lambda T-first) has to be recalled correctly at each call | Christian, 2026-08-27 |

**Bottom line:** the current order is safe only while nobody edits a call site, and the
one artifact class that would teach a newcomer the order — the module documentation —
currently teaches two different answers. A single uniform rule removes both problems.

## 2. Architecture: Why a Two-Pass, Compile-Breaking Reorder

The naive approach — edit the signatures and then hunt call sites by grep — is exactly the
approach that produced INC-114. Three `real` parameters give the compiler nothing to catch,
and two routes are invisible to a textual search: the `mFunctionRhoKohler` pointer
(§3.1 A) and any call reached through a `Material *` base pointer.

**Chosen spine: delete before you add — but it must delete the whole family at once.**
Pass 1 removes **all eight** public 3-arg overloads, keeping the 1-arg forms. Every direct
positional call then fails with `no matching function for call to 'rho(real, real, real)'`
— a hard error at every such site, including any this inventory still misses. Pass 2
reintroduces them B-first and each error is fixed at a known location.

> **REWORKED 2026-08-28 after D5/D6.** The original version deleted only four of the eight,
> which left `Metal::drhodB` and `Metal::drhodbeta` alive and still calling 3-arg `rho`
> internally — so Pass 1's error list would have looked complete while three of the four
> Kohler pointer paths and two live Calculator tangent callers sailed through untouched.
> **A partial deletion is worse than none:** it produces a confident, wrong inventory.
> The private layer (pointers, pointees, `*_table`, `UserDefined` overrides) is **not**
> reachable by this trick at all — same arity, different names — and needs the separate
> structural pass R5' below. This is exactly how INC-115 failed in July.

~~The `override` keyword does the same job for the base/derived pairs: a base signature
that moves without its overrides is a compile error, not a silent new overload.~~
**FALSE — retracted 2026-08-28 (D2, Codex, confirmed by Claude against the language rule).**
Reordering parameters that all have the same type does not change the function type at all:
`real(real,real,real) const` before and after. So `override` still binds, the member-pointer
types stay compatible, and the mangled symbols are identical. **The compiler cannot detect
this change anywhere.** R2's arity trick works — and it is the *only* mechanical protection
in the plan; everything downstream of it is structural review, not enforcement. See D2 for
what this costs the method.

Rejected alternatives:

- **Edit in place, verify by grep.** Rejected: this is the INC-114 method, and it has
  already failed once in this exact family.
- **Strong parameter types (`struct Tesla`, `struct Kelvin`).** Would make the whole class
  of error impossible, permanently. Rejected for now as a far larger change to a plugin
  API, and it conflicts with the zero-abstraction-penalty rule unless done carefully.
  Logged as O3 rather than dismissed.
- **Deprecate-and-alias (keep both orders for a release).** Rejected: two live orders is
  the failure mode, not the cure, and every argument is `real` so the compiler cannot
  route a caller to the right one.

## 3. Gap Table

| # | State / site | Needed for | Handled today? | Class | Citation |
|---|---|---|---|---|---|
| 1 | base 3-arg `rho`/`drhodT` | the reorder itself | T-first | (c) | `cl_Material.hpp:670,673` / `:2015,:2027` |
| 2 | base 3-arg `lambda`/`dlambdadT` | the reorder itself | T-first | (c) | `cl_Material.hpp:616,619` / `:1881,:1894` |
| 3 | base 4-arg HTS `lambda`/`dlambdadT` | consistency | T-first, **no caller anywhere** | (b) | `cl_Material.hpp:636,639` / `:1936,:1948` → O1 |
| 4 | `Metal` overrides | virtual dispatch | T-first | (c) | `cl_Material_Metal.hpp:183,186,204,207` |
| 5 | `Metal` definitions | — | T-first, **spelled `real T` not `const real T`** | (c) | `cl_Material_Metal.hpp:576,594,662,668` |
| 6 | `Alloy` overrides | virtual dispatch | T-first | (c) | `cl_Material_Alloy.hpp:80,83,102,105` |
| 7 | `Alloy` definitions | — | T-first, same `real T` spelling | (c) | `cl_Material_Alloy.hpp:145,157,207,213` |
| 8 | ~~`mFunctionRhoKohler` pointer type~~ **four** Kohler pointers | metal field-dependent rho + its three derivatives | all `( T, normB, angle )` | (c) | `cl_Material.hpp:387-390` (D6) |
| 9 | ~~Kohler forward call~~ **four** forwards | — | positional, **R2 cannot reach any of them** | (c) | `cl_Material_Metal.hpp:578, 596, 602, 608` (D6) |
| 9b | Kohler **pointees** + `*_table` twins | what the pointers name | `( T, B, beta )`, different names, same arity | (c) | `cl_Material.hpp:1496-1505`, `:1508-1517`; `cl_Material.cpp:560-584`; `cl_Material_Metal.hpp:435-456` |
| 9c | pointer **assignments** | wiring | 4 kohler + 4 table | (c) | `cl_Material_Metal.cpp:47-50`, `:925-928` |
| 9d | `UserDefined::rho_kohler` / `lambda_custom` | third overrider of the private 3-arg virtuals | `( T, normB, angle )` | (c) | `cl_Material_UserDefined.hpp:310, 324` (D1c) |
| 10 | ~~Calculator call sites (4)~~ **ten public positional call sites** | the callers | `( T, norm_b, mBeta )` | (c) | `cl_FEM_Calculator.hpp:2910,2929,2984,3431,**3480**,**3499**`; `src/physics/materials/main.cpp:88,94,385,402` (D1a, D5) |
| 10b | internal self-calls | Metal/Alloy quotient rules | `this->rho(T,B,beta)` etc. | (c) | `cl_Material_Metal.hpp:664,672,676,688,690,698,700`; `cl_Material_Alloy.hpp:209,217,221` |
| 11 | `MaxwellData::compute_rho/​compute_lambda` | kernel access | index-based facade, **immune** — but its *callees* `compute_rho_metal`/`compute_lambda_metal` are row 10 | (a)/(c) | `mt_maxwell_h.cpp:135`, `mt_thermal_h.cpp:35` |
| 12 | materials module docs (7 files) | teaching the order | rewritten to T-first 2026-08-25 | (c) | `dl20260825_materials_doc_overhaul.md` §1 |
| 13 | `dof_manager_usage_guide.md` | teaching the order | already B-first — **the live footgun today** | (a) | `:3118,:3154,:3205,:3269` |
| 14 | `maxwell_usage_guide.md` | teaching the order | already B-first | (a) | `:841,:952,:958` |
| 15 | `cl_Material.hpp` doc-comment | teaching the order | already B-first | (a) | `:1200` |
| 16 | ~~user-material plugin API~~ | ~~out-of-tree `.so` compatibility~~ | **path is DEAD** — 3-arg registration wires the 1-arg pointers; public 3-arg asserts `PureMetal` | (a) | D9; `cl_Material.cpp:384-394`, `cl_Material.hpp:1887,2020`, INC-116 |
| 16b | in-source plugin doc split | plugin authors | rho/lambda documented T-first, jc/n B-first, in one header | (c) | `cl_Material_UserDefined.hpp:199-200`; contract enforced `cl_Material_UserDefined.cpp:141-151` vs `:173-183` |
| 17 | `nonfree/` | — | **zero 3-arg call sites** | (a) | verified independently by Claude and Grok |
| 18 | `tests/physics/test_YBCO.cpp` | — | 1-arg `lambda(T)` only | (a) | `:99,:109,:115,**:141**` |
| 19 | `Database::evaluate( theta, log10B, angle )` | table lookup inside `Alloy::rho`, `Metal::rho_table` | 3 `real`s in **grid-axis order** — must NOT be permuted | (c) | D8 — scope guard |

### 3.1 Cross-cutting findings

**A — The private Kohler/table layer is the most dangerous part, it has already failed once,
and no mechanical device in this plan reaches it.** *(Rewritten 2026-08-28 after D6/D7; the
original text named one pointer and trusted the numeric gate, and both were wrong.)*
There are **four** pointers (`cl_Material.hpp:387-390`), **four** forwards
(`cl_Material_Metal.hpp:578, 596, 602, 608`), **eight** pointees across the `*_kohler` and
`*_table` families (`cl_Material.hpp:1496-1517`), **eight** assignments
(`cl_Material_Metal.cpp:47-50, 925-928`), and a third overrider in `UserDefined`. Every one
takes three `real`s. Consequences:
- Pass 1 cannot see any of it — these are not calls to `rho`, they are calls through a
  pointer or to a differently-named virtual of the same arity.
- "Reorder the pointer type" is not a type-level operation (D2): the before and after
  member-pointer types are identical, so nothing is checked.
- The numeric gate cannot certify it either (D3): a public-B-first / private-T-first split
  reproduces the old numbers exactly.
- **INC-115 is this exact failure, already realised.** When Christian unified the family on
  T-first in July 2026, the Kohler/table pointer path is what silently swapped. The plan's
  earlier citation of INC-286 here was wrong (D7).
Therefore the private layer needs a dedicated, enumerated, reviewed pass (R5'), not a
mechanism. Treat "the compiler is silent" as the expected state, not as evidence.

**B — The definitions are spelled differently from the declarations.** `Metal` and `Alloy`
declare `( const real T, … )` but define `( real T, … )` (`cl_Material_Metal.hpp:662,668`;
`cl_Material_Alloy.hpp:207,213`). A find-replace keyed on the declaration text silently
misses four definitions. Reorder by hand or by a pattern that matches both spellings.

**C — The plugin surface is NOT a live hazard. *(Rewritten 2026-08-28 — the original claim
was wrong in the dangerous direction, and so was Codex's correction of it.)*** Three
findings had to be reconciled: Claude said "ABI break, stale `.so` registers nothing";
Codex said "not an ABI break at all — identical mangled symbols, so a stale `.so` loads and
silently computes garbage"; Grok showed **the path is dead** and is right. The 3-arg user
callback is unreachable: 3-arg registration only `set_custom`s, which wires the **1-arg**
`mFunctionRho`/`mFunctionLambda` (`cl_Material.cpp:384-394`), while public 3-arg
`rho`/`lambda` assert `MaterialType::PureMetal` (`cl_Material.hpp:1887, 2020`) and
`UserDefined` is `MaterialType::UserDefined`. That is INC-116. Codex's mechanism is correct
C++ and unreachable here; Claude's "registers nothing" was DR-65's init-symbol issue, not a
consequence of this reorder.

What survives, and is real:
- Out-of-tree *callers* of the public 3-arg `rho` need a rebuild **and** a source rewrite —
  and the rebuild alone will not tell them, because the symbol is unchanged.
- `UserDefinedMaterial::set_user_defined_function` hard-errors on the declared dependency
  order for 3-arg rho/lambda (`cl_Material_UserDefined.cpp:141-151`: `Dependency1 == T`,
  `Dependency2 == normB`, `Dependency3 == angleBxJ`) — a loud named gate that **must** be
  updated with the reorder, or the enforced contract and the internal convention diverge.
  It guards a dead path today, but it is the contract plugin authors read.
- The same header already documents a **split** convention — rho/lambda T-first, jc/n
  B-first (`cl_Material_UserDefined.hpp:199-200`, enforced at `.cpp:141-151` vs `:173-183`)
  — and the class overview at `:51` already says `rho(B,angle,T)`. This reorder would make
  the whole file self-consistent, which is a genuine argument *for* it that the plan had
  not noticed.

**D — The gas hierarchy is a live trap for any mechanical sweep.** See the scope guard.
`cl_Gas.cpp:1295,2674,2844,3518,3524`, `cl_GT_RefGas.cpp:1260`, and both module `main.cpp`
files call `rho`/`lambda` with a different, frozen meaning.

**E — The doc work inverts.** The 2026-08-25 three-AI overhaul rewrote all seven
`src/physics/materials/doc/` files *to* `(T, B, beta)` — three days before this decision.
Those seven become wrong on the day this lands; the two fem guides and the header comment
become right. Net documentation effort is a revert of very recent work, which is worth
stating plainly so it is not mistaken for rot.

### 3.2 Audit round 1 — Codex, 2026-08-28 (plan stage, no source modified)

Thread: `tmp/ai_exchange/rho_lambda_b_first_reorder.md`. Verdict: **"not implementation-ready
— the central inventory is incomplete and the claimed compiler enforcement is materially
overstated."** Every finding below was independently re-verified by Claude against the cited
source before being recorded here.

- [ ] **D1 (HIGH) — the inventory undercounts, in three ways.** CONFIRMED.
  (a) **Eight public positional call sites, not four:** the four in `cl_FEM_Calculator.hpp`
  plus `src/physics/materials/main.cpp:88, 94, 385, 402`. *Root cause of the miss, recorded
  so it is not repeated: Claude's call-site sweep excluded `src/physics/materials/` in order
  to filter out the class definitions, and thereby excluded that module's own executable.*
  (b) **Four Kohler member pointers, not one** — `mFunctionRhoKohler`,
  `mFunctiondRhoKohlerdT`, `mFunctiondRhoKohlerdB`, `mFunctiondRhoKohlerdbeta`, all declared
  `( const real T, const real normB, const real angle )` at `cl_Material.hpp:387-390`,
  with their virtual targets (`rho_kohler`, `drhodT_kohler`, `drhodB_kohler`,
  `drhodbeta_kohler` and the table equivalents) at `:1495-1517`, base definitions at
  `cl_Material.cpp:560-609`, `Metal` declarations at `cl_Material_Metal.hpp:435-453`, and
  pointer assignments at `cl_Material_Metal.cpp:47` and `:925`.
  (c) **`UserDefined` is a third overrider** of the field-dependent family
  (`cl_Material_UserDefined.hpp:310, 324`) — §3's "only `Metal` and `Alloy` override" is
  true of the public overloads only, and false of the implementation surface.
  Gap-table rows 4, 6, 8, 9, 10 and 16 are superseded by this finding.
- [ ] **D2 (HIGH) — the change is invisible to the type system, everywhere.** CONFIRMED;
  this is a C++ language fact, not a codebase property. Consequences the plan got wrong:
  `override` gives no protection (D2 struck the §2 claim in place); the member-pointer
  "reorder" in R5 is not a type-level operation at all, only a rename; and after R4 restores
  the overloads, **stale calls compile again** — R2's window is the only moment the compiler
  helps. Internal quotient-rule calls at `cl_Material_Metal.hpp:664, 672` and
  `cl_Material_Alloy.hpp:209, 217` accept either permutation before and after.
- [ ] **D3 (HIGH) — R6's numeric gate cannot prove what §4 claims it proves.** CONFIRMED by
  construction: a compensating split — public API B-first, body still forwarding
  `( T, B, beta )` at the pointer boundary, private `rho_kohler` still T-first — reproduces
  the old numbers **exactly** while leaving two conventions alive internally. The sentence
  "this is the only artifact that can prove finding 3.1 A" is false and is withdrawn; R1/R6
  is a behavioural regression gate, and structural consistency needs a separate review pass.
  Codex adds, correctly, that `drhodT`/`dlambdadT` need their own non-trivial gate values —
  a correct value path does not exercise the derivative-pointer path.
- [ ] **D4 (MEDIUM) — the plugin claim in §3.1 C is misleading in the dangerous direction.**
  Partly confirmed. Codex: because the parameter types are unchanged, this is *not* a
  link-time ABI break; a stale plugin loads and silently feeds T-first values into B-first
  semantics. **Claude's counter-finding, verified and recorded against Codex:** the plugin
  surface is *not* silently unprotected. `UserDefinedMaterial::set_user_defined_function`
  hard-errors on the declared dependency order for exactly `rho` and `lambda`
  (`cl_Material_UserDefined.cpp:141-151`: `Dependency1 == T`, `Dependency2 == normB`,
  `Dependency3 == angleBxJ`, each its own named `BELFEM_ERROR`). If the reorder updates
  those three checks, every stale plugin fails **loudly at registration** with a named
  message — the best available failure mode, and a genuine mitigation neither the plan nor
  the audit had. It also means the reorder *must* update them, or the enforced plugin
  contract and the internal convention diverge silently. This is a user-visible contract
  change and settles O4's classification: `[F]`, freeze-relevant.

**On Appendix B, Codex dissents and Claude accepts the dissent in part.** Codex agrees the
type-system statement is correct but rejects the empirical premise: the codebase does *not*
evolve 1-arg calls into 3-arg calls. `Calculator` keeps `compute_lambda_bulk()` (`:2883`)
and `compute_lambda_metal()` (`:2895`) as separate functions, rho likewise at `:2935` and
`:2969`; `main.cpp:85-95` holds both arities side by side rather than growing one into the
other; the tests are temperature-only; and git blame shows the 3-arg calls were *reordered*
in July, never *extended* from 1-arg ones. Most of the "17 exposed sites" are normal-state
evaluations where appending a field would be physically wrong. **Verdict accepted: the
append-motion argument is a possible future pattern, not an observed one, and it is not
sufficient grounds to refuse the reorder.** Appendix B is downgraded accordingly — see the
note appended there. Note that D2 is an *independent and stronger* hazard finding than
Appendix B's: it says no permutation is protected by anything, which argues for O3
(strong types), not for either order.

### 3.3 Audit round 1 — Grok, 2026-08-28 (blind, parallel with Codex)

Verdict: *"the two-pass idea is sound for the names it deletes; the inventory is not the
whole `(T,B,beta)` family; Appendix B's '17 sites re-armed' claim fails the empirical test
Claude asked for … but the plan is not implementation-ready."* Independently confirms D1,
D2 and the Appendix B dissent. Read-only compliance checked after the round: `git status`
clean of any Grok-authored change, tool set was `read_file,grep,list_dir`. New findings:

- [ ] **D5 (HIGH) — the family has four more members, and two more live call sites.**
  CONFIRMED by Claude. `drhodB` and `drhodbeta` are `(T,B,beta)` public siblings with **live
  Calculator callers** at `cl_FEM_Calculator.hpp:3480` and `:3499` — verified by direct read.
  So the Calculator carries **six** positional 3-arg calls, not four; with
  `materials/main.cpp` the public total is **ten**, not the plan's four. `dlambdadB` /
  `dlambdadbeta` complete the family (`cl_Material.hpp:622,625`; Metal `:210,:213`;
  not overridden by Alloy). D1(a) is superseded by this count.
- [ ] **D6 (HIGH) — four pointer *forwards*, and R2 does not reach them.** The four Kohler
  pointers forward at `cl_Material_Metal.hpp:578, 596, 602, 608`, with pointees
  `rho_kohler`/`drhodT_kohler`/`drhodB_kohler`/`drhodbeta_kohler`
  (`cl_Material.hpp:1496-1505`, `cl_Material.cpp:560-584`) and a parallel `*_table` family
  (`:1508-1517`). Because R2 deletes only `rho`/`drhodT`/`lambda`/`dlambdadT`,
  `Metal::drhodB` survives Pass 1 and keeps calling 3-arg `rho` — so **Pass 1's error list
  is not the inventory the plan claims it is.**
- [ ] **D7 (HIGH) — my INC citation was wrong, and the right one is worse.** The plan cited
  "INC-286's sibling INC-114". INC-286 is the H-F moved-baseline detector — unrelated. The
  real precedents are **INC-114** (lambda order), **INC-115** (Christian's own 2026-07-11
  T-first unification) and **INC-116** (the UserDefined 3-arg callback caveat), all verified
  present. **INC-115 is the load-bearing one: the last time this exact family moved, the
  Kohler/table pointer path is what silently swapped.** That is D6's failure mode, already
  realised once. Recorded as a Claude error of the register's own standing trap #2 —
  "follow the ID you cite" — committed in the plan that cites that trap.
- [ ] **D8 (MEDIUM) — a new trap: an internal grid that must NOT be permuted.**
  `Database::evaluate( theta, log10B, angle )` inside `Alloy::rho` and `Metal::rho_table`
  takes three `real`s in a **lookup-table** order that has nothing to do with the public API.
  A sweep that permutes "every 3-real call in the materials module" corrupts the table
  lookup silently. Add to the scope guard beside the gas hierarchy.
- [ ] **D9 (MEDIUM) — the plugin hazard is smaller than either D4 or Codex's finding F
  says.** Grok, verified in part: the 3-arg user callback path is **dead** — 3-arg
  registration only `set_custom`s, which wires the **1-arg** `mFunctionRho`/`mFunctionLambda`
  (`cl_Material.cpp:384-394`), while public 3-arg `rho`/`lambda` assert `PureMetal`
  (`cl_Material.hpp:1887, 2020`) and `UserDefined` is `MaterialType::UserDefined`. That is
  INC-116. So Codex's "stale plugin silently computes garbage" is unreachable, and Claude's
  counter-finding (the loud `BELFEM_ERROR` dependency-order gate at
  `cl_Material_UserDefined.cpp:141-151`) is real but guards a dead path. **Both the risk and
  the mitigation shrink.** What survives: out-of-tree *callers* of public 3-arg `rho` need a
  rebuild and a source rewrite, and the in-source doc split at
  `cl_Material_UserDefined.hpp:199-200` (rho/lambda T-first, jc/n B-first) is not in the gap
  table.

**Grok on Appendix B — same verdict as Codex, different emphasis.** The type-system half is
"true"; the empirical half fails, with a fuller sample than Codex's (18 one-arg sites, not
17, plus ~40 `this->rho(T)` in `powerlaws.hpp` the appendix never counted, all of which are
ρ_n(T) normal-state evaluations where field dependence lives on a *different name*,
`rho_powerlaw`). Its conclusion is worth quoting because it is the sharpest statement of the
live risk either auditor made: **"the live footgun today is fem prose teaching
`mat->rho(norm_b, beta, gTbulk)` against a T-first API — a copy-paste silent swap, INC-114's
shape, without any code change. Eight guide lines close that. B-first closes those guides
and immediately opens the seven materials files just rewritten T-first on 2026-08-25."**
Grok explicitly declines to treat this as a reason to reverse Christian's decision.

**Where the two auditors disagree:** only on the plugin surface — Codex F (silent stale-`.so`
hazard) versus Grok D9 (path is dead). Grok's is the better-evidenced position and is
corroborated by INC-116; Claude's verification of the `BELFEM_ERROR` gate is consistent with
both. Resolution: the plugin surface is **not** a live hazard for this change, and §3.1 C
should be rewritten rather than merely softened.

## 4. Ordered Steps

> **REWORKED 2026-08-28 after audit round 1.** The previous R1–R9 assumed a four-function
> family, four call sites and one pointer; they are struck through where superseded. The
> ordering below fixes the two structural errors the auditors found: R2 must delete the
> whole family at once, and the private layer needs its own enumerated pass because no
> mechanism in this plan reaches it.

- [ ] **R0 — Freeze-sequencing ruling** *(blocks everything; see O4)*. This changes an
      enforced, user-visible contract (`cl_Material_UserDefined.cpp:141-151`), so it is
      `[F]`-class. Land before the freeze or accept an explicit post-freeze API break.
      **No other step starts until this is answered.**
- [ ] **R1 — Freeze the numeric baseline** *(after: R0)*. Kohler metal with genuine field
      dependence, at a state where B, angle and T are mutually distinguishable (B = 5 T,
      angle = 0.7 rad, T = 20 K — no two values close enough for a swap to hide). Record
      **all eight** functions, not two: `rho`, `drhodT`, `drhodB`, `drhodbeta`, `lambda`,
      `dlambdadT`, `dlambdadB`, `dlambdadbeta`. ~~This is the only artifact that can prove
      finding 3.1 A.~~ **False (D3)** — a public-B-first / private-T-first split reproduces
      these numbers exactly. R1/R6 is a behavioural regression gate and nothing more.
      Derivatives need their own non-trivial values: a correct value path does not exercise
      the derivative-pointer path.
- [ ] **R2 — Pass 1: delete ALL EIGHT public 3-arg overloads** *(after: R1)*. Base,
      `Metal`, `Alloy`, declarations and definitions — gap rows 1–7 in full. Keep every
      1-arg form. Compile and **record the complete error list**. Deleting only four is
      worse than deleting none: it yields a confident, wrong inventory (see §2).
- [ ] **R3 — Reconcile R2's error list** *(after: R2)*. Every compiler-named site not in the
      gap table is an inventory miss — log a `Dn` before fixing it. Every gap-table row the
      compiler does *not* name is either row 11's facade or a wrong row. **Expect the
      private layer (rows 8, 9, 9b, 9c, 9d, 10b) to be entirely absent from this list —
      that absence is the point of R5', not evidence of safety.**
- [ ] **R4 — Pass 2: reintroduce B-first** *(after: R3)*. All eight, B-first, on base and
      both overriders; fix every R2 error. Rename parameters to match the new order — no
      parameter named `T` may remain in the first slot. Note `override` validates nothing
      here (D2); it is retained for style, not safety. Watch the `real T` vs `const real T`
      spelling split (3.1 B) — it defeats find-replace on four definitions.
- [ ] **R5' — The private layer, enumerated and reviewed by hand** *(after: R4)*. **This is
      the step INC-115 failed.** Work the explicit list, ticking each: 4 pointer
      declarations (`cl_Material.hpp:387-390`) · 4 forwards (`cl_Material_Metal.hpp:578,
      596, 602, 608`) · 8 pointees `*_kohler` + `*_table` (`cl_Material.hpp:1496-1517`,
      `cl_Material.cpp:560-584`, `cl_Material_Metal.hpp:435-456`, `.cpp:727+`) · 8
      assignments (`cl_Material_Metal.cpp:47-50, 925-928`) · `UserDefined::rho_kohler` and
      `lambda_custom` (`cl_Material_UserDefined.hpp:310, 324`) · internal self-calls
      (`cl_Material_Metal.hpp:664,672,676,688,690,698,700`;
      `cl_Material_Alloy.hpp:209,217,221`). **Do not permute
      `Database::evaluate( theta, log10B, angle )`** (row 19). Second reviewer required;
      the compiler will accept every permutation of every line in this step.
- [ ] **R6 — Numeric gate** *(after: R5')*. Re-evaluate R1's eight baselines.
      Bit-identical, or a site is inconsistent. A NaN or a `PureMetal` assert means a
      material-type guard was missed. **Necessary, not sufficient** (D3).
- [ ] **R6b — Structural consistency review** *(after: R5')*. The gate R6 cannot be: walk
      the public → forward → pointee → table chain for one property end to end and confirm a
      single convention throughout. This replaces the false confidence the old R6 carried.
- [ ] **R7 — Suite gate** *(after: R6)*. `make check` and `make check-fast` green, including
      `test_YBCO` and the backend-free material gate.
- [ ] **R8 — Contract and documentation** *(after: R4)*. Update the three `BELFEM_ERROR`
      dependency-order checks (`cl_Material_UserDefined.cpp:141-151`) — **mandatory, not
      cosmetic**: leaving them makes the enforced plugin contract contradict the API. Then
      invert the seven `src/physics/materials/doc/` files back to B-first (row 12), verify
      rows 13–15 are now correct rather than editing them, and resolve the in-source split
      at `cl_Material_UserDefined.hpp:199-200` so rho/lambda and jc/n finally agree. Codex
      language sweep over the rewritten prose.
- [ ] **R9 — Record and retire** *(after: R8)*. Devlog; rewrite the Status line; mark
      `rho_lambda_argument_convention.md` HISTORICAL (its option-(c) decision is reversed);
      add an INC card for this round; tell plugin authors that out-of-tree **callers** of
      the 3-arg forms need a source rewrite and that **the rebuild will not warn them**,
      because the symbol is unchanged (D9).

## 5. Open Design Questions

- [ ] **O1 — The 4-arg HTS `lambda( T, B_par, B_perp, J )` variant** (`cl_Material.hpp:636,
      639`) has **no caller anywhere in `src/`, `tests/` or `nonfree/`** (verified
      2026-08-28). Options: (a) reorder to `( B_par, B_perp, J, T )` for consistency;
      (b) leave T-first and accept a third convention in the same header; (c) delete it as
      dead API and reintroduce it field-first when something needs it. **Claude's
      recommendation: (c)** — it is unreachable, and a dead signature in the wrong order is
      a trap for whoever first calls it. Christian's call.
      **Codex dissents (2026-08-28), and Claude now agrees with Codex:** do *not* delete it
      merely because no in-tree caller exists — it is a public `Material` method and may
      have out-of-tree callers, so deletion is a separate API-removal decision and must not
      be smuggled into a reorder. Revised recommendation: **reorder it field-first with the
      rest, or exclude it explicitly with a written rationale.** Deletion, if wanted, is its
      own row.
- [ ] **O2 — Does the powerlaw / piecewise derivative family follow?** Today it is
      `drho_powerlaw_dB( normJ, T, normB, angleNxB )` (`cl_Material.hpp:909,913`) and
      `drho_piecewise_dB` likewise (`:953,957`) — a **third** convention, with the current
      density first. Leaving it makes the tree teach two rules again, which is the problem
      this plan exists to remove; changing it enlarges an already hot diff. Not decided.
      **Both auditors: exclude it from this round.** Codex's reason is the stronger one —
      `normJ` is the *variable of the resistivity law*, not a property argument, and the
      family is a large structured overload set (`cl_Material.hpp:860-957`) whose `rho_*`,
      `d/dJ` and `d/dT` members would have to move together; reordering only `drho_*_dB`
      would break internal consistency. Grok adds the physical framing: `powerlaws.hpp`'s
      ~40 `this->rho(T)` calls are ρ_n(T), the normal-state resistivity, and field
      dependence there already lives on a **different name** (`rho_powerlaw`) — the naming
      remedy this plan's Appendix B proposes is *already in production on the HTS path*.
      **Recommendation: audit that family separately, as one unit, in its own round.**
- [ ] **O3 — Strong parameter types instead of an order convention.** Wrapping the field
      and the temperature in distinct one-member structs would make the swap a compile
      error permanently, and would have prevented INC-114 and this plan. Cost: a larger
      plugin-API change and a zero-overhead argument to make. Logged for a future session;
      **not** proposed for this round.
- [ ] **O4 — Sequencing against the design freeze. NOW THE BLOCKING QUESTION (R0).** The
      freeze is ~2026-08-29 and this changes a **contract that is enforced at runtime**, not
      merely documented: `cl_Material_UserDefined.cpp:141-151` hard-errors on the declared
      dependency order for 3-arg rho/lambda. That settles the classification — `[F]`,
      freeze-relevant — and Codex adds the decisive operational point: because the old and
      new APIs are **binary-indistinguishable**, telling plugin authors to rebuild does not
      protect them. A rebuild of unchanged source succeeds and stays wrong. So the choice is
      a real one: land before the freeze, or accept an explicit post-freeze API break with a
      written note to plugin authors. Christian's call, and **nothing starts until it is
      made** (R0).

## 6. The Interface After the Change

```cpp
// base — src/physics/materials/cl_Material.hpp
virtual real rho      ( const real B, const real angle, const real T ) const ;
virtual real drhodT   ( const real B, const real angle, const real T ) const ;
virtual real lambda   ( const real B, const real beta,  const real T ) const ;
virtual real dlambdadT( const real B, const real beta,  const real T ) const ;

// unchanged, no ordering to get wrong
real rho   ( const real T ) const ;
real lambda( const real T ) const ;
```

Rule, in one line for the module doc: **three-argument property accessors take the field
first and the temperature last; one-argument accessors take temperature.**

## 7. Definition-of-Done Checklist

- [ ] Every gap-table row maps to a step or an open question.
- [ ] Each claimed site backed by a `file:line`, not an assumption.
- [ ] R2's compiler error list reconciled against the gap table, misses logged as `Dn`.
- [ ] R6 numeric gate bit-identical on both `rho` and `lambda`.
- [ ] `make check` / `make check-fast` green.
- [ ] Seven materials docs inverted; fem guides verified correct; Codex sweep applied.
- [ ] `rho_lambda_argument_convention.md` marked HISTORICAL with a pointer here.
- [ ] Plugin ABI break recorded where plugin authors will see it.

## 8. Audit Trail

- Exchange thread: `tmp/ai_exchange/rho_lambda_b_first_reorder.md`
- Inventory performed by Claude 2026-08-28 against the working tree at `efe46051`+; every
  citation in §1 and §3 read directly, not carried from `rho_lambda_argument_convention.md`
  (whose own line references had drifted by ~80 lines).
- Plan audit dispatched to Codex and Grok 2026-08-28, blind and parallel, per the standing
  three-vendor rule. Findings to be distilled here before the thread is swept.

## Appendix A — Decision: why uniform field-first rather than per-quantity physical ordering

The rejected alternative had a real argument behind it: for resistivity the field is the
driving variable, for thermal conductivity the temperature is, so ordering each by its own
dominant dependence documents the physics in the signature — `rho( B, angle, T )` and
`lambda( T, B, beta )`.

It was rejected because the property that matters at a call site is not expressiveness but
**recallability under a compiler that cannot help.** An asymmetric rule must be recalled
correctly per function; a uniform rule is recalled once. With every parameter a `real`, the
cost of misremembering is silent — INC-114 is the recorded instance, and it survived
months and an audit round. The physics intuition survives in the documentation, where it is
free; the signature buys uniformity instead.

| | asymmetric (physical) | uniform field-first (chosen) |
|---|---|---|
| expresses dominant dependence | yes | no — doc does |
| rules to remember | one per quantity | one |
| silent-swap risk | per-function recall | single recall |
| decided | rejected 2026-08-27 | Christian, 2026-08-27 |

## Appendix B — Claude's opinion: which permutation carries the least hazard

Requested by Christian 2026-08-28. **Recorded as an opinion for the auditors to attack, not
as a decision** — the decision in the header stands until Christian changes it. Confidence:
medium-high on the mechanism, high on the two facts it rests on.

**Recommendation: keep `( T, B, beta )` and fix the eight documentation lines.** Ranking,
least hazardous first:

| rank | permutation | migration hazard | standing hazard |
|---|---|---|---|
| 1 | **`(T, B, beta)` — keep** | none | prefix-consistent; docs already fixed in 7 of 10 places |
| 2 | `(B, angle, T)` / `(B, beta, T)` | one-time, mitigable by R2 | **breaks prefix consistency** (below) |
| 3 | asymmetric per quantity | one-time | per-function recall — already rejected |

Two facts drive this, both verified against the tree on 2026-08-28.

**Fact 1 — the code is already uniform.** The decision was taken believing the tree held
`rho` and `lambda` in different orders. It does not: `rho( T, B, beta )` (`:670`) and
`lambda( T, B, beta )` (`:616`) are the same shape, and have been since the 2026-07-11
unification. The inconsistency that prompted this was **between the code and two fem
documents**, not inside the code. So B-first does not buy uniformity — the tree already has
it — and the cheapest route to one rule everywhere is eight doc lines and no code.

**Fact 2 — T-first makes the overload family prefix-consistent, and B-first does not.**
The full set is 1-arg, 3-arg and 4-arg:

```
rho( T )                              lambda( T = gTroom )
rho( T, B, beta )                     lambda( T, B, beta )
                                      lambda( T, B_par, B_perp, J )
```

Position 1 is temperature at **every arity**. The natural editing motion — "this material
is field-dependent now, add the field arguments" — is a pure append, and appending is
correct. There are **17 one-argument call sites** in the tree that are candidates for
exactly that edit.

Under B-first the same motion silently breaks. `rho( tT )` edited to `rho( tT, tNormB,
tBeta )` still compiles — it matches the 3-arg overload — but now passes temperature as
the field and the field as the angle. The compiler cannot object, because overload
resolution *succeeds*: three `real`s match a three-`real` signature. This is not a
migration risk that R2 retires; it is a permanent property of the resulting API, and it is
the same failure shape as INC-114, re-armed at 17 sites.

**Why the uniformity argument does not outweigh it.** "One rule, not one per quantity" is a
real benefit, but it is a benefit T-first already delivers. B-first would trade a hazard
that exists only in prose (two guides teach the wrong order) for a hazard that exists in
the type system (a common, natural edit compiles into a silent swap). Prose can be fixed
in five minutes and verified by grep; an overload set that punishes the obvious edit cannot
be fixed except by changing it back.

**What would change my mind:** if the append motion is not in fact how these call sites
evolve — if field-dependence is always introduced by writing a fresh call rather than
extending an existing one — then Fact 2 loses most of its force and the choice comes down
to taste, where Christian's is decisive. The auditors should test that specifically.

> **DOWNGRADED 2026-08-28 — the pre-registered condition above was met.** Codex tested it
> and refuted the empirical premise with evidence (§3.2): the codebase keeps 1-arg and
> 3-arg paths as separate functions rather than growing one into the other, and git blame
> shows the 3-arg calls were reordered in July, never extended from 1-arg ones. Fact 2 is
> therefore a possible future pattern, not an observed one. **Claude withdraws the "against
> reorder" recommendation.** Fact 1 (the tree is already uniform, so the split is in two
> documents rather than in the code) still stands and still means the *cheapest* route to
> one rule is eight doc lines — but "cheapest" is not "best", and cost was never Christian's
> criterion. The mitigation proposed at the end of this appendix — a distinct name for the
> field-dependent forms — survives the dissent intact and is now the recommendation:
> it is the only option that makes the convention unmissable rather than merely uniform,
> and D2 shows nothing else in the toolchain will catch a mistake.

**If B-first is chosen anyway**, two mitigations are worth folding into R4, and neither is
expensive:

1. **Delete the 1-arg/3-arg overload ambiguity** by giving the field-dependent forms a
   distinct name (`rho_field( B, angle, T )`). Then the append motion fails to compile
   instead of failing silently, and the uniformity goal is fully met. This is the option
   I would actually argue for if the order must change.
2. Resolve **O3** (strong parameter types) in the same round rather than deferring it —
   it is the only remedy that closes the class permanently, and the diff is already open.
