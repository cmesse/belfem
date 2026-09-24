# Materials cluster: Kim field tangent, Callaway overflow, the n = 4 knot, and two rulings measured

**Date:** 2026-09-23 (overnight into 2026-09-24)
**Purpose:** Close DR-164, DR-174 and DR-170 through the standing loop and put DR-173 and
DR-163 in front of Christian with executed numbers instead of hypotheses. Two register
diagnoses were wrong and are corrected here.
**Module:** `src/physics/materials` (`powerlaws.hpp`, `cl_JcFunction_ModifiedKim.hpp`,
`cl_JcFunction.hpp`, `debye.f90`), `tests/physics`
**AIs involved:** Claude Fable 5.1 (plan, probes, juries, close-out), Opus 5.5 (implementation),
Codex `gpt-5.6-terra`/high and Grok `grok-4.7`/high (plan jury, code jury, DR-163 physics round);
`belfem-3c` (the concurrent v0.9.2..HEAD reviewer, coordination only)
**Status:** COMPLETE — all five rows landed on `devel`: DR-174 `e3a9b6b`, DR-164 `b5c42a9`, DR-170 `9373e44` (struck on Christian's ruling), DR-173 `fcb8536`, DR-163 `31737f5` (rulings (c) and A from the 00:05 quiz; strikes owed), each behind both juries and green gates in both trees; the debug-assert leg is owed
**Exchange:** `tmp/ai_exchange/review_materials_tangent_cluster.md` (plan and code rounds),
`review_dr163_physics.md` (physics round)

## 1. Why this cluster, and how it was run

Christian asked which open register cluster to approach next; five rows share one module and
one gate shape (a `DISABLED_` pin each), so they went as one plan
(`todo/materials_tangent_cluster.md`). The session ran beside a read-only regression review of
the 89 commits since 0.9.2 (`belfem-3c`), with an explicit protocol: no edit under `src/` or
`tests/` until its literal "REVIEW DONE" (received 21:45), per-slice notes as information only,
concurrent read-only auditors allowed. Christian's overnight terms applied: plan + audit, wait
for the peer, then code; builds and one commit per row by this session; open rulings become
quiz items rather than edits.

## 2. What the probes found (executed, not inferred)

Six scratchpad probes, linked read-only against the prebuilt `build/lib/libbelfem.a`
(recipe in memory `project_scratchpad_test_probes`):

| row | register said | probe found |
|---|---|---|
| DR-173 lead cp < 0 below 3 K | "the first cp fit region is wrong below its switch" | `cp_custom` is positive everywhere (Debye cubic, then `exp(…)`); the **public 4 K spline** (`cl_Material_SplineLookupTable.cpp:131-148`) undershoots to −3.23e-2 J/(kg K) between 0.3 and 2.1 K. Tin is 48 % off `cp_custom` at 1.4 K and indium 16 % at 0.8 K — positive but wrong. Rebuilt at 1 K knots every metal is within 0.6 % of its own law |
| DR-174 YBCO λ NaN below 1 K | "the Fortran integrator overflows at θ/T ≈ 5500" (the Bose integrand) | the NaN is born in the **optical electron–phonon prefactor** (`debye.f90:326-329`): `U = ω_opt/(k_B T) = 720.8/T`, `exp(U)` overflows for T < 1.016 K, `U²V/(V−1)` = Inf/Inf. And `quadruple = selected_real_kind(r=60)` is **binary64** (kind 8, range 307) — the module has always run in double |
| DR-170 abort at the n = 4 knot | "replace the `BELFEM_ERROR` by the frozen-knot fallback `dT` uses" | the `dJ` limit at the coincident knots is exact and simple: `dρ/dJ → n_ff ρ_n / J` (the blend is `ρ_n (J/j3)^{n_ff}` there); analytic, central difference and limit agree to 1e-6 at `1 − J/j3 = 1e-7`, ratio `(1−ε)³` |
| DR-164 Kim `deval_dB` | closed form to be derived | `−α·jc·s/(B0 + s·B)` matches a central difference of `eval` to 1.8e-10 over 25 (θ, B) points |

The `e^{−z}` rewrite of `scatter_function` alone was tested first and left the NaN in place —
the register's mechanism was a masked second overflow (L-10), located only by instrumenting a
scratch copy.

## 3. Plan jury and code jury

Plan jury: Grok 21:08 (Codex failed on the ChatGPT usage limit; Christian topped up the credits
and the leg re-ran at 21:13 on the current plan). Findings folded in, all verified in source:
the Callaway site (Grok predicted the 1.016 K threshold by hand before the probe found it), the
Fortran kind, the eight `dJ` overloads that all divide by `√(b²+ac)` unguarded, the `dJ` limit
that had to be written down (probe 3 supplied it), R4's `dcpdT`/`d2cpdT2` and the pin's
log-step clause, indium's β > lead's, the wrong literature pair for DR-163 (Duron 2004, not
Messe §2.6), and — from Codex — that demoting the `a` guard would let a legal `n = 7` divide by
zero in release, so the `a` and `imag` errors stay. Reconciliation table in the exchange.

Code jury: see §6.

## 4. What landed

- **DR-164:** `JcFunctionModifiedKim::deval_dB`, the exact field derivative; the base-class
  contract comment no longer lists ModifiedKim. Kim materials now carry a field tangent in
  `drho_powerlaw_dB` and `drho_piecewise_dB`.
- **DR-174:** the optical prefactor as `U²/(1 − e^{−U})`; `scatter_function` and
  `debye_function` in the `e^{−z}` form. Scratch A/B: finite to 0.05 K, identical to the old
  code to 1e-15 from 1.5 K up. `debye_derivatives` untouched (reached at z ≈ 3–6 only).
- **DR-170:** all eight `drho_piecewise_dJ` overloads guard the discriminant before any root
  and return `n_ff ρ_n / J`; the 4-argument overload's `BELFEM_ERROR` on the discriminant is
  gone, its errors on `a` and on `imag` stay.
- Tests: `KimFieldDerivativeMatchesFiniteDifference` and `YbcoRawLambdaIsFiniteBelowOneKelvin`
  enabled; `PiecewiseDegenerateBlendStaysFinite` now calls every reachable `dJ` overload at the
  knot and checks the limit to 1e-6.

## 5. Rulings (quiz answered 2026-09-24 00:05 — DR-173: (c); DR-163: A; the three landed rows struck; the DR-170 guard band accepted as bounded)

- **DR-173:** (c) `dT = 1 K` for the cp spline only (a trailing step parameter on
  `create_spline`, callers unchanged) vs (b) 1 K for every property. Recommended (c). The lead
  pin stays `DISABLED_` until it lands.
- **DR-163:** both physics voices recommend A — the tangent follows the documented raw value
  — because the Bézier blend's first control leg has log-slope exactly `n − 1`: it is built to
  join the raw law C1 at `j1`, and B would leave a jump of `ρ1/ρn` and a kink there unless
  `ρ1`, `j2`, `j3` and the dB/dT legs are re-derived. Error size ≈ 4e-6 at `j1` for 77 K
  constants, O(1) where jc collapses. Recommended A.

## 6. Gates and commits

- **`build/` (Armadillo, MKL, release), 22:30:** `make -j16 test_physics belfem` clean (31 s
  wall). Targeted: `KimFieldDerivativeMatchesFiniteDifference`, `YbcoRawLambdaIsFiniteBelowOneKelvin`,
  `PiecewiseDegenerateBlendStaysFinite` 3/3 green — the last includes the exact-knot limit check
  at `J = j3`, so the fixture's `j3` matched the function's. Full `test_physics`: 148 ran, 147
  passed, 1 skipped by design (`magnesia has no rho`, pre-existing). The "red before" half of
  the two pins is the recorded pre-fix evidence (T-26 exchange, 2026-09-22: Kim `analytic 0`
  vs `fd 4e-15`; the value sweep's YBCO NaN at 0.05–1.0 K), not a run tonight on the unfixed
  tree; the DR-170 case had no pre-fix red — before R6 it aborted the process.
- **`cmake-build-debug/` (Blaze, release), 22:40:** `make -j16 test_physics belfem` clean
  (2 m 17 s wall); the same three targeted cases 3/3 green; full `test_physics` 147 passed + the
  same designed skip; `check-fast` 8/8 (3.8 s test time).
- **`make check-fast`, `build/` (Armadillo), 22:45:** 8/8 suites passed (4.4 s test time,
  2 m 18 s wall with the test builds).
- **Code jury (Codex terra/high 21:26, Grok 4.7/high 21:36)** — reconciled 22:55: the eight
  guards stand; Codex's P1 (the guard band admits `n` within ~1e-7 of 4, returning the coincident
  limit within ulps of `j3` where the exact left derivative is 0) is confirmed as algebra and
  bounded as a measure-zero, finite, right-magnitude tangent at a point where the tangent jumps
  for every `n` — documented in the first overload's comment and pinned by two new test cases
  (`n = 4 + 1e-6` on the root path, `4 + 1e-9` on the fallback); the "709" comment reworded
  kind-agnostic (Codex P3); the class overview names `deval_dB` (Grok). Residue recorded:
  `debye_derivatives` still on the `e^z` form (only caller at z ≈ 3–6), the one-sided `imag`
  guard, the `n = 7` window.
- **Re-gate after the post-jury edits, `build/`, 23:05:** targeted 3/3 (now including the
  near-degenerate pair), full `test_physics` 147 + 1 skip. **Blaze re-gate, 23:12:** `test_physics`
  rebuilt (1 m 07 s), targeted 3/3, full suite 147 + 1 skip.
- **Commits (23:15, one per row, private index, the shared test file split by hunk):**
  `e3a9b6b` DR-174, `b5c42a9` DR-164, `9373e44` DR-170. Each message carries the mechanism,
  the probe evidence, the jury depth and the gates. Verified at the suite rung in release on both
  backends; the debug-assert leg is owed (both trees are release-configured tonight).

## 6b. Stage 3 — the two ruled rows (2026-09-24, 00:05 ruling → implementation)

- **DR-173 (c):** `create_spline( aProperty, adYdX0, adYdX1, aStep = 4.0 )` on `Material`,
  forwarded to the function-pointer overload and its `SplineLookupTable` override; the 4 K
  literal is now the default step; only `Metal::create_cp` passes 1 K. All 26 other call sites
  keep 4 K (enumerated by the implementation agent, syntax-gated in both flag sets and with the
  asserts compiled in). `LeadCpIsPositiveBelowThreeKelvin` enabled (renamed from the "IsNegative"
  pin); the known-defect table of the value sweep is empty.
- **DR-163 (A):** the eight `drho_piecewise_dJ` power-law branches return the raw derivative;
  the dead `rhoPL` locals are gone; the family's header states the value/tangent contract.
  `PowerLawRegimeTangentMatchesRawResidual` enabled.
- **Gate, `build/` (Armadillo, release), 00:55:** rebuilt clean; targeted 10/10 (both pins, the
  lead sweeps, the degenerate and Kim cases); full `test_physics` 150 ran, 149 passed, the same
  designed skip; `check-fast` 8/8. **`cmake-build-debug/` (Blaze, release), 01:00:** rebuilt
  clean (1 m 17 s), targeted 10/10, full suite 150 ran / 149 passed / the same skip, `check-fast`
  8/8 (3.8 s). Release only in both trees; the debug-assert leg stays owed.
- **Code jury (stage 3; Codex terra/high, Grok 4.7/high), reconciled 01:20:** both found the
  unvalidated `aStep` (setup path → `BELFEM_ERROR` on finite, positive) and the six
  un-differenced sibling overloads (the raw-residual test now checks the 1-argument and both
  defect twins); Grok's comment corrections applied (the sampler comment, the value-test
  banner, the `adYdX1` rename) and the Doxygen names the zero floor. Residue recorded: the
  positional third `real` of `create_spline`, `Alloy::create_tables`' own 4 K cp grid (a
  one-component lead alloy would rebuild the old undershoot), the sweeps' 0.5 K start.
  Re-gate after these edits, 01:35: both trees rebuilt clean, targeted 10/10, `test_physics`
  149 + 1 designed skip, `check-fast` 8/8 each. **Commits 01:45:** `fcb8536` DR-173 (six files),
  `31737f5` DR-163 (two files); no shared file, private index.

## 7. Residue

- `intpoints.f90:16`: the kind named `quadruple` is double — misleading name, widen onto DR-174.
- `cl_Material.hpp:286` `mD = 2.5` is dead and the header prose "2.5 decades above Jc"
  (`powerlaws.hpp:383-384`) contradicts the formula `jc·10^{2.5/n}` (Grok, physics round).
- The `n = 7` window (`a = 0`) aborts in the 4-argument `dJ` and divides by zero elsewhere —
  pre-existing, out of DR-170's scope, recorded on its row.
- Both trees are release-configured; the debug-assert leg of the gates is owed.
- `Alloy::create_tables` samples component cp on its own 4 K grid; the third optional `real` of
  `create_spline` is positional; the value sweeps' strict positivity starts at 0.5 K (DR-173 row).
- The `( J, T )` overloads of the HTS laws abort on the factory YBCO (`jc_custom`) and remain
  untested by the derivative suite (DR-163 row).
