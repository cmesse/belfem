# Coupled h-φ/T thermal divergence: diagnosis + first fix batch

**Date:** 2026-07-09
**Purpose:** Root-cause session for the exploding thermal Picard residual in
`tape_hphiTrun` (3D tape stack, thin-shell h-φ + thermal, fully coupled);
fixes to the coupled iteration guard and the thermal matrix producers.
**Module:** fem/kernel (Controller), fem/thermal

## Symptom

Timestep 1 of `cmake-build-debug/tape_hphiTrun`: magnetic Picard converges
smoothly, thermal Picard residual pinned at the 9.0 display clamp with dB
values swinging 26-88 dB for 65+ iterations, no timestep cut.

## Diagnosis chain

1. **Display decoding.** `print_line` shows `min(eps, 9.0)` but
   `dB = 10*log10(eps)` of the unclamped value
   (`cl_FEM_Controller.cpp:1322-1341`) → thermal residual truly at 10³-10⁹.
2. **Magnetic 0.9ⁿ is healthy.** Residual = ‖A·x_relaxed − b‖/‖b‖
   (`cl_FEM_DofMgr_SolverData.cpp:2311`) leaves exactly the (1−ω) relaxation
   remainder for a linear system: 0.9ⁿ at ω = 0.1. The magnetic problem is
   essentially linear at this field level.
3. **Dead divergence guard.** `iterate_coupled` never advanced `mIteration2`,
   so `tThermalReset` (`:693`) compared `0 > mMinNumIterations2` forever —
   already tracked in `todo/iterate_refactor_plan.md`, now fixed (see below).
4. **T field diagnostic (Christian's #DIAG probe, min/max T per iteration):**
   iteration 1 already dips to **−1.9 K** (max pinned at exactly 77 K),
   iteration 3 explodes to ±3000 K. With uniform 77 K start, pure-Neumann BCs
   and f ≥ 0, T ≡ 77 solves the discrete system exactly (K·𝟙 = 0
   structurally) — so the first solve is defective, and max(T) never rising
   above 77 means something acts as a heat *sink*. Once nodes sit at T < 0,
   material lookups detonate (e.g. Wiedemann-Franz λ = L·T/ρ goes negative →
   indefinite K), explaining the ±3000 K blowup. Confirmed NOT the YBCO
   defect (tested by Christian).

## Changes

### `cl_FEM_Controller.cpp` — coupled thermal guard (with Christian)

- `iterate_coupled` now advances `mIteration2` and passes it to
  `residual()`, inside the `mKernel2 != nullptr` block (Christian's draft
  incremented unconditionally and dereferenced `mKernel2->mesh()` for the
  #DIAG print without a null check → would segfault magnetic-only runs).
  Counter lifecycle was already closed (`initialize_timestep` /
  `reset_timestep` zero it). Verified live: the run now cuts Δt at
  iteration 4 instead of grinding to 400.
- Temporary `#DIAG` min/max-T print retained for the ongoing hunt.
- Both bullets ticked in `todo/iterate_refactor_plan.md`.

### `mt_thermal_h.cpp` — audit vs the validated magnetic reference

Audited all producers against `mt_maxwell_h.cpp/.hpp`. Correct and left
alone: `j = C·q`; the `−µ0·E·q` sign (matches the framework-wide edge-dof
convention); `T = dot(Nvec, q)`; consumer sign (`aRHS += f·Δt`,
`cl_IWG_Timestep.cpp:635ff`). Three defect classes fixed:

1. **bn missing the normal projection (6 sites, all `T_h_ts_*`).** The
   producers averaged the master/slave φ-gradients but skipped the
   `bn = dot(bn,n)·n` projection that the validated `compute_bn`
   (`mt_maxwell_h.hpp:125`, "purely-normal") applies → `b = bt + bn`
   double-counted the tangential field → wrong |B| and β into Jc(B,β,T),
   `rho_powerlaw`, magnetoresistive ρ(B,β,T) and λ(T,B,β). Projection added
   at all six sites (T_h_ts_metal additionally needed the `n` declaration).
2. **`T_h_ts` used the thermal calculator for the normal field.**
   `aCalc->get_normal_calculator(phi_m, phi_s)` gathered master/slave dofs
   from the *thermal* kernel (temperatures, not φ); its five siblings use
   the Maxwell calculator. Switched to `tCalculator->...`. In this model the
   producer serves the magnesia buffer layer.
3. **Unguarded/unclamped b–j angles.** `T_h:75` / `T_h_ts:232` had bare
   `acos(dot/(|b||j|))` (0/0 → NaN at vanishing field; roundoff |cos| > 1 →
   NaN with b ∥ j, the common case); `T_h_ts_metal` / `T_h_alloy` sites were
   eps-guarded but unclamped/unfolded. All brought up to the magnetic
   reference form (zero guard + `min(|cos|,1)`, folded into [0,π/2]).
   The `acos(clamp(dot(n,b)/|b|))` HTS-beta sites were already correct.

## Open items (the −1.9 K dip is NOT yet proven fixed)

None of the fixed defects flips the sign of f *by construction*; they feed
garbage into the ρ/Jc table lookups whose out-of-domain extrapolation is the
plausible path to ρ < 0. If the dip survives this batch, remaining suspects
in order:

- **`dV(k) < 0`** on the extruded thin-shell wafer elements (cf. the CW
  orientation landmines) — sign of `w·N·dV` is the only other negative-f
  channel (mesh is linear, N ≥ 0).
- ρ/Jc lookup extrapolation behavior at out-of-range (B, β, T).
- `T_phi` still carries hardcoded placeholder properties
  (ρcp = 1231.125, λ = 0.0855, `mt_thermal_phi.cpp`).
- `LookupAlloy` in the ThinShell branch of
  `IWG_MaxwellThermal::link_to_group` maps to the volume producer
  `T_h_alloy` (hastelloy, 50 µm layer) — plausibly intentional reuse (no
  bn/anisotropy needed), worth a deliberate confirmation.
- Robustness guard worth adding once root-caused: floor T going into
  material lookups (2-4 K) so a transient undershoot cannot detonate the
  property tables.

Build/regression handed to Christian (edits are syntax-consistent with the
file's existing idioms; clangd noise in this tree is the known stale
compile-database issue).

## Addendum (same day): three-AI audit round after fixes did not cure the dip

Christian: clamped T (gTmin) and ρ (gRhoMin/gRhoMax) in the producers,
added dV/dS assert guards (no negative measures fire), tested ghosts off
(dip persists), verified computed material properties physical. Dip
survives → dispatched Codex + Grok
(`tmp/ai_exchange/thermal_cold_dip_prompt.md` /
`..._codex.md` / `..._grok.md`).

**Established this round (all three agree, citations spot-checked):**

- **No thermal-only mesh exists.** Both kernels share one mesh object:
  `KernelParameters(Kernel*)` reuses `aKernel->mesh()`
  (`cl_FEM_KernelParameters.cpp:46-52`), the second-kernel
  `distribute_mesh` path never extracts a thermal submesh
  (`cl_FEM_Kernel.cpp:624-706`; non-root `new Mesh(dim)` is an empty
  placeholder, `mesh()` returns the shared submesh).
- **q_old happy path verified correct** (shift order
  `cl_FEM_Controller.cpp:187-197`, "T0" creation + copy
  `cl_IWG_Timestep.cpp:282-314`, label resolution
  `cl_FEM_Calculator.cpp:1951-1975`, dof field_index is mesh-global
  `cl_FEM_DofMgr_FieldData.cpp:66-93`). Latent bug kept on file: `qold`'s
  length-0 branch does `set_size` WITHOUT init (uninitialized q_old) —
  real, but not reachable in this run.
- **PENTA6TS thermal path is a consistent 3D volume triple** (PENTA6
  Lagrange shapes, 3D Gauss, det(J) fallback): L4 dead.
- **No mm/m unit mismatch** (scale_mesh to SI before extrusion).
- RHS extras / writeback aliasing / scatter permutation: dead.
- **Premise correction (Codex):** the "T≡77 is exact" argument required
  f = 0, but the user-defined current is ≈ −2.4 A from t = 0
  (`lib/I_vs_t_regular_smooth.txt`), so a nonzero Joule source exists at
  the first solve.

**Surviving candidates (need runtime data, static analysis exhausted):**

1. Non-monotone consistent-FEM operator + active Joule source on
   extreme-aspect µm prisms (Codex #1, medium-high) — though the observed
   asymmetry (−79 K undershoot vs < 5e-5 K overshoot) is not a classic
   maximum-principle violation profile.
2. Conditioning / silent low-quality MUMPS solve (Grok #1) — SolverData
   checks only `isnan(mRhsNorm)` post-solve, never solution quality
   (`cl_FEM_DofMgr_SolverData.cpp:2093-2121`). Claude's estimate
   κ ~ 1e3-1e4 says this should be benign; unverified.

**Agreed next step — one-run discriminating probe** at the first thermal
solve: print (a) min/max of "T0" restricted to thermal dofs, (b) min/max
of the assembled f-contribution, (c) ‖A·77𝟙 − b‖/‖b‖, (d) min/max of the
raw MUMPS x before Picard relaxation; optionally re-solve the same A with
f = 0 (Codex's split: f-driven undershoot vs solve/RHS defect).

## RESOLVED: MUMPS symmetric-pattern mismatch

Probe added ( `SolverData::probe_cold_dip`, one-shot, thermal-only,
master-rank; `#PROBE` lines; removed once closed ). First-solve output:

```
#PROBE T0      : min 77 max 77
#PROBE A*Tb-b  : min -2.7e-07 max 1e-17 norm 2.2e-06 ( ||b|| = 3.2e-03 )
#PROBE raw x   : min -0.955 max 21.6
#PROBE |y-Tb|  : max 77.96 ( pure solver error )
```

Reads unambiguously: q_old clean (T0 ≡ 77), assembly consistent (A·77𝟙 − b
= −f·Δt, tiny), yet the raw solve returns [−0.96, 21.6] and the control
solve A·y = A·77𝟙 (exact answer 77) is off by 100%. **The defect is
strictly between the matrix and the solver, in the SYM flag.**

Root cause: `IWG_Timestep` defaulted to
`SymmetryMode::PositiveDefiniteSymmetric` (`cl_IWG_Timestep.hpp:85`),
inherited unchanged by `IWG_TransientHeatConduction` /
`IWG_MaxwellThermal`, pushed into the solver at
`cl_FEM_DofMgr_SolverData.cpp:2036`. The MUMPS wrapper hands MUMPS the
**full** CSR pattern (`create_coo_indices`, no triangle filter anywhere in
`cl_SolverMUMPS.cpp`) while setting `SYM≠0` (`:274-295`). MUMPS spec: with
SYM≠0 each off-diagonal must appear once; **duplicates are summed** → every
off-diagonal doubled → `M + ΔtK` structure destroyed (K·𝟙 = 0 broken) →
finite, silent, wrong solve; SPD mode on the now-indefinite matrix warns at
most, and the wrapper aborts only on hard errors (`:326ff`). Magnetic solve
was immune because `IWG_Maxwell` explicitly passes `Unsymmetric`
(`cl_IWG_Maxwell.cpp:40`) → SYM=0 → full pattern correct — which is exactly
why "only thermal exploded."

Fix (Christian): base-default flip
`IWG_Timestep( ... SymmetryMode::Unsymmetric ... )`
(`cl_IWG_Timestep.hpp:85`) — MUMPS now runs plain LU on the full pattern.
Covers every timestep IWG subclass at once (~2× factorization cost on a
small SPD system; correct). Latent sibling also fixed (Christian):
`IWG_MaxwellPostproc` declared `PositiveDefiniteSymmetric`
(`cl_IWG_MaxwellPostproc.cpp:30`) — if its L2 projection ever routes through
MUMPS/PARDISO it would silently solve doubled off-diagonals (subtly wrong
smoothed fields, no blow-up); now `Unsymmetric`.

**Regression to note:** the `#DIAG` / probe cleanup reverted
`cl_FEM_Controller.cpp` wholesale, which also rolled back the coupled
`mIteration2` divergence-guard fix (`iterate_coupled` is back to
`residual( mIteration )` with no `++mIteration2`, so `tThermalReset` at
`:682` compares the frozen `0 > mMinNumIterations2` again — the dead guard
is live once more). Dormant now that the MUMPS root cause is fixed, but the
two `todo/iterate_refactor_plan.md` bullets ticked earlier are stale until
the guard advance is re-applied.

## RESOLVED: PETSc KrylovMethod::AUTO not resolved (follow-on)

Switching the thermal linear solver to PETSc (`library: petsc`, no explicit
`krylov method` key → default `KrylovMethod::AUTO`,
`cl_SolverParameters.hpp:50`) crashed in `KSPSetType()` with PETSc error 86
("type name doesn't match any registered type"). Cause:
`to_string(KrylovMethod::AUTO)` returns the literal `"auto"`
(`en_SolverEnums.cpp:291-294`), handed straight to `KSPSetType`
(`cl_SolverPETSC.cpp:473-474`); the AUTO→concrete resolution the enum
comment promises ("for PETSc: GMRES") was never implemented.

The fix must be PETSc-local, NOT at the enum/`to_string`/parser level
(Codex-confirmed, two cites): STRUMPACK genuinely consumes the same AUTO
value, but as an enum — `strumpacktools.cpp:125-131` maps
`KrylovMethod::AUTO` → `strumpack::KrylovSolver::AUTO` (REFINE uncompressed
/ GMRES compressed) — so AUTO must survive as a distinct value; and the
string parser `krylov_method()` (`en_SolverEnums.cpp:305-318`) round-trips
through `to_string`, so redefining `to_string(AUTO)` to `"gmres"` would make
an explicit `"auto"` in input unparseable and alias it away.

Fix (Claude, Codex-confirmed — Codex independently re-derived the same
patch at the same site): resolve AUTO at the PETSc boundary in
`set_krylovmethod` (`cl_SolverPETSC.cpp:470-488`) — LU preconditioner →
`PREONLY` (direct), any iterative preconditioner → `GMRES`; explicit user
choices (CG/BCGS/…) pass through untouched. Default preconditioner is
JACOBI (serial) / GAMG (parallel), so AUTO → GMRES in the normal path.
Codex brief + findings: `tmp/ai_exchange/petsc_auto_krylov_codex.md`.

Same family, second instance: with krylov resolved, the PETSc run next hit
`Invalid reordering method for PETSc: automatic`. Identical shape — the
`ReorderingMethod` default is `AUTOMATIC` (`cl_SolverParameters.hpp:31`),
which `PETSC::set_matrix_ordering` (`cl_SolverPETSC.cpp:561-573`) did not
handle → `BELFEM_ERROR`, while STRUMPACK swallows it via a `default:`
pass-through to its internal default (`strumpacktools.cpp:91-95`). Fix:
add `AUTOMATIC` to the existing NATURAL/METIS/SCOTCH case →
`MATORDERINGNATURAL` (the value this wrapper already uses for every handled
method; ordering affects only factorization fill/perf, never correctness,
and the default JACOBI/GAMG PCs don't factorize). Note left in code: the
wrapper collapses METIS/SCOTCH to NATURAL too — a pre-existing limitation,
not fixed here.
