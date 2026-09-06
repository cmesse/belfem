# Coulomb Gauge Penalty for the h-φ Formulation: Theory and Effectiveness {#fem_maxwell_coulomb_gauge_penalty_theory}

**Date:** 2026-08-24
**Purpose:** Derive a dimensionally consistent Coulomb-gauge penalty for the
h-formulation, establish where it can and cannot act in BELFEM's edge-element
spaces, and record what the quench-deck conditioning numbers actually measure.
**Module:** fem/maxwell, fem/interpolation/nedelec
**Status:** Theory plus measured evidence. The gradient operator `G` is
implemented and tested on every edge-function class except `EF_LINE3`; the
penalty is implemented in the h-kernels behind the deck key
`nonlinear magnetic { coulomb gauge penalty { chi } }`. It is **opt-in** since 2026-09-01: an absent block means `chi = 0` and the term is not assembled; a deck switches it on by stating a positive `chi` (the value in `cl_IWG_Maxwell.cpp` is only the fallback for paths that never reach `Controller::set_params`). It was on by default at `chi = 1e-4` from 2026-08-27 to 2026-09-01; on the tape decks that value was invisible to the conditioning estimate, which is why the default went back to off.

---

## 1. Results in brief

A Coulomb-gauge penalty on the h-formulation is meant to regularize the curl
null space of the discrete operator. The results below say where that works and
where it cannot:

1. **A dimensionally consistent penalty exists**, and it is not the obvious one.
   Matching the K-channel forces `γ = χ·ρ*/μ²`, giving
   `P(h,v) = χ·ρ*·∫(∇·h)(∇·v) dV` for element-wise constant μ (§6). A form
   scaled as `χ·μ/Δt` acting on `μ²G′G` is not a valid SI expression with
   dimensionless χ, and its companion mass-rate term is not a divergence penalty
   at all (§6.3).
2. **On the lowest-order simplex elements (TET4, TRI3) and on the thin-shell
   kernels (PENTA6TS, QUAD4TS), every element-local penalty built from ∇h is
   structurally unable to regularize the curl null space.** Those null-space
   fields are element-wise constant, so their gradient vanishes inside every
   element, including cohomology (cut) modes (§4.3). The divergence of every
   simplex-family basis function is also identically zero.
3. On TET4 the full-gradient Gram matrix satisfies **∫∇wₐ:∇w_b = ½∫curl wₐ·curl w_b
   exactly**, so a `G′G` penalty there is precisely an artificial resistivity
   increase `Δρ = γμ²/2` on the existing curl stiffness: it perturbs the E–J
   physics and gauges nothing.
4. Channel-by-channel effectiveness elsewhere: the **divergence** penalty is
   nonzero only on non-box HEX8 (largely a mesh-shear artifact) and on the
   higher-order families (TRI6, TET10); the **full G′G** penalty is additionally
   non-vacuous on every hex mesh including axis-aligned boxes, because the hex
   curl kernel contains non-constant gradient modes (∇(xy) = (y,x,0) lies in the
   lowest-order hex Nédélec space, is curl-free, and has |∇h|² = 2).
5. **The condition numbers that motivated this work are not a pathology of the
   gauge.** Measured over 97 timesteps of a near-critical quench run, κ(A) sits
   at 1e17…7e19, is flat across the approach to a timestep collapse, and
   *anti-correlates* with solver difficulty (§5.3). The value is characteristic
   of the mixed h-φ formulation. A penalty of relative size χ cannot move it, and
   the kernel-versus-range gap it could close is a factor 1e2…1e8, not 17 decades.
6. **The G-operator is worth having regardless:** it supports divergence
   diagnostics, ∇H postprocessing, and any future weak-divergence or jump-based
   gauging. It is
   built and gated (§8).

---

## 2. Symbol table

Conventions follow Messe et al. 2023 (Eqs. 2–8) and the code in
`src/fem/maxwell/matrices/mt_maxwell_h.cpp` (anchors are greppable tokens, not
line numbers).

| symbol | meaning | code counterpart |
|---|---|---|
| h | magnetic field strength (edge dofs = circulations along edges, unit A) | dof vector `q` |
| b = μh | magnetic flux density | none |
| ρ | electric resistivity (scalar; HTS power law ρ(\|j\|,T,B,θ)) | `mx->compute_rho` |
| μ | magnetic permeability (μ0 in HTS/metal; μ(\|h\|) in ferro-conductors) | `mx->compute_mu` |
| e | electric field | boundary term only |
| j = ∇×h | current density | `mx->compute_j` |
| δh, v | test function (same edge space) | none |
| w_k | edge basis function k | column k of `E` |
| E | interpolation operator, d×n, h(x) = E·q | `Calculator::E(k)` / `EF_*::%E` |
| C | curl operator, 3×n (3D) or 1×n (2D) | `Calculator::C(k)` / `EF_*::%C` |
| G | gradient operator (§8) | `Calculator::G(k)` / `EF_*::%G` |
| M | mass matrix ∫ μ EᵀE dV | `aMatrices->M()` in `h_picard` |
| K | stiffness ∫ ρ CᵀC dV | `aMatrices->K()` |
| Δt | timestep (`mDeltaTime`) | `h` in BDF comments |
| α, βᵢ | variable-step BDF coefficients | `mAlpha`, `mBeta` |
| qhist | β-weighted dof history | `collect_qhist` |
| χ | dimensionless penalty gain | deck key `chi`, `penalty( 2 )` |
| γ | dimensional penalty coefficient (§6.2) | none |
| ρ* | reference resistivity (element, block or global) | `save_resistivity` |
| n | number of edge dofs per element | `mNumDofs` |
| d | spatial dimension | none |

---

## 3. The implemented weak form and assembly contract

The h-formulation weak form, Messe et al. 2023 Eq. (6):

    ∫_Ω δhᵀ ∂(μh)/∂t dV + ∫_Ω (∇×δh)ᵀ ρ (∇×h) dV + ∮_Γ δhᵀ (n×e) dS = 0.

As implemented (`mt_maxwell_h.cpp`, kernels `h_picard`, `h_newton_mu0`,
`h_newton_mu`):

    M += EᵀE·(μ·w·dV)        K += CᵀC·(ρ·w·dV)

Timestep assembly (`cl_IWG_Timestep.cpp`, `bdf1`…`bdf5`; `cl_TimestepMatrices.cpp`
`assemble_dJdx`): the element fixed-point system is

    A = α·M + Δt·K,          b = Δt·f + M·qhist,

and the Newton tangent adds `dJdx = α·dMdX_times_x − dMdX_times_h +
Δt·dKdX_times_x − Δt·dFdX` into a second global matrix.

Solve structure (`cl_FEM_DofMgr_SolverData.cpp`, `SolverAlgorithm::NewtonRaphson`
branch): Picard solves A·q_new = b directly with relaxation ω; Newton computes the
residual r = A·x − b at the current iterate x, solves (A + dJdx)·Δ = r, and updates
x ← x − ω·Δ. Two global matrices exist: the fixed-point matrix (used for residual
evaluation) and the tangent (used for the increment solve).

**Assembly placement consequence used below:** anything added into the element
K-position is automatically Δt-scaled into A and enters the residual through A·x;
anything in the M-position is automatically α-scaled into A and contracted with
qhist on the RHS. A penalty that is linear in h with frozen coefficients needs
no dJdx channel and no explicit RHS term.

---

## 4. Continuous theory: what a Coulomb gauge penalty is supposed to do

### 4.1 Gauss's law is already conserved, in the continuum and weakly in the discrete space

Taking the divergence of Faraday's law gives ∂(∇·b)/∂t = −∇·(∇×e) = 0: the
continuous evolution preserves ∇·b = 0 exactly if the initial condition satisfies
it (Monk 2003, §1.2 makes the same point for the Maxwell system: the divergence
conditions are consequences of the curl equations plus initial data).

Discretely: the space of gradients of the scalar (nodal) space is a subspace of
the Nédélec space, ∇S_h ⊂ X_h, and **for a simply connected domain with
connected boundary** the null space of the discrete curl is exactly
N_h(curl) = ∇S_h (Monk 2003, §7.2.1, where that hypothesis is stated). On a
multiply connected conductor the kernel is ∇S_h plus discrete harmonic
(cohomology) fields; §4.3 shows the argument below covers those as well. Testing
the discrete h-formulation with v = ∇ξ_h kills the K term (∇×∇ξ = 0), leaving
d/dt ∫ μ h·∇ξ_h dV = 0. The Galerkin scheme conserves the *weak* divergence of
b against all discrete gradients, step by step, up to solver tolerance. This is
the discrete analog of the continuum statement, and it is the reason the
h-formulation "does not require any choice of gauge" in the conductor (Arsenault
et al. 2021, introduction; the same paper notes that first-order curl elements
satisfy ∇·B = 0 only element-locally, with the deficit in the normal jumps, and
that transient runs keep the property when the initial condition is solenoidal).

**Consequence:** a penalty cannot restore a lost conservation property, because
none is lost. Its only legitimate purpose here is **spectral**: to move the
curl-null-space eigenvalues (currently set by α·M alone) closer to the curl-range
eigenvalues (set by Δt·K), reducing the condition number. That reframing matters,
because it makes the success criterion measurable: the penalty must raise the
gradient-mode diagonal without perturbing the physical (curl-range) modes.

### 4.2 The classical Coulomb-gauge penalty and its known limits

The classical construction (Bíró & Preis 1989, for the A-formulation) adds
γ(∇·u, ∇·v) to a curl-curl form to ellipticize it. Monk 2003 §7.4 ("The
ellipticized Maxwell system") analyzes exactly this and proves two things:

1. With γ large enough the penalized continuous problem has the same solution:
   the penalty is *consistent*, vanishing on the exact, solenoidal solution.
2. Its conforming discretization requires X_N = H(curl)∩H(div) conformity, which
   for piecewise polynomials forces nodal (H¹)³ elements (Monk Theorem 5.3 as
   used in §7.4), and on non-convex domains those converge to the **wrong**
   solution (Monk §3.8, Lemmas 3.55–3.56; §7.4: "we have the terrible situation
   that we compute a convergent solution … Ẽ ≠ E"). Boffi et al. 2013 §11.4
   ("Enforcing the Divergence-Free Condition by a Penalty Method") makes the same
   point for the eigenvalue problem: "this strategy can be very dangerous:
   although the method is stable, on general domains, it can converge to wrong
   values!", tracing it to Costabel's coercivity result. Monk §7.4 rules out the
   naive edge-element transfer in one sentence: "The straightforward answer is
   'no', because functions in Xh do not have a well-defined divergence."

Edge elements exist precisely to avoid that pathology (Monk, Ch. 5 introduction).
Monk §7.4 then lists the edge-element-compatible substitutes: a *weak* discrete
divergence ∇_h·, which requires a mass-matrix inversion made local by mass
lumping (Eqs. 7.43–7.44), or explicit divergence enforcement via a Lagrange multiplier.
The natural gauge for edge spaces is the tree–cotree gauge (Dular et al. 1997
§III: "natural gauge condition in the space of edge finite elements, which
requires the construction of a tree").

The operative question for BELFEM is therefore not how the Coulomb penalty is
defined, but whether an *element-local* penalty built from ∇h does anything useful in
these edge-element spaces. §6–§10 answer that quantitatively.

### 4.3 A structural no-go for element-local penalties at lowest order on simplices

Members of ∇S_h are gradients of continuous piecewise-linear scalars:
**element-wise constant vector fields**. The cohomology/cut fields that join the
kernel on multiply connected conductors are element-wise constant too: On a
single tet, every Whitney ∇w is antisymmetric. A Whitney-1 field with zero curl
therefore has no antisymmetric gradient part, so ∇h = 0 on that element and the
field is constant there. Tangential continuity glues these per-element constants
into a field that fails to be a global single-valued gradient around a
non-contractible loop, but *inside every element* ∇h ≡ 0 still holds. BELFEM's
h-φ cuts are moreover largely φ-side jumps plus multiplier rows (Messe et al.
2023, §2), untouched by any volume term.

Any element integral whose integrand is built from the pointwise gradient ∇h
(the full Jacobian, its trace ∇·h, ∇·(μh) with element-wise constant μ, or any
quadratic form in these) evaluates to zero on such fields. Hence on TET4/TRI3
and the thin-shell kernels (§12) *no element-local ∇-based penalty can raise the
curl-null-mode eigenvalues at all*, cohomology modes included. The discrete
divergence content of an edge-element field lives in the **normal jumps across
element faces** (Monk §5.5.1: "there is a singular contribution to the divergence
at the faces in the mesh"), reachable only by face/jump terms, a weak
(mass-inverted) divergence, or graph-based gauging.

**Scope caveat:** this element-wise-constant argument is a *simplex/Whitney*
property. The lowest-order **hex** kernel contains non-constant curl-free fields
(∇(xy) = (y,x,0) ∈ Q_{0,1,1}×Q_{1,0,1}×Q_{1,1,0}, curl-free, ∇h ≠ 0), so on HEX8
meshes a full G′G penalty *does* reach part of the kernel. §10 gives what that
costs, and §11 gives why it still does not move the condition number at usable χ.

---

## 5. What the conditioning numbers actually are

### 5.1 Spectral bookkeeping of A = αM + ΔtK

Order-of-magnitude entries for mesh size h_m (E ~ 1/h_m, C ~ 1/h_m², dV ~ h_m^d,
shown for d = 3):

- curl-range modes: Δt·K entries ~ Δt·ρ/h_m
- curl-null modes: α·M entries ~ α·μ·h_m

Their ratio r = Δt·ρ/(α·μ·h_m²) is the gap a gauge could close. Concretely, with
μ0 = 1.26e-6 and α = 1:

| h_m | Δt | ρ | r |
|---|---|---|---|
| 1 mm | 1e-4 s | 1e-6 Ω·m (Hastelloy-class) | ~80 |
| 1 mm | 1e-4 s | 2e-9 Ω·m (Cu at 77 K) | ~0.16 |
| 1 μm | 1e-4 s | 1e-6 Ω·m | ~8e7 |
| 1 μm | 1e-2 s | 1e-6 Ω·m | ~8e9 |

**Nothing in this range produces 17 decades.** A perfect gauge, whether
tree–cotree or an ideal weak divergence, closes at most this factor. That bounds what any
gauging work can be expected to achieve on the condition number.

Note the direction of the Δt dependence: if λ_max ~ Δt·λ_K and λ_min ~ α·λ_M,
then κ ∝ Δt. Collapsing the timestep during a quench *shrinks* this particular
split rather than widening it.

### 5.2 The other drivers, which a gauge penalty cannot touch

- **ρ contrast inside K:** superconducting zones sit at power-law ρ that is
  **unfloored** (`mRhoMin = 0.0`; the older 1e-16 floor was removed for
  Newton-tangent consistency, and ρ_PL may be IEEE 0 with the parallel
  combination staying clean). Quenched zones sit at ~1e-7…1e-6 Ω·m. The spread
  *within the curl range* is therefore unbounded from below.
- **Newton tangent stiffening:** the power-law channel `dKdX_times_x` scales like
  (n_powerlaw − 1)·ρ locally, another factor ~25–40 on quenching elements.
- **Lagrange-multiplier rows:** thin-shell/φ interfaces and cohomology cuts use
  multiplier couplings with zeros on the main diagonal (Messe et al. 2023, §2,
  the static-condensation discussion; the free cut-λ dof is legitimately free).
  These produce indefinite blocks whose conditioning is independent of any gauge.
- **Mixed-block scaling:** the H and φ blocks and the constraint rows carry very
  different scalings by construction. A mixed h-φ system reaching κ ~ 1e17 is
  characteristic of the formulation, not evidence of singularity.
- **Mesh grading:** tape-scale meshes grade h_m over 2–3 orders.

### 5.3 What was measured

A near-critical quench run (I ≈ 1.1·Ic, 97 completed steps) with the
conditioning diagnostic on in the `linear magnetic` block, κ supplied by
MUMPS from its own error analysis. (Measured before 2026-08-30; the MUMPS
numbers are now behind their own key, `mumps error analysis` — under the
current deck contract `compute conditioning` alone yields the eigenvalue
ratio instead, a different quantity.)

| quantity | value |
|---|---|
| range | 1.1e17 … 7.4e19 |
| median on steps converging in ≤ 6 iterates | 5.9e17 (n = 85) |
| median on steps needing ≥ 12 iterates | 3.0e17 (n = 5) |

Two results matter more than the magnitude:

1. **κ is flat across the approach to a timestep collapse:** there is no trend
   into the failure.
2. **κ anti-correlates with difficulty.** Hard steps had roughly half the κ of
   easy ones; the worst-conditioned step of the whole set (7.4e19) converged in
   five iterates.

The run actually stopped on a residual plateau at −90…−96 dB. Measurements
excluded timestep, linear solver, linear precision, and conditioning as causes;
loosening the nonlinear tolerance from 1e-11 to 1e-9 was what got past it. **A term that only improves κ therefore has no measured problem to
solve on that deck.** What a gauging design must answer first is whether that
residual plateau lives in the gradient null space of the curl-curl operator, a
question settled by projecting a converged residual onto the discrete-gradient /
cotree subspace, not by tuning a penalty.

κ itself needs careful interpretation: on the MUMPS path the reported number is COND1
(`RINFOG(10)`), a componentwise 1-norm condition estimate of a nonsymmetric mixed
system — computed on the **original** matrix with its actual right-hand side, not the
equilibrated one: MUMPS scales for the factorization (`ICNTL(8) = 77`), but its error
analysis hands the condition estimator an identity weight vector (MUMPS 5.9.1
`dsol_driver.F`, "Notice that D is always the identity"), so no scaling enters the
estimate. (An earlier version of this paragraph said "already-scaled"; corrected
2026-08-30.) The companion quantities from
the same factorization (COND2 and the forward and backward errors) say more
about whether a solve is actually losing the solution than COND1 does.

---

## 6. A consistent penalty

### 6.1 Derivation

Target functional (Coulomb gauge on b): Π = ½ ∫ γ |∇·(μh)|² dV, γ > 0 to be
determined. First variation against the edge test space:

    P(h, v) = ∫ γ (∇·(μh)) (∇·(μv)) dV.

Properties: (i) symmetric positive semi-definite; (ii) *consistent*: the exact
solution has ∇·(μh) = 0, so the penalty vanishes on it and adds no modeling
error in the continuum limit; (iii) linear in h for h-independent μ and γ.

It belongs in the **stationary part** of the equation (the K-channel), because the
constraint ∇·b = 0 is an algebraic condition on the state at t^{n+1}, not on the
rate. The rate variant is examined in §6.4.

### 6.2 Dimensional closure fixes γ

Bookkeeping in SI with edge dofs in A (E ~ 1/m, G ~ 1/m², dV ~ m³):

| object | entry units | ×dof (A) |
|---|---|---|
| M = ∫μEᵀE | μ·m = V·s/A | V·s |
| K = ∫ρCᵀC | ρ/m = Ω = V/A | V |
| P = ∫γμ²Gdivᵀ Gdiv | γ·μ²/m | none |

Matching P to the K-position requires γ·μ²/m = V/A, i.e.

    γ = χ · ρ* / μ²,   χ dimensionless,   ρ* a reference resistivity,

giving **P(h,v) = χ·ρ*·∫ (∇·h)(∇·v) dV** for element-wise constant μ: the μ²
cancels pointwise, so nothing is gained by carrying it. χ then has the clean
meaning *penalty resistance per physical resistance*.

**The ρ\* fork.** The choice of averaging scope is a genuine design fork, not a
detail:

- *element* ρ̄: the penalty scales covariantly with the local physics, but it
  **vanishes on superconducting elements** (ρ_PL is unfloored and may be exactly
  zero), which under the spectral reading is exactly where α·M is smallest and
  stiffening would be wanted;
- *block* or global reference ρ*: the penalty stays alive in SC zones, but there
  χ·ρ*/ρ_local is unbounded, i.e. the consistency error dominates the physical
  resistivity precisely where the HTS physics lives.

The implemented kernels use the live pointwise ρ, namely the first branch. A
floored coefficient must **not** be combined with the full G′G form: on TET4 that
is an exact resistivity injection Δρ = (χ/2)·ρ_floor on every curl mode, which
turns subcritical HTS into a poor metal and silently undoes the unfloored power
law. A floored coefficient with a divergence-only operator avoids that, but the
divergence row is identically zero on most of BELFEM's element families (§10),
so it is close to the zero operator on the decks of interest.

### 6.3 A tempting alternative form, and why it fails

A natural-looking construction adds, per element and Gauss point,
`α·[ μ²(G′G) + μ E′E dH/dt ]·w·dV` with `α = χ·μ/Δt` and χ dimensionless.
Both terms fail, independently:

**Term 1, α·μ²·G′G.** Units: (χμ/Δt)·μ²/m = χ·μ³/(Δt·m) = χ·V³s²/(A³m⁴) per
entry. The K-position needs V/A; the M-position needs V·s/A. Neither matches, so
χ would have to carry units A²m⁴/(V²s²). **With dimensionless χ this is not a
valid SI formula.**

**Term 2, α·μ·E′E·dH/dt.** Two independent problems. *Dimensions:* entry·rate =
χμ²·m·A/(Δt·s), again neither V nor consistent with Term 1. *Structure:* E′E is
the mass Gram; it contains no derivative of h at all. A term ∝ E′E·dH/dt is an
**artificial addition to the permeability**, not a divergence penalty of any
kind: it acts identically on curl-range and curl-null modes, so it cannot
selectively fix the gradient subspace, while it *does* shift the physical
eddy-current time constant.

**Where the form comes from.** Using the identity ∫|∇u|² = ∫|∇×u|² + ∫|∇·u|² +
boundary terms, one may try to build a divergence penalty as "full-gradient
penalty minus curl-curl content", then eliminate the curl-curl part through the
strong equation ∇×(ρ∇×h) = −∂(μh)/∂t, which turns it into a mass-rate term. That
produces exactly the shape above. It fails twice more: the substitution is valid
only *at* the converged solution, so inside a Newton loop it changes the tangent
and the transient; and on simplices the discrete identity degenerates, because
∫∇wₐ:∇w_b = ½∫curlₐ·curl_b holds **exactly** on TET4. "G′G minus curl content"
is then −½·C′C, a negative curl stiffness carrying no divergence information.

**Choice of channel.** In BELFEM's assembly path, the channel determines whether
a contribution belongs in the matrix or the RHS: a coefficient times E′E·dH/dt is a mass-channel term
and splits automatically (α-part into the matrix, qhist part into the RHS); a
stationary constraint G′G·h is K-channel, matrix only, Δt-scaled and
residual-consistent through r = A·x − b. Nothing should ever be added to the RHS
by hand; the BDF contract owns that split.

### 6.4 The rate form, for completeness

Penalizing ∂(∇·b)/∂t (G-terms contracted with dH/dt) enters the M-channel and
only prevents divergence *growth*: any divergence content present in the initial
state or introduced by a restart stays frozen. Combined with §4.1 (the weak
divergence is already conserved) the rate form has no target left.

---

## 7. Newton linearization and assembly placement

For P = χ·ρ*·∫(∇·h)(∇·v) dV, or the full-gradient variant, with ρ* frozen during
the nonlinear loop:

- **Matrix:** add P into the element K. Assembly then Δt-scales it into
  A = α·M + Δt·K automatically, the Picard branch sees it in the fixed-point
  matrix, and the Newton branch sees it in both the residual r = A·x − b and the
  tangent. No new code path.
- **dJdx:** zero, **iff ρ* is frozen**. If ρ* is the live power-law ρ(|j|), a
  strictly consistent tangent gains a cross term
  χ·(Gᵀg)·(Cᵀj)ᵀ·(∂ρ/∂|j|)/|j| with g = G·q, plus χ-scaled copies of the
  ρ(|B|,β) field-derivative channel. Omitting them costs Newton *rate* at
  O(χ), never correctness, because the residual is exact either way because the
  penalty is fully present in K. **Since 2026-08-27 the j-channel cross term IS
  implemented** in `h_newton_mu0` / `h_newton_mu` (`mt_maxwell_h.cpp`, inside
  the `tUseGauging` block, guarded like the power-law tangent; workspaces
  `Gq`/`GtGq`). The χ-scaled ρ(|B|,β) copies remain omitted — still an O(χ)
  rate cost, for HTS and metal kernels alike.
- **RHS:** nothing. A K-channel term contributes to the residual only through
  A·x, which is what a constraint on the end-of-step state should do.
- **Variable μ (`h_newton_mu`):** ∇·(μ(|h|)h) = μ∇·h + (dμ/d|h|)(∇|h|)·h is
  nonlinear in h. The implemented term penalizes ∇h rather than ∇·(μh), so in
  ferro-conductors it biases toward ∇·h = 0 rather than ∇·b = 0. This is a weak
  preference at small χ, and the gauge choice the implementation makes.

---

## 8. The G operator

**G is the full Jacobian of the interpolated field, not only its divergence.**

At integration point k, G(k) is the (d·d) × n matrix with

    G( i + d·j , e ) = ∂(w_e)_j / ∂x_i        i, j ∈ {0…d−1},  e ∈ {0…n−1}

(column-major vec of the d×d tensor ∇h with convention (∇h)_{ij} = ∂h_j/∂x_i).
Shapes: 9×n in 3D, 4×n in 2D. The binding layout contract, including the explicit
curl-tie rows, lives at the `mGrad` member in
`src/fem/interpolation/nedelec/cl_EF_EdgeFunction.hpp`; `Calculator::G( aIndex )`
is the assembly-side accessor.

Derived quantities:

- divergence row: ∇·h = Σ_i G(i + d·i, :)·q, the trace in one contraction;
- consistency tie: the antisymmetric part reproduces C exactly, e.g. in 3D
  (curl h)_x = [row(7) − row(5)]·q. This identity is the anchor test of every
  implementation and is gated per element family in `tests/fem/`.

Why full G and not a 1×n divergence row:

1. For TET4, TRI3 and PENTA6TS the divergence row is **identically zero**, and
   for QUAD4TS it is zero for orthogonal stacking (§10), so a divergence-only
   operator would be a literal zero matrix for the main element families,
   untestable and useless.
2. Full G supports every downstream consumer: both penalty variants, the
   antisymmetry-versus-C test, div-b diagnostics, ∇H postprocessing (forces,
   Maxwell stress), and any future jump or weak-divergence gauging.
3. The penalty needs only the trace, and taking a trace of G is one loop. A
   divergence row cannot recover full G.

Cost: everything G needs (the per-point inverse Jacobians and reference
derivatives) is already computed for E and C; G adds O(d²·n) multiplies per
Gauss point, the same order as one E evaluation. `EF_LINE3` deliberately does
not implement G: on a 1D manifold element, the ambient gradient is not defined
from the element alone. Its `G()` is a permanent stub that always fails.

---

## 9. μ dependence and where a penalty can act

In h-φ, edge dofs exist **only in conducting domains** (air/vacuum and
non-conducting iron carry φ; Messe et al. 2023 Eqs. 6–8). Therefore:

- **Air:** no edge dofs, so no penalty is applicable. None is needed, because the
  φ-Laplacian is symmetric positive definite and well conditioned.
- **HTS and normal metals:** μ = μ0, element-wise constant, so ∇·(μh) = μ0·∇·h
  exactly and the penalty reduces to χ·ρ*·(∇·h, ∇·v) with no jump term inside a
  block.
- **Ferro-conductors:** μ = μ(|h|), see §7.
- **Material interfaces:** with piecewise-constant μ the physical condition is
  the jump condition n·[μh] = 0. An element-local penalty cannot see interface
  jumps at all, since it integrates over element interiors; the weak Galerkin
  statement already controls n·[b] in the ∇S_h-weak sense (§4.1). There is nothing further to do at
  interfaces, and nothing gained.

---

## 10. Divergence of the basis functions, per element family

| family | dofs | geometry map | ∇·(basis) element-wise |
|---|---|---|---|
| TRI3 | 3 | affine | ≡ 0 (Whitney, 2D) |
| TET4 | 6 | affine | ≡ 0 (Monk §5.5.1 states it for the Whitney element) |
| TRI6 | 8 | straight or curved | ≠ 0 for the higher-order members (div ∇(quadratic) = const ≠ 0); the curved path adds mapping terms |
| TET10 | 20 | straight or curved | ≠ 0 for higher-order members; curved adds mapping terms |
| HEX8 | 12 | trilinear | box ≡ 0, sheared affine ≠ 0, general trilinear ≠ 0. **But** the full ∇h of kernel modes is ≠ 0 even on boxes (the ∇(xy) mode, §4.3) |
| PENTA6TS | 6 | flat prism, thickness τ | ≡ 0 (∇w antisymmetric-tangential; n·w = 0) |
| QUAD4TS | 2 | 2D thin shell | ≡ 0 **only for orthogonal stacking**: the basis is s·F(η)·∇ξ, so ∇·h = F′(η)·(∇η·∇ξ) ≠ 0 when the stacking is sheared. Kernel fields (equal top/bottom dofs, hence constant h) stay blind regardless |
| HEX8TB | 4 | side-connector wall, imposed orthonormal cuboid frame | ≡ 0 identically: the basis is s·F(η,ζ)·∇ξ and ∇η, ∇ζ ⊥ ∇ξ by construction |
| HEX8TS | 8 | thin-shell machinery, per-point mid-surface nablas | ≡ 0 on rectangular mid-surfaces; nonzero under sheared stacking, as for QUAD4TS |
| LINE3 | 2 | 1D manifold | divergence not meaningful on the curve |

**Mechanism for the affine cases.** Under the covariant Piola map h = J⁻ᵀĥ with
constant J, the physical divergence is ∇·h = (J⁻¹J⁻ᵀ) : ∇̂ĥ: the contraction of
the symmetric inverse metric with the reference Jacobian of the basis. Whitney
simplex functions have *antisymmetric* constant ∇̂ĥ, so the divergence is exactly
zero for every affine map. The hex reference functions have symmetric off-diagonal
∇̂ĥ content, so their divergence is zero only when the metric is diagonal (a box)
and nonzero under shear. On hex meshes the penalty magnitude is therefore partly
a measure of element *shape*, not of physics.

**Implication.** On TET4/TRI3/PENTA6TS/HEX8TB meshes the div-penalty matrix is
identically the zero matrix; on QUAD4TS and HEX8TS it is zero for orthogonal
stacking and shear-dependent otherwise. Per §4.3 the full G′G penalty cannot
touch the curl null space on any of them, so a nonzero basis divergence under
shear buys nothing. On mixed meshes the penalty switches on
and off by element type, which is itself a conditioning hazard.

---

## 11. Choosing χ, and what it can achieve

The assembled system is A = α·M + Δt·K. Two candidate scalings place the penalty
on different sides of that sum, and adaptive stepping makes the difference
operational.

**Scaling the penalty as χ·μ/Δt on a K-channel term.** The Δt from assembly
cancels the 1/Δt: the assembled penalty becomes **Δt-independent**, i.e.
proportional to the mass block. Then the penalty-to-K ratio goes as 1/Δt, so at
large Δt the penalty fades relative to the physics, and at collapsed quench
timesteps it *grows*, giving maximum distortion at maximum fragility. The gap it is
meant to close, r = Δt·ρ/(α·μ·h_m²), is itself Δt-dependent, so a Δt-free penalty
cannot track it.

**The consistent scaling γ = χ·ρ*/μ² (§6.2), K-channel.** Penalty and physical
stiffness carry the *same* Δt factor. The penalty-to-K ratio is
χ·(ρ*/ρ_local)·(div-content / curl-content), **Δt-invariant** across the adaptive
range, and the gradient-mode diagonal lifts to the same Δt·ρ*/h_m² scaling as the
curl modes.

**What χ can and cannot do.** It is easy to conflate two thresholds:

1. *Relative to K.* A K-channel term of relative size χ moves the eigenvalues it
   touches by a relative χ. At χ = 1e-4 that is a fourth-decimal effect on any
   condition number, invisible by construction.
2. *Relative to the mass block.* Lifting curl-null modes that are held by α·M
   requires χ·ρ ≳ μ0·L²/Δt. At L = 1 mm and Δt = 1e-4 s that threshold is
   1.26e-8 Ω·m, while χ·ρ = 1e-4 × 1e-6 = 1e-10 Ω·m in metal, a 0.8 % lift, and
   exactly zero in unfloored superconductor. A one-decade lift on the hex
   kernel modes would need χ ≳ 0.1 at ρ = 1e-6 Ω·m, which is far outside the
   consistent regime.

A null result at χ = 1e-4 is therefore **expected on every element family**, including
the hexes where the term is not vacuous. A K-channel penalty can visibly move a
near-1/ε condition number only if it fills *exact zeros* on the small-end
subspace, which it does nowhere on BELFEM's current meshes. Combined with §5.1
(the kernel-versus-range gap is at most 1e2…1e8) and §5.3 (κ does not predict
solver difficulty on the measured deck), χ tuning is not the lever it appears to
be.

If χ is ever needed, validate it by measuring the quantity it is supposed to
improve with and without the penalty on the same deck, not by argument alone.

---

## 12. Thin-shell elements

- **Penalty:** does not apply. PENTA6TS basis divergence ≡ 0; QUAD4TS and HEX8TS
  basis divergence is nonzero only under sheared stacking while their kernels
  stay blind; and the thin-shell h-space is a constrained trace space whose
  null-space structure is dominated by the interface and multiplier coupling
  rather than by an interior gauge freedom.
- **Aspect ratio does not rescue it.** On a thin layer the G entries grow as 1/d
  through ∇ζ = 2n̂/d, but the through-thickness curl content in CᵀC grows
  identically, so the penalty stays bounded by χ times the local resistive
  stiffness however thin the layer gets.
- **G-operator:** implemented anyway, for uniformity of the `EdgeFunction`
  interface and for diagnostics: the through-thickness derivative
  ∂h/∂n = (h_top − h_bot)/t is a physically meaningful output, carried by the
  ∇τ ⊗ w term of the PENTA6TS gradient.
- **HEX8TS scope:** its per-point nablas are a thin-shell convention rather than
  the gradients of an exact inverse map, so its G reproduces the class's own curl
  on every geometry but equals the full gradient only on affine mid-surfaces. The
  three-case scope is documented at the class declaration.

---

## 13. 2D reduction

In 2D (TRI3/TRI6; `mE` is 2×n, `mC` is 1×n scalar):

    G is 4×n:  rows (∂x h_x, ∂y h_x, ∂x h_y, ∂y h_y);
    ∇·h = row(0) + row(3);
    curl h (scalar, out-of-plane) = row(2) − row(1).

Penalty reduced form: P = χ·ρ*·∫ (∂x h_x + ∂y h_y)(∂x v_x + ∂y v_y) dA.
The effectiveness verdict is identical to the 3D simplices: TRI3 divergence ≡ 0
element-wise and the curl null space is element-wise constant, so the penalty is
vacuous on the 2D workhorse meshes. TRI6 carries nonzero divergence in its
higher-order members only.

---

## 14. Verdict and alternatives

An element-local ∇h penalty should not be used as the conditioning strategy for
BELFEM's simplicial or thin-shell quench systems: on those families it is either exactly
zero or exactly an artificial resistivity, and on the families where it is
non-vacuous it cannot reach the mass scale at any consistent χ (§11). The
implemented term is retained, defaulted off, for hex-dominant meshes and
experiments.

If conditioning or the residual plateau is to be addressed at the formulation
level, the edge-element-compatible options are, in increasing order of
invasiveness:

| option | mechanism | reference | remark |
|---|---|---|---|
| a. Measure first | attribute the number: which block owns ‖A‖, what the forward error is, whether the residual plateau lives in the gradient/cotree subspace | none | prerequisite to everything below; costs one instrumented run |
| b. Scaling / equilibration | MC64 and diagonal scaling on the assembled matrix | already load-bearing on both solver paths | zero formulation risk; on the MUMPS path the reported κ is already post-scaling |
| c. Weak-divergence penalty | γ(∇_h·(μh), ∇_h·(μv)) with a lumped-mass weak divergence | Monk §7.4, Eqs. 7.43–7.44 | acts on the true (jump) divergence; needs face-adjacent assembly, not just G. Blind to discrete harmonics: Monk's Lemma 7.26 splits E_h with ∇_h·(εẼ_{h,0}) = 0, so cohomology modes survive it |
| d. Tree–cotree gauge | eliminate the curl-null subspace in conductors | Dular et al. 1997 §III (tree gauge of the source field h_s; gauging the eddy-current unknown itself is the Albanese–Rubinacci-style use); Denis et al. 2026 | exact null-space removal, no scale threshold; interacts directly with cut/cohomology generators, and enriched meshes carry twin edges that share both end nodes, so a tree keyed on edge endpoints is wrong there |
| e. Full-G Tikhonov, hex/high-order meshes only | χρ*(∇h, ∇v) | §6.2 | inconsistent: perturbs the physics at O(χ·ρ*/ρ_local), fatal in SC zones unless support-restricted to quenched elements |

The G-operator is prerequisite for (c) and (e), useful for the diagnostics in (a)
under every option, and independently valuable for postprocessing.

---

## 15. How these results were established

The dimensional algebra, the null-space arguments, and the per-family divergence
mechanism are derivations. They are given inline above and reproducible from the
element constructions in `src/fem/interpolation/nedelec/`. Four were also checked
numerically at the mathematical level before implementation: the Whitney
divergence identity on TET4, the Gram identity ∫∇wₐ:∇w_b = ½∫curlₐ·curl_b on
TET4 (max deviation 7.1e-15), the vanishing of ∇h on reconstructed
discrete-gradient fields (max 6.9e-16), and the box-versus-sheared HEX8
divergence contrast. The compiled implementations are gated separately in
`tests/fem/test_EdgeFunctions.cpp` and `tests/fem/test_InterfaceOrientation.cpp`:
finite differences of G against E, the antisymmetric part of G against the
compiled C, the divergence trace against §10, and per-family kernel and control
modes.

The conditioning statements in §5.3 are measurements from an instrumented quench
run, not estimates. The scale bounds in §5.1 and §11 are arithmetic from those
same material and mesh parameters.

---

## 16. References

- Monk, P. (2003). *Finite Element Methods for Maxwell's Equations.* §1.2, Ch. 5
  (incl. §5.5.1), §3.8 (Lemmas 3.55–3.56), §7.2.1, §7.4 (Eqs. 7.43–7.44,
  Lemma 7.26).
- Boffi, D., Brezzi, F. & Fortin, M. (2013). *Mixed Finite Element Methods and
  Applications.* §11.4.
- Messe, C. et al. (2023). BELFEM h-φ formulation, Eqs. (2)–(8), §2
  (static condensation, cut multipliers).
- Arsenault, A. et al. (2021). h-φ implementation in COMSOL, introduction
  (no gauge required in the conductor; element-local ∇·B = 0 for first-order
  curl elements).
- Dular, P. et al. (1997). Source-field computation, §III (tree gauge in edge
  element spaces).
- Denis, A. et al. (2026). Homogenized multi-scale h-φ thin shell (co-tree gauge
  of the source field in the conductor).
- Bíró, O. & Preis, K. (1989). Classical Coulomb-gauge penalty for the
  A-formulation.

Full citations and DOIs: `doc/literature_references.md`.
