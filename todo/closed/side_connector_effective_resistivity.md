# Side Conductivity of HTS Tapes — Problem Statement and Modeling Approach

**Date:** 2026-07-09 (rev. 4 — Christian's three questions answered: conditioning argument for
extra edge DOFs hardened (§2), discretization of the resistor form spelled out (§5.3.2),
half-partition physics without/with the bridge (§5.6). Rev. 3, 2026-07-06: bridge formulation)
**Purpose:** Self-contained problem statement and modeling recommendation for copper side
conductivity of HTS coated-conductor tapes (current escape around the insulating buffer layer),
replacing three abandoned geometric/enrichment attempts.
**Module:** `fem/maxwell` (h-φ formulation, thin-shell PENTA6TS)
**Status:** CLOSED as a plan (2026-08-08) — kept as the standing r′ theory reference. The
implementation route now lives in the HEX8TB wall element campaign (`todo/hex8tb_phase2_fem_wiring.md`;
theory home `src/fem/maxwell/doc/side_connector_wall_element.md`); the bridge-as-curve-kernel
design of `sideconnector_bridge_plan.md` was superseded before execution. This file remains
the source for the wrap-resistance derivation `r′ = ρ_Cu·h_path/t_w + ΣR_ct/t_w` and the
failed-routes record. Original status: ANALYSIS COMPLETE, design note pending —
no source modified. Rev. 1–2 were cross-checked by
Codex + Grok (two rounds); the rev. 3 formulation math was **Codex-audited and confirmed** (§9; one
subclaim softened, O1 hardened to blocking). Awaiting Christian's go-ahead for the design note.
Rev. 4 additions (§2 hardening, §5.3.2, §5.6) **Codex-audited and confirmed** the same day;
precision fixes applied (§9, round 4).

> **Headline.** Do not resolve, do not enrich, do not add DOFs — **couple**. The side copper is a
> collapsed **edge-wall impedance bridge**: a curve integral `∫ r′·[[h∥]]·[[v∥]] dy` along the tape's
> lateral edge, coupling the traces on either side of the wall through the physical wrap resistance.
> In the h-formulation a resistor is `R·I(h)·I(v)` where the current `I(h)` is a **jump or
> circulation of H** — a linear functional of existing DOFs, no discrete curl evaluated, no new
> unknowns. The buffer stays φ-only; its topology is handled by the cut machinery, its physical
> bypass by the bridge. The earlier effective-ρ⊥ idea survives as an optional within-layer
> refinement (§7), rather than the primary mechanism.

---

## 1. The physical problem

BELFEM models HTS coated-conductor tapes with an **h-φ formulation**: conducting regions carry the
magnetic field intensity **h** on curl-conforming Nédélec/edge DOFs; non-conducting regions carry a
scalar potential **φ** on nodes. The dissipative term in a conductor is the curl-curl integral

```
∫_V ρ (∇×h)·(∇×h′) dV          (ρ = electrical resistivity, J = ∇×h)
```

Tapes are thin-shell **PENTA6TS** prisms — triangular in the tape plane, thin in the through-thickness
direction **n**, **one material layer per element**. REBCO is strongly nonlinear
(`E = E_c (J/J_c)^n`).

**Layer stack** (top→bottom): `Cu / Ag / Hastelloy / buffer / YBCO / Ag / Cu`. The **buffer** (an
oxide barrier required by manufacturing so the REBCO grows properly on the substrate) is essentially
insulating and is modeled **φ-only**, splitting the stack into two electrically isolated conducting
halves: upper (Cu/Ag/Hastelloy) and lower (YBCO/Ag/Cu). See §3.

**The side conductivity.** Tapes are slit from wide strips and then surround-electroplated: the
slit edge exposes the *cross-sections* of every layer (≈50 µm bare Hastelloy, ≈0.2 µm buffer,
≈1.5 µm REBCO, Ag), and the Cu plating deposits onto everything conductive there, forming a
U-shaped jacket — top plate, **edge wall**, bottom plate — one connected conductor. The physical
question is: **during a quench, can current escape around the buffer through this slim edge passage
to the other side of the tape?** This governs current sharing in CORC cables (tape-to-tape transfer
passes through the substrate side and therefore *in series* through the edge wall) and is not well
understood experimentally — which is precisely why we simulate it.

### The scale-separation obstacle

Tape width ≈ 4 mm, ~10 elements across → edge triangles ≈ 400 µm; the plated wall is ≈ 20 µm thick.
This is a 20:1 length (≈400:1 area) separation. Every attempt to *resolve* this scale failed on
conditioning or stabilization (§2). The resolution below rests on the observation that nothing at
that scale needs resolving (§4).

---

## 2. What has already failed (do not repeat)

| # | Approach | Failure mode |
|---|----------|--------------|
| 1 | Mesh the copper strip explicitly | 20 µm slivers beside 400 µm bulk → destroys conditioning |
| 2 | Extra DOFs inside the cut element | No help; ill-conditioning is intrinsic to the scale |
| 3 | Enrichment (2 edge DOFs × 2 Heavisides) | Oscillations along the cut; no ghost-penalty anchor exists (below) |

**HEX8TS "wrap" element (removed 2026-06-05/20).** A lateral wrap element shorted top and bottom
through shared/free edges. Unphysical: the wrap was not "wetted" (no `compute_bn` transport analogue),
so its binormal-H at the PENTA6TS↔HEX8TS fold remained a free parameter, which the solver drove to
~0 in low-ρ copper → ~2-order binormal-H discontinuity and suppressed side-curve `J/J_c`.
Neither geometric hanging nor `h_penalty` converged. A QUAD4TS variant was rejected
(basis-incomplete). History:
`devlog/dl20260605_remove_side_connectors.md`, `dl20260620_side_connector_removal.md`,
`dl20260417_hex8ts_thirdpass_audit.md`. **Retained for reuse:**
`ThinShellFactory::compute_binomial_vectors`, `cl_EF_HEX8TS.*`, the ghost-facet machinery.

### Why ghost-penalty / Nitsche cannot rescue the enrichment

*(Christian's reasoning, from a meeting with Mathias Schmidt, LLNL/MFEM; confirmed by both
reviewers.)* **Plain-language form:** ghost-penalty/Nitsche stabilizes a poorly cut element by
extending the solution from a **healthy, well-conditioned neighbor of the same phase** — the remedy
presupposes a healthy neighbor exists. Here the copper runs along the **entire** edge: every
perimeter element is cut identically, the whole edge row is uniformly bad, and the inboard tape is
not a healthy *copper* neighbor. Coupling bad elements to other bad elements smooths relative
oscillations along the edge but supplies none of the missing same-phase stiffness [medium-high,
both reviewers].

**BELFEM specifics (verified):** the existing stabilizer is genuinely through-stack/normal, not
in-plane: `h_ghost()` runs on ghost facets *between stacked layers*
(`src/fem/maxwell/doc/ghost_penalty_stabilization.md:26–30`); the factory header labels
`create_ghost_facets()` "ghost facets for **normal** directions"
(`src/fem/kernel/cl_ThinShellFactory.hpp:256–262`); facets couple a lower block's top face to the
upper block's bottom face (`cl_ThinShellFactory.cpp:1839–1880`,
`mt_maxwell_h.cpp:1892–1902`); ghost facets are skipped at buffer blocks
(`cl_ThinShellFactory.cpp:1849–1850`). No side-face/in-plane path exists. Even a hypothetically
stabilized enrichment would resolve fields inside a region whose exact physics is a resistor —
wrong cost/benefit even if it succeeded. The same objection applies to p-enrichment ("more edge
DOFs") — strengthened argument below.

### Why extra edge DOFs ill-condition the system — the strengthened argument

"Ill-conditioning is intrinsic to the scale" (attempt #2) deserves a sharper form. Start with a
dichotomy: any proposed extra edge DOF on the 400 µm perimeter elements either

- **(a) cannot represent the across-wall jump at all** — plain p-enrichment: smooth polynomials on
  an L ≈ 400 µm support approximate a feature governed by t_w ≈ 20 µm at a Gibbs-limited rate, so
  the added hierarchical block raises the condition number with order p while buying essentially
  no accuracy in the quantity that matters; or
- **(b) can represent it** (Heaviside partner, sub-cell sliver DOF) — and is then subject to the
  mechanisms below.

For branch (b), five mechanisms stack, and none is just an implementation artifact:

1. **Wall-resolving DOFs are energetically near-degenerate with existing ones.** A DOF whose only
   job is to let the solution differ across a volume fraction η = t_w/L ≈ 1/20 of its element is
   O(η)-distinguishable, in the energy inner product, from the standard DOF on the same edge: the
   pair spans one O(1)-energy direction and one near-null direction of energy O(η). Simultaneously,
   resolving variation *across* the band puts gradients at 1/t_w against the bulk's 1/L — a stiff
   direction at O(1/η). The local soft/stiff spread per perimeter element is therefore of order
   η⁻² ≈ 400 before any material contrast (an order-of-magnitude statement, not a proven
   eigenvalue bound — row scaling and basis normalization shift the literal values). Cut-cell analyses sharpen this: without a
   stabilizing operator, the condition number of cut/immersed discretizations diverges
   polynomially in the cut fraction, κ ~ η^(−(2p+1))-type (Burman & Hansbo's ghost-penalty papers;
   de Prenter et al. 2017 for the finite cell method — *external, from memory* [medium on the
   exponent, high on the polynomial divergence]) — 3–4 decades for η = 1/20, worsening with
   enrichment order. Mesh refinement cannot help: t_w is physical, so η is pinned unless the mesh
   descends to the sliver scale, which is attempt #1.

2. **Material contrast multiplies the geometric factor.** The new rows carry copper-scale
   coefficients (ρ_Cu ≈ 2×10⁻⁹ Ω·m at 77 K) directly coupled to quenching-REBCO rows whose
   differential resistivity (E ∝ J^n) swings decades across a quench front *within one Newton
   step*. The badly-scaled block is not even static — it migrates through the spectrum as the
   front moves, defeating one-shot equilibration/scaling [medium-high].

3. **The digit budget is already spent.** Bathe's round-off bound, s ≥ t − log₁₀ κ (Bathe,
   §8.2.6, Eq. 8.62), with double precision t ≈ 16 and BELFEM's nonlinear-residual requirement
   ε < 10⁻¹¹ (Messe et al. 2023, paper1, §2.7; `messe2023.txt:542–556`), leaves at most ~5 decades
   of condition number for the *entire* system before the linear solves can no longer feed the
   nonlinear loop accurate increments — a heuristic digit budget, not a theorem. Mechanisms 1–2 consume that headroom by themselves. The predicted symptom is
   exactly the observed one: the nonlinear iteration stalls above tolerance (checkerboarding
   regime) rather than the factorization failing outright [medium-high].

4. **The bad subspace is global, unpivotable, and uncontrolled by any operator term.** The wall
   runs along the *entire* edge, so every perimeter element is cut identically and the near-null
   directions assemble into a connected family along the side curve — smooth relative oscillations
   between the inner and outer sheets, softest at the longest wavelengths. A direct solver can
   pivot around an isolated degenerate row; it cannot rescale away a coupled O(N_edge)-dimensional
   subspace whose Schur complement against the bulk is uniformly O(η). The only known cure is
   precisely ghost penalty — import stiffness from a healthy same-phase neighbor — and this
   geometry offers none (previous subsection). Attempt #3's "oscillations along the cut" are this
   subspace made visible [high].

5. **H(curl)-specific: ad-hoc enrichment falls out of the discrete de Rham complex.** Nédélec
   spaces work because they are paired with compatible H¹/H(div)/L² spaces and commuting
   projections — subcomplex + bounded cochain projection (Monk 2003, de Rham diagram,
   `monk.txt:2416–2429`; Arnold 2018, `arnold.txt:222–230, 4288–4302`; spurious Maxwell
   eigenmodes from incompatible spaces, `arnold.txt:667–688`). A Heaviside or sub-cell enrichment
   added to H(curl) *alone* enlarges the curl space without the matching scalar space, so the
   discrete grad→curl relation is no longer guaranteed and a cluster of spurious low-energy modes
   lands on top of mechanism 4 — a structure-preservation failure unless a compatible enriched
   complex is designed and tested (none is on offer; properly paired standard p-enrichment is not
   intrinsically guilty of this — it fails by branch (a) instead). At the shell↔air boundary the situation is worse
   still: the edge trace is statically condensed onto φ-node differences (§5.4, fact ii), so an
   independent enriched trace there either cannot be condensed (zero-diagonal saddle rows — the
   structure static condensation was chosen to avoid, Messe et al. 2023, paper1) or must be tied
   back and is redundant → singular [medium-high].

**One-sentence version:** a DOF whose only job is to represent a fraction-η feature is necessarily
O(η)-degenerate with the DOFs already present; with η ≈ 1/20 fixed by physics, every added flavor
of edge DOF enlarges an uncontrolled near-null subspace instead of resolving anything — while the
wall's entire information content is one number per unit edge length, its current. Hence: couple,
don't enrich.

---

## 3. The φ-only buffer — the crux

A layer whose material lacks a resistivity property is marked `DomainType::Buffer`
(`cl_ThinShellFactory.cpp:1730–1737`), receives the air-like DOF table
(`cl_Maxwell_FieldList.cpp:279–284`), and is dispatched to the φ kernels — grad-φ Laplace, **no
curl-curl term** (`cl_IWG_Maxwell.cpp:255–269`, `mt_maxwell_phi.cpp:24–44,101–119`).
`buffer_cut_topology.md:8–24`: the buffer "electrically disconnects the upper and lower conducting
halves … no galvanic current can cross it"; the model sits in the **"insulated coil" limit**.

Consequences (established in review round 2, both reviewers [high]):

1. A φ-only element has no `∫ρ(∇×h)²` channel — **modifying ρ there is a no-op**. More
   fundamentally, `∇×∇φ ≡ 0`: a φ region *cannot carry current by construction*.
2. Each PENTA6TS is one layer thick, so an effective ρ⊥ on a conductor layer redistributes current
   only *within* that half — it never crosses the buffer.
3. Making the buffer h-active with a dummy ρ is the only way volume-ρ could act there, and it is a
   fallback at best: actually insulating against the wall path would need ρ_b ≫ ~10⁻⁴ Ω·m, so
   practically ~1 Ω·m → 10–12 decades of contrast against quenching REBCO inside every stack — the
   H-formulation disease that h-φ was built to avoid — plus an extra DOF layer per tape [medium-high].

**Resolution:** keep the buffer φ-only. Its *topology* (two disconnected conductors) is handled by
the rerouted-cut design (`buffer_cut_topology.md`); its *physical bypass* is the edge-wall bridge
(§5). The resistivity-modification instinct survives as the **constitutive law of an interface**,
not of a volume.

---

## 4. Physics reframe: the 20 µm sliver is not the interface

The picture "current squeezes through a 20 µm contact at the REBCO edge" is subtly wrong, and
correcting it is what makes the collapsed model exact rather than approximate:

- **Injection into the jacket is distributed, not edge-localized.** During quench, current exits
  the REBCO **vertically through the broad Ag/Cu top interface** (measured
  R_ct ≈ 10⁻¹²–10⁻¹⁰ Ω·m², `contact_impedance_theory.md` §10.1), spreading over the current-transfer
  length — not through the 1.5 µm REBCO edge cross-section. The wall is then reached *through the
  top plate*. No 20 µm current-crowding field structure exists that the mesh would need to resolve
  [high on structure].
- **The wall's internal physics is trivially 1D ohmic** and both electrically and
  electromagnetically thin (skin depth in Cu at 77 K ≫ 20 µm). Pristine-wall resistance
  ≈ 10⁻⁸ Ω·m of edge — nearly a perfect short [medium on numbers].
- **The genuine unknowns are interfacial**: quality of Cu plating on the slit-cut Hastelloy edge
  (native oxide), on the REBCO edge, slitting damage. Plausible slit-edge contact terms give
  10⁻⁶–10⁻⁴ Ω·m of edge — **the interfaces dominate the wall metal by 2–6 orders** [medium]. These
  are *parameters*, not fields, and belong in an impedance coefficient.
- **A second candidate path exists: buffer pinholes.** The NI-magnet community has measured
  anomalously low *through-tape* (transverse) resistance and debates edge conduction vs. buffer
  pinhole defects as the mechanism *(external literature, from memory — medium confidence,
  verification pending, e.g. transverse-resistance measurements at NHMFL/Jun Lu's group)*. A good
  model must therefore parameterize **both paths independently** and let experiment discriminate;
  the kernel family below does this with one mechanism.

---

## 5. The model: an edge-wall impedance bridge (primary)

### 5.1 How a resistor works when the DOFs are H-based

In the h-formulation we never write `V = R·I`; a resistor enters as its **dissipation added to
the energy functional**:

```
a_R(h, v) = R · I(h) · I(v)        →  a symmetric PSD contribution to K
```

where `I(·)` is the **linear functional extracting the current through the resistor from the field
DOFs**. After collapse, the apparent need to "use the curl of H" disappears — by Ampère's law the
volumetric curl degenerates into boundary expressions that Nédélec DOFs represent *exactly*:

- **Collapsed sheet:** the in-sheet current is the tangential-H jump, `K = n̂ × [[H]]` (A/m).
- **Collapsed wire (lumped 1D resistor):** the current is the circulation
  `I = ∮_γ H·dl = Σ_e ± h_e` — a signed sum of existing edge DOFs around a loop;
  `a = R·(Σ±h_e)(Σ±v_e)` is a rank-one block. *(Kinship: a transport cut imposes a circulation —
  an ideal current source; the resistor is its finite-impedance sibling; `h_ghost` penalizes a
  jump. All are `coefficient × functional × functional`.)*

No discrete curl is evaluated, no J recovery at points, no new DOFs, no loss of order. A useful
duality: the h-formulation is *current-primal*, so resistances add like spring stiffnesses
(`R·I²`); the insulating limit `R → ∞` is a hard constraint that is better imposed by the function
space itself — which is exactly what the φ-only buffer is.

### 5.2 The two TSA terms of a collapsed sheet — and the role reversal

Collapsing a conducting layer (thickness t, resistivity ρ) splits `∫ρ|∇×H|²` into exactly two
surface terms (cf. `contact_impedance_theory.md` §3.1–3.2):

| Term | Scaling | Current it represents | Is the transfer resistor when… |
|---|---|---|---|
| `(ρ/t)·[[H_t]]·[[v_t]]` | ρ/t | current **in** the sheet plane (`K = n̂×[[H]]`) | transfer runs **along** the sheet |
| `ρt·(curl_t H_t)·(curl_t v_t)` | ρt | current **through** the sheet (`J_n`) | transfer runs **across** the sheet (`R_ct = ρt`) |

**Correction flag for `contact_impedance_theory.md` [Codex-audited, high]:** for
rint/T2TCL the transfer is *across* the layer, so the contact-resistance physics lives in the
**second** term — the one §3.2/§3.4 of that note drops as "negligible," while labeling the retained
ρ/h term "the contact-impedance term." The TSA machinery Schnaubelt 2023 §III imports keeps **both**
stiffness blocks (the printed 1/h and h/6 matrices are in Alves 2022a,
`alves2022a.txt:397–436`), which is why his T2TCL reproduces radial leakage; a §3.4-style kernel on
a genuine large-R_ct T2TCL would silently lose the radial-resistance channel that makes IC2 an
output (`schnaubelt2023.txt:337–351`). For the internal Ag/YBCO rint, treating this as practically
harmless is **conditional**, not universal [medium]: the transfer length is
`λ_ct = sqrt(R_ct/(R_sq,1+R_sq,2))`, so with metal sheet resistances (Ag ≈ 2.3×10⁻³ Ω/sq,
`todo/deferred/thinshell_hphi_formulation.md:37–43`) R_ct = 10⁻¹²–10⁻¹⁰ Ω·m² gives λ_ct from a few
µm up to hundreds of µm — often sub-mesh, but the ~2 µm figure holds only against a highly
resistive sheet (e.g. normal-state YBCO). This needs to be folded into the Milestone C plan
regardless of the side-conductivity work.

**For the edge wall, the roles reverse:** the transfer current flows *vertically within* the wall
sheet, so the **first** term is the physics — and its current functional is a pure trace jump.

### 5.3 The wall functional and the bridge kernel

A horizontal Ampère loop at buffer height, crossing the wall (thickness t_w, sheet normal = outward
in-plane direction): the enclosed vertical wall current per unit tape length is

```
I′(y) = [[H_∥]] = H_∥(inside, stack edge) − H_∥(outside air)          [A/m]
```

with `∥` the direction along the edge curve. The bridge is then a **curve integral on the tape's
lateral edge**:

```
a_wall(h, v) = ∫_Γedge  r′(y) · [[h_∥]] · [[v_∥]]  dy
```

It is assembled into `K()` as Kmm/Kms/Ksm/Kss blocks over curve elements, structurally matching the
contact-impedance kernel of `contact_impedance_theory.md` §7 with dS → dy and master/slave = the
traces on either side of the wall. **Zero connector-owned DOFs.** Joule heating for the thermal
problem is directly `P′(y) = r′·|I′|²` per unit edge length. *(Sign convention: the Ampère-loop
orientation fixes the sign of `I′`; the bilinear form squares the functional, so orientation only
affects the signs of the off-diagonal assembly blocks — Codex-verified.)*

**The coefficient** is the condensed series resistance of the wall path (units Ω·m):

```
r′ = ρ_Cu·h_path/t_w  +  ( R_ct,top-edge + R_ct,bottom-edge ) / t_w
```

- Pristine metal floor: ρ_Cu(77 K)≈2×10⁻⁹ Ω·m, h_path≈100 µm, t_w≈20 µm → ≈10⁻⁸ Ω·m.
- Slit-edge interface terms R_ct/t_w: 10⁻⁶–10⁻⁴ Ω·m plausible → **interface-controlled**; r′ is the
  calibration dial, and the collapsed representation loses nothing that was ever known [medium].

**Limits behave correctly:** `r′ → ∞` forces `I′ → 0` (insulated edge — today's model recovered);
`r′ → 0` frees the jump (perfect short). Positive-semidefinite for r′ > 0; no tunable penalty — the
coefficient *is* the physics.

**Why the HEX8TS failure mode cannot recur:** the wrap demanded the wall's *field state* from a basis
and left binormal-H unwetted; the bridge needs only the wall's *current*, which Ampère determines
from existing traces. The constitutive law supplies exactly the information the wrap's basis could
not [high].

**Conditioning:** no new scale, no new DOFs; the coefficient sits at the material scale of the
copper it connects (unlike rint's ~5×10⁵ contrast against Ag). The r′ → 0 limit is the only
near-singular direction; physical interfaces keep r′ finite [high].

#### 5.3.1 Stiffness term vs. boundary term — two derivations, one kernel

*(Q&A with Christian, 2026-07-06: implement the wall (a) as an added stiffness
`∫(∇×δh)'ρ(∇×h) dV`, or (b) through the natural boundary term `∮ δh'·(n×e) dS` with Ohm's law
`e = ρ∇×h`?)* **Answer: they are the same term** — the boundary integral is the variational
residue of the volume stiffness after integration by parts. Route (b) shows *where* the wall
enters; route (a) shows *how to evaluate* it.

- **(b), done correctly:** assemble both sides' natural terms; for a thin purely resistive sheet,
  tangential e is continuous (`[[e_t]] = 0`) and tangential h jumps by the sheet current
  (`K = n×[[h]]`), so the pair collapses to `∮ [[δh_t]]'·(n×e_t) dS`. Apply the **sheet's own**
  Ohm law `e_t = ρ_s (n×[[h_t]])`, `ρ_s = ρ_w/t_w` → `a_wall = ∮ ρ_s [[δh_t]]'[[h_t]] dS`.
- **(a):** in the collapsed wall the curl degenerates: `C_wall → (E⁻ − E⁺)/t_w` (an in-plane n×
  rotation drops out of the quadratic form), so `∫C'ρC dV → ∫ ρ_s (E⁻−E⁺)'(E⁻−E⁺) dS` —
  identical, manifestly symmetric PSD, the Kmm/Kms/Ksm/Kss structure. No discrete curl is ever
  evaluated; the wall's mass term collapses to O(μ·t_w), negligible.
- **Trap [important]:** implementing (b) literally as `e = ρ·C·h` with the **adjacent tape
  elements'** C gives the *tape's* current, not the wall's — the wall current is invisible to
  neighboring C matrices (it lives entirely in the jump). Neighbor-curl evaluation is legitimate
  only for Nitsche *consistency* terms (continuity enforcement, as in `h_ghost`); for the physical
  Robin/impedance condition the pure jump–jump form is the kernel (O3 gate for consistency terms).
- **Projection:** penalize only the along-curve component `t̂·[[h_t]]` (= transfer current K_z).
  The other tangential jump component (longitudinal wall current K_y / through-stack channel) is
  physically negligible (~40× smaller cross-section than the parallel Cu plates) and is exactly
  the binormal channel of the HEX8TS failure — add no physics there.
- **Analogy:** the kernel is the electromagnetic **cohesive-zone / zero-thickness interface
  element** — `[[h_t]]` ↔ displacement jump, `n×e` ↔ traction, sheet Ohm ↔ traction–separation
  law; local stencil `r′·[+1 −1; −1 +1]` (an R, not a G, because the h-formulation is
  current-primal). The buffer needs no term by the same logic: no sheet → `e_t` continuous and
  K = 0 → the two boundary residues cancel identically; φ-only is the `r′ → ∞` limit imposed by
  the function space instead of by penalty.

#### 5.3.2 Discretization: building K from `a_R = R·I(h)·I(v)`

Both flavors of the resistor follow one recipe: express the current functional `I(·)` as a row
vector acting on element DOFs; the stiffness contribution is then `coefficient × rowᵀ·row`.

**(i) Lumped resistor (0D, circulation form).** With `q` the element/loop DOF vector and `s` the
signed incidence vector of the Ampère loop γ (`s_e ∈ {−1, 0, +1}`),

```
I(h) = ∮_γ H·dl = Σ_e s_e h_e = sᵀ q      →      K += R · s sᵀ
```

— a rank-one, symmetric PSD update; the residual contribution is `R·s·(sᵀq)`. This is the
transport cut's symmetric sibling: the cut imposes `sᵀq = I` as a hanging-node/T-matrix
*constraint*; the resistor assembles `R·s sᵀ` into `K` as *physics*.

**(ii) The distributed edge-wall bridge (what will be implemented).** The side curve is meshed by
curve elements (LINE2/LINE3 along the tape edge). On each, let `q_m`/`q_s` be the DOFs whose
traces form the two sides of the wall (their identity is O1's output, below). At integration
point k:

- `Em(k)`, `Es(k)` — 1×n row matrices: the **t̂-projected** (along-curve component only, per the
  projection rule of §5.3.1) tangential traces of the master/slave edge functions — the exact
  analogue of `Em`/`Es` in `contact_impedance_theory.md` §7;
- jump row `B(k) = [ Em(k), −Es(k) ]`, so the wall current per unit length at the point is
  `I′(k) = B(k)·[q_m; q_s]`.

```cpp
// tW( k ): quadrature weight;  tDetJ: curve Jacobian (dy);  tRp: r′( y_k )
Kmm +=  tRp * trans( Em ) * Em * tW( k ) * tDetJ ;   // master-master
Kms -=  tRp * trans( Em ) * Es * tW( k ) * tDetJ ;   // master-slave
Ksm -=  tRp * trans( Es ) * Em * tW( k ) * tDetJ ;   // slave-master
Kss +=  tRp * trans( Es ) * Es * tW( k ) * tDetJ ;   // slave-slave
```

equivalently `K_e = Σ_k w_k · detJ · r′ · B(k)ᵀ B(k)`: the quadrature rule turns the curve
integral into a **chain of rank-one resistors** — at each Gauss point, "R" is `r′·w_k·detJ` and
"I(·)" is `B(k)·(·)`. Structure and sign pattern are identical to the `h_contact_impedance`
kernel of `contact_impedance_theory.md` §7 with dS → dy.

Bookkeeping notes:

- **Goes into `K()`, not `M()`** — resistive dissipation, same slot as the ρ curl-curl blocks;
  under the time discretization it multiplies present-time h like any other resistive term. No
  mass contribution (the wall's inductance collapses to O(μ₀·t_w) — §5.3.1).
- **No h-Newton term.** `I′` is linear in h and r′ is independent of the magnetic unknowns, so
  the magnetic Jacobian block equals the stiffness block. If O4 later makes `r′ = r′(T)` inside a
  monolithic thermo-magnetic solve, `dr′/dT·I′(h)·I′(v)` belongs to the thermal cross-Jacobian,
  not the magnetic stiffness. The block is PSD, so it adds no indefinite mode; its conditioning
  effect is set by the physical magnitude of r′.
- **DOF identity is O1's output; the algebra is invariant once O1 provides a linear current
  row.** If O1 resolves to (a) distinct
  conductor edge traces, `q` are `edge_h` DOFs and `Em`/`Es` are edge-function traces. If it
  resolves to (b) the φ-jump curve form, the boundary trace is already condensed onto φ-node
  differences (§5.4, fact ii), so `B·q` becomes a signed φ-node stencil delivered by the existing
  T-matrix path (`cl_FEM_DofMgr_DofData.cpp:3555–3637`) — the `BᵀB` pattern is untouched, only
  the DOF map changes, exactly as for every thin-shell kernel that assembles through hanging DOFs.
  This is algebraic invariance (`TᵀBᵀBT`) only — it does not supply the duplicated-trace
  topology; O1 stays blocking.
- **Joule heat falls out at the same points:** `P′(y_k) = r′·(B(k)q)²` → thermal source (O4).
- **Orientation:** the ± signs in `B` follow the Ampère-loop orientation from the binormal frame
  (O5); flipping it maps `B → −B` and leaves `BᵀB` unchanged (§5.3's sign-convention note).

### 5.4 The key discrete design questions (for the design note)

The bridge only functions if the side-curve bookkeeping **allows the inside and outside traces to
differ**: wherever the wall exists, the tangential-H continuity condition at the shell edge ↔ air
interface must become this impedance (Robin-with-jump) condition. If both half traces are hard-tied
to one air trace, the jump is identically zero and the kernel is inert.

**Codex audit confirmed the trap is real and current BELFEM makes a naive side-curve jump kernel
inert or unrepresentable [high]:** (i) closed side curves are *deliberately not node-split* —
`duplicate_nodes_on_face_sidesets()` removes closed-side-curve nodes from the thin-shell duplicate
set (`cl_CutFactory.cpp:1546–1598`); (ii) air-adjacent shell boundary edges are **hung on air φ
node differences** — the boundary `edge_h` trace is not an independent H DOF but literally
`φ(node0) − φ(node1)` via the T-matrix path (`cl_MaxwellFactory.cpp:1319–1324, 1432–1477,
1492–1537`; `cl_FEM_DofMgr_DofData.cpp:3555–3637`); (iii) the conductor-air interface sideset
kernel is disabled outright ("Conductor-Air Interfaces must be disabled!",
`cl_IWG_Maxwell.cpp:398–401`) — no existing sideset weak form in which an outside H trace could live.
Cut jumps are φ-domain node objects (`cl_CutFactory.cpp:680–700`, `cl_InterfaceProcessor.cpp`),
not side-wall H traces.

- [ ] **O1 — trace bookkeeping (BLOCKING).** Introduce an explicit side-wall trace pair before any
      kernel is written. Candidates: (a) distinct conductor-layer edge traces below/above the
      buffer at the curve (rint-§7-style distinct-trace requirement — but note fact (ii): today
      those traces resolve to the *same* air φ nodes); versus (b) a carefully derived **φ-jump curve
      form** — since the boundary edges already collapse to φ-node differences, the wall current
      may be expressible as a jump between *duplicated* φ nodes along the curve
      (`I′ ∝ ∂[[φ]]/∂y`), using the duplicated-node / cut-termination machinery — which requires
      relaxing fact (i) (the closed-side-curve non-duplication rule) exactly where the wall exists.
      This is the design note's first task.
- [ ] **O2 — cut/coefficient bookkeeping.** The bridge is weak-form only — invisible to cohomology;
      the halves stay topologically separate. Verify the single rerouted cut with a free substrate
      coefficient (`buffer_cut_topology.md` §3.4, §119–126) yields the solver-determined partition,
      or whether Schnaubelt's explicit second free coefficient (IC2, `schnaubelt2023.txt:337–351`)
      is required. Testable on the 1D benchmark (§8). §5.6 spells out the physics the free
      coefficient encodes (zero-impedance terminal closure) and the floating-half caveat, which
      must become an input option, not an assumption.
- [ ] **O3 — consistency terms.** Same open question as `contact_impedance_theory.md` §3.4
      (Juntunen & Stenberg: Robin conditions are natural; unproven for this H(curl) trace pair).
      The convergence study in §8 doubles as the empirical check.
- [ ] **O4 — thermal feedback.** Route `P′ = r′|I′|²` into the thermal problem (an isothermal
      bridge understates quench propagation; broken thermal coupling was an original HEX8TS defect).
- [ ] **O5 — orientation.** Side-curve orientation/binormal frame from the retained
      `compute_binomial_vectors` (+ its side-connector-sign assert — kept in June precisely for
      this).

### 5.5 One kernel family, four physical objects

| Object | Geometry | Transfer direction | Governing term | Coefficient |
|---|---|---|---|---|
| **Edge wall** (this work) | side curve | ∥ sheet | `r′·[[H_∥]]²` on the curve | r′ ≈ 10⁻⁸…10⁻⁴ Ω·m, interface-controlled |
| **Buffer pinholes** | buffer midplane | ⊥ sheet | `ρt·(curl_t H_t)²` coupling the **adjacent conductor traces** (Hastelloy↔YBCO) — NOT the buffer φ DOFs (`∇×∇φ ≡ 0`) | R_ct,buffer large but finite |
| **Internal rint** (Ag/YBCO) | layer interface | ⊥ sheet | ditto (Milestone C; note §5.2 correction) | R_ct ≈ 10⁻¹²–10⁻¹⁰ Ω·m² |
| **CORC tape-to-tape contact** | outer-surface contact strips | ⊥ sheet | Schnaubelt T2TCL | R_ct ≈ 10⁻⁹–10⁻⁷ Ω·m² |

In a CORC winding the tape-to-tape path is REBCO → top Cu → **edge wall** → back Cu → contact spot
→ next tape: the wall is in series in the inter-tape chain. With edge impedance and buffer-face
impedance as independent dials, the simulation becomes the instrument for discriminating the two
experimentally-debated mechanisms (§4) instead of presuming one.

### 5.6 The two halves: one conductor or two? Current partition without and with the bridge

*(Q&A with Christian, 2026-07-09.)* First a correction of framing: the stack is one **cut** but
not one **conductor**. The rerouted single cut (`buffer_cut_topology.md` §3.2) imposes the
coefficient I around the **HTS-side loop only**; the substrate-side loop's net current is a free,
solver-determined circulation — Schnaubelt's IC1 (imposed) / IC2 (free) pattern folded into one
cut surface (`buffer_cut_topology.md` §3.5).

**(a) Both halves isolated (today's model, r′ → ∞).** No current crosses the buffer anywhere
along the modeled length — by construction (`∇×∇φ ≡ 0`). What then sets the substrate current?
The free coefficient's Galerkin equation is the zero-applied-EMF Faraday equation for the
substrate loop: resistive drop plus dΦ/dt around it must vanish. For the standard soldered
tape-end setup this is equivalent to closing the substrate half **through an ideal external
terminal short outside the modeled domain** — the physically right closure there, since the
termination contacts every layer, but a modeling *choice* rather than an inherent property of a
free coefficient (a finite external circuit would need a circuit equation instead).
Consequences [medium-high]:

- **Steady superconducting state:** parallel-path partition by end-to-end branch resistance,
  I_sub/I_HTS ≈ R_HTS-branch/R_sub-branch ~ 10⁻⁶ (`buffer_cut_topology.md` §3.4) — essentially
  all current in the HTS half.
- **During ramps:** the substrate half is a shorted transformer secondary — an induced net
  current decaying with τ = L_loop/R_sub — plus ordinary eddy/shielding redistribution among the
  galvanically connected layers *within* each half.
- **During a quench — where the isolated model is wrong:** current leaving a quenched REBCO
  segment can only redistribute within the lower half (Ag/Cu below) or reroute through the far
  terminations, i.e. by changing the *global* partition. There is no local escape to the upper
  copper: transfer that physically happens over λ ≈ mm–cm (§8, item 1) happens in the model over
  the entire tape length or not at all.
- **Caveat for O2:** if a configuration truly has a floating half (e.g. a coupon whose
  termination contacts only the HTS side), the closure changes: impose I_sub = 0 for the net
  transport loop (local zero-net eddy currents remain free). The free coefficient encodes a
  terminal model — it must become an input choice, not an assumption.

**(b) With the bridge — do they act as one common conductor?** Physically yes; discretely no, and
deliberately so. Finite r′ makes the tape one connected conductor (the U-jacket is real), but the
discrete model keeps **two topological conductors coupled by a weak-form impedance** — the bridge
adds no cohomology generator and touches no cut (O2). Three reasons this is the correct
structure, not a leftover [high]:

1. The buffer stays φ-only, so in the h-space the halves *are* disconnected. A literal
   one-cut-around-both model would constrain only the **sum** of the half currents and leave the
   partition mode either unrepresentable or physics-free — the very failure
   `buffer_cut_topology.md` §1 records for the pre-rerouting cut.
2. The partition — how much current has crossed the wall, and where — is the scientific output of
   the whole exercise. It must remain a solution unknown facing an impedance: exactly
   Schnaubelt's free IC2 + T2TCL resistance, with our r′ in the T2TCL role.
3. The limits close the loop: r′ → ∞ recovers (a); r′ → 0 makes the wall a near-perfect short and
   the partition follows the metals' parallel conductances over the transfer length λ — the tape
   then *behaves* as the single conductor of the pre-buffer model. "One common conductor" is the
   **limit behavior** of the coupled model, not a different topology.

With the bridge, the partition stops being set by the global terminal closure alone and becomes
the distributed transmission-line solution of §8 item 1: transfer concentrated within
λ = sqrt(r′/(R′₁+R′₂)) of a disturbance or termination — which is precisely the physics the
isolated model cannot express.

---

## 6. What the collapsed model gives up (from review round 1)

The collapsed model gives up corner current-crowding detail inside the wall, the exact 20 µm injection
footprint, copper-side longitudinal voltage gradients, local wall eddy detail, sub-element transverse
voltage variation.
Acceptable for global I_c, current-sharing trends, and total AC loss; the local hot spot in a
quench sits in the REBCO that *cannot* shed current — which the model captures — not inside the
copper wall [medium]. Not acceptable if wall-internal loss peaks ever become the quantity of
interest (unlikely).

---

## 7. Optional refinement: effective anisotropic ρ⊥ on h-active perimeter layers

After review, the original proposal is correctly scoped as follows: on h-active perimeter elements
only, a tensor resistivity lowering the through-thickness channel captures within-layer corner
shorting:

```
K += Cᵀ R C · (w·dV),   R = ρ∥·I + (ρ⊥,eff − ρ∥)·n nᵀ,
σ⊥,eff = σ⊥,material + σ_Cu·(A_Cu/A_elem)·(L_elem/L_Cu)     (parallel addition, real σ·A/L geometry)
```

Never a scalar ρ on `Cᵀ C` (it would soften in-plane modes); requires a matching nonlinear Newton
term (currents summed at common E_n); `n` from PENTA6TS (`cl_EF_PENTA6TS.cpp:69–83`), sign-free
through `n nᵀ`. The current kernels are scalar-ρ (`mt_maxwell_h.cpp:83,156,…,839–869`; `element_rho` at
`cl_IWG_Maxwell.cpp:77`; scalar Material ρ at `cl_Material.hpp:133`). **Demoted to optional:** the
layers within each half are already galvanically coupled through the stack, so this adds only
sub-element corner detail on top of the bridge. Implement only if validation (§8) shows the bridge
alone misses measured behavior [medium].

---

## 8. Validation plan

1. - [ ] **1D analytic benchmark (first).** Two parallel sheet conductors coupled by a distributed
     edge conductance: transmission-line solution with current-transfer length
     `λ = sqrt( r′ / (R′₁ + R′₂) )`, `R′ᵢ = ρᵢ/(tᵢ w)`. Back-of-envelope (Cu plates
     20 µm × 4 mm, R′₁+R′₂ ≈ 5×10⁻² Ω/m): pristine wall λ ≈ 0.5 mm; r′ = 5×10⁻⁶ → λ ≈ 1 cm;
     r′ = 5×10⁻⁴ → λ ≈ 10 cm [medium]. This validates the kernel, cut bookkeeping (O2), and
     convergence order (O3) at once — and λ *is* the physical answer to "can the current escape,
     and over what length?"
2. - [ ] **Conditioning:** with/without bridge on a representative mesh; the estimate stays near the
     no-connector baseline; no in-plane pollution.
3. - [ ] **Physics sanity:** current crosses between halves around the buffer (forbidden in the
     insulated-coil limit); side-curve `J/J_c` not spuriously suppressed (the HEX8TS failure).
4. - [ ] **Partition sweep:** upper/lower current partition moves smoothly from insulated
     (`r′ → ∞`) toward a short as r′ decreases.
5. - [ ] **Calibration:** anchor r′ (and the pinhole R_ct, if enabled) against measured transverse
     tape resistance / one fine-mesh resolved case; literature verification of the
     transverse-resistance data (§4) pending.

---

## 9. Review record

Full thread: `tmp/ai_exchange/side_connector_effective_resistivity.md`.

- **Round 1 (Codex + Grok, 2026-07-06):** effective-ρ⊥ is sound as a lumped shunt iff parallel-addition
  σ·A/L (not naive 20/400), tensor `Cᵀ R C` (never scalar), calibrated coefficients.
- **Round 2 (Codex + Grok, 2026-07-06):** no-healthy-neighbor argument confirmed; ghost machinery
  verified normal-only; φ-only buffer makes ρ⊥ a no-op; one-layer ρ⊥ cannot bridge the buffer; both
  independently recommend the lateral collapsed-impedance bridge; cut interaction flagged (IC2).
- **Round 3 (Fable session, 2026-07-06):** physics reframe (§4: distributed injection, interface-
  controlled r′, pinhole path), the resistor-in-h formulation (§5.1–5.3), the TSA role-reversal and
  the `contact_impedance_theory.md` §3.4 correction flag (§5.2), the inert-bridge trap and trace
  bookkeeping question (§5.4), λ estimates (§8). **Codex audit (same day): all four claims
  confirmed** — TSA split [high, with Alves 2022a matrices], §3.4 correction [high] with the
  transfer-length harmlessness downgraded to conditional [medium], wall functional/units/limits
  [high, sign-convention note added], inert-bridge trap [high] with the three implementation facts
  now cited in §5.4; λ formula and 0.5 mm–10 cm range verified. O1 elevated to BLOCKING.
  Devlog: `devlog/dl20260706_side_connector_bridge_formulation.md`.
- **Round 4 (Fable session, 2026-07-09, Christian's three questions):** (1) conditioning argument
  for extra edge DOFs hardened into the dichotomy + five-mechanism form (§2) — anchored on Bathe
  §8.2.6 Eq. 8.62 vs the ε < 10⁻¹¹ requirement; cut-FEM κ-scaling exponent flagged external/from
  memory. (2) Half-partition physics spelled out (§5.6): free cut coefficient = zero-impedance
  terminal closure; floating-half caveat added to O2; "one common conductor" identified as the
  r′ → 0 limit behavior, not a topology change. (3) Discretization subsection added (§5.3.2):
  `K_e = Σ w·detJ·r′·BᵀB`, `B = [Em, −Es]`, mirroring `contact_impedance_theory.md` §7 with
  dS → dy. **Codex audit (same day): all three sections confirmed** — Bathe Eq. 8.62 verified
  (`bathe.txt:27603–27611`) with the digit budget kept heuristic; spectral spread reworded to
  "order η⁻²"; cut-FEM exponent stays external [medium]; mechanism 5 re-scoped to *ad-hoc
  unpaired* enrichment (Monk/Arnold line citations added — compatible p-enrichment fails by
  branch (a), not (b)); §5.3.2 block structure/K-not-M/T-matrix collapse verified against
  `contact_impedance_theory.md:163–186` and `cl_MaxwellFactory.cpp:1319–1324,1463–1477,1524–1537`,
  with the `dr′/dT` thermal cross-Jacobian caveat and PSD phrasing fixed; §5.6 cut topology
  verified [high], free-coefficient = zero-applied-EMF loop equation [medium-high] reframed as a
  terminal-model *choice* (soldered ends ⇒ ideal short; floating half ⇒ impose I_sub = 0);
  floating-half caveat confirmed correct and important.

Rejected shortcuts, for the record: using the cohomology cut (adds no conductance — wrong
mechanism [high]); h-activating the buffer edge column (fake conductor volume, contrast disease —
fallback only [medium]); p-enrichment and XFEM+Nitsche (§2).

---

## 10. Key references

**Code (verified 2026-07-06):** curl-curl kernels `mt_maxwell_h.cpp:83,…,839–869`; ghost machinery
`ghost_penalty_stabilization.md:26–30`, `cl_ThinShellFactory.hpp:256–262`,
`cl_ThinShellFactory.cpp:1839–1880,1849–1850`; buffer classification/dispatch
`cl_ThinShellFactory.cpp:1730–1737`, `cl_Maxwell_FieldList.cpp:279–284`,
`cl_IWG_Maxwell.cpp:255–269`, `mt_maxwell_phi.cpp:24–44,101–119`; buffer topology
`buffer_cut_topology.md:8–24,65–70,109–126`; collapsed-impedance pattern + measured R_ct
`contact_impedance_theory.md` §3, §7, §9, §10; PENTA6TS frame `cl_EF_PENTA6TS.cpp:69–83,132–138,312–325`.

**Literature (`literature/papers/fem/`):** Schnaubelt 2023 — TSA §III (full matrices), parallel
effective resistivity (`schnaubelt2023.txt:353–372`), free cut coefficient (`:337–351`) — primary;
Alves 2024/2022a (`alves2024.txt:65–77`, `alves2022a.txt:177–186`); Messe 2023
(`messe2023.txt:143–147,230–235,314–320,490–497`). Measured contact resistances:
`contact_impedance_theory.md` §10 (Hayasaka & Ito 2019; Fleiter 2017; et al.). External,
unverified [medium]: NHMFL transverse-resistance / buffer-pinhole measurements — to be confirmed.

**History:** `devlog/dl20260605_remove_side_connectors.md`, `dl20260620_side_connector_removal.md`,
`dl20260417_hex8ts_thirdpass_audit.md`, `dl20260430_side_connector_plumbing_audit.md`,
`dl20260706_side_connector_bridge_formulation.md`.
