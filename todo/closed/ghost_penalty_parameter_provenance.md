# Ghost-Penalty Parameter Provenance: `eta = 4.0` and `k_reg = 1e-3`

**Date:** 2026-06-20
**Purpose:** Trace where the two hardcoded ghost-stabilization constants in `maxwell::h_ghost()` come from, whether they are literature-grounded, and how to choose better-motivated values.
**Module:** fem/maxwell
**Status:** Research note — Claude literature trace + web search, independently corroborated by Grok and Codex (both confirm: not literature-derived). No code change proposed yet.

---

## The two constants

`cl_IWG_Maxwell.cpp:50` — `mPenalty = { 4.0 , 1e-3 };`, so:

| Symbol | Source | Role |
|--------|--------|------|
| `eta`   | `penalty(0)` = **4.0** | dimensionless Nitsche / interior-penalty stabilization constant |
| `k_reg` | `penalty(1)` = **1e-3 Ω** | regularization offset added to each side's stiffness `k = ρ/h` before the harmonic mean |

The penalty coefficient (`mt_maxwell_h.cpp:1805–1847`) is
`alpha = eta · k_pen`, `k_pen = 2·km_reg·ks_reg/(km_reg+ks_reg)`, `k = ρ/h`, `k_reg = k + k_reg`.
The harmonic-mean weighting is attributed in-code to **Burman & Zunino 2006**, and `h_ghost()` is a **nonsymmetric** Nitsche scheme (penalty + consistency terms; signs verified against Ern & Guermond 2009 per `devlog/dl20260319_ghost_thinshell.md`).

## Bottom line

**Neither value is from a paper.** Both are ad-hoc engineering choices (almost certainly training-data-informed picks, as you suspected). The literature-grounded pieces are the *method* (nonsymmetric interior penalty) and the *harmonic-mean weighting* (Burman & Zunino) — not the specific multiplier or the regularization offset. The useful part: the literature gives clean *frameworks* for calibrating both.

## `eta = 4.0` — provenance and how to pick it

**Method.** Nonsymmetric interior penalty / Nitsche: Rivière–Wheeler–Girault (2001); the `eta = 0` limit is Baumann–Oden (1998). Brenner & Scott §10.5 (curated, `brenner.txt` ~13695–13925) is the cleanest reference.

**Key theory (decisive for our nonsymmetric scheme):**
- **Coercivity for *any* `eta > 0`** — Brenner & Scott (10.5.13): the nonsymmetric form `a⁺_h(v,v) = |||v|||²_h`. Unlike the *symmetric* method (which needs `eta ≥ eta*`, a trace/inverse-inequality threshold), the nonsymmetric scheme is stable for any positive penalty. **`4.0` is therefore not a stability threshold** — any positive value is stable.
- **Accuracy ↔ conditioning trade-off** — Brenner & Scott (10.5.18): the error constant scales as **`(1 + eta⁻¹)²`**. Evaluated: `eta=1 → 4×`, **`eta=4 → 1.56×`**, `eta=10 → 1.21×`, `eta=16 → 1.13×`, `eta→∞ → 1×` (but conditioning degrades with `eta`). So `eta=4` sits at a reasonable knee: most of the gain over `eta=1`, with room to tighten.

**Standard calibration (web + Codex/Grok).** The conventional Nitsche choice is **`C_pen ≈ 2·C_tr`**, where `C_tr` is the discrete trace-inequality constant for the actual element (Ern & Guermond). For the *symmetric* method, Brenner & Scott give the threshold form **`eta* = ½ + 2·C*`** (with `C*` the discrete trace estimate; brenner.txt:13904). **Shahbazi (2005)** gives an *explicit* penalty formula from element geometry for simplices (a coercive value without guessing). Concrete P1 trace constants (Warburton–Hesthaven / Shahbazi type): triangle `‖v‖²_F ≤ 3·|F|/|K|·‖v‖²_K ≈ 6/h⊥`, tetra `≈ 8/h⊥`. So `C_pen` lands roughly in the O(2)–O(20) range — `4.0` is plausibly in range, but **uncalibrated for the actual PENTA6TS/TRI3 facets**.

**Typical values in practice:** NIPG ~1–10 (any positive is coercive, but small values give poor jump control / conditioning); SIPG ~5–20, with **10 a common conservative default**. So `4.0` is low-to-moderate for our nonsymmetric scheme.

**How to choose better:**
- Cheapest: bump `eta` toward **~10–16** (error constant 1.21–1.13×) and check that conditioning / MUMPS workspace doesn't blow up — directly supported by the `(1+eta⁻¹)²` curve.
- Principled: compute `C_pen ≈ 2·C_tr` from the trace constant of the actual facet element (Shahbazi 2005 explicit formula), so it is mesh-aware rather than a magic number.
- Validate with a jump/transmission patch test (manufactured solution across a layer interface) and a contrast sweep.

## `k_reg = 1e-3` — provenance and how to pick it

**Harmonic-mean weighting** is from **Burman & Zunino (2006)** (weighted interior penalty for discontinuous diffusion; harmonic/arithmetic/geometric weights make the method *robust* w.r.t. the coefficient contrast). **But** that theory assumes coefficients bounded away from 0 and ∞ — it does **not** cover the superconductor limit (`ρ → 0`, `k → 0`), where the plain harmonic mean collapses to zero. The **regularization offset `k_reg`** that floors this collapse is **BELFEM-specific — not in Burman & Zunino**.

**Dular 2021** (curated, `dular2021.txt:1109–1111`) regularizes the power-law **resistivity** with "two limiting resistivity values" (a floor and a ceiling; refs [30],[31] for rigorous treatment) to keep the operator bounded for the inf-sup analysis. This uses the same *principle* as `k_reg` (floor the coefficient to avoid the `ρ→0` degeneracy), but it floors `ρ`, not the interface stiffness `k = ρ/h`, and gives **no numeric value**.

**The code's own rationale** (`mt_maxwell_h.cpp:1830–1837`): `k_reg` should be "a small fraction of the typical conductor stiffness... `k_conductor ~ 1e-3..1e-2 Ω`, so `k_reg = 1e-3` sits at the low end." That is internal engineering justification, not a cited value. (Codex also caught a **math error in that comment** (`:1828`): the insulator limit of the regularized harmonic mean is `2·(other_k + k_reg)`, **not** `2·min(other_k, k_reg)` — the analysis in this note uses the correct `≈ 2·other_k` limit; the comment text should be fixed.)

**Tension with the material floor (Codex):** BELFEM's HTS material model already floors resistivity at `mRhoMin = 1e-16` (`cl_Material.hpp:275`). `k_reg = 1e-3 Ω` is a *stiffness* floor (`ρ/h`) far above that — so the ghost penalty's `k_reg`, not the material floor, sets the superconducting-side stiffness. That is a deliberate (but undocumented) conditioning choice, and it is why `k_reg` directly governs stabilization strength at SC↔normal fronts.

**Why it matters at t69 (current-sharing/quench):** at a superconductor↔normal interface one side has `k → 0` (floored at `k_reg`), so `k_pen ≈ 2·k_reg` and `alpha ≈ eta·2·k_reg` *regardless of the resistive side* — i.e. the effective penalty is set almost entirely by `eta·k_reg` there. A `k_reg` chosen at the "low end" makes the stabilization weakest exactly at the moving SC/normal fronts.

**How to choose better:**
- Tie `k_reg` to **physics** (Codex + Grok agree): either `k_reg = ρ_floor / h_sc` reusing the material model's own floor policy, or `k_reg = c · median(ρ_metal / h_metal)` with `c ≈ 0.01–0.1`, documented as a conditioning knob. `ρ_floor` examples: `E_c/J_c` of the HTS, or a fraction of normal-state metal `ρ` (the contact-impedance note §10 puts the *typical* Ag/YBCO contact at `ρ ≈ 1e-4 Ω·m`).
- Or raise the current value toward the upper end of the cited range (`1e-2`) so the SC-side floor is less aggressive — the targeted knob for the t69 SC↔normal collapse.
- Either way: document the assumption and add a contrast-sweep test (`ρ_m/ρ_s → 0` and `→ ∞`).

## Recommended actions

- [ ] **Sweep `eta`** at t69 (4 → 10 → 16) and watch convergence + MUMPS workspace; the `(1+eta⁻¹)²` curve predicts diminishing accuracy gains and rising conditioning cost.
- [ ] **Sweep `k_reg`** (1e-3 → 1e-2) — the physically-targeted knob for the SC↔normal collapse.
- [ ] **(Principled `eta`)** compute `C_pen ≈ 2·C_tr` from the trace constant of the actual facet element, or implement Shahbazi (2005)'s explicit formula, so `eta` is mesh-aware.
- [ ] **(Principled `k_reg`)** replace the "fraction-of-stiffness" guess with `ρ_floor / h` tied to material physics (`E_c/J_c` or metal `ρ`).
- [ ] **Add a regression test**: manufactured transmission solution across a layer interface, swept over coefficient contrast (incl. `ρ → 0`), checking stability and convergence order — the calibration evidence that's currently missing.
- [ ] **Tidy:** clean the stray `; ;` at `mt_maxwell_h.cpp:1810`; fix the math error in the `mt_maxwell_h.cpp:1828` comment (insulator limit is `2·(other_k + k_reg)`, not `2·min(...)`); add a one-line citation/derivation comment for `eta` and `k_reg`; fix the "SIPG" wording in `contact_impedance_theory.md` (the scheme is nonsymmetric, not symmetric).

## Citations

- **Brenner & Scott**, *The Mathematical Theory of Finite Element Methods* (curated `literature/books/brenner.txt`), §10.5 — interior penalty methods; (10.5.13) nonsymmetric coercivity for any `eta>0`; (10.5.18) error `~(1+eta⁻¹)²`; (10.5.19) symmetric threshold `eta*`.
- **Rivière, Wheeler & Girault (2001)** — nonsymmetric IPG; **Baumann & Oden (1998)** — `eta=0` variant; **Arnold (1982)/Wheeler (1978)** — symmetric IPG.
- **Shahbazi (2005)**, "An explicit expression for the penalty parameter of the interior penalty method," *J. Comput. Phys.* 205(2):401–407 — explicit, geometry-based penalty for simplices.
- **Burman & Zunino (2006)**, "A domain decomposition method based on weighted interior penalties...," *SIAM J. Numer. Anal.* 44:1612–1638 — weighted (harmonic) interior penalty, robust to coefficient contrast (assumes coefficients bounded away from 0/∞).
- **Ern & Guermond**, Nitsche boundary-penalty notes — `C_pen ≈ 2·C_tr` rule of thumb; `C_tr` = discrete trace-inequality constant.
- **Dular et al. (2021)** (curated `literature/papers/fem/dular2021.txt:1109`) — regularized power-law resistivity with two limiting values (principle behind `k_reg`; no numeric value); hierarchical enrichment for interface spurious oscillations.
- **Schnaubelt et al. (2023)**, **Alves et al. (2022b/2024)** — HTS thin-shell / contact-impedance interface treatments; none parameterize a ghost-penalty `eta`/`k_reg` (searched, no mapping).

## Confidence and open items

- **High** that `4.0` and `1e-3` are not literature-derived (exhaustive search of `./literature` + web + Grok found no source value/formula).
- **High** that the nonsymmetric scheme is coercive for any `eta>0` and the error scales as `(1+eta⁻¹)²` (Brenner & Scott, direct).
- **Medium** on the *direction* to move them (depends on whether the t69 residual concentrates on a conductor ghost interface vs the buffer φ interface — still worth localizing).
- Web sources (Sources): [Shahbazi 2005](https://www.sciencedirect.com/science/article/abs/pii/S0021999104004796) · [Burman–Zunino (DG weighted penalties)](https://link.springer.com/article/10.1007/s10915-008-9219-3) · [Ern & Guermond, Nitsche penalty notes](https://people.tamu.edu/~guermond/M661_FALL_2025/chap37.pdf) · [Nitsche penalty parameter note](https://arxiv.org/pdf/1709.05832)
