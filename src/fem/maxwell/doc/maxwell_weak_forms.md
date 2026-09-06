# Maxwell Weak Forms: From the Fundamental Lemma to h-Conform and b-Conform {#fem_maxwell_maxwell_weak_forms}

**Date:** 2026-08-14
**Purpose:** The weak-form derivations behind BELFEM's electromagnetic formulations: the
fundamental lemma of variational calculus, the two divergence-theorem corollaries, the least
squares projection, a thermal warm-up example, and the b-conform and h-conform weak forms of
the quasi-magnetostatic Maxwell equations.
**Module:** src/fem/maxwell
**Origin:** Distilled from Christian Messe's quasi-magnetostatic theory notes (LaTeX, written in
the Lagrange-multiplier era of the code). Transcription errata in the notes were fixed during
extraction. Each equation block names its TeX source for traceability.
**Implementation status:** BELFEM solves the **h-phi** formulation only (`en_Maxwell_Formulations.hpp`).
The b-conform derivations (a and a-v) are included for completeness of the theory; they are not
implemented as solving formulations. The least squares projection is live in the `L2*`
postprocessing kernels (`matrices/mt_maxwell_l2_h.cpp`, `mt_maxwell_l2_b.cpp`,
`mt_maxwell_l2_phi.cpp`).

---

## 1. Fundamental Lemma of Calculus of Variations

*(source: `weakform/fundamental_lemma.tex`)*

When calculus of variations was discovered, people investigated energy levels of potential
functions. Let `phi(x)` be a potential field that describes a physics problem, and let the
functional `Pi` be a measure for the compulsion the field experiences within a domain `Omega`.
From observation we know that such systems follow Le Chatelier's principle of assuming a state
where this compulsion is minimized over the full domain:

```
Pi = integral_Omega  G( phi(x) )  dV  =  min                                       (strongform)
```

This is called the *strong form*, and `G` is the governing equation of the problem in an
integral sense. Applying the chain rule gives the field that fulfills it:

```
delta Pi = integral_Omega  delta phi(x) * dG/dphi  dV  =  0                        (dstrongform)
```

Some problems, like linear elasticity, can be formulated directly this way. For many other
problems, however, the functional `Pi` is not known or does not exist at all. Galerkin's great
discovery is that the statement generalizes to a fundamental lemma:

```
delta Pi = integral_Omega  delta phi(x,t) * g( phi(x,t) )  dV  =  0                (flcv)
```

where `g(phi) = 0` is a partial differential equation describing the conservation law of the
problem. In this context, `delta phi` loses its physical meaning as a variation of `phi` and
becomes an abstract *test function*. But that does not matter: as long as `delta phi` is smooth
and `delta Pi = 0` is fulfilled, the conservation law is automatically fulfilled within the
domain. The formulation of `(flcv)` that can be discretized in an appropriate way is called the
*weak form* of the governing equation.

### 1.1 Important Divergence Theorems

Regardless of the physical problem, the *divergence theorem* or one of its corollaries is
always required to bring the weak form into a suitable shape. Two corollaries are of
significant importance. Let `c` be a scalar field, `f` and `g` vector fields in `Omega`, and
`n` the outward normal on the boundary `dOmega`. Combining Gauss's theorem with the chain rule
gives, for the expression `div(c g)`:

```
integral_Omega  (grad c)^T g + c div g  dV  =  integral_dOmega  c n^T g  dS        (gauss_theorem)
```

while a combination of Stokes' theorem (the Stokes-Cartan theorem) with the chain rule yields,
for `div(f x g)`:

```
integral_Omega  (curl f)^T g - f^T curl g  dV  =  integral_dOmega  f^T (g x n)  dS  (stokes_theorem)
```

---

## 2. Introductory Example 1: Least Squares Mapping

*(source: `weakform/leastsquares.tex`; live in `matrices/mt_maxwell_l2_*.cpp`)*

Since the edge degrees of freedom `h_tilde` have an abstract nature, they cannot be directly
visualized. Moreover, the vast majority of visualization tools such as ParaView require a node
or element based dataset. The goal is therefore to project the edge based field onto the nodes
while minimizing the projection error. One very efficient way is a least-squares mapping, which
is one of the easiest finite element problems to formulate. We apply the fundamental lemma to
the trivial statement that the magnetic field `h` equals the magnetic field `h`:

```
delta Pi = integral_Omega  delta_h^T ( h - h )  dV  =  0                           (l2weak)
```

The test function `delta_h` introduces a set of virtual degrees of freedom. They are not
physical; their sole purpose is to span a linear system that can be solved. For the nodal
field we write `h ≈ N h_hat`; following the idea of Boris Galerkin, the same interpolation is
used for the virtual degrees of freedom, `delta_h^T ≈ delta_h_hat^T N^T`. Since the edge dofs
are known, the field can also be expressed as `h ≈ E h_tilde`. Inserting both yields the system
`A h_hat = b` that projects the edge dofs onto the nodes:

```
[ integral_Omega N^T N dV ] * h_hat  =  [ integral_Omega N^T E dV ] * h_tilde
        = A                                       = b
```

Derivative field properties such as the electric current are computed the same way with the
curl operator `C` in place of `E`. The least squares mapping is not only used to project fields
and their derivatives onto different bases and even meshes; in the Lagrange-multiplier era of
these notes it also underlay the interface conditions for mixed formulations. In today's code
the interface coupling is hanging-edge condensation, and the L2 kernels serve postprocessing.

---

## 3. Introductory Example 2: Thermal Conduction

*(source: `weakform/heat.tex`; the production thermal solver lives in `src/fem/thermal`)*

The heat conduction problem is one of the simplest problems in physics that can be solved with
finite element theory, which makes it the ideal warm-up for constructing a weak form. The
governing equation reads:

```
g = rho c T_,t + div q_dot - q_dot_v = 0
```

with density `rho`, specific heat capacity `c`, temperature `T`, heat flux `q_dot` and a
volumetric heat load `q_dot_v` (imposed, for example, by electromagnetic induction). Having
identified `T` as the degree of freedom, the fundamental lemma gives:

```
delta Pi = integral_Omega  delta_T ( rho c T_,t + div q_dot - q_dot_v )  dV = 0
```

The Gauss corollary `(gauss_theorem)` moves the divergence onto the test function:

```
integral_Omega delta_T div q_dot dV
    = integral_dOmega delta_T n^T q_dot dS - integral_Omega (grad delta_T)^T q_dot dV
```

and applying Fourier's law `q_dot = -k grad T` yields the weak form:

```
integral_Omega delta_T rho c T_,t dV + integral_Omega (grad delta_T)^T k grad T dV
    = integral_Omega delta_T q_dot_v dV - integral_dOmega delta_T n^T q_dot dS     (weakformheat)
```

The expression `-n^T q_dot` is the scalar heat load into the structure. With the interpolation
rules

```
T ≈ N T_hat,    delta_T ≈ delta_T_hat^T N^T,    grad T ≈ B T_hat,    grad delta_T ≈ delta_T_hat^T B^T
```

(the notes carried a spurious transpose on the gradient rule; `B` is `(dim x n)`, so
`grad T ≈ B T_hat`), the discretized problem becomes:

```
[ integral_Omega N^T rho c N dV ] T_hat_,t + [ integral_Omega B^T k B dV ] T_hat
        = M                                          = K
    = [ integral_Omega N^T q_dot_v dV + integral_dOmega N^T q_dot dS ]
                              = f
```

with mass matrix `M`, stiffness matrix `K` and load vector `f`. Iterative solution methods
even allow `c` and `k` to depend on the temperature. A special thought should be invested into
the temperature dependence of the density. Thermal expansion does change `rho` with
temperature; but if the problem is solved on the undeformed mesh `Omega_0`, mass conservation
requires that the reference density `rho_0 = rho(T_0)` is used at all times:

```
m = integral_Omega rho dV = integral_Omega0 rho (dV/dV_0) dV_0 = integral_Omega0 rho_0 dV_0 = const.
```

---

## 4. Maxwell Equations

*(source: `weakform/maxwell.tex`)*

The Maxwell equations are Gauss's electric law, Gauss's magnetic law, Faraday's law of
induction, and the Ampère-Maxwell equation:

```
div d  = rho_v                     (gauss)
div b  = 0                         (gaussmagnet)
curl e = - b_,t                    (faraday)
curl h = j + d_,t                  (ampere)
```

Here `d` is the displacement field, `b` the magnetic flux density, `h` the magnetic field,
`j` the electric current density and `e` the electric field. The model is closed with Ohm's
law and the constitutive equations, which in the absence of polarization and magnetization
simplify to:

```
j = sigma e                        (ohm)
d = epsilon e
h = nu b                           (transh)
```

with electric conductivity `sigma` (resistivity `rho = sigma^-1`), permittivity `epsilon` and
reluctivity `nu`. In general these are spatially dependent tensors; many metals and the vacuum
behave isotropically, in which case they collapse to a scalar times the identity.

In the *magneto-quasistatic* simplification used throughout BELFEM, the change of the
displacement field `d_,t` and the volumetric charges `rho_v` are neglected.

For thermal conduction it is obvious that the temperature is the field to solve for and that
the heat flux is the flux term, so there is only one reasonable way to construct the weak form.
This is not the case for the Maxwell equations. Two families of formulations have proven
useful: the *b-conform* formulations and the *h-conform* formulations.

---

## 5. b-Conform Formulation (a-Formulation)

*(source: `weakform/bconform.tex`; **not implemented in BELFEM**, derived for completeness)*

A vector potential `a` is postulated whose curl yields the magnetic flux density:

```
b = curl a                                                                          (vektorpot)
```

Its existence is justified by the fact that Gauss's magnetic law is always fulfilled:
`div b = div curl a = 0`. The test function `delta_a` is applied to the Ampère-Maxwell
equation (displacement currents neglected):

```
integral_Omega delta_a^T curl h dV = integral_Omega delta_a^T j dV                  (bbasic)
```

The Stokes corollary `(stokes_theorem)` is applied to the left hand side (the notes wrote the
volume integral over `dOmega`; it is over `Omega`):

```
integral_Omega delta_a^T curl h dV
    = integral_Omega (curl delta_a)^T h dV - integral_dOmega delta_a^T (h x n) dS   (btheorem)
```

Inserting `(transh)` and `(vektorpot)` into the first term on the right gives the stiffness
expression, and reassembling yields the weak form of the *a-formulation*:

```
integral_Omega (curl delta_a)^T nu curl a dV
    + integral_dOmega delta_a^T (n x h) dS
    = integral_Omega delta_a^T j dV                                                 (weaka)
```

If the a-formulation is coupled with another formulation such as the h-formulation, the
boundary integral is used to formulate the connector element. The form above assumes that all
relevant current densities `j` are known, which is useful, for example, when a current is
imposed at a coil.

### 5.1 b-Conform with Electric Voltage (a-v Formulation)

*(source: `weakform/bconform_voltage.tex`; this file was an orphan in the notes, not compiled
into the main document, and is likewise **not implemented**)*

If currents in conducting regions are to be computed, the a-formulation must be extended with
a degree of freedom for the electric voltage `v`. The continuity equation for charge
conservation in the quasistatic case reads:

```
div j = -rho_v,t = 0                                                                (coulomb)
```

(the notes wrote `div j = +rho_v,t`; charge conservation carries the minus sign, which is
irrelevant here since the right hand side is zero anyway). Applying the fundamental lemma and
the Gauss corollary:

```
integral_Omega (grad delta_v)^T j dV = integral_dOmega delta_v n^T j dS             (jweak)
```

Combining Ohm's law with Faraday's law expresses the current in terms of both potentials:

```
e = -grad v - a_,t        =>        j = -sigma ( grad v + a_,t )                    (j_av)
```

which, inserted into `(jweak)` with a volume-imposed current density, gives:

```
integral_Omega (grad delta_v)^T sigma a_,t dV + integral_Omega (grad delta_v)^T sigma grad v dV
    = - integral_Omega (grad delta_v)^T j dV
```

and inserted into `(weaka)`:

```
integral_Omega delta_a^T sigma a_,t dV + integral_Omega (curl delta_a)^T nu curl a dV
    + integral_Omega delta_a^T sigma grad v dV + integral_dOmega delta_a^T (n x h) dS = 0
```

---

## 6. h-Conform Formulation

*(source: `weakform/hconform.tex`; this is the implemented core, see `matrices/mt_maxwell_h.cpp`)*

The fundamental lemma is applied to Faraday's law of induction; the magnetic field `h` is
identified as the degree of freedom, so the test function is named `delta_h`:

```
integral_Omega delta_h^T ( b_,t + curl e ) dV = 0                                   (hbasic)
```

For the first expression, the inverse of `(transh)` is used, with the product rule picking up
a term for a time-varying permeability:

```
integral_Omega delta_h^T b_,t dV
    = integral_Omega delta_h^T mu h_,t dV + integral_Omega delta_h^T mu_,t h dV
```

The Stokes corollary is applied to the second expression:

```
integral_Omega delta_h^T curl e dV
    = integral_Omega (curl delta_h)^T e dV + integral_dOmega delta_h^T (n x e) dS   (hstokes)
```

The boundary integral in `(hstokes)` is of special interest. If a coupling with the b-conform
formulation is desired, the electric field can be expressed through the voltage (if it exists)
and the vector potential; expressing `e` through the current density instead allows current
coupling:

```
integral_dOmega delta_h^T (n x e) dS
    = integral_dOmega delta_h^T [ (grad v + a_,t) x n ] dS
    = integral_dOmega delta_h^T [ n x (rho j) ] dS
```

With Ohm's law `e = rho curl h` in the volume and `v = 0`, the weak form of the h-conform
formulation reads:

```
integral_Omega delta_h^T mu h_,t dV
    + integral_Omega delta_h^T mu_,t h dV
    + integral_Omega (curl delta_h)^T rho curl h dV
    - integral_dOmega delta_h^T (n x a_,t) dS  =  0                                 (weakh)
```

A closer look reveals that the h-formulation is well suited for discretizing a superconducting
domain, while the expression `mu_,t` can become very unhandy to compute if `mu` is not
constant. In those cases the b-conform formulation avoids that term, which is why ferromagnetic
domains are the classical b-conform candidates.

In BELFEM, `(weakh)` is discretized with edge elements in the conductors (`edge_h` dofs,
`mt_maxwell_h.cpp`) and with the scalar potential `h = -grad phi` in the non-conducting
regions (`mt_maxwell_phi.cpp`), following the magnetodynamic coupling of Arsenault et al. 2023
(read with the 2026 erratum). A sign-convention note: the notes' operator rule `h ≈ B phi_hat`
is unsigned, and the phi kernel deliberately drops the Arsenault minus on the mass term
because `integral mu |h|^2` is even in the sign (see the comment in `mt_maxwell_phi.cpp`); the
minus matters wherever the sign survives, for example in an interface term. The coupling at
the conductor boundary is performed by hanging-edge condensation, not by an interface weak
form; see the maxwell usage guide and `cl_Maxwell_TMatrix`.

---

## References

- **Messe et al. 2023** - BELFEM's h-phi implementation, static condensation, solver strategy
- **Arsenault et al. 2023** (+ **2026 erratum**) - magnetodynamic h-phi coupling
- **Dular et al. 1999** - H-formulation and circuit coupling
- **Monk 2003** - Finite Element Methods for Maxwell's Equations
- **Kuczmann & Iványi 2008**, **Meunier 2008** - general Maxwell FEM references from the notes
- Full citations: `doc/literature_references.md`

For the element-level operators (`N`, `E`, `B`, `C`) that discretize these weak forms, see
`../../interpolation/doc/nedelec_derivation.md`.
