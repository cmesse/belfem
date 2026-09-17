# Triangle derivation notes

**Date:** 2026-09-16
**Purpose:** The written-down derivations behind the TRI3 curl operator, the TRI6 edge-function ansatz, and the facet-integration remap. None of these is a generator: the C++ carries results in a different normalization or sign convention, and each header comment says which.
**Module:** `src/fem/interpolation`

| Script | Concerns | Where it differs from the C++ |
|---|---|---|
| `curl.m` | the curl of the three signed TRI3 Whitney functions, transformed with the inverse Jacobian, whose derivation `nedelec_derivation.md` §3.2 summarizes | Whitney sign convention `eta*xi_x - xi*eta_x` is the opposite of `cl_EF_TRI3.cpp` (`G*nabla_eta - H*nabla_xi`); the collapse to `C = 2/detJ [s1 s2 s3]` is left as a commented line |
| `fragment.m` | the TRI6 ansatz `theta = a*nabla_j + b*nabla_i` for the six edge functions and the two face functions, in the form that `cl_EF_TRI6.cpp` expands | the C++ polynomials carry half the coefficients (first pair checked); its normalization is unit circulation |
| `secondordertri.m` | the two-point-per-edge `(a, b)` ansatz that `fragment.m` expands | stops at the ansatz |
| `check_quad.m` | the Jacobians of the `[0,1]²` and `[-1,1]²` quad parametrizations differ by the constant factor induced by `alpha = 0.5*(1+xi)`, the remap `fn_IF_initialize_integration_points_on_facet.cpp` uses | nothing to compare |

All four run headless (`matlab -batch curl`, or `octave-cli --eval "pkg load symbolic; curl"` from this directory) and print their results, except `secondordertri.m`, which only defines the ansatz; there is nothing to assert. Octave warns that `curl.m` shadows its built-in `curl`; as a script run from here that is harmless. Symbolic Math Toolbox or the Octave `symbolic` package required.
