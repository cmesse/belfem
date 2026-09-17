# MATLAB derivation and verification scripts

**Date:** 2026-09-16
**Purpose:** The symbolic generators, zero-tests and derivation notes behind tables that the C++ in `src/fem/interpolation` carries as results. None of this is API and none of it is built; Doxygen skips these directories (`EXCLUDE_PATTERNS` in `Doxyfile.in`).
**Module:** `src/fem/interpolation`

Every script needs MATLAB with the Symbolic Math Toolbox, or Octave with the `symbolic` package, unless its README says otherwise. Under Octave the drivers run as `octave-cli --eval "pkg load symbolic; check_tet10"` from the script's directory; that is how the gate was last run. Two kinds of script live here, and each README says which is which:

- **Generators and zero-tests** reproduce a C++ table exactly. Their drivers assert every residual, so `matlab -batch <driver>` exiting 0 is the gate.
- **Derivation notes** record a derivation or a convention that the C++ states only as a result. They run, or are interactive plots, but they do not emit a C++ table, and where their conventions differ from the C++ the header comment says so.

| Directory | Concerns | Gate |
|---|---|---|
| `nedelec_tet10/` | the Nédélec edge and face tables of `cl_EF_TET10.cpp` and the Lagrange functions of `cl_IF_TET10.hpp` | `matlab -batch check_tet10` |
| `nedelec_tri/` | the TRI3 curl operator of `cl_EF_TRI3.cpp`, the TRI6 ansatz of `cl_EF_TRI6.cpp`, and the `[0,1]²` to `[-1,1]²` remap of the facet integration | notes only |
| `lagrange/` | the second-derivative table of `cl_IF_PENTA18.hpp` | `matlab -batch check_penta18` |
| `facets/` | the PENTA6 slave-face parameter maps and outward normals of `fn_IF_initialize_integration_points_on_facet.cpp` | interactive notebooks |

The zero-tests compare the generator against tables typed into the MATLAB scripts. `compare_tables.py` (needs `sympy`) checks those transcriptions against the C++ files themselves, so the two together tie the generator to the source:

```
python3 compare_tables.py        # 0 differences expected; exit 1 otherwise
```

The runtime coverage is elsewhere: `tests/fem/test_LagrangeInterpolation.cpp` checks second derivatives by central differences, `tests/fem/test_FacetIntegrationPoints.cpp` checks every slave orientation on real coordinates, and `tests/fem/test_EdgeFunctions.cpp` checks circulations. These scripts preserve what those tests cannot regenerate: the symbolic construction.
