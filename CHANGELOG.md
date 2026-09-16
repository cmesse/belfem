# Changelog

All notable changes to BELFEM are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **New example deck `examples/tapestack3d_quench`.** The `tapestack3d` geometry — eight
  soldered REBCO tapes, periodic, coupled h-φ/T — driven into a quench by a user plugin
  instead of a prescribed source: a linear 200 A/s transport ramp and a local Ic defect on
  the topmost tape, both reading one time program (`src/ramps.hpp`) so they cannot drift
  apart. It is the worked example for combining a `userdefined` current boundary condition
  with a material `defect { }` plugin in a single shared object. Nineteen decks now ship.

### Changed

- **CMake 4.0 is now the minimum**, raised from 3.11 and applied everywhere the project states
  a floor: the top-level `CMakeLists.txt`, both plugin templates
  (`src/physics/materials/User{Material,Library}Template.cmake`, previously 3.13) and the five
  example-deck plugin builds. The project is developed against a toolchain that ships CMake 4,
  and `scripts/scls_env.sh` already required `cmake` to resolve there; the declared floors were
  the last place still claiming otherwise. Every policy up to 4.0 now defaults to NEW —
  a full configure was re-run and produces no warning and no policy or deprecation message. A
  distribution CMake older than 4.0 can no longer configure BELFEM.

- **Third-party headers are system headers for diagnostics.** The include directories of the
  linked libraries (SCLS, MKL, STRUMPACK, PETSc, HDF5, …) are passed to the C++ compiler with
  `-isystem` instead of `-I`. A warning raised inside one of these headers no longer appears in
  a BELFEM build and cannot become an error under `-Wall -Werror`. This is preventative: the
  0.9.0 build failure on Ubuntu came from two BELFEM headers and was fixed in 0.9.1, but the
  same mechanism (Ubuntu's GCC injecting `-Wformat-security`) would turn the next vendor-header
  warning into a build failure. BELFEM's own headers keep ordinary `-I` and are diagnosed as
  before; a directory the compiler already searches by default also stays `-I`. The per-test
  `BEFORE` include of `$SCLS/include` was removed as redundant, so the tests now search
  BELFEM's headers before SCLS, as the library always did. Plugins built from the
  `User*Template.cmake` templates are not affected; they still receive ordinary include paths.
- **`spline::create_helpmatrix` takes its size as `const index_t` and its spacing as
  `const real`**, both by value (previously `const real &` for both). Every caller in the tree
  already passed an integer count. A separately built consumer must be recompiled.

- **Validity ceilings of copper, silver, chromium and white tin lowered** (`T_max` 1358 → 1000 K,
  1235 → 900 K, 2180 → 570 K, 505 → 400 K): above these temperatures the new elastic curves depart
  by more than 5 % in shape from Blanke's dynamic-modulus compilation, the only high-temperature
  reference available. Every property of these metals, and of the formula alloys that contain
  them, now stops there; iron (860 K) and aluminum (933 K) are unchanged.

- **The `material` report prints the Grüneisen parameter at 298.15 K**, the temperature of the
  printed density, as a diagnostic computed from the served bulk modulus, α, ρ and c_p. It is no
  longer a model input.

### Fixed

- **The remaining `-Wfloat-conversion` warnings in the library and the test headers are
  cleared.** The flag was added in 0.9.1; the implicit floating-point-to-integer conversions it
  reported — `std::ceil`/`std::round` results assigned to integers in the integration-point
  tables, `Section::get_int`, the `Alloy`, `Metal` and `SplineLookupTable` spline setups,
  `Genome::encode`, and two string helpers — are now explicit casts of the same expressions.
  Results are unchanged for the inputs the code is given in the tree. `<cmath>` is included
  directly in the three files that previously used `ceil`/`round` unqualified.
- **Input-file integers are range-checked.** `Section::get_int` rounded the stored real and
  converted it to `int` unchecked, which is undefined behaviour for a non-finite or huge value;
  it now stops with an error naming the key and section. A unit exponent (`m^2`, `s^-1`) was
  parsed with `std::stod` and truncated to `int` unchecked, so `m^1.5` silently became `m^1` and
  `m^abc` escaped as an uncaught exception; the exponent must now be an integer with magnitude
  at most 32, and anything else stops with an error naming the unit. Rounding in `get_int` is
  unchanged. Both checks run once while the input file is read and are active in every build.

- **Poisson's ratio of the pure metals fell with temperature; it now rises, as measured.** The
  elastic closure of `Metal::create_mech` held the Grüneisen parameter constant and inverted the
  Grüneisen relation for the bulk modulus, which copied the unknown temperature dependence of γ
  into K: copper's bulk modulus softened 12 % between 0 and 293 K where measurement gives 3 %, and
  ν read 0.363 → 0.344 where Ledbetter 1981 (10.1002/pssa.2210660209) gives 0.340 → 0.347. Seven of
  the nine metals were inverted. The bulk and shear moduli are now exponentials of the logarithmic
  thermal strain the expansion curve already integrates (quasi-harmonic, Garai & Laugier 2007,
  10.1063/1.2424535), fitted per metal to measured single-crystal or polycrystal constants between
  4 K and room temperature (sources with DOIs in `material_property_sources.md`), with E and ν
  derived. The served moduli are isothermal, converted from the ultrasonic (adiabatic) data with the
  material's own α, ρ and c_p. Room-temperature levels move to the measured polycrystal values:
  copper E 133 → 128 GPa, silver 84 → 81, aluminum 73 → 70, tin 55 → 49, nickel 184 → 223 (the old
  value was the demagnetized state; BELFEM's nickel is magnetically saturated, so the ΔE dip below
  the Curie point is deliberately not represented, see the class header), lead 16.5 → 16.25 at the
  static level (Blanke's E with the crystal's bulk modulus; the dynamic Hill average would be 24 GPa).
  Chromium serves a constant ν = 0.237 because its spin-density-wave anomalies collapse the bulk
  modulus by 20 % toward the Néel point and no smooth closure follows that. Iron's isothermal ν is
  flat within 1e-4. Construction refuses constants that let ν fall by more than 1e-3 anywhere.
  New regression test `tests/physics/test_MetalElastic.cpp`.

## [0.9.1] — 2026-09-13

A portability release. BELFEM 0.9.0 did not build on Ubuntu 24.04 (GCC 13.3), and once it
built, some floating-point convergence and threshold tests used the wrong `abs` overload.
Both problems are fixed. Input files and output formats are unchanged.

### Fixed

- **Floating-point `abs` resolved to the integer C `abs`.** Eleven unqualified `abs()` calls on
  `double` values could bind to `int abs(int)` from `<stdlib.h>`. This happens with libstdc++
  when no header in the include chain declares `::abs(double)`. The argument was then truncated
  to an integer, so any residual below 1 counted as zero. All eleven calls now use `std::abs`.
  Affected code:
  - `Metal::set_RRR` (`cl_Material_Metal.cpp`): the secant loop for the residual resistivity
    `rho_0` never ran, so metals using this path kept the initial guess. In the verified YBCO case
    the thermal conductivity was off by −0.94 % at 20 K. This is what failed
    `YBCOThermalConductivity.RepresentativeTemperaturesAreFinitePositive`.
  - `Metal` Debye-temperature inversion, the YBCO thermal-expansion plateau and c<sub>p</sub>
    switch searches, and the root searches in `Alloy` and `HastelloyC276`.
  - `Tmatrix` (`cl_FEM_Tmatrix.cpp`): matrix entries with magnitude below 1 were dropped from
    the sparsity pattern.
  - `fn_Mesh_ratio.cpp`: the secant loop for the mesh grading ratio.
  - `belfem` (`belfem.cpp`): the thermal/electromagnetic time-synchronisation loop could exit
    while the two clocks still differed by almost 1 s.

  Builds that already found the `double` overload produce the same results as before. On
  affected builds, material properties and coupled thermal/electromagnetic transients may
  change, because the truncation defect is removed.
- **Build failure with `-Werror=format-security`.** Ubuntu's GCC driver injects
  `-Wformat-security` by default, which turned two warnings into errors under BELFEM's
  `-Wall -Werror`:
  - `Logger` (`cl_Logger.hpp`) passed an already formatted message to `fprintf` as the format
    string. A `%` in a message, such as a file path or `50%`, was read as a conversion
    specifier, which is undefined behaviour on every platform. The message is now printed
    with `"%s"`.
  - `sprint()` (`fn_sprint.hpp`) forwards a caller-supplied format by design, and its GCC
    pragma block now also suppresses `-Wformat-security`.
- **`Matrix<T>` fill constructor took the fill value as `double`.** For integer matrices the
  value made a `double` round trip. That was exact for every in-tree caller on default builds,
  but with `USE_MKL_64BIT_API` (`index_t` = `uint64_t`) filling with the sentinel `gNoIndex`
  (2⁶⁴−1) went through a `double` that cannot hold it, which is undefined behaviour. The
  constructor now takes `T` and passes the value through unchanged.

### Changed

- GCC and Clang builds now warn on implicit floating-point to integer narrowing
  (`-Wfloat-conversion`). It is a warning only, not an error, so it cannot fail a build; it
  is there to catch the `abs` class of defect above at compile time.
- Homology: the Smith normal form no longer routes integer divisions through `floor()` (a
  no-op on `int` operands), and its integer `abs` calls are spelled `std::abs`. Reviewed by
  two independent auditors as behaviour-preserving; the coefficient arithmetic is unchanged.

## [0.9.0] — 2026-09-06

First public release, presented at ASC 2026.

[0.9.1]: https://github.com/cmesse/belfem/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/cmesse/belfem/releases/tag/v0.9.0
