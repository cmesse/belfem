# Changelog

All notable changes to BELFEM are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and version numbers follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

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
