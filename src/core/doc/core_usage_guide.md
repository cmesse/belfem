# Core Module Overview {#core_core_usage_guide}

**Date:** 2026-01-16
**Module:** src/core
**Purpose:** Comprehensive documentation of BELFEM's core utility module

**Revision History:**
- 2026-01-16: Initial documentation
- 2026-01-16: Technical corrections (Profiler, Timer, thread safety)

---

## Core Module Contracts

Before diving into details, understand these fundamental design contracts:

1. **Minimal Dependencies**: Core utilities must remain dependency-minimal (no FEM, mesh, numerics modules)
2. **Debug vs Release**: Many checks use `BELFEM_ASSERT`, which is compiled out in release (`NDEBUG`); a release build still dies on `BELFEM_ERROR`
3. **Controlled Global State**: Globals (`gLog`, `gTbulk`) are intentional and centrally managed, not accidental
4. **MPI Awareness**: Some utilities behave differently in MPI vs serial builds (random, logging, profiling)
5. **Type Safety**: BELFEM's typedefs (`real`, `index_t`, `id_t`) are the single source of truth—do not shadow or redefine

---

## Table of Contents

1. [Introduction](#introduction)
2. [Module Architecture](#module-architecture)
3. [Logging and Output](#logging-and-output)
4. [Timing and Profiling](#timing-and-profiling)
5. [Command Line Processing](#command-line-processing)
6. [Type System](#type-system)
7. [String Utilities](#string-utilities)
8. [Error Handling](#error-handling)
9. [Hashing](#hashing)
10. [Random Number Generation](#random-number-generation)
11. [Global Variables](#global-variables)
12. [Thread Safety and MPI](#thread-safety-and-mpi)
13. [Common Pitfalls](#common-pitfalls)
14. [Usage Examples](#usage-examples)

---

## Introduction

The `core` module provides fundamental utilities and infrastructure used throughout the BELFEM framework. It establishes the type system, provides logging and debugging facilities, and implements common utilities for string manipulation, timing, profiling, and command-line processing.

**Key Design Principles:**
- **Lightweight dependencies**: Core depends only on the C++ standard library, the `containers` module (`Cell`) and the `comm` layer (`gComm`, `comm_abort`) — never on FEM, mesh or numerics
- **Cross-platform compatibility**: Works on Linux and Darwin (macOS)
- **Compiler agnostic**: Supports GCC, Clang, and Intel compilers with appropriate pragmas
- **MPI-aware**: Some utilities (random, logging) adapt to MPI environments

**Location:** `src/core/`
**CMake Target:** `belfem_core`

---

## Module Architecture

### File Organization

The core module contains 26 source files organized into the following categories:

```
src/core/
├── Logging & Output
│   ├── cl_Logger.{hpp,cpp}
│   ├── cl_Progressbar.{hpp,cpp}
│   └── banner.{hpp,cpp}
│
├── Timing & Profiling
│   ├── cl_Timer.hpp
│   └── cl_Profiler.{hpp,cpp}
│
├── Command Line
│   └── cl_Arguments.{hpp,cpp}
│
├── Utilities
│   ├── cl_Hash.hpp
│   ├── stringtools.{hpp,cpp}
│   ├── fn_sprint.hpp
│   └── fn_available_memory.{hpp,cpp}
│
├── Type System
│   ├── typedefs.hpp
│   ├── constants.hpp
│   ├── units.hpp
│   ├── fn_check_unit.hpp
│   └── fn_to_enum.hpp
│
├── Error Handling
│   └── assert.{hpp,cpp}
│
├── Random Numbers
│   └── random.hpp
│
└── Globals
    └── globals.hpp
```

### Dependencies

The core module is self-contained with minimal external dependencies:

- **C++ Standard Library**: `<chrono>`, `<string>`, `<cstdio>`, `<random>`, etc.
- **Container library**: `cl_Cell.hpp` (BELFEM's dynamic array)
- **MPI** (optional): For `cl_Communicator` in random number generation

---

## Logging and Output

### Logger Class

**File:** `cl_Logger.{hpp,cpp}`

The `Logger` class provides hierarchical logging with configurable verbosity levels.

#### Info Levels

```cpp
enum class InfoLevel
{
    Silent      = 0,  // No output
    Minimal     = 1,  // Minimal output
    Default     = 2,  // Default output
    Detailed    = 3,  // More mesh information
    Verbose     = 4,  // BELFEM-specific debugging
    Everything  = 5   // Third-party library debugging
};
```

#### Key Features

- **Hierarchical verbosity**: Messages only print if `message_level <= logger_level`
- **Stream flexibility**: Output to `stdout` or file
- **Printf-style formatting**: Uses variadic templates with `sprint()`
- **Global instance**: `extern Logger gLog` for framework-wide access

#### Implementation Details

**Location:** `cl_Logger.hpp:43-108`

```cpp
class Logger
{
    uint       mInfoLevel = 0;
    std::FILE* mStream;
    bool       mWriteToAscii = false;

public:
    Logger(const InfoLevel aInfoLevel);
    Logger(const InfoLevel aInfoLevel, const std::string & aPath);
    ~Logger();

    uint info_level() const;

    template <typename ... Args>
    void message(const InfoLevel aInfoLevel,
                 const std::string & aFormat,
                 const Args ... aArgs);
};
```

**Global convenience function:**
```cpp
template <typename ... Args>
void message(const belfem::InfoLevel aInfoLevel,
             const std::string & aFormat,
             const Args ... aArgs)
{
    gLog.message(aInfoLevel, aFormat, aArgs ...);
}
```

#### Logger Lifecycle

**Critical:** The global `gLog` must be constructed **once** at application start, before any module code executes. Core and all modules assume `gLog` is already initialized. Reassigning or reinitializing `gLog` mid-run is undefined behavior.

```cpp
// In main.cpp - construct once at program start
belfem::Logger gLog(belfem::InfoLevel::Default);
```

Framework code can then use the global logger anywhere:

```cpp
// In any module
extern belfem::Logger gLog;  // Declaration (already defined in main)

// Log messages throughout the code
message(InfoLevel::Default, "Solving system with %d DOFs", num_dofs);
message(InfoLevel::Verbose, "Matrix assembly time: %.3f ms", time);
```

---

### Progressbar Class

**File:** `cl_Progressbar.{hpp,cpp}`

Visual progress indicator for long-running operations.

#### Key Features

- **Configurable steps**: Default 100 steps (1% increments)
- **Stream output**: Defaults to `stdout`, can redirect to file
- **Step control**: Manual step prescription or auto-increment

#### Implementation Details

**Location:** `cl_Progressbar.hpp`

```cpp
class Progressbar
{
    const index_t mWidth = 65;  // Bar width in characters
    const index_t mNumSteps;    // Total steps (default 100)
    FILE * mFile;
    uint mProgress = 0;         // Current progress counter
    uint mStep = 0;             // Bar width already drawn

public:
    Progressbar(const uint aNumSteps=100, FILE * aFile = stdout);

    void reset();
    void step(const uint & aProgress);  // Set specific step
    void step();                         // Auto-increment
    void finish();

private:
    void draw(const uint aStep, const uint aProgress);  // One frame
    void flush();                                       // Push the frame out
};
```

The bar only redraws when the drawn width actually grows, so the number of
frames is bounded by `mWidth` no matter how many times `step()` is called.

#### Usage Pattern

```cpp
Progressbar bar(num_iterations);
for (uint i = 0; i < num_iterations; ++i)
{
    // ... work ...
    bar.step();
}
bar.finish();
```

**⚠️ Thread Safety:** `Progressbar` is **not thread-safe**. Use only from a single thread (typically rank 0 / master in MPI applications).

**Output buffering.** The bar redraws in place with `\r` and never writes a
newline, so nothing in the C library flushes it on its own. Every frame is
therefore followed by an explicit `std::fflush`. This matters most under
`mpirun`: there `stdout` is a pipe rather than a terminal, so it is fully
buffered in 4 KiB blocks, and a short bar would otherwise reach the user in a
single burst when the program exits. Flushing per frame also emits each frame
in one `write()`, which keeps the MPI I/O forwarder from tearing it apart.
Should output still arrive in bursts, the remaining buffering is on the launcher
side; `mpirun --stream-buffering 0` turns it off.

---

### Banner System

**File:** `banner.{hpp,cpp}`

Application startup banner with seasonal variants.

#### Key Features

- **System information**: OS, CPU info via `uname` and `cpu_info()`
- **Seasonal banners**: Easter, St. Patrick's, USA, Canada, Thanksgiving, Christmas
- **Automatic selection**: Based on system date

#### API

**Location:** `banner.hpp:24-127`

```cpp
std::string exec(const std::string & aCommand);
std::string uname();        // Returns "Linux" or "Darwin"
std::string cpu_info();     // CPU information for banner
std::string version();
bool        is_built_from_git();
std::string git_commit_hash();
std::string git_commit_hash_short();
std::string git_branch();
bool        git_is_dirty();

void print_banner(const std::string aExecName = "");

namespace banners
{
    bool print_default();
    bool print_easter();
    bool print_stpatrick();
    bool print_usa();
    bool print_canada();
    bool print_thanksgiving();
    bool print_christmas();
}
```

#### Global Constants

```cpp
const std::string gLongName = "BELFEM -- The Berkeley Lab Finite Element Framework";
const std::string gURL      = "http://belfem.lbl.gov";
```

---

## Timing and Profiling

### Timer Class

**File:** `cl_Timer.hpp`

High-resolution wall-clock timer using `std::chrono::high_resolution_clock`.

#### Key Features

- **Millisecond precision**: Returns elapsed time in ms
- **Wall-clock time**: Measures actual elapsed time, not CPU time (includes MPI waits, I/O)
- **Lightweight**: Header-only implementation with inline methods
- **Three operations**: `stop()`, `next()`, `reset()`

#### Implementation Details

**Location:** `cl_Timer.hpp:21-63`

```cpp
class Timer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> mStart;

public:
    inline Timer() : mStart(std::chrono::high_resolution_clock::now()) {}

    // Returns elapsed time since construction/reset
    inline uint64_t stop() {
        return (unsigned int)(
            std::chrono::duration_cast<std::chrono::milliseconds>
            (std::chrono::high_resolution_clock::now() - mStart).count()
        );
    }

    // Returns elapsed time and resets timer
    inline uint64_t next() {
        unsigned int aTime = this->stop();  // ⚠️ Truncates to uint (49 days max)
        mStart = std::chrono::high_resolution_clock::now();
        return aTime;
    }

    // Resets timer to current time
    inline void reset() {
        mStart = std::chrono::high_resolution_clock::now();
    }
};
```

**⚠️ Limitation:** `next()` stores the result in `unsigned int` before returning `uint64_t`, so it truncates past ~49.7 days (2³² milliseconds). **`stop()` does not fix this** — it casts the same way (`cl_Timer.hpp:42-47`). Neither method is safe for intervals that long; measure them from an external clock.

#### Usage Pattern

```cpp
Timer timer;
// ... work phase 1 ...
uint64_t time1 = timer.next();  // Get time for phase 1, reset for phase 2
// ... work phase 2 ...
uint64_t time2 = timer.stop();   // Get time for phase 2

message(InfoLevel::Default, "Phase 1: %lu ms, Phase 2: %lu ms", time1, time2);
```

---

### Profiler Class

**File:** `cl_Profiler.{hpp,cpp}`

CPU profiling using Google's gperftools with Callgrind output format.

#### Key Features

- **gperftools integration**: Uses `ProfilerStart`/`ProfilerStop` from gperftools
- **Callgrind conversion**: Automatically converts profiler output to Callgrind format via `pprof --callgrind`
- **Selective profiling**: Profile only critical code sections
- **MPI-aware**: Creates separate profile files per rank in parallel runs
- **CMake control**: Enable with `-DUSE_PROFILER=ON` (generates the `BELFEM_PROFILER` define)

#### Implementation Details

**Location:** `cl_Profiler.hpp:21-46`

```cpp
class Profiler
{
    string mLogFile;
    string mCallgrindFile;

public:
    Profiler(const string aLogFile="profiler.log");
    ~Profiler() = default;

    void start();  // Begin gperftools profiling
    void stop();   // Stop profiling and convert to Callgrind format
};
```

**Implementation** (from `cl_Profiler.cpp:59-92`):
- `start()`: Calls `ProfilerStart(mLogFile.c_str())`
- `stop()`: Calls `ProfilerStop()`, then executes `pprof --callgrind` to generate `.callgrind` file

**Build requirement:** Requires gperftools library and `BELFEM_PROFILER` compile flag. Without this flag, `start()` and `stop()` are no-ops.

#### Usage Pattern

```cpp
Profiler profiler("my_analysis.log");
profiler.start();
// ... code to profile ...
profiler.stop();

// Analyze results with:
// callgrind_annotate my_analysis.callgrind
// Or use KCachegrind/QCachegrind for visualization
```

**Note:** The profiler generates **two files**: a gperftools native format log and a converted Callgrind file for compatibility with standard analysis tools.

---

## Command Line Processing

### Arguments Class

**File:** `cl_Arguments.{hpp,cpp}`

Base class for command-line argument parsing.

#### Key Features

- **Extensible design**: Virtual destructor for inheritance
- **Cell storage**: Arguments stored in `Cell<string>`
- **Index access**: Retrieve individual arguments by index
- **Built-in verbosity flag**: the constructor scans for `-v` / `--verbose`
  (GNU style) and sets the info level of the global logger

#### Built-In Flags

Every executable that constructs an `Arguments` object (or a subclass) accepts
the shared verbosity flag; no code in the executable is needed beyond the
construction itself:

| Form | Effect |
|---|---|
| `-v`, `--verbose` | `InfoLevel::Everything` (5) |
| `-v 3`, `-v3`, `--verbose 3`, `--verbose=3` | info level 3 |

Levels follow `InfoLevel` in `cl_Logger.hpp` (0 = silent … 5 = everything,
including third-party library output). Subclasses must not reuse `-v` or
`--verbose`; a version flag follows the GNU convention `-V` / `--version`
(see `gastables::Arguments`).

#### Implementation Details

**Location:** `cl_Arguments.hpp:20-57`

```cpp
class Arguments
{
protected:
    Cell<string> mArguments;  // List of arguments

public:
    Arguments(int & argc, char * argv[]);
    virtual ~Arguments() = default;

    // Returns all arguments
    const Cell<string> & data() const;

    // Access specific argument by index
    const string & data(const index_t aIndex) const;
};
```

#### Design Pattern

The `Arguments` class is designed as a base class for application-specific argument parsers. Applications derive from `Arguments` and add custom parsing logic:

```cpp
class MyAppArguments : public Arguments
{
    bool mVerbose = false;
    string mInputFile;

public:
    MyAppArguments(int & argc, char * argv[]) : Arguments(argc, argv)
    {
        // Parse mArguments for application-specific flags
    }
};
```

---

## Type System

### Type Definitions

**File:** `typedefs.hpp`

Establishes the unified type system used throughout BELFEM.

#### Fundamental Types

**Location:** `typedefs.hpp:25-37`

```cpp
typedef size_t              size_t;
typedef std::string         string;

typedef int                      sint;
typedef long int                 lsint;
typedef unsigned int             uint;
typedef short unsigned int       suint;
typedef long unsigned int        luint;
typedef long long unsigned int   lluint;

typedef double                   real;
typedef std::complex<real>       cplx;
```

#### Framework-Specific Types

**Location:** `typedefs.hpp:40-52`

```cpp
typedef unsigned int             id_t;      // Entity IDs
typedef int                      proc_t;    // Process/rank IDs
typedef long long unsigned int   key_t;     // Hash keys (64-bit)
typedef __uint128_t              key128_t;  // Hash keys (128-bit)

#ifdef BELFEM_INT64
    typedef int64_t              int_t;     // Generic integers
    typedef uint64_t             index_t;   // Array indices
#else
    typedef int32_t              int_t;     // Generic integers (32-bit default)
    typedef uint32_t             index_t;   // Array indices (32-bit default)
#endif
```

**Rationale:** The `index_t` type can be compiled in 32-bit or 64-bit mode depending on mesh size requirements. Most problems use 32-bit indices for memory efficiency.

**⚠️ Portability Note:** `key128_t` uses `__uint128_t`, a GCC/Clang extension not available in all compilers (e.g., MSVC). Code using 128-bit keys may not be portable to non-GNU toolchains.

**⚠️ Type Safety Rule:** Do not introduce `using std::size_t` or redefine fundamental typedefs in higher-level modules. BELFEM's typedefs are the single source of truth for ABI consistency.

#### Sentinel Values

**Location:** `typedefs.hpp:56-62`

```cpp
constexpr index_t gNoIndex = std::numeric_limits<index_t>::max();
constexpr id_t    gNoID    = std::numeric_limits<id_t>::max();
constexpr proc_t  gNoOwner = std::numeric_limits<proc_t>::max();

constexpr real    gTfreeze = 273.15;  // K
constexpr real    gTref    = 288.15;  // K (15°C)
constexpr real    gTroom   = 293.15;  // K (20°C)
```

#### Unit System

**Location:** `typedefs.hpp:65-68`

```cpp
// L, M, T, I, theta, N, J (Length, Mass, Time, Current, Temperature, Amount, Luminosity)
typedef std::array<real, 7> unit;

// Value with associated units
typedef std::pair<real, unit> value;
```

The `unit` type represents physical dimensions as a 7-element array following the SI base unit system.

**⚠️ Limitation:** The unit system performs **dimension checking only** via array comparison. It does **not** enforce unit consistency at compile time, nor does it track unit scaling automatically. Users must manually apply conversion factors via the constants in `units.hpp`.

#### Limits and Special Values

**Location:** `typedefs.hpp:72-85`

```cpp
#define BELFEM_SINT_MAX      std::numeric_limits<sint>::max()
#define BELFEM_UINT_MAX      std::numeric_limits<uint>::max()
#define BELFEM_REAL_MAX      std::numeric_limits<real>::max()
#define BELFEM_REAL_MIN      std::numeric_limits<real>::min()
#define BELFEM_INT_MAX       std::numeric_limits<int>::max()
#define BELFEM_LUINT_MAX     std::numeric_limits<luint>::max()
#define BELFEM_KEY_MAX       std::numeric_limits<key_t>::max()
#define BELFEM_SIGNALING_NAN std::numeric_limits<real>::signaling_NaN()
#define BELFEM_QUIET_NAN     std::numeric_limits<real>::quiet_NaN()
#define BELFEM_INFINITY      std::numeric_limits<real>::infinity()

constexpr real BELFEM_EPSILON      = 10 * std::numeric_limits<real>::epsilon();
constexpr real BELFEM_EPS          = std::numeric_limits<real>::epsilon();
constexpr real BELFEM_MESH_EPSILON = 1e-9;
```

**Usage notes:**
- `BELFEM_EPSILON`: Used for general floating-point comparisons (≈ 2.22e-15)
- `BELFEM_EPS`: Machine epsilon (≈ 2.22e-16)
- `BELFEM_MESH_EPSILON`: Geometric tolerance for mesh operations (1 nm)

**⚠️ Critical Distinction:** Use `BELFEM_EPSILON` for mathematical/algorithmic comparisons (convergence, residuals). Use `BELFEM_MESH_EPSILON` for spatial/geometric tolerances (node coincidence, element validity). Mixing these can cause topological errors in mesh operations.

---

### Physical Constants

**File:** `constants.hpp`

NIST-compliant physical and mathematical constants.

#### Mathematical Constants

**Location:** `constants.hpp:37-51`

```cpp
namespace belfem::constant
{
    const real pi  = 3.141592653589793238462643383279502884;
    const real phi = 0.5 * (1.0 + std::sqrt(5.0));  // Golden ratio
    const real deg = pi / 180.0;                     // Degree to radian
}
```

#### General Physics

**Location:** `constants.hpp:61-100`

```cpp
const real c        = 299792458.0;              // Speed of light [m/s] (exact)
const real mu0      = 1.25663706212e-6;         // Magnetic constant [V·s/(A·m)]
const real nu0      = 1.0/mu0;                  // Inverse magnetic constant
const real epsilon0 = 1.0/(mu0*std::pow(c, 2)); // Electric constant [A·s/(V·m)]
const real G        = 6.67430e-11;              // Gravitational constant [N·m²/kg²]
const real e        = 1.602176634e-19;          // Elementary charge [C]
```

All values sourced from NIST (https://physics.nist.gov/cuu/Constants/).

#### Thermodynamics

**Location:** `constants.hpp:110-178`

```cpp
const real calTh  = 4.184;                // Thermal calorie [J/cal]
const real kB     = 1.380649e-23;         // Boltzmann constant [J/K] (exact)
const real NA     = 6.02214076e23;        // Avogadro constant [1/mol] (exact)
const real u      = 0.001 / NA;           // Atomic mass unit [kg]
const real Rm     = kB * NA;              // Gas constant [J/(K·mol)] (exact)
const real Rm_cal = Rm / calTh;           // Gas constant [cal/(K·mol)]
const real h      = 6.62607015e-34;       // Planck constant [J·s] (exact)
const real sigma  = 2.0 * std::pow(pi, 5) * std::pow(kB, 4) /
                    (15.0 * std::pow(h, 3) * std::pow(c, 2));  // Stefan-Boltzmann
const real L0     = std::pow(pi * kB / e, 2) / 3.0;            // Lorenz number [V²/K²]
```

#### Earth Constants

**Location:** `constants.hpp:188-202`

```cpp
const real g0       = 9.80665;           // Standard gravity [m/s²]
const real mu_earth = 3.986004418e14;    // Earth gravity parameter [m³/s²]
const real R_earth  = std::sqrt(mu_earth / g0);  // Earth reference radius [m]
```

---

### Unit Conversions

**File:** `units.hpp`

Conversion factors from common units to SI base units.

#### Metric Units

**Location:** `units.hpp:37-41`

```cpp
namespace belfem::constant
{
    const real l   = 0.001;              // liter [m³]
    const real km  = 1000;               // kilometer [m]
    const real bar = 1.0e5;              // bar [Pa]
    const real atm = 1.01325e5;          // atmosphere [Pa]
    const real rpm = constant::pi / 30.0;  // revolutions per minute [rad/s]
}
```

#### Imperial/US Units

**Location:** `units.hpp:47-54`

```cpp
const real in  = 0.0254;                    // inch [m]
const real ft  = 12.0 * in;                 // foot [m]
const real mi  = 5280.0 * ft;               // mile [m]
const real gal = 231.0 * in * in * in;      // US gallon [m³]
const real lb  = 0.45359237;                // pound mass [kg]
const real lbf = lb * constant::g0;         // pound force [N]
const real psi = lbf / (in * in);           // pound per square inch [Pa]
const real oz  = gal / 128.0;               // US fluid ounce [m³]
```

**Note:** These are "Freedom Units" as humorously labeled in the source code.

#### Unit Checking

**File:** `fn_check_unit.hpp`

**Location:** `fn_check_unit.hpp:24-28`

```cpp
inline bool check_unit(const value & aValue, const string & aUnit)
{
    return aValue.second == unit_to_si(aUnit).second;
}
```

Validates that a `value` (pair of magnitude and unit array) matches the expected unit dimensions.

---

### Enum Conversion

**File:** `fn_to_enum.hpp`

Generic string-to-enum converter using template metaprogramming.

#### Key Features

- **Generic**: Works with any enum that defines `UNDEFINED` sentinel
- **Case-insensitive**: Uses `string_to_lower()`
- **Safe**: Returns `UNDEFINED` if no match found

#### Implementation

**Location:** `fn_to_enum.hpp:25-48`

```cpp
template <typename T>
void to_enum(const string & aString, T & aEnum)
{
    int tNumEntries = (int) T::UNDEFINED;
    aEnum = T::UNDEFINED;

    string tString = string_to_lower(aString);

    for (int k = 0; k < tNumEntries; ++k)
    {
        aEnum = (T) k;
        string tEnum = string_to_lower(to_string(aEnum));

        if (tString == tEnum)
            return;
    }

    aEnum = T::UNDEFINED;  // No match found
}
```

**Requirements:** The enum type `T` must:
1. Have an `UNDEFINED` member as the last entry
2. Provide a `to_string(T)` overload

---

## String Utilities

### String Tools

**File:** `stringtools.{hpp,cpp}`

Comprehensive string manipulation library.

#### Type Introspection

**Location:** `stringtools.hpp:62-130`

```cpp
template <typename T>
inline std::string datatype_string();   // call as datatype_string<T>()

// Specializations for:
// bool, int, unsigned int, long unsigned int, double, string, complex<double>
```

Returns human-readable type names for debugging and serialization.

#### Path Manipulation

**Location:** `stringtools.hpp:136-162`

```cpp
std::string basename(const std::string & aFilePath);   // e.g., "/path/to/file.txt" → "file.txt"
std::string dirname(const std::string & aFilePath);    // e.g., "/path/to/file.txt" → "/path/to"
std::string filetype(const std::string & aFilePath);   // e.g., "file.txt" → "txt"
std::string filename(const std::string & aFilePath);   // e.g., "/path/to/file.txt" → "file"
```

#### String Processing

**Location:** `stringtools.hpp:169-217`

```cpp
std::string clean_string(const std::string & aString);  // Trim whitespace

string first_word(const std::string & aString, const char aDelimiter = ' ');

Cell<string> string_to_words(const std::string & aString,
                               const char aDelimiter = ' ');

std::string search_and_replace(const std::string & aString,
                                const std::string & aSearch,
                                const std::string & aReplace);

std::string string_to_lower(const std::string & aString);
std::string string_to_upper(const std::string & aString);
```

#### Type Conversion

**Location:** `stringtools.hpp:220-230`

```cpp
bool string_to_bool(const std::string & aString);  // Accepts: 1, on, true, yes
real to_real(const std::string & aString);          // Returns NAN if invalid
value unit_to_si(const string & aString);           // Parse unit strings

bool is_integer(const string & aString);
```

#### Formatting

**Location:** `stringtools.hpp:234-245`

```cpp
string format_with_leading_zeros(const uint aNumber);  // e.g., 42 → "0042"

size_t utf8_character_count(const string & aString);   // Counts UTF-8 characters
```

#### Advanced Parsing

**Location:** `stringtools.hpp:249-308`

```cpp
template <typename T>
void string_to_cell(const string & aString, Cell<T> & aValues);
```

Parses complex numeric expressions into a `Cell<T>`:
- **Comma delimited** (blanks are ignored, so `"1 2 3"` reads as `123`): `"1, 2, 3"` → `{1, 2, 3}`
- **Range notation**: `"1:5"` → `{1, 2, 3, 4, 5}`
- **Reverse ranges**: `"5:1"` → `{5, 4, 3, 2, 1}`
- **Mixed format**: `"1, 3:5, 7"` → `{1, 3, 4, 5, 7}`
- **Braces ignored**: `"{1, 2, 3}"` → `{1, 2, 3}`
- **Semicolon terminator**: `"1, 2 ; ignored"` and `"1,2;ignored"` both → `{1, 2}` (the parser pads `;` into its own token)

---

### Formatted Printing

**File:** `fn_sprint.hpp`

Printf-style formatting to `std::string`.

#### Implementation

**Location:** `fn_sprint.hpp:43-57`

```cpp
template <typename ... Args>
std::string sprint(const char * aFormat, const Args ... aArgs)
{
    // Determine size of formatted string
    auto tSize = std::snprintf(nullptr, 0, aFormat, aArgs ...);

    // Allocate buffer with extra space for '\0'
    std::unique_ptr<char[]> tBuffer(new char[tSize + 1]);

    // Write formatted string into buffer
    std::snprintf(tBuffer.get(), tSize + 1, aFormat, aArgs ...);

    // Return as std::string
    return string(tBuffer.get(), tBuffer.get() + tSize);
}
```

#### Usage

```cpp
std::string msg = sprint("Iteration %d: residual = %.6e", iter, residual);
message(InfoLevel::Default, "Matrix size: %dx%d", nrows, ncols);
```

**Design note:** Compiler warnings for format strings are deliberately suppressed via pragmas (`stringtools.hpp:35-58`) to allow flexible printf-style formatting. Use with care to avoid format/argument mismatches.

---

## Error Handling

### Assertion System

**File:** `assert.{hpp,cpp}`

Custom assertion framework with formatted error messages and ASCII art dragon.

#### Macros

**Location:** `assert.hpp:243-275`

```cpp
// Debug-only assertions (disabled in release builds)
#if !defined(NDEBUG) || defined(DEBUG)
#define BELFEM_ASSERT(aCheck, ...) \
    do { \
        if (!(aCheck)) { \
            belfem::assert::belfem_assert(__FILE__, __LINE__, \
                __PRETTY_FUNCTION__, #aCheck, __VA_ARGS__); \
        } \
    } while (false)
#else
#define BELFEM_ASSERT(aCheck, ...)
#endif

// Always-on error checking
#define BELFEM_ERROR(aCheck, ...) \
    do { \
        if (!(aCheck)) { \
            belfem::assert::belfem_assert(__FILE__, __LINE__, \
                __PRETTY_FUNCTION__, #aCheck, __VA_ARGS__); \
        } \
    } while (false)
```

**Difference:**
- `BELFEM_ASSERT`: Compiled out in release builds (`NDEBUG` defined)
- `BELFEM_ERROR`: Always active, for critical runtime checks

#### Error Formatting and Reaction

**Location:** `assert.hpp` (template `error()`), `assert.cpp` (state)

```cpp
template <typename Exception>
void error(const std::string & aLocation,
           const std::string & aTask,
           const std::string & aCheck,
           const Exception & aException = Exception())
{
    // split the formatted message into lines
    ...

    print_errorbox(aLocation, aTask, aCheck, tMessage);

    // system-log copy of the message (see "System Log Integration")
    if (syslog_on_error())
        log_to_syslog(aLocation, aCheck, tMessage);

    // the reaction is a RUNTIME switch, not an #if: it initializes to
    // the build's compile-time behaviour (throw where assertions are
    // active, abort otherwise, at any rank count) and test executables
    // flip it to make BELFEM_ERROR paths catchable with EXPECT_THROW in
    // release builds
    if (throw_on_error())
        throw aException;

    error_abort();   // MPI_Abort under MPI, std::abort otherwise
}
```

#### System Log Integration (2026-08-15)

Every failed `BELFEM_ERROR` / `BELFEM_ASSERT` also writes its message to
the **system log** — identity `belfem`, priority `LOG_CRIT`, facility
`LOG_USER` — so a crashed run leaves a trace even when the terminal
scrollback or a redirected stderr is gone. Read it back with:

```bash
# Linux (systemd)
journalctl -t belfem              # all BELFEM error events
journalctl -t belfem --since -1h  # recent ones
journalctl -t belfem -f           # follow live while a run is up

# macOS (unified logging — syslog(3) is bridged into it since 10.12)
log show  --predicate 'senderImagePath CONTAINS "belfem" OR eventMessage CONTAINS "rank "' --last 1h
log stream --predicate 'eventMessage CONTAINS "rank "'   # follow live

# syslog-only systems (BSD, older Linux)
grep belfem /var/log/messages
```

`syslog(3)` is POSIX, so the hook compiles and runs on every platform
BELFEM targets — no Linux guard. Only the *retrieval* differs: macOS
routes it into the unified log rather than a text file, and `journalctl`
does not exist there. On macOS the messages are not guaranteed to reach
`/var/log/system.log`; use `log show`.

Each event carries the **MPI rank** in the payload — `LOG_PID` only
identifies local processes, and rank-local checks fire on a single rank
that is usually *not* rank 0. Errors raised before `gComm.init()` log
`rank ?`. Multi-line messages are written one syslog line per message
line, so nothing is lost to the ~1 KiB per-call truncation.

Two properties worth knowing:

- The write happens **before** the reaction branch. This is deliberate:
  debug builds *throw* by default, and a hook on the abort path only
  would never fire in exactly the builds developers run daily.
- Test binaries that exercise error paths on purpose disable it
  (`belfem::assert::set_syslog_on_error(false)`, set in every test main
  beside `set_throw_on_error(true)`), so `make check` does not spam the
  journal.

`log_to_syslog()` runs in normal control flow only — `syslog()` is not
async-signal-safe, so it must never be called from a signal handler.

#### ASCII Art Dragon

**Location:** `assert.cpp` (implementation in `hatch_dragon()`)

When an error occurs, BELFEM prints a decorative error box with an ASCII dragon character, making errors visually distinctive:

```
╔═══════════════════════════════════════════════════════╗
║                      🐉 ERROR 🐉                     ║
╠═══════════════════════════════════════════════════════╣
║ Location: my_file.cpp (line 42)                      ║
║ Task: complete call to function my_function()        ║
║ Check: Assertion index < size failed.                ║
║                                                       ║
║ Index 100 exceeds array size 50                      ║
╚═══════════════════════════════════════════════════════╝
```

#### Helper Functions

**Location:** `assert.hpp:31-58`

```cpp
std::vector<std::string> wrap_lines(std::size_t aMaxWidth, const std::string & aLine);
void hatch_dragon(std::vector<std::string> & aDragon);
void print_line(const std::vector<std::string> & aDragon, std::size_t & aCounter);
void get_lines(const string & aWhat, std::vector<std::string> & aLines);
void print_errorbox(const std::string & aLocation,
                    const std::string & aTask,
                    const std::string & aCheck,
                    const std::vector<std::string> & aMessage);
void error_abort();
std::string extract_function_name(const std::string & aPrettyFunction);
```

#### Usage

```cpp
BELFEM_ASSERT(index < size,
              "Index %lu exceeds array size %lu", index, size);

BELFEM_ERROR(file.is_open(),
             "Failed to open file: %s", filename.c_str());
```

---

## Hashing

### Hash Class

**File:** `cl_Hash.hpp`

Incremental hash computation with entropy mixing.

#### Algorithm

**Location:** `cl_Hash.hpp:68-83`

```cpp
template <typename T>
Hash & operator += (const T aValue)
{
    mValue +=
        // Compute hash of the new value
        std::hash<T>{}(aValue)
        // Add entropy constant from golden ratio for better distribution
        + 0x9e3779b97f4a7c15ULL
        // Mix in current hash with left shift (emphasizes high-order bits)
        + (mValue << 6)
        // Mix in current hash with right shift (brings in low-order bits)
        + (mValue >> 2);

    return *this;
}
```

**Entropy constant:** `0x9e3779b97f4a7c15` is derived from the golden ratio `φ = (1+√5)/2`:
```
φ × 2^64 ≈ 11400714819323198485 = 0x9e3779b97f4a7c15
```

This value is commonly used in hash functions (Knuth, boost::hash_combine) for good avalanche properties.

#### API

**Location:** `cl_Hash.hpp:30-80`

```cpp
class Hash
{
    std::size_t mValue = 0;

public:
    Hash() = default;
    ~Hash() = default;

    void reset();                              // Set to 0
    std::size_t value() const;                 // Get current hash
    void set_value(const std::size_t aValue);  // Set hash directly

    template <typename T>
    Hash & operator += (const T aValue);       // Incremental hashing
};
```

#### Usage Pattern

```cpp
Hash hash;
hash += node_id;
hash += x_coordinate;
hash += y_coordinate;
hash += z_coordinate;

std::size_t key = hash.value();  // Use as map key, cache lookup, etc.
```

**Design rationale:** The incremental interface allows computing hashes of composite objects field-by-field, useful for geometric point hashing, element signatures, etc.

**⚠️ Collision Note:** `Hash` is **not cryptographic** and does not guarantee collision-free behavior. Use only for caching, lookup heuristics, and non-security-critical map keys. For security-sensitive applications, use a proper cryptographic hash.

---

## Random Number Generation

### Random Number System

**File:** `random.hpp`

MPI-aware random number generation.

#### Key Features

- **MPI-aware**: Uses `gComm.random()` generator in MPI builds
- **Serial builds**: `random_seed()` reads a seed from `/dev/urandom` (or the clock) and seeds `std::rand()`, which `rand()` uses
- **C++ `<random>` API**: Modern random number generation

#### Seeding

**Location:** `random.hpp:26-72`

```cpp
#ifndef BELFEM_MPI
template <typename T>
void random_seed(T & aSeed)
{
    std::ifstream tStream("/dev/urandom", std::ios::binary);

    if (tStream)
    {
        char tMemblock[sizeof(T)];
        tStream.read(tMemblock, sizeof(T));
        tStream.close();
        aSeed = *reinterpret_cast<T*>(tMemblock);
    }
    else
    {
        aSeed = (T) time(NULL);  // Fallback to system clock
    }
}
#endif

inline void random_seed()
{
#ifdef BELFEM_MPI
    std::random_device rd;
    gComm.random().seed(rd());  // Seed once using random device
#else
    unsigned int tSeed;
    random_seed(tSeed);
#endif
}
```

#### Random Number Generation

**Location:** `random.hpp:80-89`

```cpp
inline real rand()
{
#ifdef BELFEM_MPI
    std::uniform_real_distribution<real> distribution(0.0, 1.0);
    return distribution(gComm.random());  // Uses MPI communicator's RNG
#else
    return ((real) std::rand()) / RAND_MAX;
#endif
}
```

#### Usage Pattern

```cpp
#include "random.hpp"

// Seed at program start
random_seed();

// Generate random numbers
for (int i = 0; i < 1000; ++i)
{
    real value = rand();  // Uniform distribution [0, 1)
}
```

**MPI considerations:** In MPI builds, each process has its own independent random stream via `gComm.random()`, ensuring reproducibility and avoiding correlation between processes.

**⚠️ Reproducibility:**
- MPI builds guarantee **independent streams per rank**, not identical sequences across ranks
- For deterministic reproducibility across runs, use explicit seeding via communicator logic with a fixed seed
- Serial builds are random run to run as well once `random_seed()` has been called; without that call `std::rand()` keeps its default seed

---

## Global Variables

### Global Variable System

**File:** `globals.hpp`

Framework for managing global variables with MPI synchronization.

#### Macro Magic

**Location:** `globals.hpp:18-22`

```cpp
#ifdef BELFEM_INITIALIZE_GLOBALS  // Only set from within Communicator
#define greal real
#define gstring string
#else
#define greal extern real
#define gstring extern string
#endif
```

**Rationale:** This pattern ensures global variables are:
1. **Declared** as `extern` in normal includes (no storage allocated)
2. **Defined** only when `Communicator::set_globals()` is called

This avoids multiple definition errors while allowing global access.

#### Current Globals

**Location:** `globals.hpp:41-53`

```cpp
namespace belfem
{
    greal   gTbulk;          // bulk temperature in K when no thermal kernel is chosen; BELFEM_QUIET_NAN until set
    greal   gRhoMin;         // minimum resistivity in Ohm*m, default 0
    greal   gRhoMax;         // maximum resistivity in Ohm*m, default 1e10
    gstring gBelfemDataPath; // data path, from $BELFEM_DATA, else the installed share directory
}
```

#### Usage Guidelines

**⚠️ CRITICAL WARNING:**
Globals should be used **only for true global physical parameters** (like bulk temperature, reference pressure) that are genuinely shared across the entire simulation. **Do not use globals for:**
- Algorithmic state or control flow
- Caches or temporary storage
- Mutable iteration counters
- Solver-specific configuration

Misuse of globals creates hidden dependencies, breaks modularity, and causes subtle MPI synchronization bugs.

From the source code comments (`globals.hpp:24-38`):

1. **Adding new parameters:**
   - Add declaration in `globals.hpp`
   - Initialize in `Communicator::set_globals()`
   - Use `BELFEM_QUIET_NAN` as default if no clear initial value
   - **Justify the need** - is this truly a framework-wide physical constant?

2. **MPI synchronization:**
   - The `BELFEM_INITIALIZE_GLOBALS` macro must **only** be defined inside `Communicator::set_globals()`
   - Defining it elsewhere causes multiple-definition linker errors
   - If value is set only on master process, it is the developer's responsibility to broadcast:
     ```cpp
     gMyVariable = compute_on_master();
     broadcast(gMyVariable);  // free function from commtools.hpp, collective, root = 0
     ```

#### Example

```cpp
// In physics module
namespace belfem
{
    extern real gTbulk;  // Declared (via macro)
}

// In Communicator::set_globals()
#define BELFEM_INITIALIZE_GLOBALS
#include "globals.hpp"

void Communicator::set_globals()
{
    gTbulk = BELFEM_QUIET_NAN;  // Initialized
}

// In application code
gTbulk = 300.0;  // Set bulk temperature
broadcast(gTbulk);  // free function from commtools.hpp, collective, root = 0
```

---

## Thread Safety and MPI

### Thread Safety Summary

Core utilities have varying thread-safety characteristics:

| Component | Thread-Safe (Read-Only) | Requires External Sync |
|-----------|------------------------|------------------------|
| `Timer` | ✓ (each thread uses own instance) | N/A |
| `Hash` | ✓ (each thread uses own instance) | N/A |
| `unit_to_si()` | ✓ (const function) | N/A |
| `string_to_bool()` | ✓ (const function) | N/A |
| `all stringtools functions` | ✓ (stateless) | N/A |
| **`Logger`** | ✗ | Needs mutex for concurrent writes |
| **`Progressbar`** | ✗ | Single-thread use only |
| **`Profiler`** | ✗ | Single-thread use only |
| **`rand()`** | ✗ | Needs per-thread RNG in OpenMP |
| **`gLog`, `gComm`, globals** | ✗ | External sync required |

### MPI Behavioral Differences

Several utilities behave differently in MPI vs serial builds:

| Function | Non-MPI Behavior | MPI Behavior |
|----------|------------------|--------------|
| `rand()` | Uses `std::rand()` | Uses `gComm.random()` per-rank generator |
| `random_seed()` | Seeds `std::rand()` from `/dev/urandom` or `time()` | Seeds from `std::random_device` |
| `error_abort()` | Calls `std::abort()` | Calls `MPI_Abort(MPI_COMM_WORLD)` |
| `print_banner()` | Always prints | Only rank 0 prints (typical usage) |
| `Progressbar` | `stdout` is a line-buffered terminal | `stdout` is a fully buffered pipe; the bar flushes each frame itself |
| `Profiler` | Single file | Separate file per rank (e.g., `profiler.4.2.log` and `profiler.4.2.callgrind` for size=4, rank=2) |

### Best Practices

1. **Logger in MPI**: Typically construct on rank 0 only, or use rank-specific log files
2. **Progressbar in MPI**: Use only on rank 0 to avoid garbled output. The bar flushes every frame, so it streams through the `mpirun` pipe instead of arriving in one burst at exit
3. **Random in OpenMP**: Create thread-local RNG instances instead of calling global `rand()`
4. **Globals**: Always broadcast after modification on master rank:
   ```cpp
   if (gComm.rank() == 0) {
       gTbulk = compute_value();
   }
   broadcast(gTbulk);  // free function from commtools.hpp, collective, root = 0
   ```

---

## Common Pitfalls

### Logger and Output

**Pitfall:** Forgetting to check `info_level()` before expensive string formatting
```cpp
// Bad: Formats even if Silent
message(InfoLevel::Verbose, "Expensive: %s", compute_expensive_string().c_str());

// Good: Check level first
if (gLog.info_level() >= static_cast<uint>(InfoLevel::Verbose)) {
    message(InfoLevel::Verbose, "Expensive: %s", compute_expensive_string().c_str());
}
```

**Pitfall:** Calling `Progressbar::step()` after `finish()` is a silent no-op — `finish()` leaves the
bar at full width, and `step()` only redraws when the width grows. Restart the bar with `reset()`
if you need it again.

**Pitfall:** Long error messages exceeding terminal width - use `wrap_lines()` handles this automatically in assert.cpp

###  Timer and Profiling

**Pitfall:** Using `Timer::next()` for long-running timers (>49 days) causes truncation

```cpp
// Bad: Truncates after 49 days
uint64_t elapsed = timer.next();

// Good: Use stop() for long durations
uint64_t elapsed = timer.stop();
```

**Pitfall:** Forgetting to enable profiler at build time - `Profiler::start()/stop()` do nothing without `-DUSE_PROFILER=ON`

### String Parsing

**Pitfall:** Range endpoints are parsed with `std::stoll` and cast to `T`
```cpp
string_to_cell("5:1", nodes);   // {5, 4, 3, 2, 1} — decreasing ranges are fine
string_to_cell("-1:1", nodes);  // a negative endpoint in an unsigned Cell<index_t> wraps silently
```

### Random Numbers

**Pitfall:** Not seeding RNG leads to deterministic sequences
```cpp
// Bad: Forgot to seed
for (int i = 0; i < 100; ++i) {
    real r = belfem::rand();  // Same sequence every run!
}

// Good: Seed at program start
belfem::random_seed();
for (int i = 0; i < 100; ++i) {
    real r = belfem::rand();  // Different each run in an MPI build; serial std::rand() stays unseeded
}
```

**Pitfall:** Assuming identical random sequences across MPI ranks - each rank has independent stream

### Type System

**Pitfall:** Mixing epsilon constants
```cpp
// Bad: Using wrong epsilon for geometry
real distance = node_a.distance_to(node_b);
if (distance < BELFEM_EPSILON) { /* coincident */ }  // Wrong! Too strict

// Good: Use mesh epsilon for geometry
if (distance < BELFEM_MESH_EPSILON) { /* coincident */ }
```

**Pitfall:** Progressbar constructor type mismatch on 64-bit builds
```cpp
// Potential issue: aNumSteps is uint, but mNumSteps is index_t
Progressbar bar(very_large_mesh_size);  // May truncate on BELFEM_INT64 builds
```

### Error Handling

**Pitfall:** Using `BELFEM_ASSERT` for runtime checks
```cpp
// Bad: Disabled in release builds!
BELFEM_ASSERT(file.is_open(), "File open failed");

// Good: Use BELFEM_ERROR for critical runtime checks
BELFEM_ERROR(file.is_open(), "File open failed");
```

**Pitfall:** Format string/argument mismatch in assertions
```cpp
// Bad: Mismatched types
BELFEM_ERROR(false, "Value: %d", some_double);  // Undefined behavior

// Good: Correct format specifiers
BELFEM_ERROR(false, "Value: %.3f", some_double);
```

### Global Variables

**Pitfall:** Defining `BELFEM_INITIALIZE_GLOBALS` outside `Communicator::set_globals()` causes linker errors

**Pitfall:** Using globals for algorithm state creates hidden dependencies and breaks modularity

**Pitfall:** Forgetting to broadcast globals after master-only computation
```cpp
// Bad: Only master has correct value
if (gComm.rank() == 0) {
    gTbulk = 300.0;
}
// Other ranks still have NaN!

// Good: Explicit broadcast
if (gComm.rank() == 0) {
    gTbulk = 300.0;
}
broadcast(gTbulk);   // free function from commtools.hpp, collective, root = 0
```

---

## Usage Examples

### Example 1: Basic Application Setup

```cpp
#include "cl_Logger.hpp"
#include "cl_Timer.hpp"
#include "banner.hpp"
#include "cl_Arguments.hpp"

belfem::Logger gLog(belfem::InfoLevel::Default);   // file scope, as in src/executables/belfem.cpp

int main(int argc, char * argv[])
{
    // Parse command line ( may raise gLog's level via -v, so the logger must exist first )
    belfem::Arguments args(argc, argv);

    // Print banner
    belfem::print_banner("MyApp");

    // Time the application
    belfem::Timer timer;

    message(belfem::InfoLevel::Default, "Starting computation...");

    // ... work ...

    uint64_t elapsed = timer.stop();
    message(belfem::InfoLevel::Default, "Completed in %lu ms", elapsed);

    return 0;
}
```

---

### Example 2: Profiling a Critical Section

```cpp
#include "cl_Profiler.hpp"
#include "cl_Timer.hpp"

void optimize_mesh()
{
    belfem::Profiler profiler("mesh_optimization.log");
    belfem::Timer timer;

    profiler.start();

    // ... expensive mesh optimization ...

    profiler.stop();

    message(belfem::InfoLevel::Verbose,
            "Mesh optimization: %lu ms", timer.stop());

    // Analyze with: callgrind_annotate mesh_optimization.callgrind
}
```

---

### Example 3: Progress Tracking

```cpp
#include "cl_Progressbar.hpp"

void assemble_system(uint num_elements)
{
    belfem::Progressbar bar(num_elements);

    message(belfem::InfoLevel::Default, "Assembling system matrix...");

    for (uint e = 0; e < num_elements; ++e)
    {
        // Assemble element contribution
        assemble_element(e);

        bar.step();
    }

    bar.finish();
}
```

---

### Example 4: Error Handling

```cpp
#include "assert.hpp"

void compute_jacobian(const Matrix & J)
{
    real det = J.determinant();

    BELFEM_ERROR(det > BELFEM_MESH_EPSILON,
                 "Singular Jacobian detected: det(J) = %.3e", det);

    real inv_det = 1.0 / det;

    // ... use inv_det ...
}
```

---

### Example 5: String Parsing

```cpp
#include "stringtools.hpp"

void parse_node_list(const string & input)
{
    // Input: "1:5, 10, 15:17"
    // Output: {1, 2, 3, 4, 5, 10, 15, 16, 17}

    Cell<index_t> nodes;
    string_to_cell(input, nodes);

    for (index_t node_id : nodes)
    {
        process_node(node_id);
    }
}
```

---

### Example 6: Hash-Based Caching

```cpp
#include "cl_Hash.hpp"

std::size_t compute_element_signature(const Element * element)
{
    belfem::Hash hash;

    // Hash element type
    hash += static_cast<uint>(element->type());

    // Hash node IDs
    for (uint i = 0; i < element->number_of_nodes(); ++i)
    {
        hash += element->node(i)->id();
    }

    return hash.value();
}

// Use as cache key
std::unordered_map<std::size_t, Matrix> shape_function_cache;

if (shape_function_cache.find(signature) == shape_function_cache.end())
{
    // Compute and cache shape functions
    shape_function_cache[signature] = compute_shape_functions();
}
```

---

### Example 7: Unit Conversion

```cpp
#include "units.hpp"
#include "constants.hpp"

void set_pressure_boundary_condition()
{
    using namespace belfem::constant;

    real pressure_psi = 14.7;  // Standard atmosphere in psi
    real pressure_si = pressure_psi * psi;  // Convert to Pa

    BELFEM_ASSERT(std::abs(pressure_si - atm) < 1.0,
                  "Conversion error: expected %f Pa, got %f Pa", atm, pressure_si);

    apply_pressure(pressure_si);
}
```

---

### Example 8: Random Mesh Perturbation

```cpp
#include "random.hpp"

void perturb_mesh(Mesh * mesh, real amplitude)
{
    belfem::random_seed();

    for (uint n = 0; n < mesh->number_of_nodes(); ++n)
    {
        Node * node = mesh->node(n);

        // Random perturbation in [-amplitude, +amplitude]
        real dx = amplitude * (2.0 * belfem::rand() - 1.0);
        real dy = amplitude * (2.0 * belfem::rand() - 1.0);
        real dz = amplitude * (2.0 * belfem::rand() - 1.0);

        node->x() += dx;
        node->y() += dy;
        node->z() += dz;
    }
}
```

---

## Summary

The `core` module provides essential infrastructure for BELFEM applications:

| Category | Components | Purpose |
|----------|-----------|---------|
| **Logging** | Logger, Progressbar, Banner | User communication and debugging |
| **Performance** | Timer, Profiler | Timing and profiling |
| **CLI** | Arguments | Command-line processing |
| **Types** | typedefs, constants, units | Type system and physical constants |
| **Strings** | stringtools, sprint | String manipulation and formatting |
| **Errors** | assert macros | Debugging and error handling |
| **Utilities** | Hash, random, globals | Common algorithmic utilities |

### Design Philosophy

1. **Minimal dependencies**: Core utilities should not depend on higher-level BELFEM modules
2. **Type safety**: Strong typing via custom typedefs (index_t, id_t, etc.)
3. **Physical correctness**: NIST-compliant constants and SI unit system
4. **Developer experience**: Informative error messages, progress feedback, profiling tools
5. **Cross-platform**: Compiler-agnostic with appropriate pragma guards

### Further Reading

- **Build system**: See `CLAUDE.md` for CMake build instructions
- **MPI integration**: See `src/comm/` for parallel communication
- **Container library**: See `src/containers/cl_Cell.hpp` for dynamic arrays
- **Mesh infrastructure**: See `src/mesh/` for mesh data structures that use core types

---

**Document version:** 1.0
**Last updated:** 2026-01-16
**Author:** Claude Code (claude.ai/code)
