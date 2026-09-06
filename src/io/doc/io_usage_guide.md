# BELFEM I/O Usage Guide {#io_io_usage_guide}

Comprehensive guide to BELFEM's file input/output system.

**Date:** 2026-01-16
**Purpose:** Usage documentation for src/io module
**Revision History:**
- 2026-01-16: Initial documentation
- 2026-01-16: Critical corrections (FileMode::NEW behavior, resource management, parallel I/O contracts)

---

## Table of Contents

1. [Overview](#overview)
2. [Common Pitfalls](#common-pitfalls)
3. [Core Classes](#core-classes)
   - [HDF5](#hdf5-class)
   - [Ascii](#ascii-class)
   - [CsvFile](#csvfile-class)
   - [InputFile and Input_Section](#inputfile-and-input_section-classes)
   - [XML](#xml-class)
4. [File Utilities](#file-utilities)
5. [Common Patterns](#common-patterns)
6. [Performance Tips](#performance-tips)
7. [MPI Considerations](#mpi-considerations)
8. [Thread Safety](#thread-safety)
9. [Debugging Tips](#debugging-tips)
10. [Enumeration Reference](#enumeration-reference)

---

## Overview

BELFEM's I/O module provides unified interfaces for reading and writing data in multiple formats:

- **HDF5**: Hierarchical data format for large-scale scientific data (requires `BELFEM_HDF5`)
- **ASCII**: Line-based text file interface with buffering
- **CSV**: Comma/delimiter-separated value files (numeric data → `Matrix<real>`)
- **InputFile**: Hierarchical configuration file parser (key-value pairs, nested sections)
- **XML**: XML document interface (requires `BELFEM_XML`)

### Build Configuration

```cmake
# Enable HDF5 support (default ON; highly recommended for large datasets)
cmake -DUSE_HDF5=ON ..

# Enable XML support via tinyxml2 (default OFF)
cmake -DUSE_TINYXML2=ON ..
```

(`BELFEM_HDF5` is the *compile definition* those options generate, not a
CMake cache entry — setting it on the command line does nothing.)

### Design Philosophy

Following BELFEM's coding philosophy:

- **Manual resource management**: Files use RAII (constructor opens, destructor closes)
- **Explicit control**: FileMode explicitly controls read/write intent
- **Performance-first**: Buffering in Ascii, direct memory copy in HDF5
- **Zero-overhead abstractions**: Template functions for type-safe I/O
- **C-style error handling**: a failed check ends the run — abort in release, throw in debug (see Error Handling Philosophy)

### Error Handling Philosophy

**IMPORTANT:** BELFEM uses **C libraries** (HDF5, BLAS, etc.) with **abort-on-error** semantics, NOT C++ exceptions.

```cpp
// Typical BELFEM error handling pattern (hdf5_tools.hpp:596-623)
T* tData = (T*) malloc(tLength * sizeof(T));     // 1. Allocate
std::memcpy(tData, aVector.data(), ...);         // 2. Use (C function, never throws)
aStatus = H5Dwrite(..., tData);                  // 3. C API (returns error code, never throws)
free(tData);                                      // 4. FREE ALWAYS EXECUTES
BELFEM_ASSERT(aStatus == 0, "...");              // 5. Abort if error (after cleanup)
```

**Why this matters:**
- HDF5 functions **return error codes** (`herr_t`), they **never throw exceptions**
- `malloc`/`free` are C functions (don't throw in C++)
- `BELFEM_ERROR` and `BELFEM_ASSERT` end the run — but **how** depends on the build: a debug
  build **throws**, a release build calls `error_abort()` (`assert.hpp:88-93,189-192`)
- So the "no unwinding" simplification holds for **release** builds only. In a debug build the
  throw does unwind, and a raw `malloc` between the allocation and the check is leaked

**Resource cleanup guarantees:**
1. **Normal path:** Destructors run when objects go out of scope (RAII)
2. **Error path:** Explicit cleanup (e.g., `free()`) runs **before** `BELFEM_ASSERT` aborts
3. **Abort/crash:** Process terminates → OS reclaims all memory, file descriptors, handles

**Implication:** in a release build, traditional C++ "exception safety" concerns (RAII wrappers,
`unique_ptr`) do not drive the design of these error paths — the process ends after the explicit
cleanup rather than unwinding through it. Do not read that as "BELFEM never throws": a debug run
does, deliberately and at every rank, so that a failure can be caught under a debugger. The
pattern above stays correct in both builds because the `free()` precedes the check.

> **Note for contributors:** If you add C++ library code that throws exceptions (e.g., `std::vector::at()`), wrap it in `try/catch` and convert to `BELFEM_ERROR`. Do NOT let exceptions propagate through C library boundaries.

---

## Common Pitfalls

### 0. Modern C++ Tool Warnings (False Positives for BELFEM)

**IMPORTANT:** Static analyzers and modern C++ code reviewers may flag BELFEM code with warnings like:

- "Use `std::unique_ptr` instead of raw pointers"
- "Potential memory leak: `malloc` without exception-safe cleanup"
- "Use RAII wrapper for file handles"
- "Exception-unsafe resource management"

**These warnings DO NOT apply to BELFEM** because:

1. **BELFEM uses C libraries** (HDF5, BLAS, etc.) that return error codes, not exceptions
2. **Error handling ends the run** (`BELFEM_ERROR`/`BELFEM_ASSERT`: abort in release, throw in debug)
3. In release there is no unwinding, so cleanup placed before the check always runs; in debug the throw unwinds, which is why the `free()` must precede the check
4. **Manual resource management is intentional** for zero-overhead performance in HPC

**Example of "false positive":**
```cpp
// Static analyzer warning: "Memory leak if H5Dwrite throws"
T* tData = (T*) malloc(tLength * sizeof(T));
aStatus = H5Dwrite(..., tData);  // C API - NEVER throws
free(tData);  // ALWAYS executes (before BELFEM_ASSERT)
BELFEM_ASSERT(aStatus == 0, "...");  // Aborts AFTER cleanup
```

This is **safe in BELFEM** because:
- `H5Dwrite()` is a C function (returns error code, never throws)
- `free()` executes before the assertion
- If assertion fails, process aborts (OS reclaims all resources)

**When to ignore modern C++ warnings:**
- ✅ `malloc`/`free` with C library calls (HDF5, BLAS, etc.)
- ✅ Raw pointers owned by RAII classes (deleted in destructor)
- ✅ Manual cleanup before `BELFEM_ERROR`/`BELFEM_ASSERT`

**When warnings ARE valid:**
- ❌ Calling C++ library functions that throw (e.g., `std::vector::at()`) between resource allocation and cleanup
- ❌ Resources not cleaned up in EITHER destructor OR before error abort
- ❌ Mixing exception-throwing C++ code with manual memory management

> **See "Error Handling Philosophy"** section above for detailed explanation of BELFEM's C-style error handling model.

### 1. Forgetting to Enable HDF5/XML at Build Time

```cpp
// WRONG: Will compile but fail at runtime if BELFEM_HDF5 not enabled
HDF5 file("data.h5", FileMode::NEW);  // May produce stubs or errors

// SOLUTION: Check CMake configuration
#ifdef BELFEM_HDF5
    HDF5 file("data.h5", FileMode::NEW);
#else
    BELFEM_ERROR(false, "BELFEM_HDF5 not enabled. Reconfigure with -DUSE_HDF5=ON");
#endif
```

### 2. Using Wrong FileMode (CRITICAL - Data Loss Hazard!)

```cpp
// DANGER: FileMode::NEW silently truncates existing files!
HDF5 file("valuable_results.h5", FileMode::NEW);  // DESTROYS existing data!

// SAFE: Check existence first for production workflows
if (file_exists("results.h5")) {
    BELFEM_ERROR(false, "File already exists. Use OPEN_RDWR or choose new name.");
}
HDF5 file("results.h5", FileMode::NEW);

// CORRECT: Use OPEN_RDWR for modifying existing files
HDF5 file("existing_results.h5", FileMode::OPEN_RDWR);

// SAFE: Use versioned or timestamped names for checkpoints
string path = "checkpoint_" + std::to_string(timestep) + ".h5";
HDF5 file(path, FileMode::NEW);
```

**FileMode options (cl_HDF5.cpp:44-97):**
- `NEW`: Create new file (**TRUNCATES if exists** — `H5F_ACC_TRUNC` — use with caution!)
- `OPEN_RDONLY`: Read-only access (`H5F_ACC_RDONLY`, fails if doesn't exist)
- `OPEN_RDONLY_PARALLEL`: **not handled by `HDF5`** — the constructor has no case for it and falls through to "unknown filemode passed". `Ascii` implements it (rank 0 reads, then broadcasts)
- `OPEN_RDWR`: Read-write access (`H5F_ACC_RDWR`, fails if doesn't exist)

> **Warning:** `FileMode::NEW` will **destroy existing results** if used with automatically generated paths.
> For restarts, always use `OPEN_RDWR`. For checkpoints, use explicit versioning.

### 3. HDF5 Group Management Errors (Resource Leaks!)

```cpp
// WRONG: Forgetting to close groups causes HDF5 handle leaks
HDF5 file("data.h5", FileMode::NEW);
file.create_group("Results");
file.save_data("GlobalMetadata", metadata);  // Writes to /Results/GlobalMetadata!

// CORRECT: Close group before operating at parent level
HDF5 file("data.h5", FileMode::NEW);
file.create_group("Results");
file.save_data("field", vector);  // Writes to /Results/field
file.close_active_group();   // Back to root /
file.save_data("GlobalMetadata", metadata);  // Writes to /GlobalMetadata

// BEST: Use RAII scope guards (destructor closes automatically)
{
    HDF5 file("data.h5", FileMode::NEW);
    file.create_group("Results");
    file.save_data("Temperature", T);
    file.save_data("Pressure", P);
    file.close_active_group();
}  // File and the active group closed by destructor
```

**HDF5 Group Navigation Invariant (cl_HDF5.cpp:78-239):**

> **Critical Rule:** HDF5 maintains exactly **one active group** at any time.
> All `save_data()` and `load_data()` calls operate **relative to the currently active group**.
>
> - `create_group(name)` opens the group (calls `H5Gcreate2`) and makes it active
> - `select_group(name)` opens an existing group (`H5Gopen2`) and makes it active
> - `close_active_group()` closes the current group (`H5Gclose`) and returns to parent
> - Destructor automatically calls `close()` → `close_active_group()` (RAII cleanup)
>
> **Resource Management:**
> - Each opened group allocates an HDF5 handle (`hid_t`)
> - Handles are **not** automatically closed when opening a new group
> - **Failure to close groups** before program exit → handle leaks (finite limit!)
> - The destructor calls `close()`, which closes only the **active** group before `H5Fclose`; `close_tree()` exists but is not called automatically. Close every group you opened, or call `close_tree()` yourself before the object dies.
>
> **Recommendation:** Always match `create_group()` / `select_group()` with `close_active_group()`.
> Rely on destructor for final cleanup, but explicitly close groups to maintain clarity.

### 4. InputFile Section Hierarchy Confusion

```cpp
// Config file: Section { subsection { key = value } }

// WRONG: Trying to access nested section directly
const input::Section* sec = input_file.section("subsection");  // BELFEM_ERROR: "Section ... does not exist"

// CORRECT: Navigate hierarchy
const input::Section* parent = input_file.section("Section");
const input::Section* child = parent->section("subsection");
real value = child->get_real("key");
```

### 5. CsvFile Assumes Numeric Data

```cpp
// data.csv contains: "Name,Age,City"
CsvFile csv("data.csv");  // WRONG: Will fail parsing non-numeric data

// SOLUTION: CsvFile only handles numeric data → Matrix<real>
// For mixed types, use Ascii and parse manually
Ascii file("data.csv", FileMode::OPEN_RDONLY);
for (index_t i = 0; i < file.length(); ++i) {
    string line = file.line(i);
    // Custom parsing logic
}
```

---

## Core Classes

### HDF5 Class

**Location**: `cl_HDF5.hpp`

Hierarchical data format interface for saving/loading scalars, vectors, matrices, and structured data.

#### Constructor

```cpp
HDF5(const string & aPath,
     const enum FileMode aMode,
     const bool aParallelMode = false);
```

**Parameters:**
- `aPath`: File path (`.h5` or `.hdf5` extension recommended)
- `aMode`: File access mode (`NEW`, `OPEN_RDONLY`, `OPEN_RDWR`; `OPEN_RDONLY_PARALLEL` is **not** handled by `HDF5` — see the FileMode table)
- `aParallelMode`: write **one file per rank** — the path is rewritten through `make_path_parallel()`. This is not PHDF5 and nothing about it is collective

#### Group Operations

```cpp
// Create and navigate hierarchical structure
hid_t create_group(const string & aLabel);  // returns the new active group
hid_t select_group(const string & aLabel);
void close_active_group();
```

**Example:**
```cpp
HDF5 file("simulation.h5", FileMode::NEW);

// Create nested structure: /TimeStep_001/Fields/Temperature
file.create_group("TimeStep_001");
file.create_group("Fields");
file.save_data("Temperature", temperature_field);  // Saved at /TimeStep_001/Fields/Temperature
file.close_active_group();  // Back to /TimeStep_001
file.close_active_group();  // Back to root /
```

#### Save Operations

```cpp
// Scalars
void save_data(const string & aLabel, const string & aValue);
void save_data(const string & aLabel, const sint & aValue);
void save_data(const string & aLabel, const uint & aValue);
void save_data(const string & aLabel, const luint & aValue);
void save_data(const string & aLabel, const lluint & aValue);
void save_data(const string & aLabel, const real & aValue);
void save_data(const string & aLabel, const bool & aValue);

// Vectors
void save_data(const string & aLabel, const Vector<sint> & aVector);
void save_data(const string & aLabel, const Vector<uint> & aVector);
void save_data(const string & aLabel, const Vector<luint> & aVector);
void save_data(const string & aLabel, const Vector<real> & aVector);

// Matrices
void save_data(const string & aLabel, const Matrix<sint> & aMatrix);
void save_data(const string & aLabel, const Matrix<uint> & aMatrix);
void save_data(const string & aLabel, const Matrix<luint> & aMatrix);
void save_data(const string & aLabel, const Matrix<real> & aMatrix);

// String arrays
void save_data(const string & aLabel, const Cell<string> & aStrings);
```

#### Load Operations

```cpp
// Scalars
void load_data(const string & aLabel, string & aValue);
void load_data(const string & aLabel, sint & aValue);
void load_data(const string & aLabel, uint & aValue);
void load_data(const string & aLabel, luint & aValue);
void load_data(const string & aLabel, lluint & aValue);
void load_data(const string & aLabel, real & aValue);
void load_data(const string & aLabel, bool & aValue);

// Vectors (auto-resizes to file data size)
void load_data(const string & aLabel, Vector<sint> & aVector);
void load_data(const string & aLabel, Vector<uint> & aVector);
void load_data(const string & aLabel, Vector<luint> & aVector);
void load_data(const string & aLabel, Vector<real> & aVector);

// Matrices (auto-resizes to file data size)
void load_data(const string & aLabel, Matrix<sint> & aMatrix);
void load_data(const string & aLabel, Matrix<uint> & aMatrix);
void load_data(const string & aLabel, Matrix<luint> & aMatrix);
void load_data(const string & aLabel, Matrix<real> & aMatrix);

// String arrays
void load_data(const string & aLabel, Cell<string> & aStrings);
```

#### Best Practices

```cpp
// GOOD: Explicit group management
{
    HDF5 file("results.h5", FileMode::NEW);

    // Organized structure
    file.create_group("Mesh");
    file.save_data("Nodes", node_coords);
    file.save_data("Elements", element_connectivity);
    file.close_active_group();

    file.create_group("Solution");
    file.save_data("Temperature", T);
    file.save_data("Pressure", P);
    file.close_active_group();
}  // File automatically closed by destructor
```

### Ascii Class

**Location**: `cl_Ascii.hpp`

Line-based ASCII file interface with an in-memory line buffer; changes must be written with `save()` explicitly.

#### Constructor

```cpp
Ascii(const string & aPath, const enum FileMode & aMode);
```

**Modes supported:**
- `NEW`: empty buffer, `save()` writes the file
- `OPEN_RDONLY`: read existing file into the buffer
- `OPEN_RDONLY_PARALLEL`: rank 0 reads, the buffer is broadcast
- `OPEN_RDWR`: not supported by `Ascii` (constructor error)

#### Interface

```cpp
// Access lines (0-indexed)
const string & line(index_t aLineNumber) const;
string & line(index_t aLineNumber);  // Mutable access

// Number of lines
index_t length() const;

// Append line to buffer (cl_Ascii.cpp:135)
void print(const string & aLine);

// Write buffer to file (cl_Ascii.cpp:76)
bool save();
```

#### Usage Patterns

**Read text file:**
```cpp
Ascii file("input.txt", FileMode::OPEN_RDONLY);

for (index_t i = 0; i < file.length(); ++i) {
    const string & line = file.line(i);  // No copy

    // Parse line
    if (line.find("PARAMETER") != string::npos) {
        // Extract parameter
    }
}
```

**Write text file (cl_Ascii.cpp:135-140):**
```cpp
Ascii output("results.txt", FileMode::NEW);

// Append lines to buffer
output.print("# Simulation results");
output.print("Time, Temperature, Pressure");
for (index_t i = 0; i < n; ++i) {
    output.print(std::to_string(time(i)) + ", " +
                 std::to_string(temp(i)) + ", " +
                 std::to_string(pres(i)));
}

// Write buffer to file
output.save();
```

**Memory considerations:**
- Entire file loaded into memory (`Cell<string>` buffer)
- Efficient for moderate-sized text files (< 10 MB)
- For large files (> 100 MB), consider streaming with standard library
- **Destructor:** raises a `BELFEM_ERROR` (always active) if the buffer changed but was never saved (`mChangedSinceLastSave`)

### CsvFile Class

**Location**: `cl_CsvFile.hpp`

CSV reader that loads numeric data into `Matrix<real>`. Extends `Ascii`.

#### Constructor

```cpp
CsvFile(const string & aPath, const char aDelimiter = ',');
```

**Parameters:**
- `aPath`: Path to CSV file
- `aDelimiter`: Column separator (default: `,`)

#### Interface

```cpp
const Matrix<real> & data() const;  // Read-only access
Matrix<real> & data();               // Mutable access
```

#### Usage

```cpp
// Load CSV with comma delimiter
CsvFile csv("data.csv");
const Matrix<real> & data = csv.data();

// data(row, col) contains numeric values
real value = data(0, 2);  // First row, third column

// Tab-separated values
CsvFile tsv("data.tsv", '\t');
```

**Assumptions:**
- All data is numeric (`real` type)
- Rectangular data (all rows same length)
- Non-numeric content will cause parsing errors

### InputFile and Input_Section Classes

**Location**: `cl_InputFile.hpp`, `cl_Input_Section.hpp`

Hierarchical configuration file parser with key-value pairs and nested sections.

#### File Format

```
// comments start with a double slash

// simple key-value pairs: key : value ;
title : my simulation ;
timesteps : 100 ;
dt : 0.001 s ;

// nested sections: the header is the line before its own-line brace
solver
{
    type : mumps ;
    tolerance : 1e-6 ;
    max iterations : 1000 ;

    // sub-sections
    preconditioner
    {
        type : ilu ;
        fill level : 2 ;
    }
}

// named sections ( type : label )
material : steel
{
    density : 7850 kg/m^3 ;
    youngs modulus : 200e9 Pa ;
}

material : aluminum
{
    density : 2700 kg/m^3 ;
    youngs modulus : 69e9 Pa ;
}
```

Every key line ends with `;`; a line without one is not a key. Keys, section types and labels
are folded to lower case by the parser.

#### InputFile Class

```cpp
// Constructor
InputFile(const string & aPath);

// Access sections
const input::Section * section(const string & aSection) const;  // By name
const input::Section * section(const index_t aIndex) const;     // By index

// Query
bool section_exists(const string & aSection) const;
index_t num_sections() const;

// Debug
void print();  // Print entire structure
```

#### Input_Section Class

```cpp
// Hierarchy navigation
const string & type() const;         // Section type (e.g., "Solver")
const string & label() const;        // Section label (e.g., "Steel" in "Material:Steel")
const string & key() const;          // Full key (type:label)

// Sub-sections
const Section * section(const string & aType) const;                      // Type only
const Section * section(const string & aType, const string & aLabel) const;  // Type + label
const Section * section(const index_t aIndex) const;                      // By index
bool section_exists(const string & aType) const;
bool section_exists(const string & aType, const string & aLabel) const;
index_t num_sections() const;

// Key-value access
bool key_exists(const string & aKey) const;
bool key_is_real(const string & aKey) const;

const string & get_string(const string & aKey) const;
bool get_bool(const string & aKey) const;
real get_real(const string & aKey) const;
value get_value(const string & aKey, const string & aUnit) const;  // With unit conversion
int get_int(const string & aKey) const;
string get_units(const string & aKey) const;

// Array parsing
void get_ids(const string & aKey, Vector<id_t> & aIDs) const;
void get_reals(const string & aKey, Vector<real> & aReals) const;
void get_id_groups(const string & aKey, Cell<Cell<id_t>> & aIDs) const;

// Keys iteration
index_t num_keys() const;
const string & key(const index_t aIndex) const;

// Hierarchy
int level() const;                   // Nesting depth
const Section * parent() const;
string tree() const;                 // Full path, joined with "->" (e.g., "input.conf->solver->preconditioner")
```

#### Usage Example

```cpp
// Load configuration
InputFile config("simulation.input");

// Access top-level key
if (config.section_exists("Solver")) {
    const input::Section* solver = config.section("Solver");

    // Get values with type conversion
    string solver_type = solver->get_string("type");        // "MUMPS"
    real tol = solver->get_real("tolerance");               // 1e-6
    int max_iter = solver->get_int("max_iterations");       // 1000

    // Navigate to sub-section
    if (solver->section_exists("Preconditioner")) {
        const input::Section* precond = solver->section("Preconditioner");
        string precond_type = precond->get_string("type");  // "ILU"
    }
}

// Access named sections
// the root stores labelled sections under "type:label"; the two-argument
// section("Material", "Steel") form exists on input::Section only
if (config.section_exists("material:steel")) {
    const input::Section* steel = config.section("material:steel");

    // Get value with unit conversion
    value rho = steel->get_value("density", "kg/m^3");       // Converted to base units
    value E = steel->get_value("youngs_modulus", "Pa");
}

// Iterate over all materials
for (index_t i = 0; i < config.num_sections(); ++i) {
    const input::Section* sec = config.section(i);
    if (sec->type() == "material") {
        string mat_name = sec->label();  // "steel", "aluminum" — the parser lower-cases types and labels
        // Process material
    }
}
```

**Unit Handling:**

The `get_value()` method automatically converts from file units to BELFEM base units.

**Supported temperature conversions (cl_Input_Section.cpp, `create_key`):**
- `K`, `°K` — Kelvin (base unit, no conversion)
- `C`, `°C` — Celsius → K: `T_K = T_C + 273.15`
- `°F` — Fahrenheit → K: `T_K = (T_F - 32) / 1.8 + 273.15` (a bare `F` is farad, not Fahrenheit)
- `R`, `°R` — Rankine → K: `T_K = T_R / 1.8`

**General units:** Uses `unit_to_si()` from physics module for automatic SI conversion:
```cpp
// Example conversions (from input file → BELFEM internal units)
density = 7.85 g/cm^3    → 7850 kg/m³
pressure = 100 MPa       → 1e8 Pa
length = 25.4 mm         → 0.0254 m
energy = 1 kJ            → 1000 J
```

**Usage:**
```cpp
// File contains: temperature = 300 C
value T = section->get_value("temperature", "C");
// T automatically converted to 573.15 K (BELFEM base unit)

// get_real() returns the converted value directly
real T_kelvin = section->get_real("temperature");  // 573.15 K

// get_units() returns the original unit string from file
string unit = section->get_units("temperature");  // "C"
```

> **Error Handling Contract:**
> - Missing keys or sections are **logic errors**, not recoverable conditions.
> - `InputFile` uses **fail-fast** behavior via `BELFEM_ERROR` assertions.
> - Always check `section_exists()` and `key_exists()` before access in production code.

### XML Class

**Location**: `cl_XML.hpp`

XML document interface using TinyXML2 backend.

> **Requires:** `BELFEM_XML` enabled at compile time and TinyXML2 library

#### Constructor

```cpp
XML(const string & aPath, const FileMode aMode = FileMode::OPEN_RDONLY);
```

#### Interface

```cpp
// Path accessor
const string & path() const;

// Navigation
void select_first_child(const string & aLabel);
void select_parent();
bool next_sibling_of_same_name();
void select_subtree(const string & aTree);  // Path like "root/section/subsection"

// Query
bool child_exists(const string & aLabel);
bool key_exists(const string & aKey);
uint number_of_children();
uint number_of_children(const string & aLabel);

// Child-element text access ( XML attributes are not read )
string get_string(const string & aKey);
int get_int(const string & aKey);
real get_real(const string & aKey);
bool get_bool(const string & aKey);
```

#### Usage Example

```xml
<!-- config.xml -->
<Configuration>
    <Solver>
        <type>MUMPS</type>
        <threads>8</threads>
        <Tolerance>1e-6</Tolerance>
        <MaxIterations>1000</MaxIterations>
    </Solver>
</Configuration>
```

```cpp
#ifdef BELFEM_XML
    XML xml("config.xml", FileMode::OPEN_RDONLY);

    // Navigate to Solver element
    xml.select_first_child("Configuration");
    xml.select_first_child("Solver");

    // Read the text of child elements ( get_* looks up a child element by
    // name; it does not read XML attributes )
    string solver_type = xml.get_string("type");     // "MUMPS"
    int threads = xml.get_int("threads");            // 8
    real tolerance = xml.get_real("Tolerance");      // 1e-6

    xml.select_parent();  // Back to Configuration
#else
    BELFEM_ERROR(false, "XML support not enabled. Reconfigure with -DUSE_TINYXML2=ON");
#endif
```

---

## File Utilities

**Location**: `filetools.hpp`, `filetools.cpp`

### FileMode Enum

```cpp
enum class FileMode
{
    NEW,                    // Create new file (TRUNCATES if exists!)
    OPEN_RDONLY,           // Read-only (error if doesn't exist)
    OPEN_RDONLY_PARALLEL,  // MPI parallel read-only
    OPEN_RDWR              // Read-write (error if doesn't exist)
};
```

### Free Functions

```cpp
// Check file existence (filetools.cpp:23)
bool file_exists(const string & aPath);

// Extract file extension — lives in src/core/stringtools.hpp, not in filetools
string filetype(const string & aPath);

// Create MPI-safe file paths (make_path_parallel, filetools.cpp)
string make_path_parallel(const string & aPath);
```

**Usage:**

```cpp
// Check before opening
if (file_exists("restart.h5")) {
    HDF5 file("restart.h5", FileMode::OPEN_RDONLY);
    // Load restart data
}

// Extract extension
string ext = filetype("data.exo");  // Returns "exo"

// MPI-parallel file naming (make_path_parallel, filetools.cpp)
string path = make_path_parallel("output.h5");
// Serial (comm_size == 1): "output.h5" (unchanged)
// Parallel with 4 ranks:
//   Rank 0: "output_4.0.h5"
//   Rank 1: "output_4.1.h5"
//   Rank 2: "output_4.2.h5"
//   Rank 3: "output_4.3.h5"
// Pattern: base_N.X.ext where N=comm_size(), X=comm_rank()

HDF5 file(path, FileMode::NEW);
```

**Implementation details:**
- `file_exists()`: `std::filesystem::exists()`
- `filetype()`: Returns substring after last '.'
- `make_path_parallel()`: reads the global communicator `gComm` for the rank and size. BELFEM is not internally thread safe (`doc/coding_philosophy.md`); this is safe only in the sense that it does not mutate `gComm`.

---

## Common Patterns

### Pattern 1: Checkpoint/Restart with HDF5

```cpp
// Save checkpoint
void save_checkpoint(const Vector<real> & aState, real aTime, int aTimestep) {
    HDF5 file("checkpoint.h5", FileMode::NEW);

    file.save_data("time", aTime);
    file.save_data("timestep", aTimestep);
    file.save_data("state", aState);
}

// Load checkpoint
void load_checkpoint(Vector<real> & aState, real & aTime, int & aTimestep) {
    if (!file_exists("checkpoint.h5")) {
        BELFEM_ERROR(false, "Checkpoint file not found");
    }

    HDF5 file("checkpoint.h5", FileMode::OPEN_RDONLY);
    file.load_data("time", aTime);
    file.load_data("timestep", aTimestep);
    file.load_data("state", aState);  // Auto-resizes aState
}
```

### Pattern 2: Configuration-Driven Simulation

```cpp
void setup_solver_from_config(const string & aConfigPath) {
    InputFile config(aConfigPath);

    // Get solver parameters
    const input::Section* solver_sec = config.section("Solver");
    string solver_type = solver_sec->get_string("type");
    real tolerance = solver_sec->get_real("tolerance");

    // Create solver based on config
    SolverType type = string_to_solver_type(solver_type);
    Solver solver(type);

    SolverParameters params;
    params.set_tolerance(tolerance);

    if (solver_sec->key_exists("max_iterations")) {
        params.set_max_iterations(solver_sec->get_int("max_iterations"));
    }

    solver.set_parameters(params);
}
```

### Pattern 3: Time-Series Data Storage

```cpp
void save_time_series(const Cell<real> & aTimes,
                      const Cell<Vector<real>> & aFields) {
    HDF5 file("timeseries.h5", FileMode::NEW);

    // Save time vector
    Vector<real> times(aTimes.size());
    for (index_t i = 0; i < aTimes.size(); ++i) {
        times(i) = aTimes(i);
    }
    file.save_data("times", times);

    // Save each timestep in separate group
    for (index_t i = 0; i < aFields.size(); ++i) {
        string group_name = "TimeStep_" + std::to_string(i);
        file.create_group(group_name);
        file.save_data("field", aFields(i));
        file.close_active_group();
    }
}
```

### Pattern 4: Parallel I/O (One File Per Rank)

```cpp
void save_distributed_data(const Vector<real> & aLocalData) {
    // Each rank writes its own file. Pass aParallelMode rather than building
    // the path by hand: the open is gated on ( aParallelMode || rank == 0 ),
    // so a per-rank path alone leaves the workers with nothing open.
    HDF5 file("distributed_output.h5", FileMode::NEW, true);
    file.save_data("local_data", aLocalData);
    file.save_data("rank", comm_rank());
    file.save_data("size", comm_size());
}

void load_distributed_data(Vector<real> & aLocalData) {
    // make_path_parallel() gives the same name the writer used, which is what
    // the existence check needs; the open still takes aParallelMode.
    string path = make_path_parallel("distributed_output.h5");

    if (!file_exists(path)) {
        BELFEM_ERROR(false, "Distributed file for rank %d not found", comm_rank());
    }

    HDF5 file("distributed_output.h5", FileMode::OPEN_RDONLY, true);
    file.load_data("local_data", aLocalData);
}
```

### Pattern 5: Reading Tabulated Data from CSV

```cpp
// Load experimental data from CSV
void load_material_curve(const string & aPath,
                         Vector<real> & aStrain,
                         Vector<real> & aStress) {
    CsvFile csv(aPath);
    const Matrix<real> & data = csv.data();

    index_t n = data.n_rows();
    aStrain.set_size(n);
    aStress.set_size(n);

    // Column 0: strain, Column 1: stress
    for (index_t i = 0; i < n; ++i) {
        aStrain(i) = data(i, 0);
        aStress(i) = data(i, 1);
    }
}
```

---

## Performance Tips

### 1. Minimize HDF5 Group Operations

```cpp
// SLOW: Open/close file multiple times
for (index_t i = 0; i < n; ++i) {
    HDF5 file("data.h5", FileMode::OPEN_RDWR);  // File open overhead!
    file.save_data("field_" + std::to_string(i), fields(i));
}

// FAST: Open once, write all data
HDF5 file("data.h5", FileMode::NEW);
for (index_t i = 0; i < n; ++i) {
    file.save_data("field_" + std::to_string(i), fields(i));
}
```

### 2. Batch HDF5 Writes

```cpp
// SLOW: Many small writes
HDF5 file("data.h5", FileMode::NEW);
for (index_t i = 0; i < 1000000; ++i) {
    file.save_data("value_" + std::to_string(i), values(i));  // 1M datasets!
}

// FAST: Write as single vector
HDF5 file("data.h5", FileMode::NEW);
file.save_data("values", values);  // Single write
```

### 3. Avoid Unnecessary String Copies in Ascii

```cpp
// INEFFICIENT: Creates string copy each access
Ascii file("large.txt", FileMode::OPEN_RDONLY);
for (index_t i = 0; i < file.length(); ++i) {
    string line_copy = file.line(i);  // Copy!
    // Process
}

// EFFICIENT: Use const reference
Ascii file("large.txt", FileMode::OPEN_RDONLY);
for (index_t i = 0; i < file.length(); ++i) {
    const string & line_ref = file.line(i);  // No copy
    // Process
}
```

### 4. Pre-Allocate for Large CSV Files

```cpp
// CsvFile loads entire file into Matrix<real>
// For very large CSVs (> 1 GB), consider streaming instead

// If CSV is too large:
Ascii file("huge_data.csv", FileMode::OPEN_RDONLY);
// Parse line-by-line, accumulate statistics, then allocate Matrix
```

### 5. One File Per Rank for Large MPI Jobs

**BELFEM has no parallel-HDF5 support.** `BELFEM_PHDF5` is defined nowhere in the tree, and
`H5Pset_fapl_mpio` is never called — every HDF5 file is opened with `H5P_DEFAULT`
(`cl_HDF5.cpp:47-51,67-71`). The per-rank pattern is not a fallback; it is the only pattern.

```cpp
// Each rank writes its own file: pass aParallelMode.
HDF5 file("output.h5", FileMode::NEW, true);   // -> output_<size>.<rank>.h5
```

**Do not hand-roll it with `make_path_parallel()` alone.** The open is gated on
`( aParallelMode || gComm.rank() == 0 )` (`cl_HDF5.cpp:40`), so passing a per-rank path without
the flag still leaves every non-root rank with no open file at all.

---

## MPI Considerations

> **MPI I/O Contract**
>
> **The only pattern:** one file per rank, via the `aParallelMode` constructor argument.
> - Each rank writes to its own file: `base_N.X.ext`
> - No communication overhead during I/O
> - Post-processing tools must gather files
>
> **There is no single-file pattern.** `aParallelMode = true` does not open a shared file — it
> rewrites the path per rank and opens that. Nothing in the HDF5 layer is collective: no MPI
> access property list is set, and `close()` is a plain `H5Fclose`. A rank may open, write and
> close entirely on its own.
>
> **Consequences:** there is no collective-close hazard to guard against, and post-processing
> always has to gather the per-rank files. If you genuinely need cooperative single-file output,
> it does not exist yet and would have to be added to the HDF5 wrapper.

### File Naming Conventions

```cpp
// Serial job (comm_size() == 1)
HDF5 file("output.h5", FileMode::NEW);

// MPI job: each rank writes its own file.
// Pass aParallelMode -- do NOT hand-roll the path, because the open is gated
// on ( aParallelMode || rank == 0 ) and a per-rank path alone leaves every
// non-root rank with no open file ( cl_HDF5.cpp:40 ).
if (comm_size() > 1) {
    HDF5 file("output.h5", FileMode::NEW, true);
    // rank 0 -> "output_4.0.h5", rank 1 -> "output_4.1.h5", ...
}
```

### Master-Rank I/O Pattern

```cpp
// Only rank 0 writes global configuration
if (comm_rank() == 0) {
    InputFile config("solver.input");
    // Parse configuration
}

// Broadcast to all ranks (using comm module)
// ... MPI_Bcast or BELFEM comm::broadcast
```

### Per-rank HDF5 output (`aParallelMode`)

> **`aParallelMode` is one file per rank, not cooperative PHDF5.** The third constructor
> argument rewrites the path through `make_path_parallel()` (`cl_HDF5.cpp:31-40`) and opens that
> per-rank file with `H5P_DEFAULT` — no MPI access property list is ever set
> (`cl_HDF5.cpp:47-51,67-71`). Nothing is collective, so there is no cooperative write and no
> collective close to synchronize.

```cpp
// Each rank opens its OWN file: global_solution_4.0.h5, _4.1.h5, ...
HDF5 file("global_solution.h5", FileMode::NEW, true);  // aParallelMode = true

file.save_data("local_field", local_data);
// Nothing here is collective. A rank may open, write and close on its own.
```

Because the files are separate, the usual parallel-I/O hazards do not apply: a rank that opens
without the others cannot deadlock, and post-processing has to stitch the per-rank files together
itself.


---

## Thread Safety

### HDF5 Thread Safety

- **Not thread-safe** by default (HDF5 library limitation)
- HDF5 output here is per-rank and process-level; there is no threaded or collective path
- **Solution:** Use OpenMP critical sections or one file per thread

```cpp
#pragma omp parallel
{
    int tid = omp_get_thread_num();

    // One file per thread
    string path = "thread_" + std::to_string(tid) + "_output.h5";
    HDF5 file(path, FileMode::NEW);

    // Thread-local writes
    file.save_data("data", thread_local_data);
}
```

### Ascii/CsvFile Thread Safety

- **Read-only safe**: Multiple threads can read same `Ascii` object
- **Write not safe**: Modifications to buffer require synchronization

```cpp
// SAFE: Concurrent reads
Ascii file("data.txt", FileMode::OPEN_RDONLY);

#pragma omp parallel for
for (index_t i = 0; i < file.length(); ++i) {
    const string & line = file.line(i);  // Safe: read-only
    // Parse line
}
```

### InputFile Thread Safety

- **Read-only after construction**: Safe for concurrent queries
- Typical usage: Parse once in serial, query in parallel

```cpp
// Parse in serial
InputFile config("params.input");
const input::Section* solver = config.section("Solver");

// Query in parallel (read-only)
#pragma omp parallel
{
    real tol = solver->get_real("tolerance");  // Safe
}
```

---

## Memory Ownership and Lifetime

All I/O classes use **RAII** (Resource Acquisition Is Initialization) — resources are automatically released in destructors.

> **Error Handling Reminder:** BELFEM ends the run on error (abort in release, throw in debug). Resources are cleaned up explicitly before abort, or reclaimed by OS on process termination. See "Error Handling Philosophy" section above for details.

| Class | Resource Owned | Ownership Model | Lifetime Rules |
|-------|----------------|-----------------|----------------|
| `HDF5` | File handle (`hid_t mFile`) | RAII | Destructor calls `close()` → `H5Fclose()` |
| | Group handles (`Cell<hid_t> mTree`) | RAII | Destructor calls `close()` → `close_active_group()` → `H5Gclose()` on the active group only; `close_tree()` must be called explicitly to close the rest |
| `Ascii` | Line buffer (`Cell<string> mBuffer`) | Owned | Destructor clears buffer; `BELFEM_ERROR` (always active) if unsaved changes |
| `CsvFile` | Matrix data (`Matrix<real> mData`) | Owned | Destructor releases matrix memory (via `Matrix` RAII) |
| `InputFile` | Section tree (`input::Section* mData`) | Owned | Destructor deletes root section (recursive delete) |
| | Section pointers returned | **Non-owning** | Valid while `InputFile` lives; **do not delete** |
| `input::Section` | Child sections (`Cell<Section*> mData`) | Owned | Destructor deletes all children recursively |
| | Map entries (`Map<string, Section*>`) | **Non-owning views** | Point to children in `mData`; no separate cleanup |
| `XML` | TinyXML2 document (`tinyxml2::XMLDocument mFile`) | RAII | Managed by tinyxml2 (auto-cleanup) |
| | Element pointers returned | **Non-owning** | Valid while `XML` object lives; **do not delete** |

### Lifetime Rules

**HDF5:**
```cpp
{
    HDF5 file("data.h5", FileMode::NEW);
    file.create_group("Results");
    // ... operations ...
}  // Destructor closes all groups + file automatically
```

**InputFile:**
```cpp
InputFile config("sim.input");
const input::Section* solver = config.section("Solver");

// solver pointer valid here (config still alive)
real tol = solver->get_real("tolerance");

// DO NOT: delete solver;  // WRONG! Not the owner
// config destructor will clean up
```

**Ascii buffer changes:**
```cpp
Ascii file("data.txt", FileMode::NEW);
file.print("Line 1");
file.print("Line 2");
file.save();       // Ascii has save(), not save_data() -- MUST call before
                   // the destructor if mChangedSinceLastSave
// Destructor raises BELFEM_ERROR (always active) if changes not saved
```

---

## Debugging Tips

### 1. Enable HDF5 Error Reporting

```cpp
// HDF5 errors printed to stderr automatically
// Check BELFEM_ERROR and BELFEM_ASSERT messages for details

// Example error:
// "Dataset 'Results/Temperature' of type Matrix<real> does already exist."
```

### 2. Check File Existence Before Opening

```cpp
if (!file_exists("restart.h5")) {
    BELFEM_ERROR(false, "Restart file 'restart.h5' not found in %s",
                 std::filesystem::current_path().c_str());
}
```

### 3. Verify HDF5 Group State

```cpp
// After complex group navigation, verify path
HDF5 file("data.h5", FileMode::NEW);
file.create_group("A");
file.create_group("B");
// Now at /A/B

// When in doubt, close all groups and start fresh
file.close_active_group();  // Back to /A
file.close_active_group();  // Back to /
```

### 4. Print InputFile Structure

```cpp
InputFile config("complex.input");
config.print();  // Prints entire hierarchical structure to stdout
```

### 5. Validate Input File Sections

```cpp
const input::Section* solver = config.section("Solver");

if (!solver->key_exists("tolerance")) {
    BELFEM_ERROR(false, "Required key 'tolerance' not found in [Solver] section");
}

if (!solver->key_is_real("tolerance")) {
    BELFEM_ERROR(false, "Key 'tolerance' must be numeric");
}
```

### 6. HDF5 File Inspection Tools

Use command-line tools to inspect HDF5 files:

```bash
# List contents
h5dump -n data.h5

# View dataset
h5dump -d /Results/Temperature data.h5

# Interactive browser (if available)
hdfview data.h5
```

### 7. Check CSV Parsing Errors

```cpp
try {
    CsvFile csv("data.csv");
} catch (...) {
    BELFEM_ERROR(false, "Failed to parse CSV. Check for non-numeric data or irregular rows.");
}
```

---

## Enumeration Reference

### FileMode

**Location**: `filetools.hpp`

```cpp
enum class FileMode
{
    NEW,                    // Create new file (TRUNCATES if exists!)
    OPEN_RDONLY,           // Open existing file read-only (fails if doesn't exist)
    OPEN_RDONLY_PARALLEL,  // handled by Ascii; NOT implemented by HDF5
    OPEN_RDWR              // Open existing file read-write (fails if doesn't exist)
};
```

**Usage Guidelines:**

| Mode | Use Case | File Must Exist | Allows Write | Behavior if Exists |
|------|----------|----------------|--------------|-------------------|
| `NEW` | Creating output files | No | Yes | **TRUNCATES** (data loss!) |
| `OPEN_RDONLY` | Reading input files | Yes | No | Opens existing |
| `OPEN_RDONLY_PARALLEL` | `Ascii` only — rank 0 reads the file and broadcasts the buffer to the others (`cl_Ascii.cpp:46-49`, `load_buffer(true)`). **`HDF5` has no case for this mode** and reaches its "unknown filemode passed" error (`cl_HDF5.cpp:42-103`) | Yes | No | Opens existing |
| `OPEN_RDWR` | Appending/modifying | Yes | Yes | Opens existing |

> **Critical:** `FileMode::NEW` uses `H5F_ACC_TRUNC` — existing files are **silently overwritten**.
> Always check `file_exists()` first or use versioned filenames for production.

---

## See Also

- [HDF5 Documentation](https://portal.hdfgroup.org/display/HDF5/HDF5) - Official HDF5 library reference
- [TinyXML2](https://github.com/leethomason/tinyxml2) - XML parser used by BELFEM
- [Mesh Module](../../mesh/doc/README.md) - Uses HDF5 for mesh I/O (Exodus format)
- [Coding Philosophy](../../../doc/coding_philosophy.md) - BELFEM design patterns
- [Communication Module](../../comm/doc/README.md) - MPI abstractions for parallel I/O
