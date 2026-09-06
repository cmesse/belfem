# I/O Module Documentation {#io_index}

Documentation for BELFEM's file input/output interfaces.

## Contents

### Guides

- [**io_usage_guide.md**](io_usage_guide.md) - Comprehensive usage guide for BELFEM I/O classes
  - HDF5 hierarchical data format interface
  - Ascii line-based text file interface
  - CsvFile CSV reader for numeric data
  - InputFile configuration file parser with hierarchical sections
  - XML document interface
  - File utilities (FileMode, file_exists, make_path_parallel)
  - MPI parallel I/O patterns
  - Performance tips and common patterns

## Quick Reference

### Core Classes

| Type | Description | File |
|------|-------------|------|
| `HDF5` | Hierarchical data format I/O | `cl_HDF5.hpp` |
| `Ascii` | Line-based ASCII file interface | `cl_Ascii.hpp` |
| `CsvFile` | CSV reader (numeric data) | `cl_CsvFile.hpp` |
| `InputFile` | Configuration file parser | `cl_InputFile.hpp` |
| `input::Section` | Hierarchical configuration section | `cl_Input_Section.hpp` |
| `XML` | XML document interface | `cl_XML.hpp` |

### Build Configuration

```cmake
# Enable HDF5 support (recommended for large datasets)
cmake -DUSE_HDF5=ON ..

# Enable XML support (requires tinyxml2)
cmake -DUSE_TINYXML2=ON ..
```

### FileMode Options

| Mode | Use Case | Must Exist | Writable | If Exists |
|------|----------|-----------|----------|-----------|
| `NEW` | Create new file | No | Yes | **TRUNCATES!** |
| `OPEN_RDONLY` | Read existing file | Yes | No | Opens |
| `OPEN_RDONLY_PARALLEL` | MPI parallel read | Yes | No | Opens |
| `OPEN_RDWR` | Modify existing file | Yes | Yes | Opens |

> **Warning:** `FileMode::NEW` silently truncates existing files (`H5F_ACC_TRUNC`). Always check `file_exists()` first!

### Basic Usage Patterns

#### HDF5 Save/Load

```cpp
// Save data
HDF5 file("output.h5", FileMode::NEW);
file.save_data("scalar", 3.14);
file.save_data("vector", my_vector);
file.save_data("matrix", my_matrix);

// Create groups
file.create_group("Results");
file.save_data("Temperature", T_field);
file.close_active_group();

// Load data (auto-resizes containers)
HDF5 input("data.h5", FileMode::OPEN_RDONLY);
Vector<real> v;
input.load_data("vector", v);  // v resized to file data
```

#### ASCII/CSV Reading

```cpp
// Line-based text file
Ascii txt("input.txt", FileMode::OPEN_RDONLY);
for (index_t i = 0; i < txt.length(); ++i) {
    const string & line = txt.line(i);
    // Parse line
}

// CSV numeric data → Matrix<real>
CsvFile csv("data.csv");
const Matrix<real> & data = csv.data();
real value = data(row, col);
```

#### Configuration Files

```cpp
// Parse hierarchical config
InputFile config("simulation.input");

// Access sections
const input::Section* solver = config.section("Solver");
string type = solver->get_string("type");
real tol = solver->get_real("tolerance");

// Navigate sub-sections
const input::Section* precond = solver->section("Preconditioner");
```

#### XML Files

```cpp
#ifdef BELFEM_XML
    XML xml("config.xml", FileMode::OPEN_RDONLY);
    xml.select_first_child("Configuration");
    xml.select_first_child("Solver");
    string solver_type = xml.get_string("type");
#endif
```

## File Utilities

| Function | Description | File |
|----------|-------------|------|
| `file_exists()` | Check if file exists | `filetools.hpp` |
| `make_path_parallel()` | Create MPI-safe paths | `filetools.hpp` |

### MPI Parallel I/O

```cpp
// One file per rank -- use the third constructor argument.
HDF5 file("output.h5", FileMode::NEW, true);   // aParallelMode
// -> each rank opens its own output_<size>.<rank>.h5

// Renaming the path by hand is NOT equivalent: the open is gated on
// ( aParallelMode || rank == 0 ), so without the flag only rank 0 opens
// anything, whatever path you hand it.
string path = make_path_parallel("output.h5");
HDF5 rank0_only(path, FileMode::NEW);
```

## Data Type Support

### HDF5 Save/Load Types

| Category | Types |
|----------|-------|
| **Scalars** | `string`, `sint`, `uint`, `luint`, `real`, `bool` |
| **Vectors** | `Vector<sint>`, `Vector<uint>`, `Vector<luint>`, `Vector<real>` |
| **Matrices** | `Matrix<sint>`, `Matrix<uint>`, `Matrix<luint>`, `Matrix<real>` |
| **Arrays** | `Cell<string>` |

### InputFile Value Types

| Method | Return Type | Use Case |
|--------|-------------|----------|
| `get_string()` | `string` | Text values |
| `get_bool()` | `bool` | true/false |
| `get_int()` | `int` | Integer values |
| `get_real()` | `real` | Floating-point |
| `get_value()` | `value` | Physical quantities; takes the expected unit as a second argument |
| `get_ids()` | `void`, out-parameter `Vector<id_t>&` | ID lists |
| `get_reals()` | `void`, out-parameter `Vector<real>&` | Real arrays |

## Common Pitfalls

1. **Modern C++ tool warnings**: Static analyzers may flag `malloc`/`free` as "unsafe" — these are **false positives** for BELFEM (uses C libraries with abort-on-error, not exceptions). See usage guide for details.
2. **HDF5 not enabled**: Check `BELFEM_HDF5` CMake flag
3. **FileMode::NEW data loss**: Silently truncates existing files! Always check `file_exists()` first
4. **Unclosed HDF5 groups**: Call `close_active_group()` after group operations (resource leaks!)
5. **InputFile hierarchy**: Navigate parent → child, don't access nested sections directly
6. **CsvFile numeric-only**: Non-numeric CSV data will fail parsing
7. **No parallel HDF5**: the wrapper never opens a shared file. Without `aParallelMode` only rank 0 opens anything; with it every rank writes its own `base_N.X.h5`. Nothing is collective, so there is no deadlock to guard — but do not expect a single output file.

## See Also

- **Source code** - header files in `src/io/`
- [Mesh module](../../mesh/doc/README.md) - Uses HDF5 for Exodus mesh I/O
- [Communication module](../../comm/doc/README.md) - MPI abstractions
- [Linear algebra module](../../linalg/doc/README.md) - Vector and Matrix types
- [Root documentation](../../../doc/README.md) - General BELFEM documentation
