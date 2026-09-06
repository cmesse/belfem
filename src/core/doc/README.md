# Core Module Documentation {#core_index}

**Module:** src/core
**Purpose:** Index of documentation for BELFEM's core utility module

---

## Overview

The `core` module provides fundamental utilities and infrastructure used throughout the BELFEM framework, including:

- Logging and output systems
- Timing and profiling tools
- Command-line argument processing
- Type system and physical constants
- String manipulation utilities
- Error handling and assertions
- Hashing and random number generation

---

## Documentation Files

### Module Documentation

- **[core_usage_guide.md](core_usage_guide.md)** - Comprehensive guide to the core module
  - Architecture and file organization
  - Logging system (Logger, Progressbar, Banner)
  - Timing and profiling (Timer, Profiler)
  - Command-line processing (Arguments)
  - Type system (typedefs, constants, units)
  - String utilities and formatting
  - Error handling and assertions
  - Hashing and random numbers
  - Usage examples

---

## Quick Reference

### Key Classes

| Class | File | Purpose |
|-------|------|---------|
| `Logger` | cl_Logger.{hpp,cpp} | Hierarchical logging with info levels |
| `Timer` | cl_Timer.hpp | High-resolution wall-clock timing |
| `Profiler` | cl_Profiler.{hpp,cpp} | gperftools CPU profiling, exported in Callgrind format |
| `Progressbar` | cl_Progressbar.{hpp,cpp} | Progress visualization |
| `Arguments` | cl_Arguments.{hpp,cpp} | Command-line argument parsing |
| `Hash` | cl_Hash.hpp | Incremental hash computation |

### Key Files

| File | Purpose |
|------|---------|
| typedefs.hpp | Framework type definitions (real, index_t, id_t, etc.) |
| constants.hpp | NIST-compliant physical and mathematical constants |
| units.hpp | Unit conversion factors (SI, metric, imperial) |
| stringtools.{hpp,cpp} | String manipulation and parsing |
| assert.{hpp,cpp} | Custom assertion system with formatted errors; failed checks also land in the system log — `journalctl -t belfem` on Linux, `log show` on macOS (see the usage guide, "System Log Integration") |
| random.hpp | MPI-aware random number generation |
| banner.{hpp,cpp} | Application startup banners |
| globals.hpp | Framework global variables |

### Common Operations

```cpp
// Logging
message(InfoLevel::Default, "Message with format: %d", value);

// Timing
Timer timer;
// ... work ...
uint64_t elapsed = timer.stop();

// Progress tracking
Progressbar bar(num_steps);
for (uint i = 0; i < num_steps; ++i) {
    bar.step();
}
bar.finish();

// Error handling — the message also lands in the system log:
//   journalctl -t belfem
BELFEM_ERROR(condition, "Error message with %d", value);

// String parsing
Cell<index_t> values;
string_to_cell("1:5, 10, 15:17", values);  // → {1,2,3,4,5,10,15,16,17}

// Hashing
Hash hash;
hash += value1;
hash += value2;
std::size_t key = hash.value();
```

---

## Source Code

**Module location:** `../../`

**Key source files:**
- Logging: `cl_Logger.{hpp,cpp}`, `cl_Progressbar.{hpp,cpp}`, `banner.{hpp,cpp}`
- Timing: `cl_Timer.hpp`, `cl_Profiler.{hpp,cpp}`
- Types: `typedefs.hpp`, `constants.hpp`, `units.hpp`
- Strings: `stringtools.{hpp,cpp}`, `fn_sprint.hpp`
- Utilities: `cl_Hash.hpp`, `random.hpp`, `assert.{hpp,cpp}`

---

## External References

### Standards and Specifications

- **NIST Physical Constants**: https://physics.nist.gov/cuu/Constants/
  - All constants in `constants.hpp` sourced from NIST
  - Values marked as "(exact)" use 2019 SI redefinition

### Related BELFEM Modules

- **Containers** (`src/containers/`): `cl_Cell.hpp` used extensively in core utilities
- **Communication** (`src/comm/`): MPI integration for Logger and random number generation
- **Mesh** (`src/mesh/`): Uses core type system (index_t, id_t, etc.)

---

## Development Notes

### Adding New Global Variables

From `globals.hpp`:
1. Add declaration in `globals.hpp`
2. Initialize in `Communicator::set_globals()`
3. Use `BELFEM_QUIET_NAN` as default if no clear initial value
4. If set only on master, synchronize using `broadcast(value)` (free function, `commtools.hpp`)

### Compiler Compatibility

The core module uses pragmas to suppress warnings across compilers:
- **GCC**: `#pragma GCC diagnostic push/pop`
- **Clang**: `#pragma clang diagnostic push/pop`
- **Intel**: `#pragma warning push/pop`

This is necessary for format string warnings in `sprint()` and `message()`.

### Type System Philosophy

- **real**: Always `double` for precision
- **index_t**: Configurable 32-bit or 64-bit (`BELFEM_INT64` flag)
- **Sentinel values**: `gNoIndex`, `gNoID`, `gNoOwner` for invalid/unset states
- **Tolerances**: `BELFEM_EPSILON` (10× machine ε), `BELFEM_MESH_EPSILON` (1 nm)

---

## See Also

- **Project README**: `../../../README.md`
- **Claude Instructions**: `../../../CLAUDE.md`
- **Documentation Guidelines**: `../../../doc/documentation_guidelines.md`
- **General Documentation**: `../../../doc/README.md`
