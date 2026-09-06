# Containers Module Documentation {#containers_index}

Documentation for BELFEM's custom container classes.

## Contents

### Guides

- [**container_usage_guide.md**](container_usage_guide.md) - Comprehensive usage guide for all BELFEM containers
  - Overview of all 10 container types
  - Detailed API reference for each container
  - Performance characteristics and selection guide
  - Common usage patterns
  - Best practices and examples

## Quick Reference

| Container | Purpose | File |
|-----------|---------|------|
| `Cell<T>` | Dynamic array (primary container) | `cl_Cell.hpp` |
| `Map<K,V>` | Hash map (unordered key-value) | `cl_Map.hpp` |
| `OrderedMap<K,V>` | Sorted map (ordered key-value) | `cl_OrderedMap.hpp` |
| `Set<T>` | Hash set with set operations | `cl_Set.hpp` |
| `Queue<T>` | FIFO queue | `cl_Queue.hpp` |
| `Bitset<N>` | Compile-time fixed bitset | `cl_Bitset.hpp` |
| `DynamicBitset` | Runtime-sized bitset | `cl_DynamicBitset.hpp` |
| `ShiftRegister<T>` | Fixed-capacity FIFO with history | `cl_ShiftRegister.hpp` |
| `Genome<B,N>` | Genetic algorithm encoding | `cl_Genome.hpp` |
| `StringList` | C-string array for I/O | `cl_StringList.hpp` |

## See Also

- **Source code** - header files in `src/containers/`
- [Root documentation](../../../doc/README.md) - General BELFEM documentation
- [Documentation guidelines](../../../doc/documentation_guidelines.md) - Documentation standards
