# Communication Module Documentation {#comm_index}

Documentation for BELFEM's MPI communication abstraction layer.

## Contents

### Guides

- [**comm_usage_guide.md**](comm_usage_guide.md) - Comprehensive usage guide for BELFEM MPI communication
  - Global communicator (`gComm`)
  - Type-safe MPI datatype mapping
  - Broadcast operations (scalars, Cells, Vectors, Matrices)
  - Point-to-point communication (send/receive)
  - Collective operations (distribute/collect)
  - Common parallel patterns
  - Performance tips and debugging

## Quick Reference

### Core Components

| Component | Purpose | File |
|-----------|---------|------|
| `Communicator` | Global MPI communicator manager | `cl_Communicator.hpp` |
| `comm_type<T>()` | Type-to-MPI datatype mapping | `commtypes.hpp` |
| Communication functions | Template wrappers for MPI operations | `commtools.hpp` |

### Utility Functions

| Function | Purpose |
|----------|---------|
| `comm_size()` | Get number of processes |
| `comm_rank()` | Get current process rank |
| `comm_barrier()` | Synchronize all processes |
| `comm_tag(src, tgt)` | Generate unique message tag |
| `comm_split(len)` | Split message into chunks |

### Communication Operations

| Operation | Types Supported |
|-----------|-----------------|
| **Broadcast** | Scalar, Cell, Vector, Matrix, Cell\<string\> |
| **Send/Receive** | Scalar, Raw array, Cell, Vector, Matrix, String |
| **Distribute** | Cell, Vector, Cell\<Vector\>, Cell\<Cell\>, Cell\<Matrix\>, Raw array |
| **Collect** | Cell, Vector, Cell\<Vector\>, Cell\<Matrix\>, Raw array |
| **Share** | Cell, Vector |

### Build Configuration

```bash
# MPI build (default if MPI available)
cmake -DUSE_MPI=ON ..

# Non-MPI build (serial fallback)
cmake -DUSE_MPI=OFF ..
```

**Open MPI is the only supported implementation.** MPICH and Intel MPI are untested — PETSc
has been observed to crash when called on that path — and the configure step refuses them.
Using the Intel compilers is fine; build Open MPI with them rather than substituting Intel
MPI. See [mpi_support.md](../../../doc/mpi_support.md) for the guard, the override, and why
the Open MPI link flags in the MUMPS and MKL configs must not be "made portable".

### Common Patterns

```cpp
// Initialize
gComm.init(argc, argv);

// Broadcast
Vector<real> v;
if (comm_rank() == 0) v = {1, 2, 3};
broadcast(v, 0);

// Send/Receive
if (comm_rank() == 0) send(data, 1);
if (comm_rank() == 1) receive(data, 0);

// Distribute + Collect: the two halves of one all-to-all exchange.
// distribute() only sends and collect() only receives, so every rank calls both.
Cell<int> send_data(comm_size(), 0);   // send_data(p) goes to rank p
distribute(send_data);
Cell<int> recv_data;                   // recv_data(p) arrives from rank p
collect(recv_data, my_value);          // my_value fills my own slot

// Finalize
return gComm.finalize();
```

## See Also

- **Source code** - header files in `src/comm/`
- [Containers module](../../containers/doc/README.md) - Cell container
- [Linear algebra module](../../linalg/doc/README.md) - Vector and Matrix types
- [Root documentation](../../../doc/README.md) - General BELFEM documentation
- [MPI Documentation](https://www.mpi-forum.org/docs/) - MPI standard reference
