# Database Module Documentation {#physics_database_index}

**Module:** src/physics/database
**Purpose:** Index of documentation for BELFEM's precomputed lookup-table module

---

## Overview

The `database` module turns an expensive property function into a cheap one. A property is
sampled once onto a structured tensor grid, projected onto a smooth basis, stored, and
thereafter read by shape-function interpolation — so the cost is paid at setup instead of at
every quadrature point of every Newton iteration.

Its principal consumer is `physics/materials`, where the field- and temperature-dependent
resistivity of normal metals is served from such a table.

---

## Documents

| Document | Read it for |
|---|---|
| [database_usage_guide.md](database_usage_guide.md) | what the module is, how to build and query a table, **why the projection step exists**, the parallel contract, and the consumer pitfalls |

---

## Quick Reference

| Class | File | Role |
|---|---|---|
| `Database` | `cl_Database.hpp` / `.cpp` | holds the finished table; `evaluate`, `evaluate_deriv{x,y,z}`, `min`, `max`, `save` |
| `database::Projector` | `cl_DatabaseProjector.hpp` / `.cpp` | build-time only: L2-projects sampled values onto a B-spline basis and returns node values |

| Question | Answer |
|---|---|
| Why project at all, instead of storing the samples? | so the **derivative** accessors return a continuous field — raw Lagrange interpolation of samples is only C⁰, and its jumps degrade Newton tangents |
| Is the build parallel? | the call is collective, the work is not — rank 0 builds, all ranks receive |
| Which solver? | compile-time preference MUMPS → PARDISO → SUPERLU → UMFPACK; STRUMPACK deliberately excluded (matrices too small) |
| Who clamps the query to the grid? | the caller, always |

---

## Development Notes

- The mesh constructor requires a **tensor mesh**; it checks `is_tensormesh()` and raises `BELFEM_ERROR` otherwise.
- Stored values may be a transformed quantity (materials store `log(rho)`) — undo the
  transform on the way out, and remember a stored derivative is the derivative of the
  transform.
- Rebuilt tables agree bit-for-bit only single-threaded; multithreaded solves carry ordinary
  reassociation noise of order 1e-7 relative.
- Saved HDF5 records carry no format-version stamp, so a caching consumer must probe for an
  expected dataset before trusting a file on disk.
