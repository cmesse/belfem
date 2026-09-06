# Optimizer Module Documentation {#numerics_opt_index}

**Date:** 2026-08-31
**Purpose:** Index for BELFEM's optimization wrapper — a thin, NLOPT-backed interface for
bound-constrained scalar minimization
**Module:** src/numerics/opt

---

## Overview

`src/numerics/opt` wraps [NLOPT](https://nlopt.readthedocs.io) behind an interface that keeps
NLOPT's headers out of the rest of the tree. A caller derives from `Objective`, passing the problem dimension to its constructor, hands it to
an `Optimizer`, and reads a `Status` back — no `nlopt.h` anywhere above this directory.

It is built by default: `src/numerics/CMakeLists.txt:5` adds the subdirectory, and `USE_NLOPT`
defaults **ON** (`CMakeLists.txt:103`). The example executable `opt` is built only under
`USE_EXAMPLES AND USE_NLOPT`.

**Nothing in the FEM path calls it today.** It is available infrastructure, not part of a solve.

---

## Quick Reference

| Class / file | Purpose |
|---|---|
| `Optimizer` (`cl_Optimizer.hpp:38`) | Owns bounds, tolerances and the algorithm choice; runs `optimize()`. Constructed from an `Objective&` and an `Algorithm` (`:71`) |
| `Objective` (`cl_Objective.hpp:31`) | Abstract base — **construct it with the dimension** (`:40`) and implement `compute_objective()` (`:62`). `dimension()` is a non-virtual accessor, not an override point |
| `Algorithm` (`en_Opt_Algorithm.hpp:27`) | Which NLOPT algorithm to use |
| `Status` (`en_Opt_Status.hpp:26`) | Why the solver stopped, without exposing NLOPT's own codes |
| `opt.cpp` | Worked example; built only with `USE_EXAMPLES` and `USE_NLOPT` |

### Algorithms

Local, **derivative-free** — `compute_objective()` may ignore its gradient argument:

| | NLOPT | Notes |
|---|---|---|
| `BOBYQA` | `LN_BOBYQA` | quadratic model, bound constrained |
| `COBYLA` | `LN_COBYLA` | linear model, supports constraints |
| `NELDERMEAD` | `LN_NELDERMEAD` | simplex |
| `SBPLX` | `LN_SBPLX` | Rowan's subplex |
| `PRAXIS` | `LN_PRAXIS` | principal-axis |

Local, **gradient-based** — the implementation **must** fill the gradient:

| | NLOPT | Notes |
|---|---|---|
| `MMA` | `LD_MMA` | method of moving asymptotes |
| `SLSQP` | `LD_SLSQP` | sequential quadratic programming |
| `LBFGS` | `LD_LBFGS` | low-storage BFGS |

---

## Usage

```cpp
#include "cl_Optimizer.hpp"
#include "cl_Objective.hpp"

class MyObjective : public opt::Objective
{
public:
    // Objective has no default constructor -- the dimension is passed up.
    MyObjective() : Objective( 2 ) {}

    real compute_objective( const Vector< real > & aX,
                                  Vector< real > & aGradient ) override
    {
        // aGradient has length zero for the derivative-free algorithms and
        // must be ignored there; for the gradient-based ones it has length
        // dimension() and must be filled with d(objective)/d(aX).
        return aX( 0 ) * aX( 0 ) + aX( 1 ) * aX( 1 );
    }
};

MyObjective tObjective;
opt::Optimizer tOptimizer( tObjective, opt::Algorithm::BOBYQA );   // takes a reference
tOptimizer.set_bounds( { -1.0, -1.0 }, { 1.0, 1.0 } );

Vector< real > tX( 2, 0.5 );
real           tValue;

opt::Status tStatus = tOptimizer.optimize( tX, tValue );
```

---

## Development Notes

- **`aX` and `aValue` mean nothing unless the run succeeded.** Check the status through
  `is_usable( status )` (`en_Opt_Status.hpp:68`) before reading either; on failure they must not be trusted
  (`cl_Optimizer.hpp:100-106`).
- **`errmsg()` carries the solver's own diagnostic** for the most recent `optimize()` — for
  instance which bound violates `lb <= ub`. It is empty when nothing was reported, and complements
  the generic `error_message( status )`.
- **The gradient argument is the trap.** Its length depends on the algorithm class: zero for the
  `LN_*` family, `dimension()` for the `LD_*` family. An implementation that always writes to it
  will write out of bounds under a derivative-free algorithm.
- **Bounds are optional.** Leaving both vectors empty is an unbounded problem
  (`cl_Optimizer.hpp:47`).
- **Defaults:** `mFtolRel` and `mXtolRel` are `1e-8`, and `mMaximize` is `false`
  (`cl_Optimizer.hpp:52,55,61`) — the wrapper minimizes unless told otherwise.
- The `Status` enum mirrors NLOPT's termination reasons deliberately, so a caller can tell *why*
  the solver stopped without linking against or including NLOPT itself.

---

**See Also:**

- [Documentation index](../../../../doc/README.md)
- [Spline](../../spline/doc/README.md) — the other documented `numerics` submodule
