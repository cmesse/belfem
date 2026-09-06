# Automatic Sub-Layering for Thick Thin-Shell Layers

**Date:** 2026-03-14
**Purpose:** Future enhancement — automatic geometric subdivision of thick layers
**Status:** Design phase (weekend thinking)
**Module:** `src/mesh/cl_ThinShellFactory.cpp`

---

## Motivation

When thin-shell edges are properly deduplicated (edge key fix), thick layers (e.g., 50 μm hastelloy) lack sufficient through-thickness DOFs. Manual sub-layering works but is tedious. Automatic subdivision with geometric grading would provide adequate resolution with minimal user input.

## Design Concept

### Bi-directional geometric grading

Refine at both interfaces, coarsen in the center. For a layer of total thickness `d` split into `N` sub-layers with growth ratio `r`:

**Half-layer (one side):** `N/2` sub-layers with geometric progression
- `d₁, d₁·r, d₁·r², ..., d₁·r^(N/2-1)`
- where `d₁ = d_half · (1-r) / (1-r^(N/2))`

Mirror for the other half. If `N` is odd, the center sub-layer is shared.

### Candidate ratios

- **φ = 1.618...** (golden ratio) — aesthetically natural transition; each sub-layer relates to the next by the most irrational number, giving the smoothest possible grading
- **Iteratively determined** — given `N` and `d`, find `r` such that the thinnest sub-layer meets a minimum thickness constraint (e.g., ≥ the adjacent material layer thickness)

### Input syntax (tentative)

```
// Option A: specify number of sub-layers (ratio = φ by default)
hastelloy : 50 mum / 10 layers ;

// Option B: specify ratio (compute N from minimum thickness)
hastelloy : 50 mum / ratio 1.618 ;

// Option C: explicit (current behavior)
hastelloy : 50 mum ;
```

### Iterative ratio determination

Given total thickness `d` and desired number of layers `N`:
1. Target: thinnest sub-layer `d_min` ≈ thickness of adjacent layer (e.g., 1 μm for YBCO neighbor)
2. Solve: `d_min · (1 + r + r² + ... + r^(N/2-1)) = d/2` for `r`
3. This is `d_min · (r^(N/2) - 1) / (r - 1) = d/2`
4. Newton iteration on `r`

Or conversely: given `r` and `d`, compute `N` such that `d_min` doesn't go below a threshold.

## Implementation Notes

- Changes confined to `ThinShellFactory::create()` — expand a single thick layer specification into multiple sub-layers before passing to `create_nodes_on_layers()`
- The layer `Material` stays the same for all sub-layers
- No changes needed downstream (hanging edges, assembly, etc.)
- Should preserve backward compatibility (single layer = no subdivision)
