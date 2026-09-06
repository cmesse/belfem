#!/usr/bin/env python3
"""
howtoread.py -- the reference reader for BELFEM material database files.

THIS FILE IS THE ROSETTA STONE for every ``*.hdf5`` material table that BELFEM
ships or generates. It reproduces, in plain numpy, exactly what the solver does
in C++ (``src/physics/database/cl_Database.hpp``), so that a table can be used,
checked, or plotted without any BELFEM installation. Only ``numpy`` and
``h5py`` are required.

File format
-----------

A material database holds one HDF5 group per tabulated property. Each group
contains five datasets that define a structured tensor-product mesh::

    order    uint32      polynomial order per axis: 1, 2 or 3
    origin   float64[d]  first node coordinate on each axis (d = 2 or 3 axes)
    points   uint32[d]   number of nodes on each axis
    step     float64[d]  node spacing on each axis
    values   float64[n]  nodal values, flattened in C order
                         (axis 0 slowest, last axis fastest)

The interpolant is piecewise-Lagrange of the given order on each axis
(order 2 in three dimensions is the familiar HEX27 element). ``values`` stores
a LOGARITHM of the physical quantity -- which one depends on the family:

============================  =========================  =====================
table family                  axes (in order)            values
============================  =========================  =====================
superconductor ``jc``, ``n``  T [K], log10(B/T),         log10(jc [A/m^2]),
(e.g. bscco-2223, sp-ap,      theta [rad] from the       log10(n)
sst-1)                        tape NORMAL
metal / alloy ``rho``         T [K], log10(B/T),         ln(rho [Ohm m])
(e.g. Copper_RRR100)          beta [rad]                 (NATURAL log!)
============================  =========================  =====================

Some files carry additional self-describing groups: ``meta`` (provenance,
DOIs, accuracy, an embedded markdown ``documentation``), ``source`` (the raw
measurement export) and ``python`` (the generator pipeline, or a copy of this
file). None of these are read by the solver. Strings written by the BELFEM
C++ side are fixed-size datasets of shape (1,); to extract one with h5py::

    text = np.ravel(f['python/howtoread.py'][()])[0].decode()

Classes
-------

``Database``
    the raw interpolant: ``evaluate`` and the partial derivatives, exactly
    like ``belfem::Database``. No clamping, no logarithm handling.
``JcFunctionDatabase``
    physical jc/n lookup: clamps T and B into the table window, wraps the
    angle with the +-pi periodicity, returns 10**(interpolated value) --
    mirrors ``src/physics/materials/cl_JcFunction_Database.hpp``.
``MetalRhoDatabase``
    physical resistivity lookup for metal/alloy tables: clamps T and B to
    the ``Tmin/Tmax/Bmin/Bmax`` scalars stored beside the group, wraps beta,
    returns exp(interpolated value) -- mirrors the ``rho_table`` family in
    ``src/physics/materials/cl_Material_Metal.hpp``.

Run this file directly to see a demonstration::

    python3 howtoread.py path/to/table.hdf5

All inputs broadcast: scalars and numpy arrays both work.
"""

import sys

import numpy as np
import h5py


# ----------------------------------------------------------------------------
# 1D Lagrange bases on the reference interval [-1, 1], equidistant nodes.
# Order 1: nodes (-1, 1); order 2: (-1, 0, 1); order 3: (-1, -1/3, 1/3, 1).
# These are the axis factors of BELFEM's QUAD4/9/16 and HEX8/27/64 shape
# functions; the tensor product over the axes reproduces them identically.
# ----------------------------------------------------------------------------

def _basis_1(x):
    return np.stack([0.5 * (1.0 - x), 0.5 * (1.0 + x)], axis=-1)


def _dbasis_1(x):
    o = np.ones_like(x)
    return np.stack([-0.5 * o, 0.5 * o], axis=-1)


def _basis_2(x):
    return np.stack([0.5 * x * (x - 1.0),
                     1.0 - x * x,
                     0.5 * x * (x + 1.0)], axis=-1)


def _dbasis_2(x):
    return np.stack([x - 0.5, -2.0 * x, x + 0.5], axis=-1)


def _basis_3(x):
    x2 = x * x
    return np.stack([0.0625 * (-9.0 * x2 * x + 9.0 * x2 + x - 1.0),
                     0.0625 * (27.0 * x2 * x - 9.0 * x2 - 27.0 * x + 9.0),
                     0.0625 * (-27.0 * x2 * x - 9.0 * x2 + 27.0 * x + 9.0),
                     0.0625 * (9.0 * x2 * x + 9.0 * x2 - x - 1.0)], axis=-1)


def _dbasis_3(x):
    x2 = x * x
    return np.stack([0.0625 * (-27.0 * x2 + 18.0 * x + 1.0),
                     0.0625 * (81.0 * x2 - 18.0 * x - 27.0),
                     0.0625 * (-81.0 * x2 - 18.0 * x + 27.0),
                     0.0625 * (27.0 * x2 + 18.0 * x - 1.0)], axis=-1)


_BASIS = {1: (_basis_1, _dbasis_1),
          2: (_basis_2, _dbasis_2),
          3: (_basis_3, _dbasis_3)}


class Database:
    """Raw tensor-mesh interpolant, the python twin of ``belfem::Database``.

    ``evaluate(x, y)`` / ``evaluate(x, y, z)`` interpolates the stored nodal
    values; ``evaluate_derivx/y/z`` return the partial derivatives with
    respect to the PHYSICAL axis coordinates (the 2/element_step scaling of
    the reference derivative is applied, as in the C++).

    Out-of-range coordinates are NOT clamped: the containing boundary element
    extrapolates its polynomial, exactly as the C++ does when a caller skips
    the clamping of the physical wrapper classes. The element INDEX is
    clamped, so the evaluation never reads outside ``values``.
    """

    def __init__(self, path_or_file, label):
        if isinstance(path_or_file, h5py.Group):
            g = path_or_file[label]
            self._load(g)
        else:
            with h5py.File(path_or_file, 'r') as f:
                self._load(f[label])
        self.label = label

    def _load(self, g):
        self.order = int(np.ravel(g['order'][()])[0])
        self.origin = np.asarray(g['origin'][()], float)
        self.step = np.asarray(g['step'][()], float)
        self.points = np.asarray(np.ravel(g['points'][()]), int)
        self.dimension = len(self.points)
        if self.order not in _BASIS:
            raise ValueError('unsupported order %d' % self.order)
        if self.dimension not in (2, 3):
            raise ValueError('unsupported dimension %d' % self.dimension)
        if ((self.points - 1) % self.order).any():
            raise ValueError('points %s incompatible with order %d'
                             % (self.points, self.order))
        # C-order reshape: axis 0 slowest -- matches Database::nidx
        self.values = np.asarray(g['values'][()], float).reshape(
            tuple(self.points))
        self.num_elems = (self.points - 1) // self.order
        if (self.num_elems < 1).any():
            raise ValueError('degenerate mesh: %s nodes at order %d'
                             % (self.points, self.order))
        self.elem_step = self.order * self.step
        # multiply by the inverse, as the C++ does: at an exact element
        # interface a divide can differ by one ulp and pick the neighbor
        # element, which matters for the (discontinuous) derivatives
        self.inv_elem_step = 1.0 / self.elem_step

    # -- coordinate window, like Database::min/max ---------------------------

    def min(self, dim):
        return float(self.origin[dim])

    def max(self, dim):
        return float(self.origin[dim]
                     + self.elem_step[dim] * self.num_elems[dim])

    # -- element location ----------------------------------------------------

    def _locate(self, dim, x):
        """Element index (clamped) and reference coordinate xi (NOT clamped).

        The C++ casts ``(x - origin) * inv_element_step`` to an unsigned
        integer before clamping; for in-range coordinates the result is
        identical to flooring. Here the FLOAT is clamped first and cast
        afterwards, which agrees with the C++ everywhere the C++ is defined
        and additionally gives coordinates below the origin by a full
        element or more the well-defined meaning "first element
        extrapolates" (the C++ cast is undefined there; its callers clamp
        their inputs beforehand).
        """
        t = (np.asarray(x, float) - self.origin[dim]) * self.inv_elem_step[dim]
        e = np.clip(t, 0, self.num_elems[dim] - 1).astype(int)
        xc = self.origin[dim] + (e + 0.5) * self.elem_step[dim]
        xi = 2.0 * (np.asarray(x, float) - xc) * self.inv_elem_step[dim]
        return e, xi

    # -- evaluation ----------------------------------------------------------

    def _eval(self, coords, deriv=None):
        basis, dbasis = _BASIS[self.order]
        p = self.order
        loc = [self._locate(d, c) for d, c in enumerate(coords)]
        fac = [(dbasis if d == deriv else basis)(xi) for d, (e, xi) in
               enumerate(loc)]
        shape = np.broadcast(*[e for e, _ in loc]).shape
        out = np.zeros(shape)
        if self.dimension == 2:
            (ei, _), (ej, _) = loc
            for a in range(p + 1):
                for b in range(p + 1):
                    out += (fac[0][..., a] * fac[1][..., b]
                            * self.values[p * ei + a, p * ej + b])
        else:
            (ei, _), (ej, _), (ek, _) = loc
            for a in range(p + 1):
                for b in range(p + 1):
                    for c in range(p + 1):
                        out += (fac[0][..., a] * fac[1][..., b]
                                * fac[2][..., c]
                                * self.values[p * ei + a, p * ej + b,
                                              p * ek + c])
        if deriv is not None:
            out *= 2.0 / self.elem_step[deriv]
        return out if out.shape else float(out)

    def evaluate(self, x, y, z=None):
        return self._eval((x, y) if z is None else (x, y, z))

    def evaluate_derivx(self, x, y, z=None):
        return self._eval((x, y) if z is None else (x, y, z), deriv=0)

    def evaluate_derivy(self, x, y, z=None):
        return self._eval((x, y) if z is None else (x, y, z), deriv=1)

    def evaluate_derivz(self, x, y, z):
        return self._eval((x, y, z), deriv=2)


_LN10 = np.log(10.0)


class JcFunctionDatabase:
    """Physical jc / n lookup on a superconductor table (log10 storage).

    Mirrors ``belfem::material::JcFunctionDatabase``: T and B are clamped
    into the table window, the angle theta (measured from the tape normal)
    is wrapped with the +-pi periodicity, values within 1e-6 rad of the
    window are snapped onto the boundary, and the returned quantity is
    ``10**(interpolated value)``. The derivative methods implement the same
    clamp-consistent tangents as the C++ (zero outside the T and B windows).
    """

    def __init__(self, path, label):
        self.db = Database(path, label)
        self.Tmin, self.Tmax = self.db.min(0), self.db.max(0)
        self.Bmin, self.Bmax = 10.0 ** self.db.min(1), 10.0 ** self.db.max(1)
        self.amin, self.amax = self.db.min(2), self.db.max(2)
        if self.amin > 1e-6 or self.amax < np.pi - 1e-6:
            raise ValueError('angular range [%g, %g] does not cover [0, pi]'
                             % (self.amin, self.amax))

    def wrap_angle(self, angle):
        # exact transcription of the C++ ternary: ONE +-pi shift for angles
        # outside the snap band, and the boundary snap (clamp) ONLY for
        # angles already inside it. An angle that is still outside the window
        # after one shift extrapolates, exactly as the solver does -- it is
        # NOT clipped (audit finding R2, 2026-08-28).
        a = np.asarray(angle, float)
        return np.where(a < self.amin - 1e-6, a + np.pi,
                        np.where(a > self.amax + 1e-6, a - np.pi,
                                 np.clip(a, self.amin, self.amax)))

    def _args(self, normB, angle, T):
        lb = np.log10(np.clip(normB, self.Bmin, self.Bmax))
        return np.clip(T, self.Tmin, self.Tmax), lb, self.wrap_angle(angle)

    def eval(self, normB, angle, T):
        return 10.0 ** self.db.evaluate(*self._args(normB, angle, T))

    def deval_dB(self, normB, angle, T):
        v = self.eval(normB, angle, T)
        t, lb, a = self._args(normB, angle, T)
        B = np.asarray(normB, float)
        inside = (B > self.Bmin) & (B < self.Bmax)
        # mask the denominator first: at B = 0 the quotient is never used,
        # but an unmasked divide would still raise under np.seterr
        d = v * self.db.evaluate_derivy(t, lb, a) / np.where(inside, B, 1.0)
        return np.where(inside, d, 0.0)

    def deval_dT(self, normB, angle, T):
        v = self.eval(normB, angle, T)
        t, lb, a = self._args(normB, angle, T)
        inside = (np.asarray(T, float) > self.Tmin) \
            & (np.asarray(T, float) < self.Tmax)
        d = v * self.db.evaluate_derivx(t, lb, a) * _LN10
        return np.where(inside, d, 0.0)

    def deval_dbeta(self, normB, angle, T):
        v = self.eval(normB, angle, T)
        t, lb, a = self._args(normB, angle, T)
        return v * self.db.evaluate_derivz(t, lb, a) * _LN10


class MetalRhoDatabase:
    """Physical resistivity lookup on a metal / alloy table (ln storage).

    Mirrors the ``rho_table`` family of ``belfem::material::Metal``: T and B
    are clamped to the ``Tmin/Tmax/Bmin/Bmax`` scalars stored in the file
    root, the field angle beta is wrapped by +-pi, and the returned quantity
    is ``exp(interpolated value)`` in Ohm m. Note the NATURAL logarithm:
    metal tables store ln(rho), not log10.
    """

    def __init__(self, path, label='rho'):
        with h5py.File(path, 'r') as f:
            self.db = Database(f, label)
            self.Tmin = float(np.ravel(f['Tmin'][()])[0])
            self.Tmax = float(np.ravel(f['Tmax'][()])[0])
            self.Bmin = float(np.ravel(f['Bmin'][()])[0])
            self.Bmax = float(np.ravel(f['Bmax'][()])[0])

    def _args(self, T, B, beta):
        t = np.clip(T, self.Tmin, self.Tmax)
        lb = np.log(np.clip(B, self.Bmin, self.Bmax)) / _LN10
        b = np.asarray(beta, float)
        a = np.where(b < 0.0, b + np.pi, np.where(b > np.pi, b - np.pi, b))
        return t, lb, a

    def rho(self, T, B, beta):
        return np.exp(self.db.evaluate(*self._args(T, B, beta)))

    def drhodT(self, T, B, beta):
        t, lb, a = self._args(T, B, beta)
        return np.exp(self.db.evaluate(t, lb, a)) \
            * self.db.evaluate_derivx(t, lb, a)

    def drhodB(self, T, B, beta):
        t, lb, a = self._args(T, B, beta)
        Bc = np.clip(B, self.Bmin, self.Bmax)
        return np.exp(self.db.evaluate(t, lb, a)) \
            * self.db.evaluate_derivy(t, lb, a) / (_LN10 * Bc)

    def drhodbeta(self, T, B, beta):
        t, lb, a = self._args(T, B, beta)
        return np.exp(self.db.evaluate(t, lb, a)) \
            * self.db.evaluate_derivz(t, lb, a)


# ----------------------------------------------------------------------------

def _demo(path):
    with h5py.File(path, 'r') as f:
        groups = [k for k in f
                  if isinstance(f[k], h5py.Group) and 'values' in f[k]
                  and k not in ('meta', 'source', 'python')]
        has_limits = all(k in f for k in ('Tmin', 'Tmax', 'Bmin', 'Bmax'))
        meta = 'meta' in f
    print('%s: property groups %s' % (path, groups))
    if meta:
        print('  self-describing: see the meta group '
              '(meta/documentation holds the construction record)')
    for grp in groups:
        if has_limits:
            m = MetalRhoDatabase(path, grp)
            T = 0.5 * (m.Tmin + m.Tmax)
            print('  %s: rho(%.0f K, 1 T, 0) = %.6e Ohm m, '
                  'drho/dT = %.3e' % (grp, T, m.rho(T, 1.0, 0.0),
                                      m.drhodT(T, 1.0, 0.0)))
        else:
            j = JcFunctionDatabase(path, grp)
            T = 0.5 * (j.Tmin + j.Tmax)
            print('  %s: value(1 T, 0 rad, %.0f K) = %.6e, '
                  'd/dT = %.3e' % (grp, T, j.eval(1.0, 0.0, T),
                                   j.deval_dT(1.0, 0.0, T)))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    _demo(sys.argv[1])
