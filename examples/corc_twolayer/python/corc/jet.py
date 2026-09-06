"""Truncated Taylor arithmetic (order 3) for exact curve derivatives.

A Jet carries the Taylor coefficients [c0, c1, c2, c3] of a function at a
point (c_k = f^(k)/k!). Arithmetic on Jets propagates them exactly, so a
curve class only implements p(t) with Jet-friendly operations and gets
machine-precision velocity, acceleration and jerk for free — no
hand-derived formulas to get wrong.
"""

import math


class Jet:

    __slots__ = ("c",)

    def __init__(self, c0, c1=0.0, c2=0.0, c3=0.0):
        self.c = [float(c0), float(c1), float(c2), float(c3)]

    # independent variable t
    @staticmethod
    def variable(t):
        return Jet(t, 1.0)

    # derivatives from coefficients
    @property
    def value(self):
        return self.c[0]

    @property
    def d1(self):
        return self.c[1]

    @property
    def d2(self):
        return 2.0 * self.c[2]

    @property
    def d3(self):
        return 6.0 * self.c[3]

    # ---------------------------------------------------------------- basic

    def _lift(self, other):
        return other if isinstance(other, Jet) else Jet(other)

    def __add__(self, other):
        o = self._lift(other)
        return Jet(*[a + b for a, b in zip(self.c, o.c)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(*[-a for a in self.c])

    def __sub__(self, other):
        o = self._lift(other)
        return Jet(*[a - b for a, b in zip(self.c, o.c)])

    def __rsub__(self, other):
        return self._lift(other).__sub__(self)

    def __mul__(self, other):
        o = self._lift(other)
        a, b = self.c, o.c
        return Jet(
            a[0] * b[0],
            a[0] * b[1] + a[1] * b[0],
            a[0] * b[2] + a[1] * b[1] + a[2] * b[0],
            a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0])

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._lift(other)
        b = o.c
        if b[0] == 0.0:
            raise ZeroDivisionError("Jet division by zero value")
        q = [0.0] * 4
        for k in range(4):
            s = self.c[k]
            for i in range(k):
                s -= q[i] * b[k - i]
            q[k] = s / b[0]
        return Jet(*q)

    def __rtruediv__(self, other):
        return self._lift(other).__truediv__(self)

    # ------------------------------------------------------------ functions


def sin(f):
    if not isinstance(f, Jet):
        return math.sin(f)
    g = [math.sin(f.c[0]), 0.0, 0.0, 0.0]   # sin
    h = [math.cos(f.c[0]), 0.0, 0.0, 0.0]   # cos
    for k in range(1, 4):
        sg = 0.0
        sh = 0.0
        for i in range(1, k + 1):
            sg += i * f.c[i] * h[k - i]
            sh += i * f.c[i] * g[k - i]
        g[k] = sg / k
        h[k] = -sh / k
    return Jet(*g)


def cos(f):
    if not isinstance(f, Jet):
        return math.cos(f)
    g = [math.sin(f.c[0]), 0.0, 0.0, 0.0]
    h = [math.cos(f.c[0]), 0.0, 0.0, 0.0]
    for k in range(1, 4):
        sg = 0.0
        sh = 0.0
        for i in range(1, k + 1):
            sg += i * f.c[i] * h[k - i]
            sh += i * f.c[i] * g[k - i]
        g[k] = sg / k
        h[k] = -sh / k
    return Jet(*h)


def power(f, m):
    """f**m for real exponent m; requires f.value > 0 (used away from
    parametrization endpoints)."""
    if not isinstance(f, Jet):
        return math.pow(f, m)
    if f.c[0] <= 0.0:
        raise ValueError("Jet power needs positive base, got {:g}".format(f.c[0]))
    p = [math.pow(f.c[0], m), 0.0, 0.0, 0.0]
    for k in range(1, 4):
        s = 0.0
        for i in range(1, k + 1):
            s += ((m + 1.0) * i - k) * f.c[i] * p[k - i]
        p[k] = s / (k * f.c[0])
    return Jet(*p)


def sqrt(f):
    return power(f, 0.5)
