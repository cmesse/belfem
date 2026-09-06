"""Light-weight exporters: VTK legacy polyline (with frame data) and STL of the tape volume."""

import struct
import numpy as np

from .basecurve import Basecurve


def frame_table(curve: Basecurve, t=None):
    """Arrays r, T, n, b (m,3) and tau, kappa_g, kappa_n (m,) at the parameters t."""
    if t is None:
        t = curve.t
    m = len(t)
    r = np.zeros((m, 3)); T = np.zeros((m, 3)); n = np.zeros((m, 3)); b = np.zeros((m, 3))
    tau = np.zeros(m); kg = np.zeros(m); kn = np.zeros(m)
    for k, tk in enumerate(t):
        r[k] = curve.r(tk)
        T[k], n[k], b[k] = curve.frame(tk)
        tau[k], kg[k], kn[k] = curve.curvatures(tk)
    return r, T, n, b, tau, kg, kn


def write_vtk_polyline(curve: Basecurve, path: str, t=None):
    """Legacy ASCII VTK polydata: the curve as a polyline with T, n, b and curvatures as point data."""
    if t is None:
        t = curve.t
    r, T, n, b, tau, kg, kn = frame_table(curve, t)
    m = len(t)
    with open(path, "w") as f:
        f.write("# vtk DataFile Version 3.0\npancake basecurve\nASCII\nDATASET POLYDATA\n")
        f.write("POINTS %d double\n" % m)
        for p in r:
            f.write("%.12e %.12e %.12e\n" % tuple(p))
        f.write("LINES 1 %d\n%d " % (m + 1, m) + " ".join(str(i) for i in range(m)) + "\n")
        f.write("POINT_DATA %d\n" % m)
        for name, arr in (("T", T), ("n", n), ("b", b)):
            f.write("VECTORS %s double\n" % name)
            for p in arr:
                f.write("%.12e %.12e %.12e\n" % tuple(p))
        for name, arr in (("s", np.asarray(t, float)), ("tau", tau), ("kappa_g", kg), ("kappa_n", kn)):
            f.write("SCALARS %s double 1\nLOOKUP_TABLE default\n" % name)
            for x in arr:
                f.write("%.12e\n" % x)


def tape_corners(curve: Basecurve, width: float, thickness: float, t=None, theta_T=0.0):
    """Corner points (m, 4, 3) of the tape cross section along the curve.

    Local cross-section coordinates: x0 along b (width), y0 along -n (thickness),
    consistent with `Basecurve.transform`.  Corner order: (-w/2,-d/2), (w/2,-d/2),
    (w/2, d/2), (-w/2, d/2).
    """
    if t is None:
        t = curve.t
    local = np.array([[-0.5 * width, -0.5 * thickness, 0.0],
                      [0.5 * width, -0.5 * thickness, 0.0],
                      [0.5 * width, 0.5 * thickness, 0.0],
                      [-0.5 * width, 0.5 * thickness, 0.0]])
    out = np.zeros((len(t), 4, 3))
    for k, tk in enumerate(t):
        R = curve.transform(tk, theta_T)
        out[k] = curve.r(tk) + local @ R.T
    return out


def write_stl_tape(curve: Basecurve, width: float, thickness: float, path: str, t=None, theta_T=0.0):
    """Binary STL of the tape as a swept box (four side faces plus two end caps)."""
    C = tape_corners(curve, width, thickness, t, theta_T)
    m = len(C)
    tris = []
    for k in range(m - 1):
        for j in range(4):
            j1 = (j + 1) % 4
            p0, p1, p2, p3 = C[k, j], C[k, j1], C[k + 1, j1], C[k + 1, j]
            tris.append((p0, p1, p2))
            tris.append((p0, p2, p3))
    tris.append((C[0, 0], C[0, 2], C[0, 1]))
    tris.append((C[0, 0], C[0, 3], C[0, 2]))
    tris.append((C[-1, 0], C[-1, 1], C[-1, 2]))
    tris.append((C[-1, 0], C[-1, 2], C[-1, 3]))
    with open(path, "wb") as f:
        f.write(b"pancake tape".ljust(80, b"\0"))
        f.write(struct.pack("<I", len(tris)))
        for p0, p1, p2 in tris:
            nrm = np.cross(p1 - p0, p2 - p0)
            ln = np.linalg.norm(nrm)
            nrm = nrm / ln if ln > 0 else nrm
            f.write(struct.pack("<3f", *nrm))
            f.write(struct.pack("<3f", *p0))
            f.write(struct.pack("<3f", *p1))
            f.write(struct.pack("<3f", *p2))
            f.write(struct.pack("<H", 0))
    return len(tris)
