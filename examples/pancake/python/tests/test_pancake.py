"""Consistency checks for the pancake geometry (run: python -m unittest discover -s tests)."""

import os
import sys
import math
import unittest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import pancake as pc  # noqa: E402

SCRATCH = os.environ.get("PANCAKE_SCRATCH", os.path.join(os.path.dirname(__file__), "_out"))


def fd(fun, s, h=1e-4):
    return (fun(s + h) - fun(s - h)) / (2.0 * h)


def relerr(x, y):
    return np.linalg.norm(x - y) / max(1.0, np.linalg.norm(y))


def make_coil(sense=+1, ramp=8.0, release=6.0):
    spiral = pc.LameSpiral(a=80.0, b=120.0, n=4, pitch=0.6, nturns=2.25, theta0=0.3, sense=sense)
    inner = [pc.Straight(25.0), pc.Twist(90.0, 40.0, ramp=ramp), pc.Bend(90.0, 30.0, ramp=ramp), pc.Straight(60.0)]
    outer = [pc.Straight(30.0), pc.Twist(-90.0, 40.0, ramp=ramp), pc.Bend(-90.0, 25.0, ramp=ramp),
             pc.Bend(30.0, 40.0, ramp=ramp, hard=True), pc.Straight(50.0)]
    return pc.PancakeCurve(spiral, inner, outer, release=release)


def jump(fun, s, h=1e-9):
    """Difference between the values just left and right of s (a jump if >> 2h |f'|)."""
    return np.linalg.norm(np.asarray(fun(s + h)) - np.asarray(fun(s - h)))


class TestSpiral(unittest.TestCase):

    def setUp(self):
        self.sp = pc.LameSpiral(a=100.0, b=150.0, n=4, pitch=0.5, nturns=3)

    def test_unit_speed_and_derivatives(self):
        sp = self.sp
        for s in np.linspace(1.0, sp.tmax - 1.0, 17):
            self.assertAlmostEqual(np.linalg.norm(sp.v(s)), 1.0, places=12)
            self.assertLess(relerr(fd(sp.r, s), sp.v(s)), 1e-7)
            self.assertLess(relerr(fd(sp.v, s), sp.a(s)), 1e-6)
            self.assertLess(relerr(fd(sp.a, s), sp.b(s)), 1e-5)

    def test_theta_inversion(self):
        sp = self.sp
        for th in np.linspace(sp.theta0, sp.theta1, 13):
            s = sp.s_of_theta(th)
            self.assertAlmostEqual(sp.theta_of_s(s), th, places=10)

    def test_frame_and_curvatures(self):
        sp = self.sp
        for s in np.linspace(0.0, sp.tmax, 11):
            T, n, b = sp.frame(s)
            self.assertLess(abs(np.dot(T, n)), 1e-13)
            self.assertLess(abs(np.dot(T, b)), 1e-13)
            self.assertLess(abs(np.dot(n, b)), 1e-13)
            self.assertLess(np.linalg.norm(np.cross(T, n) - b), 1e-13)
            tau, kg, kn = sp.curvatures(s)
            self.assertEqual(tau, 0.0)
            self.assertEqual(kg, 0.0)
            # kappa_n n must equal the curvature vector a (planar, kappa_g = 0)
            self.assertLess(np.linalg.norm(kn * n - sp.a(s)), 1e-12)
            self.assertAlmostEqual(abs(kn), sp.kappa_tau(s)[0], places=12)

    def test_curvature_vanishes_on_axes_but_frame_is_defined(self):
        sp = self.sp
        s = sp.s_of_theta(0.5 * math.pi)
        self.assertLess(abs(sp.curvatures(s)[2]), 1e-8)
        T, n, b = sp.frame(s)
        self.assertLess(np.linalg.norm(b - [0, 0, 1]), 1e-14)
        with self.assertRaises(ValueError):
            pc.Basecurve.frame(sp, s)

    def test_turn_gap_is_pitch(self):
        sp = self.sp
        for th in np.linspace(0.0, 2.0 * math.pi, 9):
            self.assertAlmostEqual(sp.turn_gap(th), sp.pitch, delta=1e-6 * sp.pitch)

    def test_sense(self):
        cw = pc.LameSpiral(a=100.0, b=150.0, n=4, pitch=0.5, nturns=1, sense=-1)
        ccw = pc.LameSpiral(a=100.0, b=150.0, n=4, pitch=0.5, nturns=1, sense=+1)
        for s in np.linspace(0.0, cw.tmax, 7):
            r1, r2 = cw.r(s), ccw.r(s)
            self.assertLess(np.linalg.norm(r1 - r2 * [1, -1, 1]), 1e-9)
            self.assertAlmostEqual(cw.curvatures(s)[2], -ccw.curvatures(s)[2], places=10)

    def test_non_even_exponent_runs(self):
        sp = pc.LameSpiral(a=100.0, b=100.0, n=2.5, pitch=0.5, nturns=1)
        s = sp.s_of_theta(0.5 * math.pi)
        for fun in (sp.r, sp.v, sp.a):
            self.assertTrue(np.all(np.isfinite(fun(s))))
        self.assertLess(relerr(fd(sp.r, 100.0), sp.v(100.0)), 1e-7)


class _LeadCurve(pc.Basecurve):
    """Adapter exposing a Lead through the Basecurve interface."""

    def __init__(self, lead):
        pc.Basecurve.__init__(self, 0.0, lead.length)
        self.lead = lead

    def r(self, s): return self.lead.derivatives(s)[0]
    def v(self, s): return self.lead.derivatives(s)[1]
    def a(self, s): return self.lead.derivatives(s)[2]
    def b(self, s): return self.lead.derivatives(s)[3]


class TestLead(unittest.TestCase):

    def test_release_fades_curvature(self):
        seg = pc.Release(10.0, kappa_n=0.1)
        self.assertEqual(seg.omega(0.0)[2], 0.1)
        self.assertEqual(seg.omega(10.0)[2], 0.0)
        self.assertAlmostEqual(seg.total_rotation[2], 0.5)

    def test_segment_factories(self):
        b = pc.Bend(90.0, 30.0, ramp=6.0)
        self.assertAlmostEqual(b.total_rotation[2], 0.5 * math.pi)
        self.assertAlmostEqual(b.length, 0.5 * math.pi * 30.0 + 6.0)
        t = pc.Twist(-45.0, 20.0, ramp=4.0)
        self.assertAlmostEqual(t.total_rotation[0], -0.25 * math.pi)
        with self.assertRaises(ValueError):
            pc.Segment(10.0, ramp=6.0)

    def test_pure_bend_is_circular_arc(self):
        lead = pc.Lead([pc.Bend(90.0, 30.0)]).build([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
        r, T, n, b = lead.end_state()
        self.assertLess(np.linalg.norm(r - [30.0, 30.0, 0.0]), 1e-12)
        self.assertLess(np.linalg.norm(T - [0.0, 1.0, 0.0]), 1e-12)
        self.assertLess(np.linalg.norm(b - [0.0, 0.0, 1.0]), 1e-12)
        rm, Tm, nm, bm = lead.state(0.5 * lead.length)
        self.assertLess(np.linalg.norm(rm - [30.0 * math.sin(math.pi / 4), 30.0 * (1 - math.cos(math.pi / 4)), 0.0]), 1e-12)
        # helix: bend + twist + hard bend at once, compare with the Frenet formulas of Sec. 3.1
        lead = pc.Lead([pc.Segment(100.0, tau=0.02, kappa_n=0.05)]).build([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
        for s in (10.0, 55.5, 100.0):
            r, T, n, b = lead.state(s)
            kappa, tau = pc.Basecurve.kappa_tau(_LeadCurve(lead), s)
            self.assertAlmostEqual(kappa, 0.05, places=12)
            self.assertAlmostEqual(tau, 0.02, places=12)

    def test_twist_rotates_frame_only(self):
        lead = pc.Lead([pc.Twist(90.0, 40.0, ramp=10.0)]).build([1, 2, 3], [0, 0, 1], [1, 0, 0], [0, 1, 0])
        r, T, n, b = lead.end_state()
        self.assertLess(np.linalg.norm(r - [1, 2, 43]), 1e-12)
        self.assertLess(np.linalg.norm(T - [0, 0, 1]), 1e-12)
        self.assertLess(np.linalg.norm(n - [0, 1, 0]), 1e-9)
        self.assertLess(np.linalg.norm(b - [-1, 0, 0]), 1e-9)

    def test_strip_equations_and_derivatives(self):
        segs = [pc.Straight(10.0), pc.Twist(90.0, 30.0, ramp=5.0), pc.Bend(90.0, 20.0, ramp=5.0),
                pc.Bend(40.0, 30.0, ramp=4.0, hard=True), pc.Segment(20.0, tau=0.02, kappa_g=0.01, kappa_n=0.03, ramp=3.0)]
        lead = pc.Lead(segs).build([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1])
        for s in np.linspace(0.5, lead.length - 0.5, 41):
            r, T, n, b = lead.state(s)
            F = np.column_stack([T, n, b])
            self.assertLess(np.linalg.norm(F.T @ F - np.eye(3)), 1e-12)
            self.assertAlmostEqual(np.linalg.det(F), 1.0, places=12)
            tau, kg, kn = lead.omega(s)
            dT = fd(lambda x: lead.state(x)[1], s)
            dn = fd(lambda x: lead.state(x)[2], s)
            db = fd(lambda x: lead.state(x)[3], s)
            self.assertLess(np.linalg.norm(dT - (kn * n - kg * b)), 1e-7)
            self.assertLess(np.linalg.norm(dn - (-kn * T + tau * b)), 1e-7)
            self.assertLess(np.linalg.norm(db - (kg * T - tau * n)), 1e-7)
            # Eqs. (19.14)-(19.16)
            self.assertAlmostEqual(np.dot(b, dn), tau, places=7)
            self.assertAlmostEqual(np.dot(T, db), kg, places=7)
            self.assertAlmostEqual(np.dot(n, dT), kn, places=7)
            # derivatives of the position
            rr, v, a, j = lead.derivatives(s)
            self.assertLess(relerr(fd(lambda x: lead.state(x)[0], s), v), 1e-8)
            self.assertLess(relerr(fd(lambda x: lead.derivatives(x)[1], s), a), 1e-7)
            self.assertLess(relerr(fd(lambda x: lead.derivatives(x)[2], s), j), 1e-6)


class TestPancake(unittest.TestCase):

    def setUp(self):
        self.coil = make_coil()

    def test_layout(self):
        c = self.coil
        self.assertAlmostEqual(c.tmax, c.L_inner + c.L_spiral + c.L_outer)
        self.assertEqual(c.piece(0.5 * c.L_inner), "inner")
        self.assertEqual(c.piece(c.s_spiral0 + 1.0), "spiral")
        self.assertEqual(c.piece(c.tmax - 1.0), "outer")

    def test_derivatives_everywhere(self):
        c = self.coil
        for s in np.linspace(0.5, c.tmax - 0.5, 60):
            self.assertAlmostEqual(np.linalg.norm(c.v(s)), 1.0, places=10)
            self.assertLess(relerr(fd(c.r, s), c.v(s)), 1e-7, msg="v at s=%g" % s)
            self.assertLess(relerr(fd(c.v, s), c.a(s)), 1e-6, msg="a at s=%g" % s)
            self.assertLess(relerr(fd(c.a, s), c.b(s)), 1e-5, msg="b at s=%g" % s)

    def test_frame_consistency_everywhere(self):
        c = self.coil
        for s in np.linspace(0.5, c.tmax - 0.5, 60):
            T, n, b = c.frame(s)
            F = np.column_stack([T, n, b])
            self.assertLess(np.linalg.norm(F.T @ F - np.eye(3)), 1e-12)
            self.assertAlmostEqual(np.linalg.det(F), 1.0, places=12)
            self.assertLess(np.linalg.norm(T - c.v(s)), 1e-12)
            tau, kg, kn = c.curvatures(s)
            # a = T' = kappa_n n - kappa_g b   (Eq. 19.8) must hold for the reported curvatures
            self.assertLess(np.linalg.norm(c.a(s) - (kn * n - kg * b)), 1e-9, msg="s=%g" % s)
            dn = fd(lambda x: c.frame(x)[1], s)
            db = fd(lambda x: c.frame(x)[2], s)
            self.assertAlmostEqual(np.dot(b, dn), tau, places=7, msg="tau at s=%g" % s)
            self.assertAlmostEqual(np.dot(T, db), kg, places=7, msg="kappa_g at s=%g" % s)
            R = c.transform(s)
            self.assertAlmostEqual(np.linalg.det(R), 1.0, places=12)
            self.assertLess(np.linalg.norm(R.T @ R - np.eye(3)), 1e-12)

    def test_continuity_at_breakpoints(self):
        c = self.coil
        self.assertEqual(len(c.breakpoints), len(c.inner.segments) + 2 + len(c.outer.segments))
        for s in c.breakpoints:
            if s <= c.tmin or s >= c.tmax:
                continue
            for fun, tol, name in ((c.r, 1e-8, "r"), (c.v, 1e-8, "v"), (c.a, 1e-7, "a"), (c.b, 1e-6, "jerk"),
                                   (lambda x: c.frame(x)[1], 1e-8, "n"), (lambda x: c.frame(x)[2], 1e-8, "b"),
                                   (lambda x: c.curvatures(x), 1e-7, "curv")):
                if name == "jerk" and s in (c.s_spiral0, c.s_spiral1):
                    continue   # winding meets the release with a different d(kappa)/ds: C^2 only
                self.assertLess(jump(fun, s), tol, msg="%s jumps at s=%g" % (name, s))

    def test_without_release_only_curvature_jumps(self):
        c = make_coil(release=0.0)
        for s in (c.s_spiral0, c.s_spiral1):
            for fun in (c.r, c.v, lambda x: c.frame(x)[1], lambda x: c.frame(x)[2]):
                self.assertLess(jump(fun, s), 1e-8)
            self.assertGreater(jump(c.curvatures, s), 1e-4)

    def test_leads_are_attached_to_the_spiral(self):
        c = self.coil
        r0 = c.spiral.r(c.spiral.tmin)
        r1 = c.spiral.r(c.spiral.tmax)
        self.assertLess(np.linalg.norm(c.r(c.s_spiral0) - r0), 1e-12)
        self.assertLess(np.linalg.norm(c.r(c.s_spiral1) - r1), 1e-12)
        # inner lead: after straight + twist + 90 degree bend the tangent points along +z
        s_after_bend = c.L_inner - c.inner.offsets[4]
        T = c.frame(s_after_bend)[0]
        self.assertLess(np.linalg.norm(-T - [0, 0, 1]), 1e-9)   # leaving direction is -T
        # outer lead with twist -90 / bend -90 also goes up
        s_after_bend = c.s_spiral1 + c.outer.offsets[4]
        T = c.frame(s_after_bend)[0]
        self.assertLess(np.linalg.norm(T - [0, 0, 1]), 1e-9)

    def test_clockwise_coil(self):
        c = make_coil(sense=-1)
        for s in np.linspace(0.5, c.tmax - 0.5, 30):
            self.assertLess(relerr(fd(c.r, s), c.v(s)), 1e-7)
            self.assertLess(relerr(fd(c.v, s), c.a(s)), 1e-6)
            T, n, b = c.frame(s)
            tau, kg, kn = c.curvatures(s)
            self.assertLess(np.linalg.norm(c.a(s) - (kn * n - kg * b)), 1e-9)

    def test_sampling_and_export(self):
        c = self.coil
        t = c.sample_spacing(5.0)
        self.assertTrue(np.all(np.diff(t) > 0))
        for s in c.breakpoints:
            self.assertTrue(np.any(np.abs(t - s) < 1e-12))
        os.makedirs(SCRATCH, exist_ok=True)
        pc.write_vtk_polyline(c, os.path.join(SCRATCH, "coil.vtk"))
        ntri = pc.write_stl_tape(c, 4.0, 0.3, os.path.join(SCRATCH, "coil.stl"))
        self.assertGreater(ntri, 0)
        geo = pc.Geometry(c, pc.CrossSection(2, 4.0, 0.3))
        geo.save(os.path.join(SCRATCH, "coil.geo"))
        with open(os.path.join(SCRATCH, "coil.geo")) as f:
            txt = f.read()
        self.assertIn("Volume(1)", txt)

    def test_rigid_transform(self):
        c = self.coil
        R = pc.rotation_z(0.7) @ pc.rotation_x(0.3)
        d = np.array([1.0, -2.0, 3.0])
        w = pc.RigidTransform(c, R, d)
        for s in (0.3 * c.tmax, 0.8 * c.tmax):
            self.assertLess(np.linalg.norm(w.r(s) - (R @ c.r(s) + d)), 1e-12)
            self.assertLess(relerr(fd(w.r, s), w.v(s)), 1e-7)
            self.assertEqual(w.curvatures(s), c.curvatures(s))


if __name__ == "__main__":
    unittest.main()
