"""Pancake coil in a 30 degree sector of a 12-coil machine, meshed for BELFEM.

The machine axis is parallel to y at x = -R0, z = 0.  The sector is bounded by
the planes at +-15 degrees about that axis (periodic pair), the cylinders
rho = R_in and rho = R_out and the planes y = +-Y.  Both leads leave radially
outward and end on the outer cylinder, where the stack end faces are the
current terminals.

    python3 main.py [workdir]
"""

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import pancake as pc

WORKDIR = sys.argv[1] if len(sys.argv) > 1 else "/tmp"

# --- tape stack
numTapes = 8
tapeWidth = 4.0
tapeDistance = 0.25
stack = tapeDistance * (numTapes - 1)

# --- winding: 1.75 turns so that the outer end sits at the bottom of the
#     coil (theta = 3 pi / 2) with the tangent pointing away from the axis
turnGap = 1.0
spiral = pc.LameSpiral(a=25.0, b=35.0, n=4, pitch=stack + turnGap, nturns=1.75)

# --- sector of the machine
wedge = pc.WedgeDomain(R0=60.0, R_in=10.0, R_out=120.0, Y=55.0, half_angle_deg=15.0)

# --- leads (leaving direction)
#     outer: straight on, radially outward
#     inner: rise hard-way out of the winding plane, then bend easy-way to run
#            radially outward above the winding to the outer cylinder
inner = [pc.Bend(-90.0, 5.0, ramp=2.0, hard=True), pc.Straight(2.0), pc.Bend(90.0, 5.0, ramp=2.0)]
outer = []
C = pc.leads_to_cylinder(spiral, inner, outer, wedge, release=3.0)
print("curve length %.1f mm (inner lead %.1f, winding %.1f, outer lead %.1f)"
      % (C.tmax, C.L_inner, C.L_spiral, C.L_outer))

# --- clearance and a picture
t = C.sample_spacing(1.0)
S = pc.TapeStack(C, numTapes, tapeWidth, tapeDistance)
d, i, j = S.clearance()
print("stack thickness %.2f mm, smallest clearance %.3f mm (s = %.1f %s / s = %.1f %s)"
      % (S.thickness, d, t[i], C.piece(t[i]), t[j], C.piece(t[j])))
if d < 0.5 * turnGap:
    raise SystemExit("stack too close to itself")

# distance of the stack envelope from the periodic planes
env = S.envelope().reshape(-1, 3)
rho = np.hypot(env[:, 0] + wedge.R0, env[:, 2])
phi = np.arctan2(env[:, 2], env[:, 0] + wedge.R0)
gap_planes = (rho * np.sin(wedge.alpha - np.abs(phi))).min()
print("smallest distance of the stack to the periodic planes: %.2f mm" % gap_planes)

ax = plt.figure(figsize=(11, 8)).add_subplot(111, projection="3d")
r = np.array([C.r(s) for s in t])
k1, k2 = np.searchsorted(t, C.s_spiral0), np.searchsorted(t, C.s_spiral1)
ax.plot(r[:k1 + 1, 0], r[:k1 + 1, 1], r[:k1 + 1, 2], "-r")
ax.plot(r[k1:k2 + 1, 0], r[k1:k2 + 1, 1], r[k1:k2 + 1, 2], "-b")
ax.plot(r[k2:, 0], r[k2:, 1], r[k2:, 2], "-r")
S.plot(ax, tapes=(0, numTapes - 1), sections_every=40)
# wedge outline
c, s_ = np.cos(wedge.alpha), np.sin(wedge.alpha)
for sgn in (-1, 1):
    for y in (-wedge.Y, wedge.Y):
        ax.plot([-wedge.R0 + wedge.R_in * c, -wedge.R0 + wedge.R_out * c], [y, y],
                [sgn * wedge.R_in * s_, sgn * wedge.R_out * s_], color="gray", lw=0.8)
    for rr in (wedge.R_in, wedge.R_out):
        ax.plot([-wedge.R0 + rr * c] * 2, [-wedge.Y, wedge.Y], [sgn * rr * s_] * 2, color="gray", lw=0.8)
ph = np.linspace(-wedge.alpha, wedge.alpha, 30)
for y in (-wedge.Y, wedge.Y):
    for rr in (wedge.R_in, wedge.R_out):
        ax.plot(-wedge.R0 + rr * np.cos(ph), y * np.ones_like(ph), rr * np.sin(ph), color="gray", lw=0.8)
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
ax.set_xlim(-60, 70); ax.set_ylim(-65, 65); ax.set_zlim(-40, 40); ax.set_box_aspect((130, 130, 80))
plt.savefig("pancake_wedge.png", dpi=120)

# --- FEM mesh
FM = pc.PancakeMesh(C, wedge, numTapes, tapeWidth, tapeDistance, ds=0.5, n_across=12,
                    domain_resolution=2.0, workdir=WORKDIR, name="pancake_wedge")
FM.build()
FM.print()
FM.save("pancake.msh")
FM.write_topology("pancake_topology.conf")
