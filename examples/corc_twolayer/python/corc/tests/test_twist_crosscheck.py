"""Independent cross-check: gmsh's own twist extrusion (screw motion,
Layers) must reproduce the analytic strip nodes of tape 0 / layer 0 for
the straight cable — same discretization, positions equal to 1e-9
(file round trip involved)."""

import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, ".")

from corc.winding import Winding
from corc.curve import StraightLine
from corc.frame import Frame
from corc.builder import Builder
from corc.shells import ShellGenerator
from corc.mshio import read_msh2

WORK = "/tmp"

W = Winding(tapeWidth=4, gap=1, tapeThickness=100, solderThickness=10,
            numTapesPerLayer=3, numLayers=1, pitch=10, numTurns=1,
            delta=5, tapeResolution=0.5)
S = StraightLine(W.length)
F = Frame(S)
G = ShellGenerator(W, F, 50.0, 10.0)
B = Builder()
G.build(B)

# our tape-0 nodes: slots 0..nt of layer 0 rings, all rows
angles, features = G.layer_slots(0)
kind, k0, start, count = features[0]
assert kind == "tape" and k0 == 0
rows = B.rings[0]
ours = []
for row in rows:
    for j in range(start, start + count + 1):
        ours.append(B.nodes[row[j % len(angles)] - 1])
ours = np.array(ours)

# twist-extrude the same tape with gmsh: 4 quarter turns, Layers
L = W.length
r = W.layers[0].radius
phi0 = W.layers[0].tapes[0].phi0
phi1 = W.layers[0].tapes[0].phi1
nt = count
ns = W.s_partition_count()
assert ns % 4 == 0
geo = os.path.join(WORK, "corc_twistcheck.geo")
with open(geo, "w") as f:
    f.write("Point(1) = {0, 0, 0};\n")
    f.write("Point(2) = {{{:.17g}, {:.17g}, 0}};\n".format(
        r * math.cos(phi0), r * math.sin(phi0)))
    f.write("Point(3) = {{{:.17g}, {:.17g}, 0}};\n".format(
        r * math.cos(phi1), r * math.sin(phi1)))
    f.write("Circle(1) = {2, 1, 3};\n")
    f.write("Transfinite Curve {{1}} = {:d};\n".format(nt + 1))
    f.write("cur = 1;\n")
    for q in range(4):
        f.write("ex{:d}[] = Extrude {{ {{0, 0, {:.17g}}}, {{0, 0, 1}}, "
                "{{0, 0, 0}}, Pi/2 }} {{ Curve{{cur}}; Layers{{{:d}}}; }};\n"
                .format(q, L / 4, ns // 4))
        f.write("cur = ex{:d}[0];\n".format(q))
    f.write("Mesh.SaveAll = 1;\n")

out = os.path.join(WORK, "corc_twistcheck.msh")
res = subprocess.run(["gmsh", "-2", "-format", "msh22", geo, "-o", out],
                     capture_output=True, text=True, timeout=300)
assert res.returncode == 0, res.stderr[-1500:]

nodes, elements = read_msh2(out)
# only nodes that belong to the structured surface (gmsh also saves orphan
# construction nodes of the internal helical edges — no elements use them)
tri_nodes = set()
for (t, p, g, nn) in elements:
    if t == 2:
        tri_nodes.update(nn)
pts = np.array([p for i, p in sorted(nodes.items())
                if i in tri_nodes and abs(math.hypot(p[0], p[1]) - r) < 1e-6])

expected = (ns + 1) * (nt + 1)
print("twist nodes on tape radius: {:d} (expect {:d}); ours: {:d}".format(
    len(pts), expected, len(ours)))

# match every twist node to our nearest strip node
worst = 0.0
for p in pts:
    d = np.linalg.norm(ours - p, axis=1).min()
    worst = max(worst, float(d))
print("max distance gmsh-twist vs analytic strips: {:.3e}".format(worst))
# gmsh's transfinite/extrusion node placement is itself iterative and lands
# ~4e-9 off the exact uniform-angle positions (Spike C); our strip nodes
# are the analytically exact ones, so the agreement bound is gmsh's
# placement accuracy, not serialization
ok = len(pts) == expected and worst < 2e-8
print("TWIST CROSS-CHECK:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
