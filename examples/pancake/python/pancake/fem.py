"""FEM mesh of a tape stack along a pancake base curve, BELFEM conventions.

Follows the CORC pipeline (old/python/corc): every surface that gmsh must
respect is generated as a structured triangulation in Python, the air volume
is filled by ``gmsh -3`` over the merged shell, and the result is written as
msh2 with named physical groups plus a matching ``topology`` section.

Structure (compare tapestack3d.geo):

    tape k  (k = 1..N)       thin shell, physical Tape_k
    slab k  (k = 1..N-1)     conductor volume between tape k and k+1,
                             physical Volume_{k+1}, meshed as structured tets
    air                      Volume_1, tets from gmsh
    DomainBoundary           domain faces without terminals or periodicity
    <terminal face>          domain face carrying the two stack end faces as holes
                             (DomainTop for the box, DomainOuter for the wedge)
    DomainLower/Upper        periodic pair of the wedge (master / slave)
    SolderInner_k/Outer_k    end faces of slab k (terminals)
    WallLeft_k/Right_k       tip walls of slab k (conductor-air interface)

The domain is supplied as a `Domain` object (BoxDomain, WedgeDomain).  It
meshes its faces with gmsh -2, reuses the stack nodes on the hole outlines
and, for periodic faces, provides the node pairs and the vertex triples
BELFEM derives the affine map from.
"""

import os
import math
import subprocess

import numpy as np
from scipy.spatial import cKDTree

import mesh as meshlib


class FemError(Exception):
    pass


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def read_msh2(path):
    with open(path) as f:
        lines = f.read().splitlines()
    i = lines.index("$Nodes")
    n = int(lines[i + 1])
    nodes = {}
    for line in lines[i + 2:i + 2 + n]:
        d = line.split()
        nodes[int(d[0])] = np.array([float(d[1]), float(d[2]), float(d[3])])
    i = lines.index("$Elements")
    n = int(lines[i + 1])
    elements = []
    for line in lines[i + 2:i + 2 + n]:
        d = [int(x) for x in line.split()]
        ntags = d[2]
        elements.append((d[1], d[3] if ntags > 0 else 0, d[4] if ntags > 1 else 0, d[3 + ntags:]))
    return nodes, elements


def run_gmsh(args, cwd, timeout, verbose=False):
    """Run the gmsh binary.  verbose=True streams gmsh's output to the terminal,
    otherwise it is captured and shown only when gmsh fails."""
    if verbose:
        print("gmsh " + " ".join(args), flush=True)
        r = subprocess.run(["gmsh"] + args, cwd=cwd, text=True, timeout=timeout)
        if r.returncode != 0:
            raise FemError("gmsh failed (see output above)")
        return r
    r = subprocess.run(["gmsh"] + args, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise FemError("gmsh failed:\n" + r.stderr[-3000:] + "\n" + r.stdout[-2000:])
    return r


def _quad_tris(p00, p10, p11, p01):
    """Two triangles of a grid quad, split along the (0,0)-(1,1) diagonal (Kuhn rule)."""
    return [(p00, p10, p11), (p00, p11, p01)]


_KUHN = [((1, 0, 0), (1, 1, 0)), ((1, 0, 0), (1, 0, 1)),
         ((0, 1, 0), (1, 1, 0)), ((0, 1, 0), (0, 1, 1)),
         ((0, 0, 1), (1, 0, 1)), ((0, 0, 1), (0, 1, 1))]


class GeoWriter:
    """Tiny helper to write planar .geo files with several disjoint surfaces."""

    def __init__(self):
        self.lines = []
        self.np = 0
        self.nc = 0
        self.nl = 0
        self.ns = 0

    def point(self, x, y, res):
        self.np += 1
        self.lines.append("Point(%d) = {%.17g, %.17g, 0, %.17g};" % (self.np, x, y, res))
        return self.np

    def line(self, a, b, n=None):
        self.nc += 1
        self.lines.append("Line(%d) = {%d, %d};" % (self.nc, a, b))
        if n:
            self.lines.append("Transfinite Curve {%d} = %d;" % (self.nc, n))
        return self.nc

    def circle(self, a, c, b, n=None):
        self.nc += 1
        self.lines.append("Circle(%d) = {%d, %d, %d};" % (self.nc, a, c, b))
        if n:
            self.lines.append("Transfinite Curve {%d} = %d;" % (self.nc, n))
        return self.nc

    def loop(self, curves):
        self.nl += 1
        self.lines.append("Curve Loop(%d) = {%s};" % (self.nl, ", ".join(str(c) for c in curves)))
        return self.nl

    def plane(self, loops):
        self.ns += 1
        self.lines.append("Plane Surface(%d) = {%s};" % (self.ns, ", ".join(str(l) for l in loops)))
        return self.ns

    def rectangle(self, x0, x1, y0, y1, res, nx=None, ny=None):
        p = [self.point(x0, y0, res), self.point(x1, y0, res), self.point(x1, y1, res), self.point(x0, y1, res)]
        c = [self.line(p[0], p[1], nx), self.line(p[1], p[2], ny), self.line(p[2], p[3], nx), self.line(p[3], p[0], ny)]
        return self.loop(c)

    def save(self, path, size_max):
        with open(path, "w") as f:
            f.write("\n".join(self.lines) + "\nMesh.MeshSizeMax = %.17g;\nMesh.SaveAll = 1;\n" % size_max)


# ----------------------------------------------------------------------
# domains
# ----------------------------------------------------------------------
class Domain:
    """Interface: check_ends(fm), build(fm).  build adds surfaces with roles
    'domain', 'terminal', 'periodic_master', 'periodic_slave' to fm and fills
    fm.periodic_pairs [(slave, master)] and fm.vertex_nodes."""

    terminal_name = "DomainTerminal"

    def check_ends(self, fm):
        raise NotImplementedError()

    def build(self, fm):
        raise NotImplementedError()

    # ---- shared machinery ------------------------------------------------
    @staticmethod
    def _hole_uv(fm, i, to_uv):
        """Corner points and outline nodes of the stack end face at station i in face coordinates."""
        g = fm.grid
        cids = [g[0, i, 0], g[0, i, -1], g[-1, i, -1], g[-1, i, 0]]
        outline = fm.hole_outline(i)
        return cids, [to_uv(fm.nodes[c - 1]) for c in cids], outline

    @staticmethod
    def _add_hole(G, uv, J, N, res):
        p = [G.point(u, v, res) for (u, v) in uv]
        c = [G.line(p[0], p[1], J), G.line(p[1], p[2], N), G.line(p[2], p[3], J), G.line(p[3], p[0], N)]
        return G.loop(c)

    @staticmethod
    def _merge_faces(fm, faces, hole_outline, hole_outline_uv, hole_tol):
        """faces: list of (role, meta, points (n,3), tris (m,3) local indices, uv (n,2) or None).
        Nodes of a face with uv are matched against the hole outline in the face's
        parameter plane and take the ids of the stack nodes; all other nodes are
        deduplicated across faces by their 3D position.  Returns per-face id arrays."""
        otree = cKDTree(np.asarray(hole_outline_uv)) if hole_outline else None
        prepared = []
        for role, meta, pts, tris, uv in faces:
            ids = np.zeros(len(pts), dtype=int)
            if uv is not None and otree is not None:
                claimed = set()
                for k in range(len(pts)):
                    dist, m = otree.query(uv[k])
                    if dist < hole_tol:
                        ids[k] = hole_outline[m]
                        claimed.add(m)
                if len(claimed) != len(hole_outline):
                    raise FemError("terminal face: only %d of %d hole outline nodes matched"
                                   % (len(claimed), len(hole_outline)))
            prepared.append((role, meta, pts, tris, ids))
        pool = np.array([pts[k] for _, _, pts, _, ids in prepared for k in range(len(pts)) if ids[k] == 0])
        parent = list(range(len(pool)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for a, b in cKDTree(pool).query_pairs(1e-6):
            parent[find(a)] = find(b)
        rep_id = {}
        pool_ids = []
        for k in range(len(pool)):
            r = find(k)
            if r not in rep_id:
                rep_id[r] = fm._add_node(pool[r])
            pool_ids.append(rep_id[r])
        c = 0
        out = []
        for role, meta, pts, tris, ids in prepared:
            for k in range(len(pts)):
                if ids[k] == 0:
                    ids[k] = pool_ids[c]
                    c += 1
            tag = fm._new_surface(role, **meta)
            for (a, b, d) in tris:
                fm.tris[tag].append((int(ids[a]), int(ids[b]), int(ids[d])))
            out.append(ids)
        return out


def _classify_2d(gnodes, gelems):
    """Split a planar gmsh mesh into per-surface (points, local tris)."""
    ids = sorted(gnodes)
    index = {g: k for k, g in enumerate(ids)}
    pts = np.array([gnodes[g] for g in ids])
    per = {}
    for (etype, phys, geo, nn) in gelems:
        if etype == 2:
            per.setdefault(geo, []).append([index[n] for n in nn])
    return pts, per


class BoxDomain(Domain):
    """Axis-aligned box; both leads end vertically on the top face z = z_top."""

    terminal_name = "DomainTop"

    def __init__(self, margin=20.0, margin_bottom=15.0):
        self.margin = float(margin)
        self.margin_bottom = float(margin_bottom)

    def check_ends(self, fm):
        C = fm.curve
        zs = []
        for s, name in ((C.tmin, "inner"), (C.tmax, "outer")):
            T = C.frame(s)[0]
            if abs(abs(T[2]) - 1.0) > 1e-9:
                raise FemError("%s lead does not end vertically (T = %s)" % (name, T))
            zs.append(C.r(s)[2])
        if abs(zs[0] - zs[1]) > 1e-9:
            raise FemError("leads end at different heights z = %g and %g" % tuple(zs))
        self.z_top = float(zs[0])
        fm.snap_end_faces(lambda p: np.array([p[0], p[1], self.z_top]))

    def build(self, fm):
        pts = np.array(fm.nodes)
        lo = pts.min(axis=0) - self.margin
        hi = pts.max(axis=0) + self.margin
        z0 = pts[:, 2].min() - self.margin_bottom
        z1 = self.z_top
        h = fm.domain_resolution
        G = GeoWriter()
        # five plain faces in their own planes, then the top face with holes
        faces_meta = []
        # bottom (x, y), sides: y=lo (x, z), x=hi (y, z), y=hi (x, z), x=lo (y, z)
        specs = [("bottom", (lo[0], hi[0], lo[1], hi[1]), lambda u, v: (u, v, z0)),
                 ("ymin", (lo[0], hi[0], z0, z1), lambda u, v: (u, lo[1], v)),
                 ("xmax", (lo[1], hi[1], z0, z1), lambda u, v: (hi[0], u, v)),
                 ("ymax", (lo[0], hi[0], z0, z1), lambda u, v: (u, hi[1], v)),
                 ("xmin", (lo[1], hi[1], z0, z1), lambda u, v: (lo[0], u, v))]
        off = 0.0
        for name, (u0, u1, v0, v1), fmap in specs:
            nu = int(math.ceil((u1 - u0) / h)) + 1
            nv = int(math.ceil((v1 - v0) / h)) + 1
            loop = G.rectangle(u0 + off, u1 + off, v0, v1, h, nu, nv)
            faces_meta.append((G.plane([loop]), off, fmap, False))
            off += (u1 - u0) + 10.0 * h
        nu = int(math.ceil((hi[0] - lo[0]) / h)) + 1
        nv = int(math.ceil((hi[1] - lo[1]) / h)) + 1
        loop = G.rectangle(lo[0] + off, hi[0] + off, lo[1], hi[1], h, nu, nv)
        holes = []
        outline = []
        res = min(fm.ds, fm.w / (fm.J - 1))
        for i in (0, fm.M - 1):
            cids, uv, ol = self._hole_uv(fm, i, lambda p: (p[0] + off, p[1]))
            holes.append(self._add_hole(G, uv, fm.J, fm.N, res))
            outline += ol
        outline_uv = [(fm.nodes[n - 1][0], fm.nodes[n - 1][1]) for n in outline]
        faces_meta.append((G.plane([loop] + holes), off, lambda u, v: (u, v, z1), True))
        geo = os.path.join(fm.workdir, fm.name + "_domain.geo")
        out = os.path.join(fm.workdir, fm.name + "_domain.msh")
        G.save(geo, h)
        run_gmsh(["-2", "-format", "msh22", geo, "-o", out], fm.workdir, 600, fm.verbose)
        pts2, per = _classify_2d(*read_msh2(out))
        faces = []
        for sid, off, fmap, is_term in faces_meta:
            tris = per.get(sid, [])
            if not tris:
                raise FemError("domain face %d has no triangles" % sid)
            used = sorted({k for t in tris for k in t})
            loc = {k: i for i, k in enumerate(used)}
            P = np.array([fmap(pts2[k][0] - off, pts2[k][1]) for k in used])
            T = [[loc[a] for a in t] for t in tris]
            uv = np.array([[pts2[k][0] - off, pts2[k][1]] for k in used]) if is_term else None
            faces.append(("terminal" if is_term else "domain", {}, P, T, uv))
        self._merge_faces(fm, faces, outline, outline_uv, 0.2 * min(fm.d, fm.w / (fm.J - 1)))
        fm.periodic_pairs = []
        fm.vertex_nodes = []


class WedgeDomain(Domain):
    """One 1/n sector of a machine with the axis parallel to y at x = -R0, z = 0.

    Bounded by the planes at +-half_angle about the axis (periodic pair), the
    cylinders rho = R_in and rho = R_out about the axis, and the planes
    y = +-Y.  Both leads must end on the outer cylinder pointing radially
    outward; the stack end faces become holes in the outer cylinder.
    """

    terminal_name = "DomainOuter"

    def __init__(self, R0, R_in, R_out, Y, half_angle_deg=15.0):
        self.R0 = float(R0)
        self.R_in = float(R_in)
        self.R_out = float(R_out)
        self.Y = float(Y)
        self.alpha = math.radians(half_angle_deg)

    def radial(self, p):
        """(rho, phi) of a point about the axis."""
        return math.hypot(p[0] + self.R0, p[2]), math.atan2(p[2], p[0] + self.R0)

    def check_ends(self, fm):
        C = fm.curve
        for s, name, sign in ((C.tmin, "inner", -1.0), (C.tmax, "outer", 1.0)):
            r = C.r(s)
            T = C.frame(s)[0] * sign
            rho, phi = self.radial(r)
            e_rho = np.array([math.cos(phi), 0.0, math.sin(phi)])
            if abs(rho - self.R_out) > 1e-6:
                raise FemError("%s lead ends at rho = %.6f, not on R_out = %g" % (name, rho, self.R_out))
            if np.dot(T, e_rho) < math.cos(math.radians(10.0)):
                raise FemError("%s lead does not leave radially (T = %s)" % (name, T))
            if abs(phi) > self.alpha or abs(r[1]) > self.Y:
                raise FemError("%s lead end lies outside the wedge" % name)

    def build(self, fm):
        R0, Ri, Ro, Y, al = self.R0, self.R_in, self.R_out, self.Y, self.alpha
        h = fm.domain_resolution
        c, s = math.cos(al), math.sin(al)
        n_y = int(math.ceil(2 * Y / h)) + 1
        n_rho = int(math.ceil((Ro - Ri) / h)) + 1
        n_po = int(math.ceil(2 * al * Ro / h)) + 1
        n_pi = max(3, int(math.ceil(2 * al * Ri / h)) + 1)

        G = GeoWriter()
        meta = []
        off = 0.0

        # 1. lower periodic plane in (rho, y); the upper plane is its rotated image
        loop = G.rectangle(Ri + off, Ro + off, -Y, Y, h, n_rho, n_y)
        meta.append(("plane", G.plane([loop]), off))
        off += (Ro - Ri) + 10 * h

        # 2. outer cylinder in (Ro*phi, y) with the two terminal holes
        u0, u1 = -Ro * al, Ro * al
        loop = G.rectangle(u0 + off, u1 + off, -Y, Y, h, n_po, n_y)
        holes = []
        outline = []
        res = min(fm.ds, fm.w / (fm.J - 1))
        to_uv = lambda p: (Ro * self.radial(p)[1], p[1])
        for i in (0, fm.M - 1):
            cids, uv, ol = self._hole_uv(fm, i, lambda p: (to_uv(p)[0] + off, to_uv(p)[1]))
            holes.append(self._add_hole(G, uv, fm.J, fm.N, res))
            outline += ol
        outline_uv = [to_uv(fm.nodes[n - 1]) for n in outline]
        meta.append(("outer", G.plane([loop] + holes), off))
        off += (u1 - u0) + 10 * h

        # 3. inner cylinder in (Ri*phi, y)
        loop = G.rectangle(-Ri * al + off, Ri * al + off, -Y, Y, h, n_pi, n_y)
        meta.append(("inner", G.plane([loop]), off))
        off += 2 * Ri * al + 10 * h

        # 4./5. end faces y = -Y and y = +Y: annular sectors drawn in the (X, Z) plane
        for name in ("ymin", "ymax"):
            cen = G.point(off, 0.0, h)
            p = [G.point(Ri * c + off, -Ri * s, h), G.point(Ro * c + off, -Ro * s, h),
                 G.point(Ro * c + off, Ro * s, h), G.point(Ri * c + off, Ri * s, h)]
            cv = [G.line(p[0], p[1], n_rho), G.circle(p[1], cen, p[2], n_po),
                  G.line(p[2], p[3], n_rho), G.circle(p[3], cen, p[0], n_pi)]
            sid = G.plane([G.loop(cv)])
            meta.append((name, sid, off))
            off += Ro + 10 * h

        geo = os.path.join(fm.workdir, fm.name + "_domain.geo")
        out = os.path.join(fm.workdir, fm.name + "_domain.msh")
        G.save(geo, h)
        run_gmsh(["-2", "-format", "msh22", geo, "-o", out], fm.workdir, 600, fm.verbose)
        pts2, per = _classify_2d(*read_msh2(out))

        def face(sid, off):
            tris = per.get(sid, [])
            if not tris:
                raise FemError("domain face %d has no triangles" % sid)
            used = sorted({k for t in tris for k in t})
            loc = {k: i for i, k in enumerate(used)}
            uv = np.array([[pts2[k][0] - off, pts2[k][1]] for k in used])
            return uv, [[loc[a] for a in t] for t in tris]

        faces = []
        plane_uv = None
        for name, sid, off in meta:
            uv, T = face(sid, off)
            if name == "plane":
                plane_uv = uv
                lower = np.array([[-R0 + r * c, y, -r * s] for r, y in uv])
                upper = np.array([[-R0 + r * c, y, r * s] for r, y in uv])
                faces.append(("periodic_master", {}, lower, T, None))
                faces.append(("periodic_slave", {}, upper, T, None))
            elif name == "outer":
                P = np.array([[-R0 + Ro * math.cos(u / Ro), y, Ro * math.sin(u / Ro)] for u, y in uv])
                faces.append(("terminal", {}, P, T, uv))
            elif name == "inner":
                P = np.array([[-R0 + Ri * math.cos(u / Ri), y, Ri * math.sin(u / Ri)] for u, y in uv])
                faces.append(("domain", {"face": "inner"}, P, T, None))
            else:
                yy = -Y if name == "ymin" else Y
                P = np.array([[X - R0, yy, Z] for X, Z in uv])
                faces.append(("domain", {"face": name}, P, T, None))
        ids = self._merge_faces(fm, faces, outline, outline_uv, 0.2 * min(fm.d, fm.w / (fm.J - 1)))
        lower_ids, upper_ids = ids[0], ids[1]
        fm.periodic_pairs = [(int(u), int(l)) for u, l in zip(upper_ids, lower_ids)]
        # vertex triples: three corners of the lower plane and their images
        corners = [(Ri, -Y), (Ro, -Y), (Ro, Y)]
        vtx_lower, vtx_upper = [], []
        for cr in corners:
            k = int(np.argmin(np.linalg.norm(plane_uv - np.array(cr), axis=1)))
            vtx_lower.append(int(lower_ids[k]))
            vtx_upper.append(int(upper_ids[k]))
        fm.vertex_nodes = vtx_lower + vtx_upper


# ----------------------------------------------------------------------
class PancakeMesh:

    def __init__(self, curve, domain: Domain, numtapes: int, tapewidth: float, tapedistance: float,
                 ds: float = 1.0, n_across: int = 5, domain_resolution: float = 8.0,
                 workdir: str = "/tmp", name: str = "pancake", verbose: bool = True):
        self.curve = curve
        self.verbose = bool(verbose)
        self.domain = domain
        self.N = int(numtapes)
        self.w = float(tapewidth)
        self.d = float(tapedistance)
        self.ds = float(ds)
        self.J = int(n_across)
        self.domain_resolution = float(domain_resolution)
        self.workdir = workdir
        self.name = name
        if self.N < 2 or self.J < 2:
            raise FemError("need at least two tapes and two nodes across the width")
        self.nodes = []
        self.tris = {}
        self.surfaces = {}
        self.slab_tets = []
        self.air_tets = []
        self.periodic_pairs = []
        self.vertex_nodes = []
        self._mesh = None
        self.physicals = {}

    # ------------------------------------------------------------------
    def _add_node(self, p):
        self.nodes.append(np.asarray(p, float))
        return len(self.nodes)

    def _new_surface(self, role, **meta):
        tag = len(self.surfaces) + 1
        self.surfaces[tag] = dict(role=role, **meta)
        self.tris[tag] = []
        return tag

    def _tags(self, role, **match):
        return [t for t, m in self.surfaces.items()
                if m["role"] == role and all(m.get(k) == v for k, v in match.items())]

    def build(self):
        self._build_stack()
        self.domain.check_ends(self)
        self.domain.build(self)
        self._mesh_air()
        self._assemble()
        return self

    # ------------------------------------------------------------ stack
    def _build_stack(self):
        C = self.curve
        t = C.sample_spacing(self.ds)
        M = len(t)
        self.t = t
        x = np.linspace(-0.5 * self.w, 0.5 * self.w, self.J)
        h = self.d * (self.N - 1)
        y = np.linspace(-0.5 * h, 0.5 * h, self.N)
        self.stack_thickness = h
        grid = np.zeros((self.N, M, self.J), dtype=int)
        self._station_R = []
        self._station_r = []
        for i, ti in enumerate(t):
            r = C.r(ti)
            R = C.transform(ti)
            self._station_r.append(r)
            self._station_R.append(R)
            for k in range(self.N):
                for j in range(self.J):
                    grid[k, i, j] = self._add_node(r + R @ np.array([x[j], y[k], 0.0]))
        self.grid = grid
        self.M = M
        self._station_tree = cKDTree(np.array(self._station_r))

        for k in range(self.N):
            tag = self._new_surface("tape", tape=k)
            for i in range(M - 1):
                for j in range(self.J - 1):
                    self.tris[tag] += _quad_tris(grid[k, i, j], grid[k, i + 1, j],
                                                 grid[k, i + 1, j + 1], grid[k, i, j + 1])
        for k in range(self.N - 1):
            left = self._new_surface("wall", slab=k, side="left")
            right = self._new_surface("wall", slab=k, side="right")
            for i in range(M - 1):
                self.tris[left] += _quad_tris(grid[k, i, 0], grid[k, i + 1, 0],
                                              grid[k + 1, i + 1, 0], grid[k + 1, i, 0])
                self.tris[right] += _quad_tris(grid[k, i, -1], grid[k, i + 1, -1],
                                               grid[k + 1, i + 1, -1], grid[k + 1, i, -1])
            inner = self._new_surface("end", slab=k, end="inner")
            outer = self._new_surface("end", slab=k, end="outer")
            for j in range(self.J - 1):
                self.tris[inner] += _quad_tris(grid[k, 0, j], grid[k, 0, j + 1],
                                               grid[k + 1, 0, j + 1], grid[k + 1, 0, j])
                self.tris[outer] += _quad_tris(grid[k, -1, j], grid[k, -1, j + 1],
                                               grid[k + 1, -1, j + 1], grid[k + 1, -1, j])
        for k in range(self.N - 1):
            for i in range(M - 1):
                for j in range(self.J - 1):
                    def v(a, b, c):
                        return grid[k + c, i + a, j + b]
                    v0, v7 = v(0, 0, 0), v(1, 1, 1)
                    for e1, e2 in _KUHN:
                        tet = [v0, v(*e1), v(*e2), v7]
                        if self._volume(tet) < 0.0:
                            tet[1], tet[2] = tet[2], tet[1]
                        self.slab_tets.append((k, tet))
        self._fix_orientation()

    def snap_end_faces(self, project):
        """Move the nodes of both end faces with a projection function (e.g. onto a plane)."""
        for i in (0, self.M - 1):
            for n in self.grid[:, i, :].ravel():
                q = project(self.nodes[n - 1])
                if np.linalg.norm(q - self.nodes[n - 1]) > 1e-6:
                    raise FemError("end face at station %d is far from the terminal surface" % i)
                self.nodes[n - 1] = q

    def hole_outline(self, i):
        """Node ids around the stack end face at station i (closed loop, first node not repeated)."""
        g = self.grid
        return (list(g[0, i, :]) + list(g[1:, i, -1]) + list(g[-1, i, ::-1][1:]) + list(g[-2::-1, i, 0][:-1]))

    def _volume(self, tet):
        a, b, c, d = (self.nodes[n - 1] for n in tet)
        return float(np.dot(np.cross(b - a, c - a), d - a)) / 6.0

    def _frame_dir(self, p, col):
        i = int(self._station_tree.query(p)[1])
        return self._station_R[i][:, col]

    def _fix_orientation(self):
        for tag, meta in self.surfaces.items():
            role = meta["role"]
            fixed = []
            for (a, b, c) in self.tris[tag]:
                pa, pb, pc = self.nodes[a - 1], self.nodes[b - 1], self.nodes[c - 1]
                nrm = np.cross(pb - pa, pc - pa)
                cen = (pa + pb + pc) / 3.0
                if role == "tape":
                    ref = self._frame_dir(cen, 1)
                elif role == "wall":
                    ref = -self._frame_dir(cen, 0) if meta["side"] == "left" else self._frame_dir(cen, 0)
                else:
                    ref = -self._frame_dir(cen, 2) if meta["end"] == "inner" else self._frame_dir(cen, 2)
                if float(np.dot(nrm, ref)) < 0.0:
                    a, b = b, a
                fixed.append((a, b, c))
            self.tris[tag] = fixed

    # ------------------------------------------------------------ air
    def _mesh_air(self):
        loop = (self._tags("domain") + self._tags("terminal") + self._tags("periodic_master")
                + self._tags("periodic_slave") + self._tags("tape", tape=0)
                + self._tags("tape", tape=self.N - 1) + self._tags("wall"))
        shell = os.path.join(self.workdir, self.name + "_shell.msh")
        geo = os.path.join(self.workdir, self.name + "_air.geo")
        out = os.path.join(self.workdir, self.name + "_air.msh")
        used = sorted({n for tag in loop for tri in self.tris[tag] for n in tri})
        with open(shell, "w") as f:
            f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
            f.write("$Nodes\n%d\n" % len(used))
            for n in used:
                p = self.nodes[n - 1]
                f.write("%d %.17g %.17g %.17g\n" % (n, p[0], p[1], p[2]))
            f.write("$EndNodes\n$Elements\n%d\n" % sum(len(self.tris[tag]) for tag in loop))
            e = 0
            for tag in loop:
                for (a, b, c) in self.tris[tag]:
                    e += 1
                    f.write("%d 2 2 %d %d %d %d %d\n" % (e, tag, tag, a, b, c))
            f.write("$EndElements\n")
        with open(geo, "w") as f:
            f.write('Merge "%s";\n' % os.path.basename(shell))
            f.write("Surface Loop(1) = {%s};\n" % ", ".join(str(t) for t in loop))
            f.write("Volume(1) = {1};\n")
            f.write("Mesh.MeshSizeExtendFromBoundary = 1;\n")
            f.write("Mesh.MeshSizeMax = %.17g;\n" % self.domain_resolution)
            f.write("Mesh.OptimizeNetgen = 1;\n")
            f.write("Mesh.SaveAll = 1;\n")
        run_gmsh(["-3", "-format", "msh22", os.path.basename(geo), "-o", os.path.basename(out)],
                 self.workdir, 3600, self.verbose)

        gnodes, gelems = read_msh2(out)
        gids = sorted(gnodes)
        tree = cKDTree(np.array([gnodes[g] for g in gids]))
        gid_to_id = {}
        for n in used:
            dist, k = tree.query(self.nodes[n - 1])
            if dist > 1e-9:
                raise FemError("shell node %d was moved by the volume run (%.3e)" % (n, dist))
            gid_to_id[gids[k]] = n
        if len(set(gid_to_id.values())) != len(used):
            raise FemError("shell node matching is not injective")
        for g in gids:
            if g not in gid_to_id:
                gid_to_id[g] = self._add_node(gnodes[g])
        want = {tag: {frozenset(tri) for tri in self.tris[tag]} for tag in loop}
        got = {}
        for (etype, phys, geo_tag, nn) in gelems:
            if etype == 2:
                got.setdefault(geo_tag, set()).add(frozenset(gid_to_id[n] for n in nn))
        bad = [tag for tag in loop if got.get(tag) != want[tag]]
        if bad:
            raise FemError("shell connectivity changed on surfaces %s" % bad)
        self.air_tets = []
        for (etype, phys, geo_tag, nn) in gelems:
            if etype == 4:
                tet = [gid_to_id[n] for n in nn]
                if self._volume(tet) < 0.0:
                    tet[1], tet[2] = tet[2], tet[1]
                self.air_tets.append(tet)
        if not self.air_tets:
            raise FemError("no air tets")

    # ------------------------------------------------------------ final
    def _assemble(self):
        M = meshlib.Mesh()
        node_obj = [meshlib.Node(p[0], p[1], p[2]) for p in self.nodes]

        def make(cls, ids, tag):
            E = cls()
            E.nodes = [node_obj[i - 1] for i in ids]
            E.geometryTag = tag
            E.physicalTag = tag
            return E

        elements = []
        # vertex elements "k 15 2 0 k k" on the periodic triples (nodes 1..6 after reordering)
        for k, n in enumerate(self.vertex_nodes):
            V = meshlib.Vertex()
            V.nodes = [node_obj[n - 1]]
            V.geometryTag = k + 1
            V.physicalTag = 0
            elements.append(V)

        M.add_physical(3, 1, "Volume_1")
        for tet in self.air_tets:
            elements.append(make(meshlib.Tet4, tet, 1))
        for k in range(self.N - 1):
            M.add_physical(3, k + 2, "Volume_%d" % (k + 2))
        for k, tet in self.slab_tets:
            elements.append(make(meshlib.Tet4, tet, k + 2))

        tag = self.N
        phys = {}

        def emit(tags, name):
            nonlocal tag
            tag += 1
            M.add_physical(2, tag, name)
            for t in tags:
                for tri in self.tris[t]:
                    elements.append(make(meshlib.Tri3, tri, tag))
            return tag

        phys["tapes"] = [emit(self._tags("tape", tape=k), "Tape_%d" % (k + 1)) for k in range(self.N)]
        phys["boundary"] = emit(self._tags("domain"), "DomainBoundary")
        phys["terminal"] = emit(self._tags("terminal"), self.domain.terminal_name)
        if self.periodic_pairs:
            phys["master"] = emit(self._tags("periodic_master"), "DomainLower")
            phys["slave"] = emit(self._tags("periodic_slave"), "DomainUpper")
        phys["solder_inner"] = [emit(self._tags("end", slab=k, end="inner"), "SolderInner_%d" % (k + 1))
                                for k in range(self.N - 1)]
        phys["solder_outer"] = [emit(self._tags("end", slab=k, end="outer"), "SolderOuter_%d" % (k + 1))
                                for k in range(self.N - 1)]
        phys["wall_left"] = [emit(self._tags("wall", slab=k, side="left"), "WallLeft_%d" % (k + 1))
                             for k in range(self.N - 1)]
        phys["wall_right"] = [emit(self._tags("wall", slab=k, side="right"), "WallRight_%d" % (k + 1))
                              for k in range(self.N - 1)]
        M.elements = elements

        # node order: vertex nodes first, then everything else
        first = [node_obj[n - 1] for n in self.vertex_nodes]
        fs = set(self.vertex_nodes)
        M.nodes = first + [node_obj[i] for i in range(len(node_obj)) if (i + 1) not in fs]
        for i, Nn in enumerate(M.nodes):
            Nn.id = i + 1
        if self.periodic_pairs:
            pairs = [(node_obj[s - 1].id, node_obj[m - 1].id) for s, m in self.periodic_pairs]
            M.periodic = [meshlib.PeriodicLink(2, phys["slave"], phys["master"], None, pairs)]
        self._mesh = M
        self.physicals = phys
        self._validate()

    def _validate(self):
        worst = math.inf
        for group, tets in (("slab", [t for _, t in self.slab_tets]), ("air", self.air_tets)):
            gam = []
            for tet in tets:
                a, b, c, d = (self.nodes[n - 1] for n in tet)
                v = float(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
                if v <= 0.0:
                    raise FemError("inverted %s tet" % group)
                area = sum(0.5 * np.linalg.norm(np.cross(f[1] - f[0], f[2] - f[0]))
                           for f in ((a, b, c), (a, b, d), (a, c, d), (b, c, d)))
                rin = 3.0 * v / area
                aa = np.linalg.norm(b - a) * np.linalg.norm(c - d)
                bb = np.linalg.norm(c - a) * np.linalg.norm(b - d)
                cc = np.linalg.norm(d - a) * np.linalg.norm(b - c)
                s = 0.5 * (aa + bb + cc)
                rc = math.sqrt(max(s * (s - aa) * (s - bb) * (s - cc), 0.0)) / (6.0 * v)
                gam.append(3.0 * rin / rc if rc > 0 else 0.0)
            gam = np.array(gam)
            print("%s tets: %d, gamma min %.3f / p5 %.3f / median %.3f"
                  % (group, len(gam), gam.min(), np.percentile(gam, 5), np.median(gam)))
            worst = min(worst, gam.min())
        if worst < 1e-3:
            raise FemError("degenerate tets")
        if self.periodic_pairs:
            # the pairs must be images under the rotation defined by the vertex triple
            P = [self.nodes[n - 1] for n in self.vertex_nodes]
            U = np.column_stack([P[1] - P[0], P[2] - P[0], np.cross(P[1] - P[0], P[2] - P[0])])
            V = np.column_stack([P[4] - P[3], P[5] - P[3], np.cross(P[4] - P[3], P[5] - P[3])])
            R = V @ np.linalg.inv(U)
            if np.abs(R.T @ R - np.eye(3)).max() > 1e-6:
                raise FemError("vertex-defined periodic map is not a rotation")
            worst = max(np.linalg.norm(R @ (self.nodes[m - 1] - P[0]) + P[3] - self.nodes[s - 1])
                        for s, m in self.periodic_pairs)
            if worst > 1e-6:
                raise FemError("periodic pair residual %.3e" % worst)
            print("periodic: %d node pairs, max residual %.2e" % (len(self.periodic_pairs), worst))
        print("mesh: %d nodes, %d elements" % (len(self.nodes), len(self._mesh.elements)))

    # ------------------------------------------------------------ output
    def save(self, path):
        if self._mesh is None:
            raise RuntimeError("call build() first")
        self._mesh.save(path)

    def print(self):
        p = self.physicals
        print("Volumes: 1 (air), 2..%d (solder slabs)" % self.N)
        print("Tapes:", p["tapes"])
        print("DomainBoundary:", p["boundary"], " %s:" % self.domain.terminal_name, p["terminal"])
        if "master" in p:
            print("DomainLower (master):", p["master"], " DomainUpper (slave):", p["slave"])
        print("SolderInner:", p["solder_inner"], " SolderOuter:", p["solder_outer"])
        print("WallLeft:", p["wall_left"], " WallRight:", p["wall_right"])

    def topology_section(self, amplitude="60 A", period="10 s", offset="-0.01 s"):
        p = self.physicals
        L = []
        a = L.append
        a("topology")
        a("{")
        a("\tthinshell : tape")
        a("\t{")
        # BELFEM puts the master of an air-adjacent tape on its solder side and leaves the
        # inner tapes with the master towards the smaller element ids, i.e. the slab above
        # (slabs are numbered from tape 1 downwards, after the air).  Tape N then agrees
        # with the inner tapes and tape 1 comes out mirrored, so its sideset is signed:
        # a negative id flips the facet orientation and with it the layer stacking.
        a("\t\t// tape 1 is signed: BELFEM puts the master of an air-adjacent tape on its solder side,")
        a("\t\t// which mirrors tape 1 against the other %d (see the generator)" % (self.N - 1))
        if self.N > 1:
            a("\t\tsidesets : -%d, %d:%d ;" % (p["tapes"][0], p["tapes"][1], p["tapes"][-1]))
        else:
            a("\t\tsidesets : -%d ;" % p["tapes"][0])
        a("\t}")
        a("")
        a("\tair")
        a("\t{")
        a("\t\tblock : 1 ;")
        a("\t}")
        a("")
        a("\tconductor : solder")
        a("\t{")
        a("\t\tblocks : 2:%d ;" % self.N if self.N > 2 else "\t\tblocks : 2 ;")
        a("\t\tmaterial : solder ;")
        a("\t}")
        a("")
        a("\t// curves: tape edge on the end face of the solder slab below it. Tape 1 borders the air,")
        a("\t// but both of its ends lie on the same domain face (%s), so it uses the slab" % self.domain.terminal_name)
        a("\t// shared with tape 2: an intersection with two disconnected chains is not a curve.")
        a("\tcurves")
        a("\t{")
        cid = 0
        for k in range(self.N):
            cid += 1
            face = p["solder_inner"][max(k - 1, 0)]
            a("\t\t%2d : %2d @ %2d ;   // tape %d, inner end" % (cid, face, p["tapes"][k], k + 1))
        for k in range(self.N):
            cid += 1
            face = p["solder_outer"][max(k - 1, 0)]
            a("\t\t%2d : %2d @ %2d ;   // tape %d, outer end" % (cid, face, p["tapes"][k], k + 1))
        a("\t}")
        if self.periodic_pairs:
            a("")
            a("\t// periodicity: three corners of the lower plane (vertex elements 1, 2, 3)")
            a("\t// map to their images on the upper plane (4, 5, 6)")
            a("\tperiodic")
            a("\t{")
            a("\t\tsource : 1, 2, 3 ;")
            a("\t\ttarget : 4, 5, 6 ;")
            a("\t}")
        a("}")
        a("")
        a("boundary conditions")
        a("{")
        a("\tcurrent")
        a("\t{")
        a("\t\tinput terminals  : [%d:%d] ;   // solder end faces at the inner connector"
          % (p["solder_inner"][0], p["solder_inner"][-1]))
        a("\t\toutput terminals : [%d:%d] ;   // solder end faces at the outer connector"
          % (p["solder_outer"][0], p["solder_outer"][-1]))
        a("\t\ttype : ramp ;")
        a("\t\tamplitude : %s ;" % amplitude)
        a("\t\tperiod : %s ;" % period)
        a("\t\toffset : %s ;" % offset)
        a("\t}")
        a("}")
        return "\n".join(L) + "\n"

    def write_topology(self, path, **kw):
        with open(path, "w") as f:
            f.write(self.topology_section(**kw))


# ----------------------------------------------------------------------
# lead planning helpers
# ----------------------------------------------------------------------
def _probe(spiral, inner, outer, release):
    from .pancake import PancakeCurve
    from .lead import Straight
    return PancakeCurve(spiral, list(inner) + [Straight(1.0)], list(outer) + [Straight(1.0)], release=release)


def leads_to_plane(spiral, inner, outer, z_top: float, release: float = 0.0):
    """PancakeCurve whose leads end with a vertical straight at z = z_top (box domain)."""
    from .pancake import PancakeCurve
    from .lead import Straight
    probe = _probe(spiral, inner, outer, release)
    out = []
    for s_end, sign in ((probe.tmin, -1.0), (probe.tmax, +1.0)):
        T = probe.frame(s_end)[0] * sign
        if abs(T[2] - 1.0) > 1e-9:
            raise FemError("lead does not point upwards after its last segment (T = %s)" % T)
        z_end = probe.r(s_end)[2] - 1.0
        if z_top <= z_end:
            raise FemError("z_top = %g lies below the end of the lead (z = %g)" % (z_top, z_end))
        out.append(z_top - z_end)
    return PancakeCurve(spiral, list(inner) + [Straight(out[0])], list(outer) + [Straight(out[1])], release=release)


def leads_to_cylinder(spiral, inner, outer, wedge: WedgeDomain, release: float = 0.0):
    """PancakeCurve whose leads end with a straight on the outer cylinder of the wedge.

    The last user segment of each lead must leave (approximately) radially, i.e.
    in the +x direction away from the axis at x = -R0."""
    from .pancake import PancakeCurve
    from .lead import Straight
    probe = _probe(spiral, inner, outer, release)
    out = []
    for s_end, sign in ((probe.tmin, -1.0), (probe.tmax, +1.0)):
        T = probe.frame(s_end)[0] * sign
        r0 = probe.r(s_end) - 1.0 * T
        # solve |(r0 + L T) - axis|_{xz} = R_out for L > 0
        ax = np.array([r0[0] + wedge.R0, r0[2]])
        tx = np.array([T[0], T[2]])
        A = np.dot(tx, tx)
        B = 2.0 * np.dot(ax, tx)
        Cc = np.dot(ax, ax) - wedge.R_out ** 2
        disc = B * B - 4 * A * Cc
        if A < 1e-12 or disc < 0.0:
            raise FemError("lead direction %s never reaches the outer cylinder" % T)
        L = (-B + math.sqrt(disc)) / (2 * A)
        if L <= 0.0:
            raise FemError("lead already outside the outer cylinder")
        out.append(L)
    return PancakeCurve(spiral, list(inner) + [Straight(out[0])], list(outer) + [Straight(out[1])], release=release)
