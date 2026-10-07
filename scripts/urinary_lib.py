"""Helpers for build_urinary: pressure-grown hollow organs, neighbours that give way (make_room), shape-key stages
driven by one control, a slit cut into the skin. Arrays are numpy (n, 3) in one frame (atlas space unless said)."""
import bpy, bmesh, os, json, math
import numpy as np
from mathutils import Vector, Matrix, kdtree
from mathutils.bvhtree import BVHTree


def save(path, out, report):
    os.makedirs(out, exist_ok=True)
    json.dump(report, open(os.path.join(out, "report.json"), "w"), indent=1, default=float)
    bpy.ops.wm.save_as_mainfile(filepath=path, relative_remap=False)


def ss(a, b, x):
    t = np.clip((np.asarray(x, float) - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


# ---------------------------------------------------------------- mesh <-> arrays
def bm_in(ob, M=None, weld=False):
    """bmesh of ob in the frame M @ matrix_world (M=None: world)"""
    bm = bmesh.new(); bm.from_mesh(ob.data)
    bm.transform((M @ ob.matrix_world) if M is not None else ob.matrix_world)
    if weld:
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7)
    bm.verts.ensure_lookup_table(); bm.faces.ensure_lookup_table(); bm.normal_update()
    return bm


def tree_in(ob, M=None, cap=True):
    """BVH of ob's shell in frame M, seams welded, open rims capped (inside tests / signed distance)"""
    bm = bm_in(ob, M, weld=True)
    if cap:
        bnd = [e for e in bm.edges if e.is_boundary]
        if bnd:
            bmesh.ops.holes_fill(bm, edges=bnd, sides=0)
    bm.normal_update(); t = BVHTree.FromBMesh(bm); bm.free(); return t


def tree_np(X, F):
    return BVHTree.FromPolygons([tuple(map(float, p)) for p in X], [tuple(map(int, f)) for f in F])


def vnormals(X, F):
    F = np.asarray(F); tri = X[F]
    fn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n = np.zeros_like(X)
    for k in range(3):
        np.add.at(n, F[:, k], fn)
    return n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-18)


def volume(X, F):
    F = np.asarray(F); tri = X[F]
    return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)


def neighbours(n, F):
    nb = [set() for _ in range(n)]
    for f in F:
        for a in range(len(f)):
            for b in range(len(f)):
                if a != b:
                    nb[f[a]].add(f[b])
    return [np.array(sorted(s), int) for s in nb]


def laplacian_matrix(nb):
    """sparse-free umbrella operator as (rows, cols, w) for np.add.at use"""
    rows = np.concatenate([np.full(len(s), i) for i, s in enumerate(nb)]); cols = np.concatenate(nb)
    w = np.concatenate([np.full(len(s), 1.0 / max(len(s), 1)) for s in nb])
    return rows, cols, w


def umbrella(X, L):
    rows, cols, w = L
    avg = np.zeros_like(X); np.add.at(avg, rows, X[cols] * w[:, None])
    return avg - X


# ---------------------------------------------------------------- signed distances
DIRS = [Vector(d).normalized() for d in ((0.31, 0.52, 0.79), (-0.6, 0.2, -0.77), (0.1, -0.95, 0.3),
                                         (0.83, -0.29, -0.47), (-0.44, -0.71, 0.55))]


def inside(tree, p):
    """ray-parity inside test (majority of five rays: three can be fooled by a ray grazing an edge), robust to
    capped rims and collapsed walls"""
    votes = 0
    for d in DIRS:
        c = 0; q = p.copy()
        while True:
            h = tree.ray_cast(q, d, 5.0)
            if h[0] is None:
                break
            c += 1; q = h[0] + d * 1e-7
        votes += c % 2
    return votes >= 3


def signed_dist(tree, P, band=0.006):
    """(d, q, out) per point: distance to the shell (- inside, by ray parity, tested only within `band`), the nearest
    point, and the unit direction that leads out (away from the shell, or out of it when inside)"""
    d = np.empty(len(P)); q = np.empty_like(P); out = np.empty_like(P)
    for i, p in enumerate(P):
        v = Vector(p); hit = tree.find_nearest(v)
        if hit[0] is None:
            d[i] = 1e9; q[i] = p; out[i] = (0, 0, 1); continue
        qq = hit[0]; dv = v - qq; dl = dv.length
        if dl > band:
            d[i] = dl; q[i] = qq; out[i] = dv / dl; continue
        ins = inside(tree, v)
        if dl < 1e-9:
            dv = hit[1].copy(); dl_ = 0.0
        else:
            dl_ = dl
        dirv = (dv.normalized() if dl_ > 0 else dv) * (-1.0 if ins else 1.0)
        d[i] = -dl_ if ins else dl_; q[i] = qq; out[i] = dirv
    return d, q, out


class Obstacles:
    """hard limits for a growing organ. solid: [(tree, gap)] - stay `gap` outside these closed shells (bones);
    container: [(tree, wall)] - stay `wall` inside (the skin: the abdominal wall); planes: [(point, normal, mask
    fn or None)] - stay on the normal's side (mask(X) -> bool per vertex selects the vertices it applies to)"""

    def __init__(self, solid=(), container=(), planes=()):
        self.solid = list(solid); self.container = list(container); self.planes = list(planes)

    def project(self, X, movable, near=0.02, debug=None):
        X = X.copy(); hits = np.zeros(len(X), bool)
        for k, (tree, gap) in enumerate(self.solid):
            d, q, n = signed_dist(tree, X)
            bad = (d < gap) & movable
            if debug is not None and bad.any():
                debug.append(("solid", k, int(bad.sum()), float(np.linalg.norm(q[bad] + n[bad] * gap - X[bad], axis=1).max())))
            X[bad] = q[bad] + n[bad] * gap; hits |= bad
        for tree, wall in self.container:            # the skin is open (a body section): side by its normal
            d, q, n = signed_dist(tree, X, band=0.0)
            for i in np.nonzero(movable & (np.abs(d) < near + wall))[0]:
                nn = tree.find_nearest(Vector(X[i]))[1]
                s_ = float((X[i] - q[i]) @ np.array(nn))
                if s_ > -wall:
                    X[i] = q[i] - np.array(nn) * wall; hits[i] = True
        for p0, nrm, mask in self.planes:
            p0 = np.asarray(p0); nrm = np.asarray(nrm) / np.linalg.norm(nrm)
            s = (X - p0) @ nrm
            m = (s < 0) & movable
            if mask is not None:
                m &= mask(X)
            X[m] -= np.outer(s[m], nrm); hits |= m
        return X, hits


def icosphere(subdiv=5):
    """unit icosphere as (X, F) arrays"""
    bm = bmesh.new(); bmesh.ops.create_icosphere(bm, subdivisions=subdiv, radius=1.0)
    X = np.array([tuple(v.co) for v in bm.verts]); F = np.array([[v.index for v in f.verts] for f in bm.faces])
    bm.free(); return X / np.linalg.norm(X, axis=1, keepdims=True), F


def grow(X, F, mob, V_target, obst, iters=1500, gain=0.25, dmax=6e-4, smooth=0.25, fair=0.06, tol=0.003,
         check_every=3, log=None):
    """pressure-grow (or shrink) the closed shell X toward volume V_target: each pass pushes every vertex along its
    normal (proportional + integral control on the volume error, at most dmax per pass, scaled by mob), relaxes the
    membrane (tangential umbrella `smooth` for even triangles, full umbrella `fair` for a fair surface - a soap
    bubble pressed against its surroundings) and projects back out of the obstacles. mob 0 = pinned."""
    X = X.copy(); ref = X.copy(); L = laplacian_matrix(neighbours(len(X), F)); movable = mob > 1e-6
    integ = 0.0; R = (3 * abs(V_target) / (4 * math.pi)) ** (1 / 3)
    for it in range(iters):
        V = volume(X, F); err = (V_target - V) / V_target
        if abs(err) < tol and it > 30:
            break
        integ = float(np.clip(integ + 0.02 * err, -0.5, 0.5))
        n = vnormals(X, F)
        dn = float(np.clip(gain * (err + integ) * R, -dmax, dmax))
        X += n * (dn * mob)[:, None]
        lap = umbrella(X, L); ln = np.einsum("ij,ij->i", lap, n)[:, None] * n
        X += ((lap - ln) * smooth + lap * fair) * np.sqrt(mob)[:, None]   # (smoothing reaches further into the
                                                                           # still zone's rim than the pressure)
        if it % check_every == 0 or it > iters - 5:
            X, _ = obst.project(X, movable)
        X[~movable] = ref[~movable]
        if log is not None and it % 100 == 0:
            log.append((it, round(V * 1e6, 1), round(dn * 1e3, 3)))
    for _ in range(4):                           # the limits are applied in turn: repeat until none moves a vertex
        X, hit = obst.project(X, movable); X[~movable] = ref[~movable]
        if not hit.any():
            break
    return X, volume(X, F)


# ---------------------------------------------------------------- rings, loops, triangulated annuli
def slit_ring(K, centre, e_lat, e_ap, half_w, gap, bend=0.0):
    """K points round a closed slit (lips `gap` apart, `half_w` from the middle to each tip) in the plane e_lat /
    e_ap, bent into a crescent `bend` toward +e_ap at the middle. Point k sits at angle 2 pi k / K (the open ring
    uses the same angles, so vertex k always means the same spot)"""
    c = np.array(centre); a = np.array(e_lat); b = np.array(e_ap); out = []
    for k in range(K):
        ph = 2 * math.pi * k / K; x = half_w * math.cos(ph)
        y = gap * math.sin(ph) + bend * (1 - (x / half_w) ** 2)
        out.append(c + a * x + b * y)
    return np.array(out)


def round_ring(K, centre, e_lat, e_ap, r, rx=None):
    c = np.array(centre); a = np.array(e_lat); b = np.array(e_ap)
    return np.array([c + a * (rx or r) * math.cos(2 * math.pi * k / K) + b * r * math.sin(2 * math.pi * k / K)
                     for k in range(K)])


def boundary_loops(bm):
    """ordered boundary vertex loops of a bmesh"""
    E = {e for e in bm.edges if e.is_boundary}; seen = set(); loops = []
    for e0 in E:
        if e0 in seen:
            continue
        v0 = e0.verts[0]; loop = [v0]; prev = None; cur = v0; e = e0
        while True:
            seen.add(e); nxt = e.other_vert(cur)
            if nxt is v0:
                break
            loop.append(nxt); cur = nxt
            cand = [g for g in cur.link_edges if g in E and g not in seen]
            if not cand:
                break
            e = cand[0]
        loops.append(loop)
    return loops


def fill_annulus(bm, outer, inner_pts, origin, ax_u, ax_v, flip=False):
    """triangulate between the existing boundary loop `outer` (BMVerts) and new verts at inner_pts (a closed loop,
    the hole that stays open) by a constrained Delaunay triangulation in the plane (origin, ax_u, ax_v). Returns the
    new inner BMVerts in order."""
    from mathutils.geometry import delaunay_2d_cdt
    o = Vector(origin); u = Vector(ax_u); v = Vector(ax_v)
    inner = [bm.verts.new(Vector(p)) for p in inner_pts]
    allv = list(outer) + inner
    co2 = [Vector(((p.co - o).dot(u), (p.co - o).dot(v))) for p in allv]
    no, ni = len(outer), len(inner)
    edges = [(i, (i + 1) % no) for i in range(no)] + [(no + i, no + (i + 1) % ni) for i in range(ni)]
    faces = [list(range(no)), list(range(no, no + ni))]
    vo, eo, fo, ov, oe, of = delaunay_2d_cdt(co2, edges, faces, 2, 1e-12, True)
    vid = []
    for k, ids in enumerate(ov):
        if not ids:
            raise RuntimeError("fill_annulus: the triangulation added a vertex (loops cross in the plane)")
        vid.append(allv[ids[0]])
    n_new = 0
    po = [tuple(c) for c in co2[:no]]; pi = [tuple(c) for c in co2[no:]]
    for f in fo:
        cen = sum((vo[i] for i in f), Vector((0, 0))) / len(f)
        if not point_in_poly(cen, po) or point_in_poly(cen, pi):
            continue                             # only the ring between the loops
        vs = [vid[i] for i in f]
        if flip:
            vs = vs[::-1]
        try:
            bm.faces.new(vs); n_new += 1
        except ValueError:
            pass
    bm.normal_update()
    return inner, n_new


# ---------------------------------------------------------------- a chart of the skin round a crest (the meatus)
class CrestChart:
    """2D coordinates (u, v) on the skin round point M: u along t_ap (the slit's direction), v the arc length
    along the skin's cross-section (the plane square to t_ap) from the crest point, + toward e_lat. On a flat patch
    this is the ordinary tangent plane; in the vestibule's slot (the labial walls meet in a knife-edge crest) it
    unfolds the two walls on either side of the crest, which no projection can. Sections are cached per u."""

    def __init__(self, bm, M, t_ap, e_lat, radius=0.008):
        self.M = Vector(M); self.t = Vector(t_ap).normalized(); self.e = Vector(e_lat).normalized()
        self.tris = []                           # ((co, co, co), (key, key, key)): keys identify the vertices
        for f in bm.faces:
            if (f.calc_center_median() - self.M).length < radius:
                vs = list(f.verts)
                for k in range(1, len(vs) - 1):
                    t_ = (vs[0], vs[k], vs[k + 1])
                    self.tris.append((tuple(v.co.copy() for v in t_), tuple(v.index for v in t_)))
        self.cache = {}

    def section(self, u):
        key = round(u, 7)
        if key in self.cache:
            return self.cache[key]
        P0 = self.M + self.t * u; segs = []; skeys = []
        for tri, ids in self.tris:
            d = [(p - P0).dot(self.t) for p in tri]; pts = []; ks = []
            for a in range(3):
                b = (a + 1) % 3
                if (d[a] < 0) != (d[b] < 0):
                    lo, hi = (a, b) if ids[a] < ids[b] else (b, a)    # same point from both triangles
                    pts.append(tri[lo].lerp(tri[hi], d[lo] / (d[lo] - d[hi]))); ks.append((ids[lo], ids[hi]))
            if len(pts) == 2:
                segs.append(pts); skeys.append(ks)
        # chain the segments into polylines: crossings on the same mesh edge are the same point
        adj = {}
        for i, (ka, kb) in enumerate(skeys):
            adj.setdefault(ka, []).append((i, 0)); adj.setdefault(kb, []).append((i, 1))
        key_of = {}
        for i, (ka, kb) in enumerate(skeys):
            key_of[(i, 0)] = ka; key_of[(i, 1)] = kb
        used = set(); lines = []
        for i0 in range(len(segs)):
            if i0 in used:
                continue
            used.add(i0); line = [segs[i0][0], segs[i0][1]]; ends = [key_of[(i0, 0)], key_of[(i0, 1)]]
            for direction in (1, -1):
                while True:
                    endk = ends[-1] if direction == 1 else ends[0]
                    nxt = [(i, e) for i, e in adj.get(endk, []) if i not in used]
                    if not nxt:
                        break
                    i, e = nxt[0]; used.add(i); p = segs[i][1 - e]; k2 = key_of[(i, 1 - e)]
                    if direction == 1:
                        line.append(p); ends.append(k2)
                    else:
                        line.insert(0, p); ends.insert(0, k2)
            lines.append(line)
        # close T-junction gaps (one wall's vertex on the middle of the other's edge): join ends within 20 um
        merged = True
        while merged and len(lines) > 1:
            merged = False
            for i in range(len(lines)):
                for j in range(len(lines)):
                    if i == j:
                        continue
                    A, B = lines[i], lines[j]
                    for a_end, b_end in ((-1, 0), (-1, -1), (0, 0), (0, -1)):
                        if (A[a_end] - B[b_end]).length < 2e-5:
                            A2 = A if a_end == -1 else A[::-1]
                            B2 = B if b_end == 0 else B[::-1]
                            lines[i] = A2 + B2[1:]; lines.pop(j); merged = True
                            break
                    if merged:
                        break
                if merged:
                    break
        # the polyline through the crest: the one passing nearest the chart's axis point
        best = None
        for line in lines:
            for p in line:
                dd = (p - P0).length
                if best is None or dd < best[0]:
                    best = (dd, line)
        line = best[1]
        s = [0.0]
        for a, b in zip(line, line[1:]):
            s.append(s[-1] + (b - a).length)
        # the crest: highest point (along e_lat x t_ap, into the body) near the axis point - the apex of a knife
        # edge; on a flat patch simply the nearest point
        up = self.e.cross(self.t)
        k = max(range(len(line)), key=lambda i: line[i].dot(up) - (line[i] - P0).length)
        sc = s[k]
        # + toward e_lat: compare the points a little way along either side
        def at(t_):
            t_ = max(0.0, min(s[-1], t_))
            for i in range(len(s) - 1):
                if s[i] <= t_ <= s[i + 1]:
                    return line[i].lerp(line[i + 1], 0.0 if s[i + 1] == s[i] else (t_ - s[i]) / (s[i + 1] - s[i]))
            return line[-1]
        sign = 1.0 if (at(sc + 3e-4) - at(sc - 3e-4)).dot(self.e) >= 0 else -1.0
        res = (line, s, sc, sign); self.cache[key] = res
        return res

    def to3d(self, u, v):
        line, s, sc, sign = self.section(u); target = sc + sign * v
        for k in range(len(s) - 1):
            if s[k] <= target <= s[k + 1] or k == len(s) - 2:
                t = 0.0 if s[k + 1] == s[k] else (target - s[k]) / (s[k + 1] - s[k])
                return line[k].lerp(line[k + 1], max(0.0, min(1.0, t)))

    def to2d(self, p):
        p = Vector(p); u = (p - self.M).dot(self.t); line, s, sc, sign = self.section(u); best = None
        for k in range(len(line) - 1):
            a, b = line[k], line[k + 1]; ab = b - a; L2 = ab.dot(ab)
            t = 0.0 if L2 == 0 else max(0.0, min(1.0, (p - a).dot(ab) / L2))
            dd = (a + ab * t - p).length
            if best is None or dd < best[0]:
                best = (dd, s[k] + (s[k + 1] - s[k]) * t)
        if best[0] > 2e-5:                       # not on the crest's section (a broken piece): straight distance
            c = self.to3d(u, 0.0); dv = p - c
            return u, math.copysign(dv.length, dv.dot(self.e)), best[0]
        return u, sign * (best[1] - sc), best[0]


def point_in_poly(pt, poly):
    x, y = pt; inside_ = False; n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside_ = not inside_
    return inside_


def cut_meatus(skin, M, t_ap, e_lat, K, half, lip, reg_u, reg_v, step=0.00022):
    """open a slit of K vertices (half-length `half` along t_ap, lips `lip` apart, vertex k at angle 2 pi k / K:
    u = half cos, v = -lip/2 sin) in the skin round M; the skin within |u| < reg_u, |v| < reg_v is re-triangulated
    (constrained Delaunay in the chart, extra points every `step`), new vertices on the original surface, UVs
    sampled on the same side of the midline seam, every existing shape key carried over (= Basis there).
    Returns (slit vertex indices in order, chart, {new vertex index: (u, v)}, report)."""
    from mathutils.geometry import delaunay_2d_cdt
    me = skin.data; rep = {}
    bm = bmesh.new(); bm.from_mesh(me)
    # refine the skin round M first (the region may sit on coarse faces): split long edges near M, a few passes;
    # bmesh interpolates UVs and shape keys on the new vertices, the surface itself does not change
    R_ref = 2.0 * max(reg_u, reg_v); n_split = 0
    for _ in range(4):
        long_ = [e for e in bm.edges if e.calc_length() > 2.2 * step and
                 min((v.co - Vector(M)).length for v in e.verts) < R_ref]
        if not long_:
            break
        bmesh.ops.subdivide_edges(bm, edges=long_, cuts=1, use_grid_fill=True); n_split += len(long_)
        bmesh.ops.triangulate(bm, faces=[f for f in bm.faces if len(f.verts) > 4])
    rep["edges_split"] = n_split
    # the crest carries an unwelded seam (coincident duplicates): weld it round M so the walls are one surface there
    nearv = [v for v in bm.verts if (v.co - Vector(M)).length < R_ref]
    n0 = len(bm.verts); bmesh.ops.remove_doubles(bm, verts=nearv, dist=5e-6); rep["welded"] = n0 - len(bm.verts)
    # ...and zero-thickness flaps standing on it in the midline plane (faces with every vertex on x = 0): hidden in
    # the closed slot, they make the crest branch
    fins = [f for f in bm.faces if all(abs(v.co.x) < 1e-7 for v in f.verts)
            and (f.calc_center_median() - Vector(M)).length < R_ref]
    bmesh.ops.delete(bm, geom=fins, context='FACES')
    bmesh.ops.delete(bm, geom=[v for v in bm.verts if not v.link_faces], context='VERTS')
    rep["fins_removed"] = len(fins)
    bm.verts.ensure_lookup_table(); bm.faces.ensure_lookup_table(); bm.normal_update(); bm.verts.index_update()
    chart = CrestChart(bm, M, t_ap, e_lat, radius=max(reg_u, reg_v) * 3)
    # originals for sampling: one BVH per side of the midline (x = 0 is a UV seam)
    uv_layers = list(bm.loops.layers.uv.values())
    side_faces = {1: [f for f in bm.faces if f.calc_center_median().x >= 0], -1: [f for f in bm.faces if f.calc_center_median().x <= 0]}
    samp = {}
    for sd, fs in side_faces.items():
        near = [f for f in fs if (f.calc_center_median() - Vector(M)).length < 4 * max(reg_u, reg_v)]
        tris = []
        for f in near:
            loops = list(f.loops)
            for k in range(1, len(loops) - 1):
                t_ = (loops[0], loops[k], loops[k + 1])
                tris.append(([l.vert.co.copy() for l in t_], [[l[lay].uv.copy() for lay in uv_layers] for l in t_]))
        tree = BVHTree.FromPolygons([tuple(c) for t_ in tris for c in t_[0]], [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(len(tris))])
        samp[sd] = (tree, tris)

    def sample_uv(p, sd):
        from mathutils.geometry import barycentric_transform
        tree, tris = samp[sd]; q, n, i, d = tree.find_nearest(p)
        (a, b, c), uvs = tris[i]
        w = barycentric_transform(q, a, b, c, Vector((1, 0, 0)), Vector((0, 1, 0)), Vector((0, 0, 1)))
        return [sum((uvs[j][li] * w[j] for j in range(3)), Vector((0, 0))) for li in range(len(uv_layers))]

    # region: faces grown outward from the one under M while their centre lies inside the (u, v) box on the chart
    near_f = [f for f in bm.faces if (f.calc_center_median() - Vector(M)).length < 2.5 * max(reg_u, reg_v)]
    fuv = {}
    for f in near_f:
        u, vv, dd = chart.to2d(f.calc_center_median())
        fuv[f] = (u, vv, dd)
    ok_ = lambda f: f in fuv and abs(fuv[f][0]) < reg_u and abs(fuv[f][1]) < reg_v and fuv[f][2] < 2e-5
    f0 = min(near_f, key=lambda f: (f.calc_center_median() - Vector(M)).length)
    kill = set(); st = [f0]
    while st:
        f = st.pop()
        if f in kill or not ok_(f):
            continue
        kill.add(f)
        for e in f.edges:
            st += [g for g in e.link_faces if g not in kill]
    # islands: kept faces enclosed by the region (the flood skipped them) join it
    rest = [f for f in near_f if f not in kill]; seen = set()
    for f0_ in rest:
        if f0_ in seen:
            continue
        comp = []; st = [f0_]; touches_out = False
        while st:
            f = st.pop()
            if f in seen:
                continue
            seen.add(f); comp.append(f)
            for e in f.edges:
                for g in e.link_faces:
                    if g in kill or g in seen:
                        continue
                    if g not in fuv:
                        touches_out = True
                    else:
                        st.append(g)
        if not touches_out and len(comp) < 40:
            kill |= set(comp)
    kill = list(kill)
    uv_of = {}
    for f in kill:
        for v in f.verts:
            if v not in uv_of:
                u, vv, dd = chart.to2d(v.co); uv_of[v] = (u, vv)
    inreg = {v for v in uv_of if all(g in kill for g in v.link_faces)}
    rep["faces_removed"] = len(kill)
    normals_ref = {f: f.normal.copy() for f in bm.faces}
    bmesh.ops.delete(bm, geom=kill, context='FACES_ONLY')
    bmesh.ops.delete(bm, geom=[v for v in inreg if v.is_valid and not v.link_faces], context='VERTS')
    allloops = boundary_loops(bm)
    loops = [lp for lp in allloops if all(v in uv_of for v in lp) and len(lp) > 2]
    if len(loops) != 1:
        info = [(len(lp), sum(1 for v in lp if v in uv_of), tuple(round(x * 1000, 2) for x in sum((v.co for v in lp), Vector()) / len(lp))) for lp in allloops]
        raise RuntimeError(f"meatus region: {len(loops)} boundary loops: {info}")
    outer = loops[0]; o2 = [uv_of[v] for v in outer]
    # slit + inner points
    slit2 = [(half * math.cos(2 * math.pi * k / K), -lip / 2 * math.sin(2 * math.pi * k / K)) for k in range(K)]
    pts = []
    nu, nv = int(reg_u / step), int(reg_v / step)
    allp = o2 + slit2
    for i in range(-nu, nu + 1):
        for j in range(-nv, nv + 1):
            p = (i * step, j * step)
            if not point_in_poly(p, o2) or (abs(p[0]) <= half + step * 0.6 and abs(p[1]) < step * 0.6):
                continue
            if min(math.hypot(p[0] - a, p[1] - b) for a, b in allp) < step * 0.55:
                continue
            pts.append(p)
    no, ni = len(o2), K
    co2 = [Vector(p) for p in o2 + slit2 + pts]
    edges = [(i, (i + 1) % no) for i in range(no)] + [(no + i, no + (i + 1) % ni) for i in range(ni)]
    vo, eo, fo, ov, oe, of = delaunay_2d_cdt(co2, edges, [list(range(no)), list(range(no, no + ni))], 2, 1e-12, True)
    newv = []
    for k, ids in enumerate(ov):
        if not ids:
            json.dump({"outer": o2, "slit": slit2, "pts": pts, "outer3d": [tuple(v.co) for v in outer]},
                      open(os.environ.get("MEATUS_DEBUG", "/tmp/meatus_debug.json"), "w"))
            raise RuntimeError("cut_meatus: triangulation added a vertex (debug: MEATUS_DEBUG json)")
        i0 = ids[0]
        if i0 < no:
            newv.append(outer[i0])
        else:
            u, vv = co2[i0]; nvx = bm.verts.new(chart.to3d(u, vv)); newv.append(nvx); uv_of[nvx] = (u, vv)
    slit_v = [newv[[k for k, ids in enumerate(ov) if ids and ids[0] == no + j][0]] for j in range(K)]
    created = []
    po = [tuple(p) for p in o2]
    for f in fo:
        cen = sum((vo[i] for i in f), Vector((0, 0))) / len(f)
        if not point_in_poly(cen, po) or point_in_poly(cen, slit2):
            continue
        try:
            created.append(bm.faces.new([newv[i] for i in f]))
        except ValueError:
            pass
    bm.normal_update()
    # winding: agree with the removed faces' neighbours (the surrounding skin)
    agree = 0
    for f in created:
        for e in f.edges:
            for g in e.link_faces:
                if g is not f and g in normals_ref:
                    agree += 1 if g.normal.dot(f.normal) > 0 else -1
    if agree < 0:
        bmesh.ops.reverse_faces(bm, faces=created)
    rep["faces_added"] = len(created); rep["winding_votes"] = agree
    for f in created:
        sd = 1 if f.calc_center_median().x >= 0 else -1
        for l in f.loops:
            for lay, uvv in zip(uv_layers, sample_uv(l.vert.co, sd)):
                l[lay].uv = uvv
    # shape keys: the new vertices take their Basis position in every key (the region does not move in them)
    for lay in bm.verts.layers.shape.values():
        for f in created:
            for v in f.verts:
                if v not in outer:
                    v[lay] = v.co.copy()
    bm.verts.index_update()
    bm.to_mesh(me); me.update()
    slit_idx = [v.index for v in slit_v]
    new_uv = {v.index: uv_of[v] for f in created for v in f.verts if v not in outer}
    bm.free()
    return slit_idx, chart, new_uv, rep


# ---------------------------------------------------------------- neighbours that give way
def make_room(objs, M_to, M_from, expander_trees, pinned, rigid, gap=0.0015, decay=0.9992, iters=700, rounds=6,
              band=0.06, log=None, expanders=None, carry_from=None, carry_sigma=0.010, carry_reach=0.015, start=None):
    """corrective positions for soft organs `objs` while an organ grows through stages. At each stage, per round:
      - every soft vertex inside (or within `gap` of) the grown shell is pushed just outside it;
      - where the grown shell bulges through the middle of a soft organ's face (coarse faces: no soft vertex is
        inside, yet the shell is), that face's vertices are pushed past the shell's vertex (expanders: per stage the
        shell's vertices + outward normals, needed for this);
      - soft organs pushed into each other (a vertex newly inside another soft organ) step back out of it;
    the pushes spread over each organ's surface (a decaying average: a dent that eases out over a few cm, shared by
    organs welded together, so their seams stay closed), then the `rigid` Obstacles are applied. Repeated until clear.
    M_to: world -> work frame, M_from: work frame -> world. pinned(name, P) -> bool array (never moves).
    Returns {name: [positions (local coords) per stage]} and the worst remaining depth (mm) per stage."""
    P_all = []; owner = []; F_all = []; off = 0
    for ob in objs:
        Mw_ = ob.matrix_world
        src = [Vector(p) for p in start[ob.name]] if start and ob.name in start else [v.co for v in ob.data.vertices]
        P = np.array([tuple(M_to @ (Mw_ @ c)) for c in src])   # (start: local positions to begin from)
        P_all.append(P); owner.append((ob, off, len(P)))
        F_all.append([[off + i for i in p.vertices] for p in ob.data.polygons]); off += len(P)
    P = np.concatenate(P_all)
    kd = kdtree.KDTree(len(P))
    for i, p in enumerate(P):
        kd.insert(Vector(p), i)
    kd.balance()
    rep_idx = np.array([min(k_ for _, k_, _ in kd.find_range(Vector(p), 1e-7)) for p in P])
    uniq, inv = np.unique(rep_idx, return_inverse=True)
    X = P[uniq]
    Fo = [[[int(inv[i]) for i in f] for f in fl] for fl in F_all]      # faces per organ, union indices
    F = [f for fl in Fo for f in fl]
    nb = neighbours(len(X), F); L = laplacian_matrix(nb)
    pin = np.zeros(len(X), bool); org_v = []
    for ob, o_, n_ in owner:
        ids = inv[o_:o_ + n_]; pin[ids] |= pinned(ob.name, P[o_:o_ + n_]); org_v.append(np.unique(ids))

    def organ_trees(X_):
        out = []
        for fl in Fo:
            bm = bmesh.new(); vm = {}; fmap = []
            for fi, f in enumerate(fl):
                vs = []
                for i in f:
                    if i not in vm:
                        vm[i] = bm.verts.new(Vector(X_[i]))
                    vs.append(vm[i])
                try:
                    bm.faces.new(vs); fmap.append(fi)
                except ValueError:
                    pass
            bnd = [e for e in bm.edges if e.is_boundary]
            if bnd:
                bmesh.ops.holes_fill(bm, edges=bnd, sides=0)
            t_ = BVHTree.FromBMesh(bm)
            out.append((t_, fmap)); bm.free()
        return out

    def inside_other(X_, trees_):
        """(vertex, organ) pairs: vertex of one soft organ inside another's shell"""
        res = []
        for a, ids in enumerate(org_v):
            for b, tb in enumerate(trees_):
                if a == b:
                    continue
                for i in ids:
                    p = Vector(X_[i]); h = tb.find_nearest(p)
                    if h[0] is not None and h[3] < 0.004 and h[3] > 2e-5 and inside(tb, p):
                        res.append((i, b))
        return res
    base_in = set(inside_other(X, [t_ for t_, _ in organ_trees(X)]))   # contacts that exist at rest (seams, atlas) are not ours
    stages = []; resid = []
    for k, tree in enumerate(expander_trees):
        if carry_from is not None and expanders is not None:
            # carried along: soft tissue next to the growing wall moves the way the wall next to it moves (a
            # weighted average of the wall's motion since the last stage, fading with distance from the wall) -
            # the bowel rides up with the dome as one piece instead of being pushed every which way
            prevP = carry_from if k == 0 else expanders[k - 1][0]; curP = expanders[k][0]; dW = curP - prevP
            kdw = kdtree.KDTree(len(prevP))
            for i, p in enumerate(prevP):
                kdw.insert(Vector(p), i)
            kdw.balance()
            Dc = np.zeros_like(X)
            for i in np.nonzero(~pin)[0]:
                near = kdw.find_range(Vector(X[i]), 2.5 * carry_sigma + carry_reach)
                if not near:
                    continue
                dmin = min(d_ for _, _, d_ in near)
                w = np.array([math.exp(-((d_ - dmin) / carry_sigma) ** 2) for _, _, d_ in near])
                Dc[i] = (w[:, None] * dW[[j for _, j, _ in near]]).sum(0) / w.sum() * math.exp(-(dmin / carry_reach) ** 2)
            X = X + Dc
            X, _ = rigid.project(X, ~pin)
        for r in range(rounds):
            D = np.zeros_like(X); C = np.zeros(len(X), bool)
            d, q, out = signed_dist(tree, X, band=band)
            m = (d < gap) & ~pin; D[m] = q[m] + out[m] * gap - X[m]; C |= m
            tf_ = organ_trees(X); trees_ = [t_ for t_, _ in tf_]
            n_bulge = 0
            if expanders is not None:            # the shell bulging through soft faces
                EP, EN = expanders[k]
                for ob_i, (tb, fmap) in enumerate(tf_):
                    for p, nrm in zip(EP, EN):
                        pv = Vector(p); h = tb.find_nearest(pv)
                        if h[0] is None or h[3] > 0.004 or not inside(tb, pv):
                            continue
                        fl = Fo[ob_i][fmap[h[2]]] if h[2] < len(fmap) else None
                        if fl is None:
                            continue
                        push = np.array(pv - h[0]) + np.array(nrm) * gap
                        for i in fl:
                            if not pin[i] and np.dot(push, push) > np.dot(D[i], D[i]):
                                D[i] = push; C[i] = True
                        n_bulge += 1
            n_soft = 0
            for i, b in inside_other(X, trees_):   # both give way, half each
                if (i, b) in base_in or pin[i]:
                    continue
                h = trees_[b].find_nearest(Vector(X[i])); dv = np.array(h[0]) - X[i]; ln = np.linalg.norm(dv)
                if ln < 1e-12:
                    continue
                full = dv + dv / ln * gap * 0.5
                cand = [(i, full * 0.5)]
                fb = tf_[b][1]
                if h[2] < len(fb):
                    cand += [(j, -full * 0.5) for j in Fo[b][fb[h[2]]]]
                for j, push in cand:
                    if pin[j]:
                        continue
                    if np.dot(push, push) > np.dot(D[j], D[j]):
                        D[j] = push; C[j] = True
                n_soft += 1
            if log is not None:
                log.append((k, r, int(m.sum()), n_bulge, n_soft))
            if not C.any():
                break
            fixed = C | pin; free = ~fixed
            rows, cols, w = L
            for it in range(iters):
                avg = np.zeros_like(D); np.add.at(avg, rows, D[cols] * w[:, None])
                D[free] = decay * avg[free]
            X = X + D
            X, _ = rigid.project(X, ~pin)
        d, _, _ = signed_dist(tree, X, band=band)
        resid.append(round(float(max(0.0, -d.min())) * 1000, 2))
        stages.append(X.copy())
    out_ = {}
    for ob, o_, n_ in owner:
        Mi_ = ob.matrix_world.inverted()
        out_[ob.name] = [np.array([tuple(Mi_ @ (M_from @ Vector(p))) for p in Xk[inv[o_:o_ + n_]]]) for Xk in stages]
    return out_, resid


# ---------------------------------------------------------------- shape keys + drivers
def add_key(ob, name, positions=None, relative=None):
    """shape key `name` on ob (Basis made if missing); positions: (n, 3) local coords or {index: Vector}"""
    if ob.data.shape_keys is None:
        ob.shape_key_add(name="Basis", from_mix=False)
    kb = ob.data.shape_keys.key_blocks.get(name) or ob.shape_key_add(name=name, from_mix=False)
    if relative is not None:
        kb.relative_key = ob.data.shape_keys.key_blocks[relative]
    if positions is not None:
        if isinstance(positions, dict):
            for i, p in positions.items():
                kb.data[i].co = Vector(p)
        else:
            kb.data.foreach_set("co", np.asarray(positions, dtype=np.float32).ravel())
    kb.slider_min, kb.slider_max = 0.0, 1.0
    return kb


def drive(kb, ctl, prop, expr):
    """kb.value = expr of the control's custom property (variable f): a simple expression, no Python needed"""
    fc = kb.driver_add("value"); dr = fc.driver; dr.type = 'SCRIPTED'
    while dr.variables:
        dr.variables.remove(dr.variables[0])
    v = dr.variables.new(); v.name = "f"; v.type = 'SINGLE_PROP'
    v.targets[0].id = ctl; v.targets[0].data_path = f'["{prop}"]'
    dr.expression = expr
    return fc


def settle(X, faces, fixed, trees, gap=0.0006, decay=0.985, iters=300, rounds=8, band=0.02):
    """the grown organ itself gives way where it still presses into its neighbours (an organ trapped against bone
    cannot move further): vertices inside or within `gap` of any shell in `trees` move back out, the correction eases
    out over the organ's surface (fixed vertices stay). Returns (X, worst depth before, after) in metres."""
    X = np.array(X, float); nb = neighbours(len(X), faces); L = laplacian_matrix(nb); worst0 = None
    for r in range(rounds):
        D = np.zeros_like(X); C = np.zeros(len(X), bool); worst = 0.0
        for tree in trees:
            d, q, out = signed_dist(tree, X, band=band)
            m = (d < gap) & ~fixed
            worst = max(worst, float(max(0.0, -d[m].min())) if m.any() else 0.0)
            push = q + out * gap - X
            upd = m & (np.einsum("ij,ij->i", push, push) > np.einsum("ij,ij->i", D, D))
            D[upd] = push[upd]; C |= m
        if worst0 is None:
            worst0 = worst
        if not C.any():
            break
        free = ~(C | fixed); rows, cols, w = L
        for it in range(iters):
            avg = np.zeros_like(D); np.add.at(avg, rows, D[cols] * w[:, None])
            D[free] = decay * avg[free]
        X = X + D
    after = 0.0                                  # (movable vertices: the still ones never move)
    for tree in trees:
        d, _, _ = signed_dist(tree, X, band=band); after = max(after, float(max(0.0, -d[~fixed].min())))
    return X, worst0 or 0.0, after


# ---------------------------------------------------------------- pelvis: measurements, female reshaping
def pelvis_metrics(bones_P):
    """bones_P: {name: (n,3) atlas points}. Subpubic angle (deg, from the medial edges of the inferior pubic rami in
    the coronal plane), intertuberous and brim (inlet transverse) widths (mm), sacrum height (mm)"""
    H = np.concatenate([bones_P["AN_HipBone_L"], bones_P["AN_HipBone_R"]])
    zb = H[:, 2].min()
    out = {}
    lines = {}
    for s, P in (("L", bones_P["AN_HipBone_L"]), ("R", bones_P["AN_HipBone_R"])):
        F = P[P[:, 1] < -0.05]                    # the front: pubis and its rami
        pts = []
        for z in np.linspace(zb + 0.004, zb + 0.022, 7):
            sl = F[np.abs(F[:, 2] - z) < 0.0015]
            if len(sl):
                pts.append((np.abs(sl[:, 0]).min(), z))
        pts = np.array(pts); a, b = np.polyfit(pts[:, 1], pts[:, 0], 1)   # |x| = a z + b
        lines[s] = a
    out["subpubic_angle_deg"] = round(math.degrees(2 * math.atan(abs((lines["L"] + lines["R"]) / 2))), 1)
    tub = []
    for P in (bones_P["AN_HipBone_L"], bones_P["AN_HipBone_R"]):
        q = P[P[:, 1] > -0.09]; tub.append(q[q[:, 2].argmin()])
    out["intertuberous_mm"] = round(float(abs(tub[0][0] - tub[1][0])) * 1000, 1)
    brim = []
    for z in np.linspace(0.835, 0.875, 9):
        sl = H[(np.abs(H[:, 2] - z) < 0.002) & (H[:, 1] > -0.10) & (H[:, 1] < -0.04)]
        if len(sl):
            L = sl[sl[:, 0] > 0][:, 0]; R = sl[sl[:, 0] < 0][:, 0]
            if len(L) and len(R):
                brim.append(L.min() - R.max())
    out["brim_width_mm"] = round(float(max(brim)) * 1000, 1) if brim else None
    S = bones_P["AN_Sacrum"]; out["sacrum_height_mm"] = round(float(S[:, 2].max() - S[:, 2].min()) * 1000, 1)
    out["sacrum_width_mm"] = round(float(S[:, 0].max() - S[:, 0].min()) * 1000, 1)
    return out


def female_pelvis_field(P, amt=1.0, sym_x=0.0):
    """smooth displacement (atlas space, real size) taking a male pelvis toward the female form; one field for the
    hip bones and the sacrum together, so the symphysis and sacroiliac joints stay matched. amt scales all of it.
      - outlet: below the acetabula the inferior pubic / ischial rami and the tuberosities spread outward (up to
        OUT_MM a side), nothing at the symphysis (|x| < 6 mm) -> wider subpubic arch, wider intertuberous distance;
      - mid-pelvis: the ischial spines (posterior, mid height) move out SPINE_MM;
      - inlet: the brim and the sacrum's wings widen BRIM_MM a side (oval instead of heart-shaped);
      - the lower sacrum and coccyx straighten backward (COCCYX_MM at the tip), the iliac wings sit a little lower."""
    OUT_MM, SPINE_MM, BRIM_MM, COCCYX_MM, ILIUM_LOWER = 0.012, 0.006, 0.005, 0.006, 0.08
    x, y, z = P[:, 0] - sym_x, P[:, 1], P[:, 2]; ax = np.abs(x); sg = np.sign(x)
    g = ss(0.006, 0.030, ax)                     # nothing at the midline (symphysis, sacral middle)
    w_out = ss(0.805, 0.765, z)                  # below the acetabula
    bump = lambda c, h: np.exp(-((z - c) / h) ** 2)
    post = ss(-0.090, -0.060, y)
    dx = sg * g * (OUT_MM * w_out + SPINE_MM * bump(0.795, 0.02) * post + BRIM_MM * bump(0.86, 0.03))
    dy = COCCYX_MM * ss(0.86, 0.79, z) * ss(-0.07, -0.03, y)
    dz = -ILIUM_LOWER * np.maximum(0.0, z - 0.875) * ss(0.035, 0.06, ax)
    return np.stack([dx, dy, dz], 1) * amt


def drive2(kb, ctl, props, expr):
    """kb.value = expr of several control properties (variables a, b, ...), a simple expression"""
    fc = kb.driver_add("value"); dr = fc.driver; dr.type = 'SCRIPTED'
    while dr.variables:
        dr.variables.remove(dr.variables[0])
    for name, prop in zip("abcdefgh", props):
        v = dr.variables.new(); v.name = name; v.type = 'SINGLE_PROP'
        v.targets[0].id = ctl; v.targets[0].data_path = f'["{prop}"]'
    dr.expression = expr
    return fc


def weld_union(objs, M_to, tol=1e-6):
    """one welded mesh (work frame) from several objects that share seams: X, faces, and per object the union index
    of each of its vertices"""
    P_all = []; F = []; owner = []; off = 0
    for ob in objs:
        P = np.array([tuple(M_to @ (ob.matrix_world @ v.co)) for v in ob.data.vertices])
        P_all.append(P); F += [[off + i for i in p.vertices] for p in ob.data.polygons]; owner.append((ob.name, off, len(P))); off += len(P)
    P = np.concatenate(P_all); kd = kdtree.KDTree(len(P))
    for i, p in enumerate(P):
        kd.insert(Vector(p), i)
    kd.balance()
    rep = np.array([min(k_ for _, k_, _ in kd.find_range(Vector(p), tol)) for p in P])
    uniq, inv = np.unique(rep, return_inverse=True)
    X = P[uniq]; Fu = []
    for f in F:
        g = [int(inv[i]) for i in f]
        if len(set(g)) == len(g):
            Fu.append(g)
    return X, Fu, {n: inv[o:o + c] for n, o, c in owner}


def cap_holes(X, F):
    """close every boundary loop of (X, F) with a fan round its centroid (new vertices appended); returns X, F
    (triangles), the loops (vertex index lists) and the indices of the cap centres"""
    edges = {}
    for f in F:
        for a, b in zip(f, f[1:] + f[:1]):
            edges[(a, b)] = edges.get((a, b), 0) + 1
    bnd = {(a, b) for (a, b) in edges if (b, a) not in edges}
    nxt = {a: b for a, b in bnd}; loops = []; seen = set()
    for a in list(nxt):
        if a in seen:
            continue
        lp = [a]; seen.add(a); c = nxt[a]
        while c != a and c not in seen:
            lp.append(c); seen.add(c); c = nxt.get(c, a)
        loops.append(lp)
    X = list(map(tuple, X)); T = []
    for f in F:
        for k in range(1, len(f) - 1):
            T.append([f[0], f[k], f[k + 1]])
    centres = []
    for lp in loops:
        c = np.mean([X[i] for i in lp], axis=0); X.append(tuple(c)); ci = len(X) - 1; centres.append(ci)
        for a, b in zip(lp, lp[1:] + lp[:1]):
            T.append([b, a, ci])                 # the boundary edge a->b belongs to a face; the cap runs b->a
    return np.array(X), T, loops, centres
