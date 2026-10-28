"""The vagina's lumen (mucosal surface): a swept tube whose cross-sections morph, without twisting, through the
shapes of the collapsed vagina (`lumen_ring`), plus helpers for joining it to the skin and the cervix.

Cross-section ("ring") model, in the ring's plane (x lateral, y front - the anterior wall's side): N points in four
fixed blocks of N/4, one per ARM - right (+x), front (+y), left (-x), back (-y). Each block walks round its arm
from the core corner on its clockwise side to the one on its counter-clockwise side: along one wall out to the
arm's end, round a T-bar there (the lateral recesses of an "H"), back along the other wall. The walls are `g` from
the arm's midline (collapsed: g = 0.01 mm, the walls touching). Arms never vanish (>= 3 g), so the point blocks keep
their meaning in every ring and neighbouring rings can be joined index by index:
  - introitus: a sagittal slit = long front / back arms, stub side arms;
  - lower / middle vagina: an H = long side arms with T-bars, stub front / back arms;
  - upper vagina: a transverse slit (the T-bars shrink away);
  - vault: each block blends into its quarter of an ellipse round the cervix (an open space; `ellipse_ring`)."""
import math
import numpy as np

SEGS = (7, 1, 2, 4, 2, 1, 7)                       # points per walk segment (sum = points per arm)


def _poly(pts, counts):
    """points along a polyline's segments, `counts[i]` evenly spaced on segment i (its end excluded)"""
    out = []
    for (a, b), n in zip(zip(pts, pts[1:]), counts):
        for k in range(n):
            out.append(a + (b - a) * (k / n))
    return out


def arm_walk(L, V, g, segs=SEGS):
    """one arm along +a: walls at b = -g / +g from the core (a = g) to the T-bar (inner side a = L - g), the bar
    (half height V, half thickness g) round its end; returns (a, b) points, the first at the core corner (g, -g)"""
    L = max(L, 3 * g); V = max(V, 2 * g)
    P = [np.array(p, float) for p in ((g, -g), (L - g, -g), (L - g, -V), (L + g, -V), (L + g, V), (L - g, V),
                                      (L - g, g), (g, g))]
    return _poly(P, segs)


def lumen_ring(arms, g):
    """arms: [(L, V)] for right, front, left, back. Returns an (N, 2) array (x, y), counter-clockwise"""
    out = []
    for i, (L, V) in enumerate(arms):
        d = np.array([math.cos(i * math.pi / 2), math.sin(i * math.pi / 2)]); n = np.array([-d[1], d[0]])
        out += [d * a + n * b for a, b in arm_walk(L, V, g)]
    return np.array(out)


def ellipse_ring(rx, ry_front, ry_back, n_per=sum(SEGS)):
    """an ellipse-like ring (half-width rx, reaching ry_front forward and ry_back back), sampled in the same four
    blocks: block i covers angles [i pi/2 - pi/4, i pi/2 + pi/4), counter-clockwise"""
    out = []
    for i in range(4):
        for k in range(n_per):
            th = i * math.pi / 2 - math.pi / 4 + (math.pi / 2) * k / n_per
            c, s = math.cos(th), math.sin(th)
            out.append((rx * c, (ry_front if s >= 0 else ry_back) * s))
    return np.array(out)


def plane_section(X, F, origin, normal, xdir, ydir):
    """the mesh (vertices X (n,3), triangles F (m,3)) cut by a plane: its segments as (k, 2, 2) points in the
    plane's (xdir, ydir) coordinates"""
    o = np.asarray(origin, float); nrm = np.asarray(normal, float)
    d = (X - o) @ nrm; T = d[F]; segs = []
    s = np.sign(T); cross = (s.max(1) > 0) & (s.min(1) < 0)
    for f, t in zip(F[cross], T[cross]):
        pts = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (t[a] > 0) != (t[b] > 0):
                u = t[a] / (t[a] - t[b]); p = X[f[a]] + (X[f[b]] - X[f[a]]) * u - o
                pts.append((p @ np.asarray(xdir, float), p @ np.asarray(ydir, float)))
        if len(pts) == 2:
            segs.append(pts)
    return np.array(segs).reshape(-1, 2, 2)


def near(segs, rmax):
    """the section's points within rmax of the ring's centre (a plane also cuts far parts of a curved organ)"""
    P = segs.reshape(-1, 2)
    P = P[np.linalg.norm(P, axis=1) < rmax]
    return P if len(P) else None


def extents(P, centre=(0.0, 0.0)):
    """reach of section points P from `centre` along +x, -x, +y, -y (only points near each axis count)"""
    if P is None or len(P) == 0:
        return None
    P = P - np.asarray(centre)
    out = []
    for ax, sg in ((0, 1), (0, -1), (1, 1), (1, -1)):
        other = P[:, 1 - ax]; along = P[:, ax] * sg
        band = np.abs(other) < 0.0006
        out.append(float(along[band & (along > 0)].max()) if np.any(band & (along > 0)) else 0.0)
    return out                                      # [+x, -x, +y, -y]


def sweep(centre_line, frames, rings2d):
    """vertices (rings x N, 3) of a tube: ring j = centre_line[j] + x * frames[j][0] + y * frames[j][1]"""
    out = []
    for c, (ex, ey), P in zip(centre_line, frames, rings2d):
        out += [c + ex * x + ey * y for x, y in P]
    return np.array(out)


def tangents(C):
    C = np.asarray(C, float); out = []
    for j in range(len(C)):
        t = C[min(j + 1, len(C) - 1)] - C[max(j - 1, 0)]; out.append(t / np.linalg.norm(t))
    return np.array(out)


def frames_along(C, xref, normals=None):
    """(x, y) frames of ring planes along a polyline C: x = xref made square to the plane's normal (default: the
    tangent), y = x cross normal"""
    C = np.asarray(C, float); out = []; TT = tangents(C) if normals is None else np.asarray(normals, float)
    for j in range(len(C)):
        t = TT[j] / np.linalg.norm(TT[j])
        x = np.asarray(xref, float) - t * (np.asarray(xref, float) @ t); x /= np.linalg.norm(x)
        y = np.cross(x, t); out.append((x, y / np.linalg.norm(y)))
    return out


def keyed_snapshot(ob):
    """an object's shape keys as displacements from Basis, with a BVH of its Basis surface (local coords), so the
    keys can be rebuilt after the mesh is re-cut (keys_restore)"""
    import bmesh
    from mathutils.bvhtree import BVHTree
    me = ob.data
    if not me.shape_keys:
        return None
    X = np.array([tuple(v.co) for v in me.vertices]); F = []
    for p in me.polygons:
        vs = list(p.vertices)
        for k in range(1, len(vs) - 1):
            F.append((vs[0], vs[k], vs[k + 1]))
    F = np.array(F)
    D = {kb.name: np.array([tuple(d.co) for d in kb.data]) - X for kb in me.shape_keys.key_blocks}
    tree = BVHTree.FromPolygons([tuple(p) for p in X], [tuple(f) for f in F])
    return X, F, D, tree


def keys_restore(ob, snap):
    """after the mesh changed: every vertex's key position = its position + the displacement interpolated
    (barycentric) over the old Basis triangle it lies on"""
    from mathutils import Vector
    from mathutils.geometry import barycentric_transform
    if snap is None:
        return 0
    X, F, D, tree = snap; me = ob.data; n_new = 0
    P = [v.co.copy() for v in me.vertices]; W = []
    for p in P:
        q, _, fi, dist = tree.find_nearest(p)
        a, b, c = (Vector(X[i]) for i in F[fi])
        # barycentric weights of q in triangle abc
        v0, v1, v2 = b - a, c - a, q - a
        d00, d01, d11, d20, d21 = v0.dot(v0), v0.dot(v1), v1.dot(v1), v2.dot(v0), v2.dot(v1)
        den = d00 * d11 - d01 * d01 or 1e-30
        wb = (d11 * d20 - d01 * d21) / den; wc = (d00 * d21 - d01 * d20) / den
        W.append((F[fi], (1 - wb - wc, wb, wc)))
        n_new += dist > 2e-6
    for kb in me.shape_keys.key_blocks:
        if kb.name not in D:
            continue
        Dk = D[kb.name]
        co = [tuple(Vector(p) + Vector(sum(w * Dk[i] for i, w in zip(tri, ws)))) for p, (tri, ws) in zip(P, W)]
        kb.data.foreach_set("co", np.array(co, np.float32).ravel())
    me.update()
    return n_new


def reach_y_at(P, x0, centre=(0.0, 0.0), half=0.0006):
    """the section's reach up (+y) and down (-y) from the centre's height, among points within `half` of x = x0"""
    if P is None or len(P) == 0:
        return None
    Q = P - np.asarray(centre); Q = Q[np.abs(Q[:, 0] - x0) < half]
    if len(Q) == 0 or Q[:, 1].max() <= 0 or Q[:, 1].min() >= 0:
        return None
    return float(Q[:, 1].max()), float(-Q[:, 1].min())
