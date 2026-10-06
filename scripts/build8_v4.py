"""v4: like v3, but the inner pucker is an even, fine triangle mesh (no grid lines converging at the slit tips).
Detail pass on the anus of Hips.blend (v8: squeezed pucker, randomised 3D creases):
skin crease smoothing, cleft walls pressed in around the anus, high-detail pucker rebuilt as an AP slit with
creases fanning out of it, laid out by arc length over the skin so it rides up the cheek walls,
pigment as an editable shader mask, internal fit raised 1 cm, angled anal canal + lower rectal ampulla,
sphincter re-seated."""
import bpy, bmesh, math, sys, json, os, random
import numpy as np
from mathutils import Vector, Matrix

OUT = sys.argv[-1]          # folder for report / optional save path
SAVE_AS = sys.argv[-2]      # "INPLACE" or a .blend path
REF = os.environ.get("ANATOMY_REF", r"C:\Users\Parker\Anatomy Project\anatomy_ref.blend")
report = {}


def ss(e0, e1, x):
    t = min(max((x - e0) / (e1 - e0), 0.0), 1.0)
    return t * t * (3 - 2 * t)


def r4(v):
    return [round(float(x), 4) for x in v]


hips = bpy.data.objects["Hips"]; me = hips.data
assert hips.matrix_world == Matrix.Identity(4)

# ======================= 1. local frame at the anus =======================
O = Vector((0.0, 0.0565, 0.790))                 # on the cleft floor, centre of the new anus
bm = bmesh.new(); bm.from_mesh(me)
pig_layer = bm.verts.layers.float.get("anus_pigment") or bm.verts.layers.float.new("anus_pigment")  # before taking refs
cre_layer = bm.verts.layers.float.get("anus_crease") or bm.verts.layers.float.new("anus_crease")
fs_layer = bm.verts.layers.float.get("anus_fold_shade") or bm.verts.layers.float.new("anus_fold_shade")
bn_layer = bm.verts.layers.float_vector.get("anus_base_normal") or bm.verts.layers.float_vector.new("anus_base_normal")
bm.verts.ensure_lookup_table(); bm.faces.ensure_lookup_table()
bm.normal_update()
ann = np.array([tuple(v.co) for v in bm.verts if 0.010 < (v.co - O).length < 0.016])
c = ann.mean(0); _, V = np.linalg.eigh((ann - c).T @ (ann - c))
n = Vector(V[:, 0]); n.x = 0.0; n.normalize()
if n.z > 0:
    n = -n                                        # outward: down / back
u = (Vector((0, 1, 0)) - n * n.y).normalized()    # along the cleft, towards the back
bv = n.cross(u)                                   # lateral
if bv.x < 0:
    bv = -bv


def to_local(p):
    d = p - O
    return d.dot(u), d.dot(bv), d.dot(n)


# ======================= 2. Taubin-smooth the cleft crease =======================
def r_eff(p):
    a, b, h = to_local(p)
    return math.sqrt((a * (2.2 if a < 0 else 1.0)) ** 2 + b * b + h * h)


region = [v for v in bm.verts if not v.is_boundary and r_eff(v.co) < 0.035]
wt = {v: 1.0 - ss(0.022, 0.035, r_eff(v.co)) for v in region}
for _ in range(20):
    for f in (0.5, -0.53):
        new = {}
        for v in region:
            nb = [e.other_vert(v).co for e in v.link_edges]
            avg = sum(nb, Vector()) / len(nb)
            new[v] = v.co + (avg - v.co) * (f * wt[v])
        for v, p in new.items():
            v.co = p
for v in bm.verts:                                # keep the mesh exactly mirror-symmetric on the midline
    if abs(v.co.x) < 1e-5:
        v.co.x = 0.0
bm.normal_update()

# ======================= 2b. cheeks pressed together: narrow the cleft around the anus =======================
SQ = 0.35                        # max fraction of lateral distance removed (0 = cleft untouched)


def squeeze_w(p):
    a, b, h = (x * 1000 for x in to_local(p))
    fa = 1 - (ss(15, 40, a) if a > 0 else ss(8, 25, -a))     # towards the coccyx the cleft stays closed longer
    fb = 1 - ss(10, 22, abs(b))
    fh = ss(-25, -10, h) * (1 - ss(12, 22, h))                 # cleft floor and walls, not the outer cheeks
    return fa * fb * fh


for v in bm.verts:               # x' = x * (1 - SQ*w): monotonic in x, so nothing folds or crosses the midline
    if not v.is_boundary:
        w = squeeze_w(v.co)
        if w > 0:
            v.co.x *= 1 - SQ * w
bm.normal_update()
report["cleft_squeeze"] = SQ

# ======================= 2c. refine the cleft =======================
# The body is a low-poly game mesh: across its big triangles the toon shader's hard light/shadow cut comes out as a
# zigzag along the cheek walls, and the strip joining the dense pucker to it gets long thin triangles (white slivers).
# Split every cleft edge into 4 and put the new vertices on the Phong surface of the original triangles (curved by the
# vertex normals): old vertices stay put, so the silhouette is unchanged.
from mathutils.bvhtree import BVHTree
from mathutils.interpolate import poly_3d_calc
REFINE_CUTS = 3
A_FRONT, A_BACK = 45, 80         # mm along the cleft: towards the perineum / up towards the coccyx

# the cleft floor curves away from the local frame further up (19 mm off at 45 mm, 56 mm at 75 mm), so the region
# follows a sampled midline floor profile instead of a fixed height band
_t = BVHTree.FromBMesh(bm); _fa, _fh = [], []
for a_mm in range(-A_FRONT - 10, A_BACK + 11, 5):
    P = O + u * (a_mm / 1000)
    hit = _t.ray_cast(P - n * 0.08, n, 0.2)[0]
    if hit is not None:
        _fa.append(a_mm); _fh.append((hit - O).dot(n) * 1000)


def cleft_floor_h(a_mm):
    return float(np.interp(a_mm, _fa, _fh))


def cleft_coords(p):
    a, b, h = (x * 1000 for x in to_local(p))
    return a, b, h - cleft_floor_h(a)


def in_refine(p):
    a, b, hr = cleft_coords(p)
    return -A_FRONT < a < A_BACK and abs(b) < 28 and -15 < hr < 30


old = bm.copy(); bmesh.ops.triangulate(old, faces=old.faces[:]); old.faces.ensure_lookup_table(); old.normal_update()
old_tree = BVHTree.FromBMesh(old)


def phong_old(p, alpha=0.75):
    loc, _, fi, _ = old_tree.find_nearest(p)
    vs = [l.vert for l in old.faces[fi].loops]
    ws = poly_3d_calc([v.co for v in vs], loc)
    q = Vector()
    for w_, v in zip(ws, vs):
        q += (loc - v.normal * (loc - v.co).dot(v.normal)) * w_
    return loc.lerp(q, alpha)


before = set(bm.verts)
ref_edges = [e for e in bm.edges if not e.is_boundary and all(in_refine(v.co) for v in e.verts)]
bmesh.ops.subdivide_edges(bm, edges=ref_edges, cuts=REFINE_CUTS, use_grid_fill=True)
refined_new = [v for v in bm.verts if v not in before]
for v in refined_new:
    v.co = phong_old(v.co)
    if abs(v.co.x) < 2e-5:
        v.co.x = 0.0
old.free()
bm.verts.ensure_lookup_table(); bm.faces.ensure_lookup_table(); bm.normal_update()
report["refined_edges"] = len(ref_edges); report["refine_new_verts"] = len(refined_new)

# ======================= 3. region to rebuild + smooth base surface (old pucker faired away) =======================
# The hole is cut in the same (a, s) surface coordinates the pucker is laid out in (a along the cleft, s = arc
# length across it), a fixed margin outside the outer ring, so the bridge strip is narrow everywhere - a projected
# ellipse would reach far up the steep cheek walls and leave tall sliver triangles there.
H0 = 0.007                       # ray origin height above the cleft floor (in the air of the cleft)
D_AP, D_LAT = 8.8 + 1.2, 0.85 * math.sqrt(8.8 ** 2 - 3.0 ** 2) + 1.2     # outer ring semi-axes + margin, mm
_cur = BVHTree.FromBMesh(bm)


def surf_coords(p, tree=_cur):
    """(a, s) in mm of a point on the skin: s = arc length across the cleft from the midline (signed by side)"""
    a = to_local(p)[0]
    C = O + u * a + n * H0; q = p - C
    x_, y_ = q.dot(bv), -q.dot(n)
    sg = 1.0 if x_ >= 0 else -1.0; phi_v = math.atan2(abs(x_), y_)
    step = math.radians(0.5); phi = 0.0; prev = None; acc = 0.0
    while True:
        last = phi >= phi_v
        ph = min(phi, phi_v)
        hit = tree.ray_cast(C, -n * math.cos(ph) + bv * (sg * math.sin(ph)), 0.05)[0]
        if hit is None:
            return a * 1000, sg * 1e9
        if prev is not None:
            acc += (hit - prev).length
        if last:
            return a * 1000, sg * acc * 1000
        prev = hit; phi += step


def in_D(p):
    a, b, h = to_local(p)
    if abs(a) > 0.0115 or abs(b) > 0.012 or not -0.008 < h < 0.02:
        return False
    a_mm, s_mm = surf_coords(p)
    return (a_mm / D_AP) ** 2 + (s_mm / D_LAT) ** 2 < 1.0


seed = min(bm.verts, key=lambda v: (v.co - O).length)
D = set(); stack = [seed]
while stack:
    v = stack.pop()
    if v in D or v.is_boundary or not in_D(v.co):
        continue
    D.add(v); stack += [e.other_vert(v) for e in v.link_edges]

# base = copy of the skin where D is replaced by a biharmonic (C1-smooth) fill of the cleft
base = bm.copy(); base.verts.ensure_lookup_table(); base.faces.ensure_lookup_table()
Didx = sorted(v.index for v in D); col = {i: k for k, i in enumerate(Didx)}
rows_idx = set(Didx) | {e.other_vert(bm.verts[i]).index for i in Didx for e in bm.verts[i].link_edges}
Am = []; Bm = []
for i in rows_idx:
    vi = base.verts[i]; nb = [e.other_vert(vi) for e in vi.link_edges]
    row = np.zeros(len(Didx)); rhs = np.zeros(3)
    if i in col: row[col[i]] += 1.0
    else: rhs -= np.array(vi.co)
    for w_ in nb:
        if w_.index in col: row[col[w_.index]] -= 1.0 / len(nb)
        else: rhs += np.array(w_.co) / len(nb)
    Am.append(row); Bm.append(rhs)
# minimise |Lx|^2 over the rows -> stack Laplacian rows and solve the normal equations of L^T L (biharmonic)
Am = np.array(Am); Bm = np.array(Bm)
Xs, *_ = np.linalg.lstsq(Am, Bm, rcond=None)
for i in Didx:
    base.verts[i].co = Vector(Xs[col[i]])
for v in base.verts:
    if abs(v.co.x) < 2e-5: v.co.x = 0.0
base.normal_update()
# the base is still the coarse body mesh, so draping onto its flat triangles would copy their facets into the dense
# pucker. Phong tessellation curves each triangle with the vertex (shading) normals: it passes through the original
# vertices, so it meets the surviving skin at the hole rim, and matches how that skin looks when smooth-shaded.
btri = base.copy(); bmesh.ops.triangulate(btri, faces=btri.faces[:]); btri.faces.ensure_lookup_table()
btri.normal_update()
base_tree = BVHTree.FromBMesh(btri)


def phong(p, fi, alpha=0.75):
    vs = [l.vert for l in btri.faces[fi].loops]
    ws = poly_3d_calc([v.co for v in vs], p)
    q = Vector()
    for w_, v in zip(ws, vs):
        q += (p - v.normal * (p - v.co).dot(v.normal)) * w_
    return p.lerp(q, alpha)


def cast(C, dr):
    hit, nrm, fi, _ = base_tree.ray_cast(C, dr, 0.05)
    return (None, nrm) if hit is None else (phong(hit, fi), nrm)

# side-aware UV sampler (the UV layout has a seam on the midline)
side_trees = {}
for s in (1, -1):
    faces = [f for f in base.faces if s * f.calc_center_median().x >= -1e-6]
    side_trees[s] = (BVHTree.FromPolygons([v.co for v in base.verts], [[v.index for v in f.verts] for f in faces]), faces)
base_uv = {lay.name: base.loops.layers.uv[lay.name] for lay in bm.loops.layers.uv.values()}


def sample_uv(p, side, lname):
    tree_, faces = side_trees[side]
    _, _, idx, _ = tree_.find_nearest(p)
    f = faces[idx]; ws = poly_3d_calc([v.co for v in f.verts], p)
    uv = Vector((0.0, 0.0))
    for w_, l in zip(ws, f.loops):
        uv += l[base_uv[lname]].uv * w_
    return uv


# remember surviving UVs of the rim (per side) before deleting
uv_layers = list(bm.loops.layers.uv.values())
bmesh.ops.delete(bm, geom=list(D), context='VERTS')
bm.verts.ensure_lookup_table()
hole = [v for v in bm.verts if v.is_boundary and (v.co - O).length < 0.03]
start = hole[0]; L = [start]; prev = None; cur = start
while True:
    nxt = [e.other_vert(cur) for e in cur.link_edges if e.is_boundary and e.other_vert(cur) != prev][0]
    if nxt == start:
        break
    prev, cur = cur, nxt; L.append(cur)
report["deleted_verts"] = len(D); report["hole_loop"] = len(L)

# ======================= 4/5. new pucker: AP slit, confocal rings laid out over the skin, creases fanning from the slit ===
# (a, s) surface coordinates: a = distance along the cleft (u), s = arc length over the base surface across the cleft,
# measured in the cross-section plane at that a. Laying the rings out in arc length lets the lateral sides climb the
# cheek walls instead of being projected onto the cleft floor.
SLIT = 1.8                       # mm, half-length of the closed slit (short, so the creases converge like a star)
K_LAT = 0.85                     # lateral (arc length) / AP for the confocal rings -> AP-elongated
RINGS = [1.82, 7.2, 7.6, 8.2, 8.8]   # v4: lip ring + outer rings only; inside 7.2 mm is an even triangle mesh
# Creases: each one sits exactly on a straight radial grid line (a "spoke"), so it is one clean line of constant
# width running from the slit's lip outwards, and it ends in the centre crease. Spokes every 10 deg exist in every
# ring out to the 36-vertex rings, so no crease breaks at a ring step. Randomness = which spokes carry a crease
# (picked separately for each side, never mirrored), how long each runs and how deep it is. No wobble, no smoothing.
CREASE_SEED = 11
_rng = random.Random(CREASE_SEED)
STEP = 2 * math.pi / 36
# front (perineum) half, approved as is: a crease every 20 deg fanning from the front tip (deg, depth, runs to r)
_FRONT = [(120, 0.97, 5.5), (140, 0.93, 6.4), (160, 0.92, 6.6), (180, 0.93, 6.5),
          (-160, 0.91, 5.6), (-140, 0.85, 6.6), (-120, 0.88, 6.3)]
# back half and sides at the same 20 deg density; the two sides are offset by 10 deg, so the back tip gets an
# uneven V of creases instead of one running straight up the cleft, and the sides are not mirror images
_BACK = [20, 40, 60, 80, 100, -30, -50, -70, -90]
CREASES = [(math.radians(d) % (2 * math.pi), w, r1) for d, w, r1 in _FRONT]
CREASES += [(math.radians(d) % (2 * math.pi), _rng.uniform(0.85, 1.0), _rng.uniform(6.2, 7.0)) for d in _BACK]
CREASES.sort()
_ang = [c[0] for c in CREASES]
CREASE_HALF = []                                                             # half the gap to the nearer neighbour
for i, t0 in enumerate(_ang):
    g = min((t0 - _ang[i - 1]) % (2 * math.pi), (_ang[(i + 1) % len(_ang)] - t0) % (2 * math.pi))
    CREASE_HALF.append(min(max(g / 2, math.radians(5)), math.radians(14)))
FUNNEL = 0.3; R_FOLD = 7.4; FOLD_DEPTH = 0.45; SLIT_DEPTH = 0.5; RIDGE = 0.0
# pigment outline: a few low harmonics with random phases (+-~4 %) so it isn't a perfect, mirrored oval
PIG_WOBBLE = [(k, _rng.uniform(0.008, 0.016), _rng.uniform(0, 2 * math.pi)) for k in (2, 3, 5)]


def pig_r(r, th):
    return r * (1 + sum(a * math.sin(k * th + ph) for k, a, ph in PIG_WOBBLE))


def crease(theta, r):
    """0..1 groove at grid point (r, theta): a narrow V centred on each crease's spoke, flat pad elsewhere"""
    best = 0.0
    for (t0, w, r1), hg in zip(CREASES, CREASE_HALF):
        d = abs((theta - t0 + math.pi) % (2 * math.pi) - math.pi)
        if d < hg:
            best = max(best, w * (1 - d / hg) ** 2.5 * (1 - ss(r1 - 1.2, r1, r)))
    # the lip edge itself stays level so both lips meet (closed); full depth is reached 0.08 mm out, at the next ring
    return best * ss(RINGS[0], RINGS[0] + 0.08, r)


def fold_amp(r):                 # full depth from the lip outwards; each crease's own length fades it
    return FOLD_DEPTH


def paint_amp(r):                # painted crease lines run a little further than the relief
    return ss(RINGS[0], RINGS[0] + 0.4, r) * (1 - ss(5.4, 6.8, r))


def surf_point(a_mm, s_mm):
    """point on the (Phong-curved) base at AP coordinate a and signed lateral arc length s (mm)"""
    P0 = O + u * (a_mm / 1000); C = P0 + n * H0
    sg = 1.0 if s_mm >= 0 else -1.0; target = abs(s_mm) / 1000
    step = math.radians(0.25); phi = 0.0; prev = None; acc = 0.0
    while True:
        dr = -n * math.cos(phi) + bv * (sg * math.sin(phi))
        hit, nrm = cast(C, dr)
        assert hit is not None and nrm.dot(dr) < 0, (a_mm, s_mm, math.degrees(phi))
        if prev is None:
            if target == 0.0:
                hit.x = 0.0
                return hit
        else:
            seg = (hit - prev).length
            if acc + seg >= target:
                ph = phi - step * (1 - (target - acc) / seg)
                return cast(C, -n * math.cos(ph) + bv * (sg * math.sin(ph)))[0]
            acc += seg
        prev = hit; phi += step
        assert phi < math.radians(150), (a_mm, s_mm)


K_SPOKE = 0.97                   # lateral growth per mm of AP growth (outer ring ~8.8 x 6.8 mm)


def ring_xy(r, th):
    """ring of AP semi-axis r around the slit. Lateral offset grows linearly with r, so every grid spoke (constant
    th) is a straight line from its point on the slit outwards - the creases sit on these lines. The innermost ring
    is ~0.06 mm from the slit: the lips touch (closed passage)."""
    return r * math.cos(th), K_SPOKE * (r - SLIT) * math.sin(th)


S_OUT = 18 if len(L) < 27 else (36 if len(L) < 54 else 72)
S_OUT = max(S_OUT, 72)
report["outer_ring_verts"] = S_OUT
MESH_H = 0.07                    # mm, edge length of the even inner mesh


def layout_rt(x, y):
    """inverse of ring_xy: (r, th) of a layout point"""
    lo, hi = max(SLIT + 1e-7, abs(x) + 1e-9), 40.0
    for _ in range(60):
        mid = (lo + hi) / 2
        f = (x / mid) ** 2 + (y / (K_SPOKE * (mid - SLIT))) ** 2 - 1
        lo, hi = (mid, hi) if f > 0 else (lo, mid)
    r = (lo + hi) / 2
    return r, math.atan2(y / (K_SPOKE * (r - SLIT)), x / r) % (2 * math.pi)


def surf_column(x, s_list):
    """surface points at AP coordinate x for many signed lateral arc lengths: one march per side, then one exact cast
    per point"""
    out = {}
    P0 = O + u * (x / 1000); C = P0 + n * H0
    for sg in (1.0, -1.0):
        mine = [t_ for t_ in s_list if (t_ >= 0) == (sg > 0)]
        if not mine:
            continue
        targets = sorted(abs(t_) / 1000 for t_ in mine)
        step = math.radians(0.25); phi = 0.0; prev = None; acc = 0.0; table = [(0.0, 0.0)]
        while acc < targets[-1] + 1e-6:
            dr = -n * math.cos(phi) + bv * (sg * math.sin(phi))
            hit, nrm = cast(C, dr)
            assert hit is not None and nrm.dot(dr) < 0, (x, math.degrees(phi))
            if prev is not None:
                acc += (hit - prev).length; table.append((acc, phi))
            prev = hit; phi += step
            assert phi < math.radians(150), x
        ta = np.array([t_[0] for t_ in table]); tp = np.array([t_[1] for t_ in table])
        for t_ in mine:
            ph = float(np.interp(abs(t_) / 1000, ta, tp))
            out[t_] = cast(C, -n * math.cos(ph) + bv * (sg * math.sin(ph)))[0]
    return out


rings = []; new_faces = []; info = {}; layout = {}
# --- lip ring (two lips + tips), vertices evenly spaced along the slit, and the zipper midline
N0 = 2 * int(round(2 * RINGS[0] / MESH_H)); H = N0 // 2
LIP_X = [RINGS[0] * (1 - 2 * k / H) for k in range(H + 1)]
LIP_X += [LIP_X[N0 - k] for k in range(H + 1, N0)]   # lower lip mirrors the upper one
R0 = []
for k in range(N0):
    x = LIP_X[k]; sg = 1.0 if k <= H else -1.0
    y = sg * K_SPOKE * (RINGS[0] - SLIT) * math.sqrt(max(1 - (x / RINGS[0]) ** 2, 0.0))
    if k in (0, H):
        y = 0.0
    v = bm.verts.new(surf_point(x, y)); r_, th_ = layout_rt(x, y)
    info[v] = (RINGS[0], th_); layout[v] = (x, y); R0.append(v)
rings.append(R0)
M = {0: R0[0], H: R0[H]}; slit_verts = {}
for k in range(1, H):
    x = LIP_X[k]
    v = bm.verts.new(surf_point(x, 0.0)); info[v] = (0.0, math.acos(max(-1, min(1, x / RINGS[0]))))
    slit_verts[v] = info[v][1]; M[k] = v; layout[v] = (x, 0.0)
# --- outer polar rings, stepping down 2:1 to the join with the skin
OUTER_S = [576, 288, 144, S_OUT]
assert all(a_ in (b_, 2 * b_) for a_, b_ in zip(OUTER_S, OUTER_S[1:])), OUTER_S
for r, S in zip(RINGS[1:], OUTER_S):
    ring = []
    for k in range(S):
        th = 2 * math.pi * k / S
        x, y = ring_xy(r, th)
        v = bm.verts.new(surf_point(x, y)); info[v] = (r, th); layout[v] = (x, y); ring.append(v)
    rings.append(ring)
# --- even triangle lattice between the lips and the 7.2 mm ring
cols = {}
dx = MESH_H * math.sqrt(3) / 2
for i in range(-int(9.0 / dx), int(9.0 / dx) + 1):
    x = i * dx
    for j in range(-int(7.0 / MESH_H), int(7.0 / MESH_H) + 1):
        y = (j + 0.5 * (i % 2)) * MESH_H
        r_, th_ = layout_rt(x, y)
        if r_ > RINGS[1] - 0.6 * MESH_H:
            continue
        if r_ < RINGS[0] + 0.6 * MESH_H:                       # keep clear of the lips
            continue
        cols.setdefault(x, []).append((y, r_, th_))
lat = []
for x, pts in cols.items():
    sp = surf_column(x, [p_[0] for p_ in pts])
    for y, r_, th_ in pts:
        v = bm.verts.new(sp[y]); info[v] = (r_, th_); layout[v] = (x, y); lat.append(v)
report["inner_mesh_verts"] = len(lat)
# --- triangulate the ring band between the lips and the 7.2 mm ring (constrained Delaunay, lips = hole)
from mathutils.geometry import delaunay_2d_cdt
from mathutils import Vector as _V2
band = R0 + rings[1] + lat
idx = {v: i for i, v in enumerate(band)}
co2 = [_V2(layout[v]) for v in band]
outer_poly = [idx[v] for v in rings[1]]; lip_poly = [idx[v] for v in R0]
cdt = delaunay_2d_cdt(co2, [], [outer_poly, lip_poly], 2, 1e-7, True)
vo, eo, fo, orig_v = cdt[0], cdt[1], cdt[2], cdt[3]
assert len(vo) == len(band), ("cdt added vertices", len(vo), len(band))
omap = [band[ov[0]] for ov in orig_v]
for f in fo:
    vs = [omap[i] for i in f]
    if len(set(vs)) == 3 and not all(i < N0 for i in (idx[w] for w in vs)):   # skip any face between the lips
        new_faces.append(bm.faces.new(vs))
report["inner_mesh_faces"] = len(new_faces)


def add_face(vs):
    vs = [x for i, x in enumerate(vs) if x is not vs[i - 1]]
    new_faces.append(bm.faces.new(vs))


for k in range(H):
    add_face((M[k], R0[k], R0[k + 1], M[k + 1]))                       # lip with s > 0
    add_face((M[k + 1], R0[N0 - 1 - k], R0[(N0 - k) % N0], M[k]))      # lip with s < 0
for ri in range(1, len(rings) - 1):
    A_, B_ = rings[ri], rings[ri + 1]
    if len(A_) == len(B_):
        S = len(A_)
        for k in range(S):
            new_faces.append(bm.faces.new((A_[k], B_[k], B_[(k + 1) % S], A_[(k + 1) % S])))
    else:
        nB = len(B_); nA = len(A_)
        for k in range(nB):
            new_faces.append(bm.faces.new((A_[2 * k], B_[k], A_[2 * k + 1])))
            new_faces.append(bm.faces.new((A_[2 * k + 1], B_[k], B_[(k + 1) % nB], A_[(2 * k + 2) % nA])))
bmesh.ops.recalc_face_normals(bm, faces=new_faces)
bm.normal_update()
if sum(f.normal.dot(n) for f in new_faces) < 0:
    for f in new_faces:
        f.normal_flip()
bm.normal_update()
# relax the inner surface into a smooth membrane (normal-only steps, numpy), fading out over the outer rings
_vl = list(info); _vi = {v: i for i, v in enumerate(_vl)}
P = np.array([tuple(v.co) for v in _vl])
W = np.array([1.0 if v in slit_verts else 1 - ss(6.4, RINGS[-1], info[v][0]) for v in _vl])
E_ = np.array([(_vi[e.verts[0]], _vi[e.verts[1]]) for f in new_faces for e in f.edges
               if e.verts[0] in _vi and e.verts[1] in _vi])
E_ = np.unique(np.sort(E_, axis=1), axis=0)
T_ = np.array([(_vi[f.verts[0]], _vi[f.verts[i]], _vi[f.verts[i + 1]]) for f in new_faces for i in range(1, len(f.verts) - 1)])
deg = np.bincount(E_.ravel(), minlength=len(_vl)).astype(float)
N_ROUND = 800                   # numpy: relax fully into a smooth membrane (fills the old dip, no lift needed)
for _ in range(N_ROUND):
    fn = np.cross(P[T_[:, 1]] - P[T_[:, 0]], P[T_[:, 2]] - P[T_[:, 0]])
    vn = np.zeros_like(P)
    for c_ in range(3):
        np.add.at(vn, T_[:, c_], fn)
    vn /= np.maximum(np.linalg.norm(vn, axis=1, keepdims=True), 1e-20)
    acc_ = np.zeros_like(P)
    np.add.at(acc_, E_[:, 0], P[E_[:, 1]]); np.add.at(acc_, E_[:, 1], P[E_[:, 0]])
    d_ = acc_ / np.maximum(deg, 1)[:, None] - P
    P += vn * (0.9 * W * np.einsum('ij,ij->i', d_, vn))[:, None]
for v, p_ in zip(_vl, P):
    v.co = Vector(p_)
for v in info:
    if abs(v.co.x) < 2e-5:
        v.co.x = 0.0
bm.normal_update()
base_pos = {v: v.co.copy() for v in info}


def base_normal(p):
    loc, _, fi, _ = base_tree.find_nearest(p)
    vs = [l.vert for l in btri.faces[fi].loops]
    ws = poly_3d_calc([v.co for v in vs], loc)
    return sum((v.normal * w_ for w_, v in zip(ws, vs)), Vector()).normalized()


bnorm = {v: base_normal(base_pos[v]) for v in info}
# The filled base still sags towards the original model's deep anus dimple (~1.3 mm low in the middle). Lift it so
# the centre sits level with the skin at the front / back ends of the pucker; fades to nothing at the outer ring.
_ends = [v for r_, ring in zip(RINGS, rings) if r_ == RINGS[-1] for v in ring
         if abs(math.sin(info[v][1])) < 0.2]
_end_h = sum(to_local(base_pos[v])[2] for v in _ends) / len(_ends)
_ctr_h = to_local(base_pos[M[H // 2]])[2]
LIFT = max(0.0, _end_h - _ctr_h)
report["base_lift_mm"] = round(LIFT * 1000, 2)
for v, (r, th) in info.items():
    pass
# v3: the leftover dip is lopsided (deeper towards the back), so a lift that falls off evenly in every direction left
# a U-shaped ridge around the back half. Lift by the dip's own profile instead: along the cleft, how far the
# midline sits below the straight line between the pucker's front and back ends; faded out sideways.
_mid = sorted((ring_xy(r, th)[0], to_local(base_pos[v])[2]) for v, (r, th) in info.items()
              if abs(ring_xy(r, th)[1]) < 0.25 and r > 0)
_mx = np.array([m[0] for m in _mid]); _mh = np.array([m[1] for m in _mid]) * 1000
_ref = np.interp(_mx, [_mx[0], _mx[-1]], [_mh[0], _mh[-1]])
_def = np.maximum(_ref - _mh, 0.0)
_def = np.convolve(np.pad(_def, 6, mode="edge"), np.ones(13) / 13, mode="same")[6:-6]   # smooth the profile
report["base_lift_mm"] = round(float(_def.max()), 2)
for v, (r, th) in info.items():
    x_, y_ = ring_xy(max(r, RINGS[0]), th)
    lift_mm = float(np.interp(x_, _mx, _def)) * (1 - ss(0.0, 2.5, abs(y_))) * (1 - ss(RINGS[-1] - 1.6, RINGS[-1] - 0.4, r))
    base_pos[v] = base_pos[v] + bnorm[v] * (lift_mm / 1000)
# relief along the (smooth) base normal: shallow funnel, creases fanning from the slit, closed slit groove
# Profile across the slit (from your sketch): the skin rises gently towards the slit, peaks right at its lips, then
# drops into a narrow V - no wide valley. The rise goes straight out of the cleft (along n): on the steep cheek walls
# the surface normal points sideways and would pinch the pucker.
PUFF = 0.35                      # mm, height of the rise at the slit's lips
SLIT_V = 0.7                     # mm, depth of the V below the lips
LIP_ROLL = 0.9                   # mm the lips curl down into the open fissure
relief_of = {}; puff_of = {}
for v, (r, th) in info.items():
    puff_of[v] = n * (PUFF / 1000 * (1 - ss(RINGS[0], RINGS[0] + 2.8, max(r, RINGS[0]))))   # gentle climb over the last ~2.8 mm, peak at the lips
    if v in slit_verts:
        sw = math.sin(th)
        relief = -SLIT_V * sw ** 0.6
        v[pig_layer] = 1.0; v[cre_layer] = 0.9 * sw ** 0.5
    else:
        _x, _y = ring_xy(r, th)                    # creases converge at the slit tips: calm them there
        _tip = min(math.hypot(_x - SLIT, _y), math.hypot(_x + SLIT, _y))
        relief = -fold_amp(r) * crease(th, r)  - LIP_ROLL * math.exp(-abs(ring_xy(r, th)[1]) / 0.12) * (1 - ss(SLIT - 0.7, SLIT, abs(ring_xy(r, th)[0])))   # curl in real sideways mm, tapering out before the tips
        v[pig_layer] = 1 - ss(5.6, 7.8, pig_r(r, th))
        v[cre_layer] = crease(th, r)
    relief_of[v] = relief / 1000
# round off the creases: near the slit (and at its ends) a crease's angular width covers very little skin, so the
# raw relief makes knife-edge grooves there. A few averaging passes over the grid give every groove a natural
# minimum width without moving the pattern.
RELIEF_SMOOTH = 2                # v4: on the fine mesh this only rounds the groove bottoms (~0.1 mm)
_lips = set()
for _ in range(RELIEF_SMOOTH):
    relief_of = {v: d if (v in slit_verts or v in _lips) else           # the closed slit keeps its groove
                 0.5 * d + 0.5 * sum(relief_of.get(e.other_vert(v), d) for e in v.link_edges) / len(v.link_edges)
                 for v, d in relief_of.items()}
for v in info:
    v.co = base_pos[v] + bnorm[v] * relief_of[v] + puff_of.get(v, Vector())
for v, (r, th) in info.items():
    v[fs_layer] = 1 - ss(6.5, 8.6, r)              # where the fold shading acts (0 at the outer ring)
bm.normal_update()
for v in bm.verts:                                 # smooth (un-creased) normal, for the fold shading's reference light
    v[bn_layer] = bnorm.get(v, v.normal)

outer = rings[-1]
outer_edges = [bm.edges.get((outer[k], outer[(k + 1) % len(outer)])) for k in range(len(outer))]
L_edges = [bm.edges.get((L[k], L[(k + 1) % len(L)])) for k in range(len(L))]
bridge_faces = bmesh.ops.bridge_loops(bm, edges=outer_edges + L_edges)["faces"]
for f in bridge_faces:                          # wind consistently with the surviving skin
    for e in f.edges:
        others = [g for g in e.link_faces if g is not f and g not in bridge_faces]
        if others:
            g = others[0]
            fa = [l.vert for l in f.loops]; gv = [l.vert for l in g.loops]
            same = ((fa.index(e.verts[1]) - fa.index(e.verts[0])) % len(fa) == 1) == \
                   ((gv.index(e.verts[1]) - gv.index(e.verts[0])) % len(gv) == 1)
            if same:
                f.normal_flip()
            break
for f in new_faces + bridge_faces:
    f.smooth = True; f.material_index = 0
bm.normal_update()

# UVs: sample the faired base on the same side of the midline seam as the face
for f in new_faces + bridge_faces:
    side = 1 if f.calc_center_median().x >= 0 else -1
    for l in f.loops:
        v = l.vert
        if v in info:
            p = base_pos[v]
        else:                                      # rim vertex: keep its own surviving UV from the same side
            same = [ll for ll in v.link_loops if ll.face not in new_faces and ll.face not in bridge_faces
                    and (1 if ll.face.calc_center_median().x >= 0 else -1) == side]
            if same:
                for lay in uv_layers:
                    l[lay].uv = same[0][lay].uv
                continue
            p = v.co
        for lay in uv_layers:
            l[lay].uv = sample_uv(p, side, lay.name)

# light blend of the bridge strip
strip = (set(L) | set(outer)) - {v for v in L if v.is_boundary and v not in L}
for _ in range(3):
    for fct in (0.5, -0.53):
        new = {v: v.co + (sum((e.other_vert(v).co for e in v.link_edges), Vector()) / len(v.link_edges) - v.co) * fct
               for v in strip}
        for v, p in new.items():
            v.co = p
for v in bm.verts:
    if abs(v.co.x) < 2e-5:
        v.co.x = 0.0
bm.normal_update()
ring0 = [v.co.copy() for v in rings[0]]
A = sum(ring0, Vector()) / len(ring0)
widths = [abs(v.co.x) for v in rings[-1]]
report["new_verts"] = len(info); report["degenerate_faces"] = sum(1 for f in bm.faces if f.calc_area() < 1e-12)
report["anus_centre_A"] = r4(A); report["boundary_edges_total"] = sum(e.is_boundary for e in bm.edges)
report["pucker_outer_lateral_halfwidth_mm"] = round(max(widths) * 1000, 1)
report["slit_length_mm"] = round((R0[0].co - R0[H].co).length * 1000, 1)
pr = rings[-2]; q = len(pr) // 4         # ring where the pigment has faded out
floor_h = to_local(base_pos[M[18]])[2]
report["pigment_outline_mm"] = {
    "AP": round((pr[0].co - pr[2 * q].co).length * 1000, 1),
    "lateral_over_surface": round(2 * K_SPOKE * (8.0 - SLIT), 1),
    "lateral_projected": round((pr[q].co.x - pr[3 * q].co.x) * 1000, 1),
    "climb_up_cheek_wall": round((to_local(pr[q].co)[2] - floor_h) * 1000, 1)}
# The slit's lips bunch up at its two tips (elliptic spacing), leaving near-zero-width triangles whose shading
# normals spike. Weld pucker vertices closer than 0.06 mm there.
TIP_CO = [R0[0].co.copy(), R0[H].co.copy()]          # the slit's real tip positions
_tips = [v for v in info if v.is_valid and (v.co - R0[0].co).length < 0.0006 or (v.co - R0[H].co).length < 0.0006]
_n0 = len(bm.verts)                                  # (v3: no proximity weld - the lips are ~0.02 mm apart along
                                                     # the whole slit, so it collapsed both lips into tangled fans)
bm.verts.ensure_lookup_table(); bm.normal_update()
# soften the little pits where creases meet the slit's lips: normal-only smoothing on the skin right around the slit
_near = [v for v in bm.verts if v.is_valid and (min(abs(to_local(v.co)[0]) - SLIT / 1000, 0) ** 2 + 0) >= 0
         and abs(to_local(v.co)[0]) < (SLIT + 1.6) / 1000 and abs(to_local(v.co)[1]) < 0.0006 and not v.is_boundary
         and not (v in slit_verts) and v not in set(rings[0])]   # the slit line itself stays, so the groove keeps its depth
for _ in range(0):                                   # off: any averaging erases the narrow grooves
    bm.normal_update()
    _new = {}
    for v in _near:
        d_ = sum((e.other_vert(v).co for e in v.link_edges), Vector()) / len(v.link_edges) - v.co
        _new[v] = v.co + v.normal * (0.5 * d_.dot(v.normal))
    for v, p_ in _new.items():
        v.co = p_
bm.normal_update()
_skin_tree = BVHTree.FromBMesh(bm)                # the canal's bottom ring follows the final skin
ring0 = [_skin_tree.find_nearest(p_)[0] for p_ in ring0]
# open the fissure: drop the skin that zipped the slit shut, so it leads straight into the anal canal (the canal's
# bottom ring sits exactly on the lips)
_n0 = len(bm.verts)
OPEN_MARGIN = 0.35                                   # mm at each end of the slit that stays zipped shut, so the
bmesh.ops.delete(bm, geom=[v for v in slit_verts if v.is_valid], context='VERTS')    # lips joined pairwise below
_tmap = {}
for k in range(1, H):                                # lip vertex k and its partner N0-k on the other lip
    if abs(LIP_X[k]) >= SLIT - OPEN_MARGIN:
        a_, b_ = R0[k], R0[N0 - k]
        a_.co = (a_.co + b_.co) / 2
        _tmap[b_] = a_
bmesh.ops.weld_verts(bm, targetmap=_tmap)            # a clean seam: every seam edge has one face on each side
report["tip_seam_pairs"] = len(_tmap)
report["fissure_open_removed_verts"] = _n0 - len(bm.verts)
bmesh.ops.dissolve_degenerate(bm, dist=0.00002, edges=bm.edges[:])
bm.normal_update()
# the lips are pressed almost together, so at each tip several rings bunch into one spot and fold over: merge
# everything within 0.4 mm of the tip into a single point (the smoothing pass below rounds the patch)
# (v3: no tip merge - it made non-manifold fans; the tips stay zipped instead)
bm.verts.ensure_lookup_table(); bm.normal_update()
# tidy the tips after the weld / zipper removal: smooth a small patch around each end of the fissure
_tipc = [O + u * (SLIT / 1000), O - u * (SLIT / 1000)]
_tipv = [v for v in bm.verts if v.is_valid and not v.is_boundary and min((v.co - t_).length for t_ in TIP_CO) < 0.0006]
for _ in range(0):                                   # v3: off - its disc edge cut across the creases
    bm.normal_update()
    _new = {v: v.co + v.normal * (0.5 * (sum((e.other_vert(v).co for e in v.link_edges), Vector())
                                          / len(v.link_edges) - v.co).dot(v.normal)) for v in _tipv}
    for v, p_ in _new.items():
        v.co = p_
_inner = [v for v in _tipv if min((v.co - t_).length for t_ in TIP_CO) < 0.0006]
for _ in range(0):                                   # v3: off (no fold-overs left to untangle)
    _new = {v: v.co.lerp(sum((e.other_vert(v).co for e in v.link_edges), Vector()) / len(v.link_edges), 0.5)
            for v in _inner}
    for v, p_ in _new.items():
        v.co = p_
bm.normal_update()
report["tip_tidied_verts"] = len(_tipv)
bm.verts.ensure_lookup_table(); bm.normal_update()
_skin_tree = BVHTree.FromBMesh(bm)                    # canal opening follows the FINAL skin (after the tip merge),
ring0 = [_skin_tree.find_nearest(p_)[0] - n * 0.00008 for p_ in ring0]   # tucked 0.08 mm under it so none shows
# the canal's opening IS the skin's final slit edge (after the tip merge): same points, so nothing shows between them
_slit_edges = {e for e in bm.edges if e.is_boundary and all((v.co - A).length < 0.006 for v in e.verts)}
_loops, _seen = [], set()
for e0 in _slit_edges:
    if e0 in _seen:
        continue
    comp, st = [], [e0]
    while st:
        e = st.pop()
        if e in _seen:
            continue
        _seen.add(e); comp.append(e)
        st += [g for v in e.verts for g in v.link_edges if g in _slit_edges and g not in _seen]
    _loops.append(comp)
_loop = max(_loops, key=len); _lv = {v for e in _loop for v in e.verts}
_start = next(iter(_lv)); _order = [_start]; _prev = None; _cur = _start
while True:
    _nx = [e.other_vert(_cur) for e in _cur.link_edges if e in _loop and e.other_vert(_cur) is not _prev]
    if not _nx or _nx[0] is _start:
        break
    _prev, _cur = _cur, _nx[0]; _order.append(_cur)
ring0 = [v.co.copy() for v in _order]
report["canal_opening_points"] = len(ring0); report["slit_edge_loops"] = len(_loops)
refined_ids = None
bm.verts.ensure_lookup_table()
refined_ids = [v.index for v in bm.verts if in_refine(v.co)]
bm.to_mesh(me); me.update(); bm.free(); base.free(); btri.free()
if me.has_custom_normals:
    keep = [tuple(c.vector) for c in me.corner_normals]
    reg = set(refined_ids)
    me.normals_split_custom_set([(0.0, 0.0, 0.0) if l.vertex_index in reg else keep[i] for i, l in enumerate(me.loops)])
    me.update()


def cleft_nm_off(p):             # 1 over the refined cleft, fading to 0 before the edge of the refined region
    a, b, hr = cleft_coords(p)
    fa = 1 - (ss(A_BACK - 18, A_BACK - 3, a) if a > 0 else ss(A_FRONT - 17, A_FRONT - 3, -a))
    return fa * (1 - ss(16, 26, abs(b))) * ss(-14, -8, hr) * (1 - ss(22, 29, hr))


if "cleft_nm_off" in me.attributes:
    me.attributes.remove(me.attributes["cleft_nm_off"])
_nm = me.attributes.new("cleft_nm_off", 'FLOAT', 'POINT')
_nm.data.foreach_set("value", [cleft_nm_off(v.co) for v in me.vertices])

# ======================= 6. remove the old pigment dot from the textures =======================
SPOT = (0.1566, 0.5109); RUV = 0.0098
for name in ["cf_m_body_MT_CT.png", "cf_m_body_MT_DT.png", "cf_m_body_CM.png", "cf_m_body_DM.png",
             "cf_m_body_LM.png", "cf_m_body_NMP_CNV.png", "cf_m_body_NMPD_CNV.png"]:
    im = bpy.data.images.get(name) or next(i for i in bpy.data.images  # datablock may be named without the cf_m_ prefix
                                             if i.filepath.replace("\\", "/").split("/")[-1] == name)
    w, h = im.size
    px = np.array(im.pixels[:]).reshape(h, w, 4)
    cx, cy, R = SPOT[0] * w, SPOT[1] * h, RUV * w
    yy, xx = np.mgrid[0:h, 0:w]
    dist = np.hypot(xx + 0.5 - cx, yy + 0.5 - cy)
    inside = dist < R
    ring_ = (dist > R) & (dist < R + max(2, R * 0.3))
    diff = np.abs(px[inside, :3].mean(0) - px[ring_, :3].mean(0)).max()
    if diff < 0.01:
        report.setdefault("textures_untouched", []).append(name); continue
    out = px.copy()
    for x in range(int(cx - R) - 1, int(cx + R) + 2):
        dx = x + 0.5 - cx
        if abs(dx) >= R:
            continue
        hh = math.sqrt(R * R - dx * dx)
        y0, y1 = int(math.floor(cy - hh - 1)), int(math.ceil(cy + hh + 1))
        for y in range(y0 + 1, y1):
            t = (y - y0) / (y1 - y0)
            fill = px[y0, x] * (1 - t) + px[y1, x] * t
            k = ss(0.8, 1.0, math.hypot(dx, y + 0.5 - cy) / R)
            out[y, x] = fill * (1 - k) + px[y, x] * k
    im.pixels = out.ravel(); im.update(); im.pack()
    report.setdefault("textures_inpainted", []).append([name, round(float(diff), 3)])

# ======================= 7. pigment as an editable shader layer =======================
def srgb2lin(c):
    return tuple(x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in c)


mat = bpy.data.materials["Body"]; nt = mat.node_tree
for nd in [x for x in nt.nodes if x.name.startswith("AnusPig")]:
    nt.nodes.remove(nd)
gen, shd = nt.nodes["Gentex"], nt.nodes["Shader"]
frame = nt.nodes.new("NodeFrame"); frame.name = "AnusPig_Frame"; frame.label = "Anus pigment (edit colours / strength here)"
attr = nt.nodes.new("ShaderNodeAttribute"); attr.name = "AnusPig_Mask"; attr.attribute_name = "anus_pigment"
attr.attribute_type = 'GEOMETRY'
stren = nt.nodes.new("ShaderNodeMath"); stren.name = "AnusPig_Strength"; stren.label = "Strength"
stren.operation = 'MULTIPLY'; stren.inputs[1].default_value = 1.0; stren.use_clamp = True
nt.links.new(attr.outputs["Fac"], stren.inputs[0])
cattr = nt.nodes.new("ShaderNodeAttribute"); cattr.name = "AnusPig_CreaseMask"; cattr.attribute_name = "anus_crease"
cattr.attribute_type = 'GEOMETRY'
cstren = nt.nodes.new("ShaderNodeMath"); cstren.name = "AnusPig_CreaseStrength"; cstren.label = "Crease line strength"
cstren.operation = 'MULTIPLY'; cstren.inputs[1].default_value = 0.0; cstren.use_clamp = True   # 0 = no painted lines
nt.links.new(cattr.outputs["Fac"], cstren.inputs[0])
cols = {"maintex": ("Lit colour", (0.98, 0.66, 0.71), "Crease colour (lit)", (0.82, 0.47, 0.55)),
        "darktex": ("Shadow colour", (0.84, 0.55, 0.61), "Crease colour (shadow)", (0.66, 0.37, 0.45))}
x0, y0 = shd.location.x - 700, shd.location.y + 420
for i, (sock, (lab, rgb, clab, crgb)) in enumerate(cols.items()):
    out_s = gen.outputs[f"Body {sock}"]; in_s = shd.inputs[f"Body {sock}"]
    for lk in list(out_s.links):
        if lk.to_socket == in_s:
            nt.links.remove(lk)
    col = nt.nodes.new("ShaderNodeRGB"); col.name = f"AnusPig_{sock}_Colour"; col.label = lab
    col.outputs[0].default_value = (*srgb2lin(rgb), 1.0)
    mix = nt.nodes.new("ShaderNodeMix"); mix.name = f"AnusPig_{sock}_Mix"; mix.data_type = 'RGBA'
    assert mix.inputs[6].type == 'RGBA' and mix.inputs[7].type == 'RGBA'
    nt.links.new(stren.outputs[0], mix.inputs[0]); nt.links.new(out_s, mix.inputs[6])
    nt.links.new(col.outputs[0], mix.inputs[7])
    ccol = nt.nodes.new("ShaderNodeRGB"); ccol.name = f"AnusPig_{sock}_CreaseColour"; ccol.label = clab
    ccol.outputs[0].default_value = (*srgb2lin(crgb), 1.0)
    cmix = nt.nodes.new("ShaderNodeMix"); cmix.name = f"AnusPig_{sock}_CreaseMix"; cmix.data_type = 'RGBA'
    nt.links.new(cstren.outputs[0], cmix.inputs[0]); nt.links.new(mix.outputs[2], cmix.inputs[6])
    nt.links.new(ccol.outputs[0], cmix.inputs[7]); nt.links.new(cmix.outputs[2], in_s)
    col.location = (x0, y0 - 260 * i); mix.location = (x0 + 200, y0 - 260 * i)
    ccol.location = (x0 + 200, y0 - 260 * i - 130); cmix.location = (x0 + 400, y0 - 260 * i)
    for nd in (col, mix, ccol, cmix):
        nd.parent = frame
attr.location = (x0 - 420, y0); stren.location = (x0 - 210, y0)
cattr.location = (x0 - 420, y0 - 300); cstren.location = (x0 - 210, y0 - 300)
for nd in (attr, stren, cattr, cstren):
    nd.parent = frame

# 3D creases in a two-tone shader: the hard cut hides small relief, so inside the anus the toon result is darkened
# by how much the creased surface is lit relative to the smooth surface under it (light ratio from two Diffuse
# evaluations). It follows the scene light and uses the shader's own light / shadow colours.
for nd in [x for x in nt.nodes if x.name.startswith("AnusFold")]:
    nt.nodes.remove(nd)
raw_out = nt.nodes["RawShade"].outputs["Raw shading"]; raw_in = nt.nodes["Shader"].inputs["Raw Shading"]
for lk in list(raw_out.links):
    if lk.to_socket == raw_in:
        nt.links.remove(lk)
fframe = nt.nodes.new("NodeFrame"); fframe.name = "AnusFold_Frame"; fframe.label = "Anus crease shading (3D)"


def fnode(kind, name, label=None, op=None, vals=None, clamp=False):
    nd = nt.nodes.new(kind); nd.name = "AnusFold_" + name; nd.parent = fframe
    if label: nd.label = label
    if op: nd.operation = op
    for i, val in (vals or {}).items():
        nd.inputs[i].default_value = val
    if clamp: nd.use_clamp = True
    return nd


fmask = fnode("ShaderNodeAttribute", "Mask"); fmask.attribute_name = "anus_fold_shade"; fmask.attribute_type = 'GEOMETRY'
fbn = fnode("ShaderNodeAttribute", "BaseNormal"); fbn.attribute_name = "anus_base_normal"; fbn.attribute_type = 'GEOMETRY'
fnorm = fnode("ShaderNodeVectorMath", "BaseNormalUnit", op='NORMALIZE')
dgeo = fnode("ShaderNodeBsdfDiffuse", "LightCreased"); dbase = fnode("ShaderNodeBsdfDiffuse", "LightSmooth")
for d_ in (dgeo, dbase):
    d_.inputs["Color"].default_value = (1, 1, 1, 1)
rgeo = fnode("ShaderNodeShaderToRGB", "RGBCreased"); rbase = fnode("ShaderNodeShaderToRGB", "RGBSmooth")
bgeo = fnode("ShaderNodeRGBToBW", "BWCreased"); bbase = fnode("ShaderNodeRGBToBW", "BWSmooth")
eps = fnode("ShaderNodeMath", "Eps", op='ADD', vals={1: 1e-4})
ratio = fnode("ShaderNodeMath", "Ratio", op='DIVIDE', clamp=True)
gam = fnode("ShaderNodeMath", "Contrast", label="Contrast (power)", op='POWER', vals={1: 1.6}, clamp=True)
inv = fnode("ShaderNodeMath", "Darkening", op='SUBTRACT', vals={0: 1.0})
fstr = fnode("ShaderNodeMath", "Strength", label="Crease shading strength", op='MULTIPLY', vals={1: 0.8})
fmul = fnode("ShaderNodeMath", "Masked", op='MULTIPLY')
keep = fnode("ShaderNodeMath", "Keep", op='SUBTRACT', vals={0: 1.0})
fout = fnode("ShaderNodeMath", "Apply", op='MULTIPLY', clamp=True)
L = nt.links.new
L(fbn.outputs["Vector"], fnorm.inputs[0]); L(fnorm.outputs["Vector"], dbase.inputs["Normal"])
L(dgeo.outputs[0], rgeo.inputs[0]); L(dbase.outputs[0], rbase.inputs[0])
L(rgeo.outputs["Color"], bgeo.inputs[0]); L(rbase.outputs["Color"], bbase.inputs[0])
L(bbase.outputs[0], eps.inputs[0]); L(bgeo.outputs[0], ratio.inputs[0]); L(eps.outputs[0], ratio.inputs[1])
L(ratio.outputs[0], gam.inputs[0]); L(gam.outputs[0], inv.inputs[1])
L(inv.outputs[0], fstr.inputs[0]); L(fstr.outputs[0], fmul.inputs[0]); L(fmask.outputs["Fac"], fmul.inputs[1])
L(fmul.outputs[0], keep.inputs[1]); L(raw_out, fout.inputs[0]); L(keep.outputs[0], fout.inputs[1])
L(fout.outputs[0], raw_in)
_x, _y = nt.nodes["Shader"].location.x - 1500, nt.nodes["Shader"].location.y - 700
for i, nd in enumerate((fmask, fbn, fnorm, dgeo, dbase, rgeo, rbase, bgeo, bbase, eps, ratio, gam, inv, fstr, fmul, keep, fout)):
    nd.location = (_x + 180 * (i // 2), _y - 160 * (i % 2))

# body normal map off over the rebuilt cleft: the 1024 px map is magnified ~15x here and still carries the old
# cleft/anus shape, so under the toon cut it draws texel-sized spikes (sawtooth pigment edge, white sliver above the
# anus). The refined geometry (2c) carries the shape instead. Drives RawShade's own "Use Normals?" input.
for nd in [x for x in nt.nodes if x.name.startswith("CleftNM")]:
    nt.nodes.remove(nd)
rs = nt.nodes["RawShade"]; use_nm = rs.inputs["Use Normals?"]
nframe = nt.nodes.new("NodeFrame"); nframe.name = "CleftNM_Frame"; nframe.label = "Body normal map off in the cleft"
nma = nt.nodes.new("ShaderNodeAttribute"); nma.name = "CleftNM_Mask"; nma.attribute_name = "cleft_nm_off"
nma.attribute_type = 'GEOMETRY'
nms = nt.nodes.new("ShaderNodeMath"); nms.name = "CleftNM_Strength"; nms.label = "Strength (1 = normal map fully off)"
nms.operation = 'MULTIPLY'; nms.inputs[1].default_value = 1.0; nms.use_clamp = True
nmi = nt.nodes.new("ShaderNodeMath"); nmi.name = "CleftNM_UseNormals"; nmi.label = "Use Normals?"
nmi.operation = 'MULTIPLY_ADD'; nmi.inputs[1].default_value = -float(use_nm.default_value)
nmi.inputs[2].default_value = float(use_nm.default_value)          # = original * (1 - mask)
nt.links.new(nma.outputs["Fac"], nms.inputs[0]); nt.links.new(nms.outputs[0], nmi.inputs[0])
nt.links.new(nmi.outputs[0], use_nm)
nma.location = (rs.location.x - 640, rs.location.y - 260); nms.location = (rs.location.x - 430, rs.location.y - 260)
nmi.location = (rs.location.x - 220, rs.location.y - 260)
for nd in (nma, nms, nmi):
    nd.parent = nframe

# ======================= 8. internals: raise fit 1 cm, re-seat sphincter =======================
E = bpy.data.objects["Internal_Fit_Xform"]
E.location = (0.0, 0.0736, 0.252)
col_fit = bpy.data.collections["Internal_Fit"]
for nm in ["AN_AnalSphincter", "AN_AnalCanal", "AN_Rectum_LowerAmpulla"]:
    ob = bpy.data.objects.get(nm)
    if ob:
        mesh = ob.data; bpy.data.objects.remove(ob)
        if mesh and mesh.users == 0:
            bpy.data.meshes.remove(mesh)
with bpy.data.libraries.load(REF) as (src, dst):
    dst.objects = ["AN_AnalSphincter"]
sph = dst.objects[0]; col_fit.objects.link(sph); sph.parent = E
bpy.context.view_layer.update()

# ======================= 9. anal canal + lower ampulla =======================
P = np.array([tuple(v.co) for v in me.vertices])
navel_c = P[(np.abs(P[:, 0]) < 0.002) & (P[:, 1] < 0) & (P[:, 2] > 0.90) & (P[:, 2] < 0.95)]
navel = Vector(navel_c[navel_c[:, 1].argmax()])
d = (navel - A).normalized()
L_CANAL = 0.025
J = A + d * L_CANAL

rect = bpy.data.objects["AN_Rectum"]; Mr = rect.matrix_world.copy()
sac = bpy.data.objects["AN_Sacrum"]
_sb = bmesh.new(); _sb.from_mesh(sac.data); _sb.transform(sac.matrix_world); bone_tree = BVHTree.FromBMesh(_sb)
GAP = 0.005                                     # clearance between rectum and sacrum/coccyx (mesorectal fat)


def band_push(points, z_of, fixed_w):
    """forward (-y) shift per 3 mm height band so every band clears the bone by GAP; returns per-point shift"""
    gaps = [bone_tree.find_nearest(p)[3] for p in points]
    zs = np.array([z_of(p) for p in points]); bins = np.floor(zs / 0.003).astype(int)
    need = {}
    for bi, g in zip(bins, gaps):
        need[bi] = max(need.get(bi, 0.0), GAP - g)
    keys = np.array(sorted(need)); vals = np.array([max(0.0, need[k]) for k in keys])
    sm = np.array([max(vals[i], (vals * np.exp(-((keys - keys[i]) * 3.0) ** 2 / (2 * 4.0 ** 2))).max())
                   for i in range(len(keys))])      # spread the push over ~4 mm of height
    lut = dict(zip(keys, sm))
    return [lut[b] * fixed_w(p) for b, p in zip(bins, points)]


# 1) lower rectum (atlas) clears the sacrum/coccyx; fades out before the junction with the colon
rmesh = rect.data; Mi_r = Mr.inverted()
for _ in range(3):
    W = [Mr @ v.co for v in rmesh.vertices]
    sh_ = band_push(W, lambda p: p.z, lambda p: 1 - ss(0.862, 0.876, p.z))
    for v, p, s in zip(rmesh.vertices, W, sh_):
        v.co = Mi_r @ (p - Vector((0, s, 0)))
rmesh.update()
_rg = min(bone_tree.find_nearest(Mr @ v.co)[3] for v in rmesh.vertices if (Mr @ v.co).z < 0.86)
report["rectum_gap_to_bone_below_860mm_mm"] = round(_rg * 1000, 1)

rb = bmesh.new(); rb.from_mesh(rect.data); rb.transform(Mr)
loops = []; seen = set()
for e in rb.edges:
    if not e.is_boundary or e in seen:
        continue
    st = [e]; comp = []
    while st:
        x = st.pop()
        if x in seen:
            continue
        seen.add(x); comp.append(x)
        for vv in x.verts:
            st += [y for y in vv.link_edges if y.is_boundary and y not in seen]
    loops.append(comp)
low = max(loops, key=len)                          # the clean cut rim at the bottom of the rectum
vs = {v for e in low for v in e.verts}
s0 = next(iter(vs)); order = [s0]; prev = None; cur = s0
while True:
    nx = [e.other_vert(cur) for e in cur.link_edges if e in low and e.other_vert(cur) != prev]
    if not nx or nx[0] == s0:
        break
    prev, cur = cur, nx[0]; order.append(cur)
rect_ring = [v.co.copy() for v in order]
Rc = sum(rect_ring, Vector()) / len(rect_ring)
Pr = np.array([tuple(p) for p in rect_ring]); _, Vr = np.linalg.eigh((Pr - Pr.mean(0)).T @ (Pr - Pr.mean(0)))
nr = Vector(Vr[:, 0]); nr = nr if nr.z > 0 else -nr
rb.free()

# bezier for the ampulla from J (anorectal junction) up into the rectum rim
L2 = (Rc - J).length
T0 = Vector((0.0, math.sin(math.radians(35)), math.cos(math.radians(35))))
B0, B1, B2, B3 = J, J + T0 * L2 * 0.4, Rc - nr * L2 * 0.35, Rc


def bez(t):
    s = 1 - t
    return B0 * s ** 3 + B1 * 3 * s * s * t + B2 * 3 * s * t * t + B3 * t ** 3


def bez_t(t):
    s = 1 - t
    return ((B1 - B0) * 3 * s * s + (B2 - B1) * 6 * s * t + (B3 - B2) * 3 * t * t).normalized()


X = Vector((1, 0, 0))


def frame_at(T):
    e1 = (X - T * X.dot(T)).normalized(); return e1, T.cross(e1)


def ring_fn(pts, ctr, T):
    e1, e2 = frame_at(T); arr = []
    for p in pts:
        q = p - ctr; q = q - T * q.dot(T)
        arr.append((math.atan2(q.dot(e2), q.dot(e1)), q.length))
    arr.sort(); th = np.array([x for x, _ in arr]); rr = np.array([x for _, x in arr])
    th = np.concatenate([th - 2 * np.pi, th, th + 2 * np.pi]); rr = np.concatenate([rr, rr, rr])
    return lambda t: float(np.interp(t, th, rr))


f_bot = ring_fn(ring0, A, d); f_top = ring_fn(rect_ring, Rc, nr)
NS = 36


def canal_r(th, t):   # collapsed canal: slit at the verge, ~3.2 mm with anal columns higher up, 4 mm at the junction
    rc = (3.2 + 0.8 * ss(0.5, 1.0, t)) * (1 + 0.15 * ss(0.4, 0.8, t) * math.cos(8 * th)) / 1000
    return f_bot(th) * (1 - ss(0.08, 0.35, t)) + rc * ss(0.08, 0.35, t)   # collapsed (slit-thin) for the first ~2 mm


def build_tube(name, ring_pts_list, exact_bot, exact_top):
    tb = bmesh.new()
    loops_v = [[tb.verts.new(p) for p in exact_bot]] + [[tb.verts.new(p) for p in rp] for rp in ring_pts_list] + \
              [[tb.verts.new(p) for p in exact_top]]
    for lv in loops_v:
        for k in range(len(lv)):
            tb.edges.new((lv[k], lv[(k + 1) % len(lv)]))
    for a_, b_ in zip(loops_v[:-1], loops_v[1:]):
        ea = [tb.edges.get((a_[k], a_[(k + 1) % len(a_)])) for k in range(len(a_))]
        eb = [tb.edges.get((b_[k], b_[(k + 1) % len(b_)])) for k in range(len(b_))]
        bmesh.ops.bridge_loops(tb, edges=ea + eb)
    bmesh.ops.recalc_face_normals(tb, faces=tb.faces)
    tb.normal_update()
    f = tb.faces[len(tb.faces) // 2]; cc = f.calc_center_median()
    axis_pts = [sum(lv_, Vector()) / len(lv_) for lv_ in [[v.co for v in lv] for lv in loops_v]]
    near_c = min(axis_pts, key=lambda q: (q - cc).length)
    if f.normal.dot(cc - near_c) < 0:
        bmesh.ops.reverse_faces(tb, faces=tb.faces)
    for f in tb.faces:
        f.smooth = True
    mesh = bpy.data.meshes.new(name); tb.to_mesh(mesh); tb.free()
    ob = bpy.data.objects.new(name, mesh); col_fit.objects.link(ob)
    mesh.materials.append(rect.data.materials[0])
    return ob


canal_rings = []
for i in range(1, 9):
    t = i / 9; ctr = A.lerp(J, t); e1, e2 = frame_at(d)
    canal_rings.append([ctr + (e1 * math.cos(th) + e2 * math.sin(th)) * canal_r(th, t)
                        for th in (-math.pi + 2 * math.pi * k / NS for k in range(NS))])
e1, e2 = frame_at(d)
J_ring = [J + (e1 * math.cos(th) + e2 * math.sin(th)) * canal_r(th, 1.0)
          for th in (-math.pi + 2 * math.pi * k / NS for k in range(NS))]
canal = build_tube("AN_AnalCanal", canal_rings, ring0, J_ring)

amp_rings = []
for i in range(1, 8):
    t = i / 8; ctr = bez(t); T = bez_t(t); e1, e2 = frame_at(T); w = ss(0.0, 1.0, t)
    ring_ = []
    for k in range(NS):
        th = -math.pi + 2 * math.pi * k / NS
        rr = canal_r(th, 1.0) * (1 - w) + f_top(th) * w
        ring_.append(ctr + (e1 * math.cos(th) + e2 * math.sin(th)) * rr)
    amp_rings.append(ring_)
# 2) intermediate ampulla rings also clear the coccyx (junction ring and rectum rim stay put)
for _ in range(3):
    need = []
    for ring_ in amp_rings:
        need.append(max(0.0, max(GAP - bone_tree.find_nearest(p)[3] for p in ring_)))
    need = [max(need[i], max(need[j] * math.exp(-((i - j) ** 2) / 2.0) for j in range(len(need))))
            for i in range(len(need))]
    amp_rings = [[p - Vector((0, s * math.sin(math.pi * (i + 1) / 8) ** 0.5, 0)) for p in ring_]
                 for i, (ring_, s) in enumerate(zip(amp_rings, need))]
amp = build_tube("AN_Rectum_LowerAmpulla", amp_rings, J_ring, rect_ring)

# sphincter (external anal sphincter ring) around the lower canal
Ps = np.array([tuple(sph.matrix_world @ v.co) for v in sph.data.vertices]); sc_ = Vector(Ps.mean(0))
_, Vs = np.linalg.eigh((Ps - Ps.mean(0)).T @ (Ps - Ps.mean(0))); sax = Vector(Vs[:, 0])
sax = sax if sax.dot(d) > 0 else -sax
rot = sax.rotation_difference(d).to_matrix().to_4x4()
tgt = A + d * (0.4 * L_CANAL)
Mnew = Matrix.Translation(tgt) @ rot @ Matrix.Translation(-sc_) @ sph.matrix_world
sph.data.transform(sph.matrix_world.inverted() @ Mnew); sph.data.update()

# ======================= 10. anatomical checks =======================
sac = bpy.data.objects["AN_Sacrum"]
Sw = np.array([tuple(sac.matrix_world @ v.co) for v in sac.data.vertices]); coccyx = Vector(Sw[Sw[:, 2].argmin()])
hipL = bpy.data.objects["AN_HipBone_L"]; hipR = bpy.data.objects["AN_HipBone_R"]
Hw = np.vstack([[tuple(o.matrix_world @ v.co) for v in o.data.vertices] for o in (hipL, hipR)])
med = Hw[np.abs(Hw[:, 0]) < 0.004]; symph_low = Vector(med[med[:, 2].argmin()])
tub = Hw[Hw[:, 2].argmin()]
pcl_z_at_J = symph_low.z + (coccyx.z - symph_low.z) * (J.y - symph_low.y) / (coccyx.y - symph_low.y)
amp_dir = bez_t(0.0)
sb = bmesh.new(); sb.from_mesh(sac.data); sb.transform(sac.matrix_world); stree = BVHTree.FromBMesh(sb)
report["min_gap_ampulla_to_sacrum_coccyx_mm"] = round(min(stree.find_nearest(amp.matrix_world @ v.co)[3]
                                                         for v in amp.data.vertices) * 1000, 1)
report["min_gap_canal_to_coccyx_mm"] = round(min(stree.find_nearest(canal.matrix_world @ v.co)[3]
                                                 for v in canal.data.vertices) * 1000, 1)
sb.free()
report.update({
    "navel": r4(navel), "canal_axis": r4(d), "canal_length_mm": round(L_CANAL * 1000, 1),
    "canal_angle_from_vertical_deg": round(math.degrees(math.acos(d.z)), 1),
    "anorectal_angle_deg": round(math.degrees((A - J).angle(amp_dir)), 1),
    "junction_J": r4(J), "junction_to_coccyx_tip_mm": round((J - coccyx).length * 1000, 1),
    "junction_vs_pubococcygeal_line_mm": round((J.z - pcl_z_at_J) * 1000, 1),
    "anus_vs_ischial_tuberosity_height_mm": round((A.z - tub[2]) * 1000, 1),
    "coccyx_tip": r4(coccyx), "symphysis_low": r4(symph_low), "ampulla_len_mm": round(L2 * 1000, 1),
})
print("REPORT", json.dumps(report))
json.dump(report, open(os.path.join(OUT, "report.json"), "w"), indent=1)
if SAVE_AS == "INPLACE":
    bpy.ops.wm.save_mainfile()
else:
    bpy.ops.wm.save_as_mainfile(filepath=SAVE_AS, relative_remap=False)
