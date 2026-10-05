"""Detail pass on the anus of Hips.blend (v8: squeezed pucker):
skin crease smoothing, cleft walls pressed in around the anus, high-detail pucker rebuilt as an AP slit with
creases fanning out of it, laid out by arc length over the skin so it rides up the cheek walls,
pigment as an editable shader mask, internal fit raised 1 cm, angled anal canal + lower rectal ampulla,
sphincter re-seated."""
import bpy, bmesh, math, sys, json, os
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
SLIT = 3.0                       # mm, half-length of the closed AP slit
K_LAT = 0.85                     # lateral (arc length) / AP for the confocal rings -> AP-elongated
RINGS = [3.06, 3.5, 3.9, 4.4, 5.0, 5.7, 6.4, 7.2, 8.0, 8.8]     # AP semi-axis of each ring, mm
# folds on one lip (theta in [0, pi]); mirrored to the other lip so the skin stays exactly mirror-symmetric.
# (angle deg, weight, width rad, AP radius where the fold dies out) - uneven so they don't read as a ladder
FOLDS = [(0, 0.8, 0.07, 6.0), (18, 0.6, 0.055, 5.4), (37, 0.95, 0.065, 6.4), (57, 0.75, 0.06, 5.7),
         (78, 1.05, 0.07, 6.6), (97, 0.85, 0.06, 6.1), (117, 1.0, 0.065, 6.5), (139, 0.7, 0.055, 5.6),
         (160, 0.9, 0.06, 6.2), (180, 0.8, 0.07, 6.0)]
FOLDS += [(-a, w, sg, re) for a, w, sg, re in FOLDS if 0 < a < 180]
FUNNEL = 1.2; R_FOLD = 7.4; FOLD_DEPTH = 0.14; SLIT_DEPTH = 0.6
def crease(theta, r=None):
    s_ = 0.0
    for a, w, sg, re in FOLDS:
        dt = (theta - math.radians(a) + math.pi) % (2 * math.pi) - math.pi
        fade = 1.0 if r is None else 1 - ss(re - 1.2, re, r)
        s_ += w * fade * math.exp(-dt * dt / (2 * sg * sg))
    return s_


def fold_amp(r):                 # wrinkles are deepest just outside the slit and die out before the cheek walls
    return FOLD_DEPTH * ss(RINGS[0], RINGS[0] + 0.5, r) * (1 - ss(4.5, 6.4, r))


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


def ring_xy(r, th):
    """confocal ellipse with foci at the slit ends; r = AP semi-axis"""
    return r * math.cos(th), K_LAT * math.sqrt(max(r * r - SLIT * SLIT, 0.0)) * math.sin(th)


S_OUT = 18 if len(L) < 27 else (36 if len(L) < 54 else 72)
report["outer_ring_verts"] = S_OUT
rings = []; new_faces = []; info = {}
for r in RINGS:
    S = 144 if r < 6.0 else (72 if r < 7.0 else (max(S_OUT, 36) if r < 8.4 else S_OUT))
    ring = []
    for k in range(S):
        th = 2 * math.pi * k / S
        v = bm.verts.new(surf_point(*ring_xy(r, th))); info[v] = (r, th); ring.append(v)
    rings.append(ring)
# closed slit: a line of vertices on the midline zips the two lips of the innermost ring together
R0 = rings[0]; N0 = len(R0); H = N0 // 2; M = {0: R0[0], H: R0[H]}; slit_verts = {}
for k in range(1, H):
    th = 2 * math.pi * k / N0
    v = bm.verts.new(surf_point(RINGS[0] * math.cos(th), 0.0)); info[v] = (0.0, th); slit_verts[v] = th; M[k] = v


def add_face(vs):
    vs = [x for i, x in enumerate(vs) if x is not vs[i - 1]]
    new_faces.append(bm.faces.new(vs))


for k in range(H):
    add_face((M[k], R0[k], R0[k + 1], M[k + 1]))                       # lip with s > 0
    add_face((M[k + 1], R0[N0 - 1 - k], R0[(N0 - k) % N0], M[k]))      # lip with s < 0
for ri in range(len(rings) - 1):
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
# round the tight floor/wall corner of the squeezed cleft inside the pucker. Normal-only Laplacian steps, so nothing
# slides along the surface (ring / pigment layout stays where it was laid out); fades out over the outer rings so
# the rim still sits on the skin.
N_ROUND = 25
round_w = {v: 1 - ss(len(rings) - 4, len(rings) - 1, i) for i, ring in enumerate(rings) for v in ring}
round_w.update({v: 1.0 for v in slit_verts})
for _ in range(N_ROUND):
    bm.normal_update()
    new = {}
    for v, w_ in round_w.items():
        if w_ > 0:
            d_ = sum((e.other_vert(v).co for e in v.link_edges), Vector()) / len(v.link_edges) - v.co
            new[v] = v.co + v.normal * (0.5 * w_ * d_.dot(v.normal))
    for v, p_ in new.items():
        v.co = p_
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
# relief along the (smooth) base normal: shallow funnel, creases fanning from the slit, closed slit groove
for v, (r, th) in info.items():
    if v in slit_verts:
        sw = math.sin(th)
        relief = -FUNNEL - SLIT_DEPTH * sw ** 1.5
        v[pig_layer] = 1.0; v[cre_layer] = 0.9 * sw ** 0.5
    else:
        relief = -FUNNEL * (1 - ss(RINGS[0], R_FOLD, r)) + fold_amp(r) * (0.15 - crease(th, r))
        v[pig_layer] = 1 - ss(5.0, 8.0, r)
        v[cre_layer] = min(1.0, crease(th, r)) * paint_amp(r)
    v.co = base_pos[v] + bnorm[v] * (relief / 1000)

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
pr = rings[RINGS.index(8.0)]; q = len(pr) // 4         # ring where the pigment has faded out
floor_h = to_local(base_pos[M[18]])[2]
report["pigment_outline_mm"] = {
    "AP": round((pr[0].co - pr[2 * q].co).length * 1000, 1),
    "lateral_over_surface": round(2 * K_LAT * math.sqrt(8.0 ** 2 - SLIT ** 2), 1),
    "lateral_projected": round((pr[q].co.x - pr[3 * q].co.x) * 1000, 1),
    "climb_up_cheek_wall": round((to_local(pr[q].co)[2] - floor_h) * 1000, 1)}
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
cstren.operation = 'MULTIPLY'; cstren.inputs[1].default_value = 0.6; cstren.use_clamp = True
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
    return f_bot(th) * (1 - ss(0, 0.3, t)) + rc * ss(0, 0.3, t)


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
