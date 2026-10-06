"""build5_passage: build5_rise with a real passage - the skin is open at the centre and its edge loop IS the anal
canal's bottom ring (same vertices, same positions), so skin and canal form one continuous tube. At rest the opening
is a collapsed slit (lips LIP_GAP apart), so no hole shows; the canal stays collapsed for its first CANAL_COLL mm.
build5_rise: build5_adapted with the funnel turned into a rise - the surface climbs towards the centre (same
smoothstep profile, straight out of the body along n), and no sink at the very centre.
build5_adapted: the user's build5, unchanged except for three fixes (each marked "adapted:"):
  1. the region rebuilt is cut narrow, in surface coordinates, so the fill no longer flattens the cheek walls;
  2. that also stops the fill bridging the cleft, which made the anus stand up like a lump;
  3. the rings are laid out in arc length over the skin, so the edges ride up the cheek walls.
Plus portability only: ANATOMY_REF env var, report path, texture lookup by file name.
Detail pass on the anus of Hips.blend:
skin crease smoothing, high-detail pucker rebuild, pigment moved from texture to an editable shader mask,
internal fit raised 1 cm, anatomically angled anal canal + lower rectal ampulla, sphincter re-seated."""
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
# build5_rise: behind the anus this pass only raised the cleft floor (~1.3 mm), which tilted the back half of the anus
# up ~19 deg; it is needed in front (perineum side), where it hides the original's low-poly facets and crease. Fade it
# out behind the anus between TAUBIN_BACK[0] and [1] mm (ANUS_TAUBIN_BACK="a0,a1"; "0,0" = off: as build5)
TAUBIN_BACK = [float(x) for x in os.environ.get("ANUS_TAUBIN_BACK", "0,0").split(",")]   # off: faded, the original's facets show behind the anus
if TAUBIN_BACK[1] > 0:
    wt = {v: w_ * (1 - ss(TAUBIN_BACK[0], TAUBIN_BACK[1], to_local(v.co)[0] * 1000)) for v, w_ in wt.items()}
# build5_rise: this pass (from build7) raises the cleft floor behind the anus ~1.3 mm above the original, which
# left a steep climb ("cliff") at the back edge of the rebuilt anus. ANUS_TAUBIN = number of passes (build5: 20)
TAUBIN_ITERS = int(os.environ.get("ANUS_TAUBIN", "20"))
for _ in range(TAUBIN_ITERS):
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

# (build5_rise: ported from build8 - the low-poly cheek walls made the toon light/shadow cut zigzag)
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
from mathutils.bvhtree import BVHTree
from mathutils.interpolate import poly_3d_calc

# adapted: the region is cut in surface coordinates (a along the cleft, s = arc length across it), a fixed margin
# outside the outer ring. The projected ellipse (9.8 x 8 mm) reached far up the steep cheek walls: the fill flattened
# them and bridged the cleft, so the anus stood up like a lump.
H0 = 0.007                       # ray origin height above the cleft floor (in the air of the cleft)
D_AP, D_LAT = 8.6 + 1.2, 0.75 * 8.6 + 1.2     # outer ring semi-axes + margin, mm
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
base_tree = BVHTree.FromBMesh(base)
# adapted: the rings now lie on the real cheek walls (coarse game-mesh triangles), so cast onto a Phong-curved copy
# of the base: it passes through the original vertices and matches how that skin looks smooth-shaded
btri = base.copy(); bmesh.ops.triangulate(btri, faces=btri.faces[:]); btri.faces.ensure_lookup_table()
btri.normal_update()
btri_tree = BVHTree.FromBMesh(btri)


def phong(p, fi, alpha=0.75):
    vs = [l.vert for l in btri.faces[fi].loops]
    ws = poly_3d_calc([v.co for v in vs], p)
    q = Vector()
    for w_, v in zip(ws, vs):
        q += (p - v.normal * (p - v.co).dot(v.normal)) * w_
    return p.lerp(q, alpha)


def cast(C, dr):
    hit, nrm, fi, _ = btri_tree.ray_cast(C, dr, 0.05)
    return (None, nrm) if hit is None else (phong(hit, fi), nrm)


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

# ======================= 4/5. new pucker draped on the base: squeezed by the cheeks, creases, AP slit =======================
RINGS = [0.9, 1.4, 2.0, 2.7, 3.5, 4.4, 5.4, 6.5, 7.6, 8.6]          # mm (AP); lateral is squeezed
# rings inside the entrance, so the dip into the centre can curve smoothly (build5 had only the pole there, which can
# only make a straight cone); the creases don't reach inside 0.9 mm, so they are unaffected
RINGS = sorted([0.06, 0.15, 0.3, 0.45, 0.6, 0.75] + RINGS + [1.15, 1.65, 2.35])   # + finer steps over the shoulder,
# which otherwise bent only at build5's 0.9 / 1.4 / 2.0 mm rings and showed corners
# extra rings over the outer creases: with 288 spokes but build5's ~1 mm ring steps the faces there were long thin
# slivers, and the toon light/shadow cut followed them (crease ends split into streaks). RING_STEP mm, 0 = off
# build5_passage: the pole and the rings inside R_HOLE are gone; the R_HOLE ring is the skin's edge loop = the canal's
# mouth. Inside R_COLL the rings are squeezed sideways until, at R_HOLE, the two lips are LIP_GAP apart (a closed slit
# 2 * R_HOLE long at the bottom of the plunge)
R_HOLE = float(os.environ.get("ANUS_HOLE", "0.3"))       # mm
R_COLL = float(os.environ.get("ANUS_COLL", "0.75"))      # mm where the sideways collapse starts
LIP_GAP = 0.02                                            # mm between the lips at rest
RINGS = [r_ for r_ in RINGS if r_ >= R_HOLE - 1e-9]
assert abs(RINGS[0] - R_HOLE) < 1e-9, "ANUS_HOLE must be one of the ring radii"
RING_STEP = float(os.environ.get("ANUS_RING_STEP", "0"))
if RING_STEP > 0:
    _extra = [round(2.35 + RING_STEP * k, 3) for k in range(1, int((6.9 - 2.35) / RING_STEP) + 1)]
    RINGS = sorted(set(RINGS) | {r_ for r_ in _extra if all(abs(r_ - q) > RING_STEP / 3 for q in RINGS)})
SQUEEZE = 0.75
FOLD_W = [1.2, 0.8, 1.0, 0.75, 1.05, 0.85, 1.2, 0.85, 1.05, 0.75, 1.0, 0.8]
SIG = float(os.environ.get("ANUS_SIG", "0.06")); R_FOLD = 7.2   # rad, a crease's angular half-width (build5: 0.09; user picked 0.06)
# crease width: build5's fixed angular width made each crease's real width grow with r, so the grooves took ~2/3 of
# the circumference all the way out. Past R_W the angular width now shrinks as (R_W / r)^K, so the creases still
# narrow into the centre but stay slim further out (K = 0: build5's behaviour)
R_W = 1.5
CREASE_TAPER = float(os.environ.get("ANUS_CREASE_TAPER", "0.0"))
N_SPOKE = 288                    # build5: 72 (5 deg apart) - the slimmer creases need finer spokes to stay clean
# rise: build5's 2 mm funnel profile turned upside down - the surface climbs towards the centre and fades to nothing
# at the outer ring
RISE = float(os.environ.get("ANUS_RISE", "0.2"))   # mm (0.5 and up read as an unnatural dome)
CREASE_SCALE = float(os.environ.get("ANUS_CREASE_SCALE", "2.0"))   # x build5's crease depth (user's pick; 1x read too faint
                                 # once the anus follows the curved cleft and cheek walls)
ENTRANCE_DIP = 0.0               # mm the pole drops below the rise (build5: 0.6 - on a rise it read as an abrupt sink)
# the user's sketch: gentle shoulders rising towards the centre that roll over and curve down into a narrow plunge.
# D * (1 - sqrt(r / R_DIP))^2 starts with zero slope at R_DIP (no corner), steepens inwards and is near vertical at
# the centre
DIP = float(os.environ.get("ANUS_DIP", "1.5"))      # mm at the centre
R_DIP = 2.5                      # mm where the roll-over begins


def crease(theta, r=0.0):
    sig = SIG * min(1.0, R_W / max(r, 1e-9)) ** CREASE_TAPER
    s = 0.0
    for i, w in enumerate(FOLD_W):
        dt = (theta - 2 * math.pi * i / 12 + math.pi) % (2 * math.pi) - math.pi
        s += w * math.exp(-dt * dt / (2 * sig * sig))
    return s


def fold_amp(r):
    t = min(max((r - 0.9) / (R_FOLD - 0.9), 0.0), 1.0)
    return 0.42 * math.sin(math.pi * t) ** 0.8


def lat(r, theta):
    """lateral offset (mm) of the layout; build5_passage: collapses to a closed slit towards the R_HOLE edge loop"""
    b = r * math.sin(theta) * SQUEEZE * (0.45 + 0.55 * ss(0.9, 5.0, r))
    w = ss(R_HOLE, R_COLL, r)
    return b * w + (LIP_GAP / 2) * math.sin(theta) * (1 - w)


def base_point(r, theta):
    # adapted: same a and (squeezed) lateral offsets as before, but the lateral one is now arc length over the skin,
    # so the sides climb the cheek walls instead of being projected flat across the cleft
    return surf_point(r * math.cos(theta), lat(r, theta))


rings = []; new_faces = []; info = {}
for r in RINGS:
    S = N_SPOKE if r < 7.0 else {7.6: N_SPOKE // 2, 8.6: N_SPOKE // 4}[r]   # 2:1 steps to the join with the skin
    ring = []
    for k in range(S):
        th = 2 * math.pi * k / S
        v = bm.verts.new(base_point(r, th)); info[v] = (r, th); ring.append(v)
    rings.append(ring)
mouth = rings[0]                                  # build5_passage: open - no pole, no fan
for ri in range(len(rings) - 1):
    A_, B_ = rings[ri], rings[ri + 1]
    if len(A_) == len(B_):
        S = len(A_)
        for k in range(S):
            new_faces.append(bm.faces.new((A_[k], B_[k], B_[(k + 1) % S], A_[(k + 1) % S])))
    else:
        nA, nB = len(A_), len(B_)
        for k in range(nB):
            new_faces.append(bm.faces.new((A_[2 * k], B_[k], A_[2 * k + 1])))
            new_faces.append(bm.faces.new((A_[2 * k + 1], B_[k], B_[(k + 1) % nB], A_[(2 * k + 2) % nA])))
bm.normal_update()
if sum(f.normal.dot(n) for f in new_faces) < 0:
    for f in new_faces:
        f.normal_flip()
bm.normal_update()
base_pos = {v: v.co.copy() for v in info}
# level along the cleft (the "cliff" fix): the narrow fill sits low near the original pit and then climbs steeply in a
# narrow band at the back edge to meet the cleft floor rising toward the coccyx. Smooth the heights ALONG the cleft
# only (edge weights (da / length)^2, so nothing diffuses across the cleft and the walls keep their shape), on the cleft
# floor (fading out up the walls), with the outer ring held fixed: the climb spreads over the whole anus.
# ANUS_LEVEL = smoothing iterations (0 = off)
LEVEL_ITERS = int(os.environ.get("ANUS_LEVEL", "3000"))
# build5_passage: the levelling runs on build5_rise's closed layout - the removed inner rings and pole come back as
# virtual points, and edge weights use the uncollapsed lateral offsets. Without the pole the lips' only links are
# along the dense mouth loop, which converges far slower: the middle of each lip lagged ~1 mm below the tips.
def lat_closed(r, theta):
    return r * math.sin(theta) * SQUEEZE * (0.45 + 0.55 * ss(0.9, 5.0, r))


_vl = list(info); _vi = {v: i for i, v in enumerate(_vl)}
_lay = [(r * math.cos(th), lat_closed(r, th)) for v, (r, th) in ((v, info[v]) for v in _vl)]
_h0 = [to_local(base_pos[v])[2] * 1000 for v in _vl]
_Eset = {tuple(sorted((_vi[e.verts[0]], _vi[e.verts[1]]))) for f in new_faces for e in f.edges
         if e.verts[0] in _vi and e.verts[1] in _vi}
_vring = []                                      # virtual rings inside R_HOLE (build5_rise: 0.06, 0.15 mm) + pole
for r in [0.06, 0.15]:
    if r < R_HOLE - 1e-9:
        _vring.append([])
        for k in range(N_SPOKE):
            th = 2 * math.pi * k / N_SPOKE
            _vring[-1].append(len(_lay)); _lay.append((r * math.cos(th), lat_closed(r, th)))
            _h0.append(to_local(surf_point(r * math.cos(th), lat_closed(r, th)))[2] * 1000)
_vring.append([_vi[v] for v in mouth])
for ri, rg in enumerate(_vring):
    for k in range(N_SPOKE):
        if ri < len(_vring) - 1:
            _Eset.add(tuple(sorted((rg[k], rg[(k + 1) % N_SPOKE]))))
            _Eset.add(tuple(sorted((rg[k], _vring[ri + 1][k]))))
_pole = len(_lay); _lay.append((0.0, 0.0)); _h0.append(to_local(surf_point(0.0, 0.0))[2] * 1000)
for k in range(N_SPOKE):
    _Eset.add(tuple(sorted((_vring[0][k], _pole))))
_E = np.array(sorted(_Eset)); _L = np.array(_lay)
_d = _L[_E[:, 1]] - _L[_E[:, 0]]
_w = (_d[:, 0] ** 2) / np.maximum((_d ** 2).sum(1), 1e-12)          # 1 for an edge along the cleft, 0 across it
_h0 = np.array(_h0); _h = _h0.copy()
_mob = np.array([(1 - ss(3.0, 6.0, abs(_lay[i][1]))) * (info[_vl[i]][0] < 8.6 if i < len(_vl) else 1.0)
                 for i in range(len(_lay))])
_W = np.zeros(len(_lay)); np.add.at(_W, _E[:, 0], _w); np.add.at(_W, _E[:, 1], _w)
for _ in range(LEVEL_ITERS):
    _acc = np.zeros(len(_lay))
    np.add.at(_acc, _E[:, 0], _w * (_h[_E[:, 1]] - _h[_E[:, 0]])); np.add.at(_acc, _E[:, 1], _w * (_h[_E[:, 0]] - _h[_E[:, 1]]))
    _h += 0.5 * _mob * _acc / np.maximum(_W, 1e-9)
_h = _h[:len(_vl)]; _h0 = _h0[:len(_vl)]
report["level_change_mm"] = [round(float((_h - _h0).min()), 2), round(float((_h - _h0).max()), 2)]
for v in _vl:
    base_pos[v] = base_pos[v] + n * ((_h[_vi[v]] - _h0[_vi[v]]) / 1000)
# relief along the (smooth) base normal: radial creases; the closed centre keeps rising with the surface (no sink at
# the end). The rise goes straight out of the body (along n): along the surface normal it would push the steep
# cheek walls sideways into the cleft.
for v, (r, th) in info.items():
    relief = CREASE_SCALE * fold_amp(r) * (0.15 - crease(th, r))
    rise = RISE * (1 - ss(0, R_FOLD, r)) - DIP * (1 - min(r / R_DIP, 1.0) ** 0.5) ** 2
    v.co = base_pos[v] + v.normal * (relief / 1000) + n * (rise / 1000)
    v[pig_layer] = min(1.0, (1 - ss(3.0, 8.5, r)) * (1 + 0.1 * crease(th, r) * fold_amp(r) / 0.42))
relief_off = {v: v.co - base_pos[v] for v in info}   # corridor (5b): the relief is kept out of the smoothing

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
for v in bm.verts:                                 # only the midline spokes of the new rings sit on the mirror plane
    if (abs(math.sin(info[v][1])) < 1e-9) if v in info else abs(v.co.x) < 2e-5:   # (by distance, the dense
        v.co.x = 0.0                                   # centre rings would fold)


# ======================= 5b. corridor: spread the climb toward the coccyx along the cleft floor =======================
# build5_rise_corridor: behind the anus the cleft floor climbs steeply (the original rises ~27-30 deg toward the
# coccyx, and the Taubin pass leaves it ~1.3 mm higher), so the back half of the anus tilted up and met a "cliff".
# Here the finished mesh (anus + bridge strip + cleft skin) gets one along-cleft smoothing of its heights (along n)
# over a corridor from CORR_A[0] to CORR_A[1] mm, so the climb spreads over the whole corridor:
#   * energy: integral of (dh/du_t)^2 over the surface (linear FEM on the real triangles, u_t = the cleft direction in
#     each triangle's plane): smoothing ALONG the cleft only, independent of the mesh density (the dense anus and the
#     coarse skin get the same slope per mm);
#   * + CORR_LAT x integral of (d(change)/d lateral)^2: the change is laterally coherent, so the cross-section of the
#     cleft (the groove, the walls) keeps its shape and nothing diffuses across the cleft;
#   * + anchoring by mobility (1 = free, 0 = fixed): full on the floor for |b| < CORR_B[0], none from |b| = CORR_B[1]
#     (the cheek walls stay put) and fading to zero over CORR_FADE mm at both corridor ends (no seams);
#   * the anus relief (creases, rise, plunge) is taken off before and put back after, so its look is unchanged.
# ANUS_CORRIDOR = max solver (PCG) iterations, 0 = off (as build5_rise)
CORR_ITERS = int(os.environ.get("ANUS_CORRIDOR", "3000"))   # on: the user's "cliff" behind the anus
CORR_A = [float(x) for x in os.environ.get("ANUS_CORR_A", "-10,28").split(",")]        # mm along the cleft (u)
CORR_FADE = [float(x) for x in os.environ.get("ANUS_CORR_FADE", "6,8").split(",")]     # mm: front, back end fade
CORR_B = [float(x) for x in os.environ.get("ANUS_CORR_B", "3,7").split(",")]           # mm |b|: full .. none
CORR_LAT = float(os.environ.get("ANUS_CORR_LAT", "0.5"))      # lateral coherence of the change
CORR_ANCHOR = float(os.environ.get("ANUS_CORR_ANCHOR", "0.03"))    # 1/mm^2: anchoring strength scale
CORR_STIFF = float(os.environ.get("ANUS_CORR_STIFF", "1"))  # along-cleft stiffness of the anus faces (>1: flatter)
CORR_ALONGP = float(os.environ.get("ANUS_CORR_ALONGP", "0"))   # along term x mobility^p (0: everywhere in the set)
CORR_HR = [float(x) for x in os.environ.get("ANUS_CORR_HR", "0,0").split(",")]   # mm above the floor: full .. none
CORR_ALONG_HR = [float(x) for x in os.environ.get("ANUS_CORR_ALONG_HR", "0,0").split(",")]  # along term: floor only
CORR_MODE = os.environ.get("ANUS_CORR_MODE", "floor")      # "floor" (1D floor profile) or "fem" (2D surface energy)
CORR_ORDER = int(os.environ.get("ANUS_CORR_ORDER", "2"))   # floor mode: 1 = straight ramp, 2 = least bending
CORR_FB = [float(x) for x in os.environ.get("ANUS_CORR_FB", "-2,-1,0,1,2").split(",")]   # floor sample lines, mm
CORR_LATSM = float(os.environ.get("ANUS_CORR_LATSM", "0.75"))   # grid mode: blur of the change across the lines, mm


# --- CORRIDOR_SOLVER_BEGIN (pure numpy; also exec'd by the offline tuning script) ---
def corridor_delta(Pb, tris, tri_anus, fixed_extra, prm):
    """Pb: (N,3) local coords in mm (a, b, h) of the relief-free surface; tris: (T,3) loop triangles of the whole mesh;
    tri_anus: (T,) bool; fixed_extra: (N,) bool (boundary verts). Returns (delta_h (N,) mm, info dict)."""
    import numpy as np
    def _ss(e0, e1, x):
        t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
        return t * t * (3 - 2 * t)
    a0, a1 = prm["A"]; f0, f1 = prm["FADE"]; b0, b1 = prm["B"]
    Nv = len(Pb); a, b, h = Pb[:, 0], Pb[:, 1], Pb[:, 2]
    in_set = (a > a0 - 1.5) & (a < a1 + 1.5) & (np.abs(b) < b1 + 1.5) & (h > -12) & (h < 25) & ~fixed_extra
    mob = _ss(a0, a0 + f0, a) * (1 - _ss(a1 - f1, a1, a)) * (1 - _ss(b0, b1, np.abs(b)))
    # height above the cleft floor (the floor = lowest point near the midline, per mm of a)
    mid = in_set & (np.abs(b) < 1.5)
    bins = np.arange(np.floor(a0) - 2, np.ceil(a1) + 3)
    fl = np.array([h[mid & (np.abs(a - c) < 0.75)].min() if (mid & (np.abs(a - c) < 0.75)).any() else np.nan
                   for c in bins])
    okb = ~np.isnan(fl); fl = np.interp(bins, bins[okb], fl[okb])
    fl = np.convolve(np.pad(fl, 1, mode='edge'), np.ones(3) / 3, 'valid')
    hr = h - np.interp(a, bins, fl)
    if prm.get("HR", [0, 0])[1] > 0:              # also fade out by the height above the cleft floor (the walls)
        mob = mob * (1 - _ss(prm["HR"][0], prm["HR"][1], hr))
    mob = np.where(in_set, mob, 0.0)
    tri_ok = in_set[tris].all(1)
    T = tris[tri_ok]; ta = tri_anus[tri_ok]
    # a vertex of the set that touches a triangle outside it is held fixed (keeps the outer edge continuous)
    touch_out = np.zeros(Nv, bool); np.logical_or.at(touch_out, tris[~tri_ok].ravel(), True)
    p0, p1, p2 = Pb[T[:, 0]], Pb[T[:, 1]], Pb[T[:, 2]]
    Nr = np.cross(p1 - p0, p2 - p0); A2 = np.linalg.norm(Nr, axis=1)
    ok = A2 > 1e-10
    T, ta, p0, p1, p2, Nr, A2 = T[ok], ta[ok], p0[ok], p1[ok], p2[ok], Nr[ok], A2[ok]
    Nh = Nr / A2[:, None]
    G = np.stack([np.cross(Nh, e) / A2[:, None] for e in (p2 - p1, p0 - p2, p1 - p0)], 1)  # (T,3,3) hat gradients
    ut = np.array([1.0, 0.0, 0.0])[None, :] - Nh * Nh[:, :1]
    ul = np.linalg.norm(ut, axis=1); good = ul > 0.05
    ut = ut / np.maximum(ul, 1e-9)[:, None]; vt = np.cross(Nh, ut)
    gu = np.einsum('tkd,td->tk', G, ut) * good[:, None]
    gv = np.einsum('tkd,td->tk', G, vt)
    Ar = A2 / 2; st = np.where(ta, prm["STIFF"], 1.0)
    if prm.get("ALONGP", 0) > 0:                  # the along-cleft term only where the surface may move (the floor):
        st = st * mob[T].mean(1) ** prm["ALONGP"]  # the converging cheek walls behind the anus don't drive it
    if prm.get("ALONG_HR", [0, 0])[1] > 0:        # ... or only on the floor itself (up to ALONG_HR mm above it)
        ma = _ss(a0, a0 + f0, a) * (1 - _ss(a1 - f1, a1, a))
        st = st * ((1 - _ss(prm["ALONG_HR"][0], prm["ALONG_HR"][1], hr)) * ma)[T].mean(1)
    Wa = (Ar * st)[:, None, None] * gu[:, :, None] * gu[:, None, :]
    Wl = Ar[:, None, None] * (gv[:, :, None] * gv[:, None, :] +
                              (~good)[:, None, None] * np.einsum('tkd,tjd->tkj', G, G))
    R = np.repeat(T[:, :, None], 3, 2).ravel(); C = np.repeat(T[:, None, :], 3, 1).ravel()
    Va = Wa.ravel(); Vk = (Wa + prm["LAT"] * Wl).ravel()
    M = np.bincount(T.ravel(), weights=np.repeat(Ar / 3, 3), minlength=Nv)
    used = M > 0
    fixed = ~in_set | touch_out | ~used | (mob < 1e-3)
    kap = prm["ANCHOR"] * (1 - mob) / np.maximum(mob, 1e-3)
    Dg = kap * M + 1e-9 * M
    def mv(x, vals=Vk):
        return np.bincount(R, weights=vals * x[C], minlength=Nv)
    rhs = -mv(h, Va); rhs[fixed] = 0.0
    diag = np.bincount(R, weights=Vk * (R == C), minlength=Nv) + Dg
    pre = np.where(fixed, 0.0, 1.0 / np.maximum(diag, 1e-30))
    x = np.zeros(Nv); r = rhs.copy(); z = pre * r; p = z.copy(); rz = r @ z; r0 = np.sqrt(r @ r) + 1e-30
    it = 0
    for it in range(1, prm["ITERS"] + 1):
        Ap = mv(p) + Dg * p; Ap[fixed] = 0.0
        al = rz / (p @ Ap); x += al * p; r -= al * Ap
        if np.sqrt(r @ r) < 1e-9 * r0:
            break
        z = pre * r; rz2 = r @ z; p = z + (rz2 / rz) * p; rz = rz2
    x[fixed] = 0.0
    return x, {"iters": it, "rel_res": float(np.sqrt(r @ r) / r0), "n_free": int((~fixed).sum()),
               "n_tris": int(len(T)), "dh_min": float(x.min()), "dh_max": float(x.max())}


def corridor_floor_delta(Pb, Fa, Fh, prm):
    """mode "floor": smooth the cleft floor's height profile Fh(Fa) (mm, sampled on the relief-free surface near the
    midline) along the cleft over the corridor, then carry the change to every vertex, fading out up the walls.
    The whole floor moves together at each a, so the cross-section of the cleft (groove, walls) keeps its shape.
    1D energy: ORDER 1 = slope^2 (straight ramp), ORDER 2 = curvature^2 (no bends), + anchoring by mobility."""
    import numpy as np
    def _ss(e0, e1, x):
        t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
        return t * t * (3 - 2 * t)
    a0, a1 = prm["A"]; f0, f1 = prm["FADE"]; b0, b1 = prm["B"]
    dk = 0.25; ak = np.arange(a0 - 4, a1 + 4 + 1e-9, dk); Fk = np.interp(ak, Fa, Fh); K = len(ak)
    mk = _ss(a0, a0 + f0, ak) * (1 - _ss(a1 - f1, a1, ak))
    fixed = mk < 1e-3
    order = int(prm.get("ORDER", 2))
    kap = prm["ANCHOR"] * (1 - mk) / np.maximum(mk, 1e-3)
    Dm = np.diff(np.eye(K), order, axis=0) / dk ** order      # finite differences
    Amat = Dm.T @ Dm * dk + np.diag(kap * dk)
    rhs = np.diag(kap * dk) @ Fk
    free = ~fixed
    x = Fk.copy()
    x[free] = np.linalg.solve(Amat[np.ix_(free, free)], rhs[free] - Amat[np.ix_(free, fixed)] @ Fk[fixed])
    dk_ = x - Fk
    a, b, h = Pb[:, 0], Pb[:, 1], Pb[:, 2]
    w = (1 - _ss(b0, b1, np.abs(b))) * ((h > -12) & (h < 25))
    if prm.get("HR", [0, 0])[1] > 0:              # fade by the height above the floor as well (the walls)
        w = w * (1 - _ss(prm["HR"][0], prm["HR"][1], h - np.interp(a, Fa, Fh)))
    dh = np.interp(a, ak, dk_, left=0.0, right=0.0) * w
    return dh, {"dh_min": float(dh.min()), "dh_max": float(dh.max()),
                "floor_change_mm": [[float(round(q, 1)), float(round(np.interp(q, ak, dk_), 2))]
                                    for q in range(int(a0), int(a1) + 1, 2)]}


def corridor_grid_delta(Ga, Gb, Gh, prm):
    """mode "grid": the along-cleft smoothing of build5_rise's levelling, done line by line on a regular (a, b) grid
    of the relief-free surface (Gh[i, j] = height at a = Ga[i], b = Gb[j], NaN = no surface), so it does not depend
    on the mesh density (dense anus / coarse skin). Each line b: 1D energy (ORDER 1: slope^2 -> straight ramps,
    ORDER 2: curvature^2 -> least bending) + anchoring sum kappa (h - h0)^2, kappa = ANCHOR (1 - m) / m from the
    mobility m = (corridor ends fade) x (lateral fade |b|: B) x (fade with the height above the floor: HR, the walls).
    The change is then blurred across the lines (LATSM mm) so neighbouring lines stay together. Returns dH[i, j]."""
    import numpy as np
    def _ss(e0, e1, x):
        t = np.clip((x - e0) / (e1 - e0), 0.0, 1.0)
        return t * t * (3 - 2 * t)
    a0, a1 = prm["A"]; f0, f1 = prm["FADE"]; b0, b1 = prm["B"]
    da = Ga[1] - Ga[0]; NA = len(Ga)
    order = int(prm.get("ORDER", 1))
    Dm = np.diff(np.eye(NA), order, axis=0) / da ** order
    S = Dm.T @ Dm * da
    mid = np.abs(Gb) <= 1.5
    floor = np.nanmin(np.where(mid[None, :], Gh, np.nan), axis=1)
    ok = ~np.isnan(floor); floor = np.interp(Ga, Ga[ok], floor[ok])
    ma = _ss(a0, a0 + f0, Ga) * (1 - _ss(a1 - f1, a1, Ga))
    dH = np.zeros_like(Gh)
    for j, bj in enumerate(Gb):
        h0 = Gh[:, j]; have = ~np.isnan(h0)
        m = ma * (1 - _ss(b0, b1, abs(bj)))
        if prm.get("HR", [0, 0])[1] > 0:
            m = m * (1 - _ss(prm["HR"][0], prm["HR"][1], np.where(have, h0 - floor, 1e9)))
        m = np.where(have, m, 0.0)
        fixed = m < 1e-3
        if fixed.all():
            continue
        hh = np.where(have, h0, 0.0)
        kap = prm["ANCHOR"] * (1 - m) / np.maximum(m, 1e-3)
        Am = S + np.diag(kap * da); rhs = kap * da * hh
        fr = ~fixed
        x = hh.copy()
        x[fr] = np.linalg.solve(Am[np.ix_(fr, fr)], rhs[fr] - Am[np.ix_(fr, fixed)] @ hh[fixed])
        dH[:, j] = np.where(have, x - hh, 0.0)
    sig = prm.get("LATSM", 0.75)
    if sig > 0:                                   # blur across the lines
        db = Gb[1] - Gb[0]; k = np.arange(-int(3 * sig / db), int(3 * sig / db) + 1) * db
        g = np.exp(-k * k / (2 * sig * sig)); g /= g.sum()
        dH = np.apply_along_axis(lambda r: np.convolve(np.pad(r, len(k) // 2, mode='edge'), g, 'valid'), 1, dH)
    return dH


def grid_sample(Ga, Gb, G, a, b):
    """bilinear interpolation of G (on the Ga x Gb grid) at points (a, b); 0 outside"""
    import numpy as np
    fa = (a - Ga[0]) / (Ga[1] - Ga[0]); fb = (b - Gb[0]) / (Gb[1] - Gb[0])
    inside = (fa >= 0) & (fa <= len(Ga) - 1) & (fb >= 0) & (fb <= len(Gb) - 1)
    ia = np.clip(np.floor(fa).astype(int), 0, len(Ga) - 2); ib = np.clip(np.floor(fb).astype(int), 0, len(Gb) - 2)
    ta = np.clip(fa - ia, 0, 1); tb = np.clip(fb - ib, 0, 1)
    v = (G[ia, ib] * (1 - ta) * (1 - tb) + G[ia + 1, ib] * ta * (1 - tb) + G[ia, ib + 1] * (1 - ta) * tb
         + G[ia + 1, ib + 1] * ta * tb)
    return np.where(inside, v, 0.0)
# --- CORRIDOR_SOLVER_END ---


if CORR_ITERS > 0 or os.environ.get("ANUS_CORR_DUMP"):
    bm.verts.ensure_lookup_table(); bm.verts.index_update()
    _P = np.array([tuple(v.co - relief_off[v]) if v in relief_off else tuple(v.co) for v in bm.verts]) - np.array(O)
    _Pb = np.stack([_P @ np.array(u), _P @ np.array(bv), _P @ np.array(n)], 1) * 1000
    _lt = bm.calc_loop_triangles()
    _tris = np.array([[l.vert.index for l in tr] for tr in _lt])
    _nf = set(new_faces)
    _tan = np.array([tr[0].face in _nf for tr in _lt])
    _mset = set(mouth)                             # build5_passage: the mouth loop moves with the floor (not held)
    _bnd = np.array([v.is_boundary and v not in _mset for v in bm.verts])
    _prm = {"A": CORR_A, "FADE": CORR_FADE, "B": CORR_B, "LAT": CORR_LAT, "ANCHOR": CORR_ANCHOR,
            "STIFF": CORR_STIFF, "ITERS": CORR_ITERS, "ALONGP": CORR_ALONGP, "HR": CORR_HR,
            "ALONG_HR": CORR_ALONG_HR, "MODE": CORR_MODE, "ORDER": CORR_ORDER, "FB": CORR_FB,
            "LATSM": CORR_LATSM}
    if os.environ.get("ANUS_CORR_DUMP"):          # for offline tuning: the relief-free surface + final positions
        _F = np.array([tuple(v.co) for v in bm.verts]) - np.array(O)
        np.savez(os.environ["ANUS_CORR_DUMP"], Pb=_Pb, Pf=np.stack([_F @ np.array(u), _F @ np.array(bv),
                 _F @ np.array(n)], 1) * 1000, tris=_tris, tri_anus=_tan, bnd=_bnd)
    if CORR_ITERS > 0 and CORR_MODE == "floor":
        # floor profile on the relief-free surface: mean height of the lines b = -FB .. FB mm
        _bt = BVHTree.FromPolygons([Vector(p_) for p_ in _Pb], _tris.tolist())
        _Fa = np.arange(CORR_A[0] - 6, CORR_A[1] + 6.01, 0.25); _Fh = []
        for a_ in _Fa:
            hs = [_bt.ray_cast(Vector((a_, b_, -30.0)), Vector((0.0, 0.0, 1.0)), 60.0)[0] for b_ in CORR_FB]
            hs = [q.z for q in hs if q is not None]
            _Fh.append(np.mean(hs) if hs else np.nan)
        _Fh = np.array(_Fh); _ok = ~np.isnan(_Fh); _Fh = np.interp(_Fa, _Fa[_ok], _Fh[_ok])
        _dh, _ci = corridor_floor_delta(_Pb, _Fa, _Fh, _prm)
    elif CORR_ITERS > 0 and CORR_MODE == "grid":
        _bt = BVHTree.FromPolygons([Vector(p_) for p_ in _Pb], _tris.tolist())
        _Ga = np.arange(CORR_A[0] - 2, CORR_A[1] + 2.01, 0.25); _Gb = np.arange(-CORR_B[1] - 1, CORR_B[1] + 1.01, 0.25)
        _Gh = np.full((len(_Ga), len(_Gb)), np.nan)
        for i_, a_ in enumerate(_Ga):
            for j_, b_ in enumerate(_Gb):
                q = _bt.ray_cast(Vector((a_, b_, -30.0)), Vector((0.0, 0.0, 1.0)), 60.0)[0]
                if q is not None:
                    _Gh[i_, j_] = q.z
        _dH = corridor_grid_delta(_Ga, _Gb, _Gh, _prm)
        # a vertex takes the change of the surface it lies on (the first surface along n: not an overhang behind it)
        _hs = grid_sample(_Ga, _Gb, np.nan_to_num(_Gh, nan=-99.0), _Pb[:, 0], _Pb[:, 1])
        _on = np.abs(_hs - _Pb[:, 2]) < 0.6
        _dh = grid_sample(_Ga, _Gb, _dH, _Pb[:, 0], _Pb[:, 1]) * _on
        _ci = {"dh_min": float(_dh.min()), "dh_max": float(_dh.max()), "verts_moved": int((np.abs(_dh) > 1e-4).sum()),
               "off_surface_skipped": int(((~_on) & (np.abs(grid_sample(_Ga, _Gb, _dH, _Pb[:, 0], _Pb[:, 1])) > 1e-3)).sum())}
    elif CORR_ITERS > 0:
        _dh, _ci = corridor_delta(_Pb, _tris, _tan, _bnd, _prm)
    if CORR_ITERS > 0:
        for v in bm.verts:
            if _dh[v.index] != 0.0:
                v.co += n * (_dh[v.index] / 1000)
        report["corridor"] = {k: (round(x_, 4) if isinstance(x_, float) else x_) for k, x_ in _ci.items()}
        report["corridor_params"] = _prm
bm.normal_update()
ring0 = [v.co.copy() for v in rings[RINGS.index(0.9)]]   # A (and so the canal axis, J) as before: the 0.9 mm ring
A = sum(ring0, Vector()) / len(ring0)
hole_loop = [v.co.copy() for v in mouth]                 # build5_passage: the canal's mouth = the skin's edge loop
H = sum(hole_loop, Vector()) / len(hole_loop)
widths = [abs(v.co.x) for v in rings[-1]]
report["new_verts"] = len(info); report["degenerate_faces"] = sum(1 for f in bm.faces if f.calc_area() < 1e-12)
report["anus_centre_A"] = r4(A); report["boundary_edges_total"] = sum(e.is_boundary for e in bm.edges)
report["pucker_outer_lateral_halfwidth_mm"] = round(max(widths) * 1000, 1)
report["pucker_surface_span_mm"] = {
    "AP": round((rings[-1][0].co - rings[-1][len(rings[-1]) // 2].co).length * 1000, 1),
    "lateral_over_surface": round(sum((rings[-1][k].co - rings[-1][k + 1].co).length
                                      for k in range(len(rings[-1]) // 4, 3 * len(rings[-1]) // 4)) * 1000, 1)}
bm.verts.ensure_lookup_table()
refined_ids = [v.index for v in bm.verts if in_refine(v.co)]
bm.to_mesh(me); me.update(); bm.free(); base.free()
# (build5_rise: ported from build8) the body's custom split normals don't fit the refined / rebuilt cleft: automatic
# normals there
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
cols = {"maintex": ("Lit colour", (0.98, 0.66, 0.71)), "darktex": ("Shadow colour", (0.84, 0.55, 0.61))}
x0, y0 = shd.location.x - 520, shd.location.y + 420
for i, (sock, (lab, rgb)) in enumerate(cols.items()):
    out_s = gen.outputs[f"Body {sock}"]; in_s = shd.inputs[f"Body {sock}"]
    for lk in list(out_s.links):
        if lk.to_socket == in_s:
            nt.links.remove(lk)
    col = nt.nodes.new("ShaderNodeRGB"); col.name = f"AnusPig_{sock}_Colour"; col.label = lab
    col.outputs[0].default_value = (*srgb2lin(rgb), 1.0)
    mix = nt.nodes.new("ShaderNodeMix"); mix.name = f"AnusPig_{sock}_Mix"; mix.data_type = 'RGBA'
    assert mix.inputs[6].type == 'RGBA' and mix.inputs[7].type == 'RGBA'
    nt.links.new(stren.outputs[0], mix.inputs[0]); nt.links.new(out_s, mix.inputs[6])
    nt.links.new(col.outputs[0], mix.inputs[7]); nt.links.new(mix.outputs[2], in_s)
    col.location = (x0, y0 - 220 * i); mix.location = (x0 + 220, y0 - 220 * i)
    col.parent = frame; mix.parent = frame
attr.location = (x0 - 420, y0); stren.location = (x0 - 210, y0); attr.parent = frame; stren.parent = frame

# (build5_rise: ported from build8; diagnosed on this build - with the map off the spiky shadow edge and the
# streaky crease ends both disappear)
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


# build5_passage: the canal starts ON the skin's edge loop (hole_loop: the same N_SPOKE points, in the same places)
# and is collapsed at rest along its whole length (item 3): the closed slit is carried straight up for CANAL_COLL mm,
# then the lumen is a closed STAR - an AP slit that lengthens, then lateral and diagonal arms (8 in all) growing in
# between, the tissue between neighbouring arms being the anal columns. Walls LUMEN_HW mm from the arm's midline.
# Over the top 40 % it eases out to the round junction ring the ampulla starts from. Every arm owns a fixed block of
# N_SPOKE / 8 points (tip in the middle), so points never slide across an arm between rings, and point k always
# sits near angle -pi + 2 pi k / N_SPOKE (the Open shape key can spread them round a circle).
CANAL_COLL = float(os.environ.get("ANUS_CANAL_COLL", "2.0"))   # mm the skin's slit runs straight up
LUMEN_HW = 0.01                                                  # mm, half the gap between touching walls
N_ARM = 8; SEC = N_SPOKE // N_ARM
Lc = (J - H).length; _Lmm = Lc * 1000
e1, e2 = frame_at(d)
_psi = [math.atan2((p_ - H).dot(e2), (p_ - H).dot(e1)) for p_ in hole_loop]
_turn = sum((_psi[(k + 1) % N_SPOKE] - _psi[k] + math.pi) % (2 * math.pi) - math.pi for k in range(N_SPOKE))
_sg = 1 if _turn > 0 else -1                     # does the skin's spoke order run with the canal frame's angle?
_k0 = round((_psi[0] + math.pi) / (2 * math.pi / N_SPOKE)) % N_SPOKE
_spoke = [(_sg * (kk - _k0)) % N_SPOKE for kk in range(N_SPOKE)]   # canal point k' sits on skin spoke _spoke[k']


def canal_rc(th, t):   # junction ring: 4 mm with a light 8-fold column ripple (as build5_rise)
    return (3.2 + 0.8 * ss(0.5, 1.0, t)) * (1 + 0.15 * ss(0.4, 0.8, t) * math.cos(8 * th)) / 1000


def arm_lengths(s_mm):
    """mm; arm i points at angle -pi + i pi / 4 (i = 2, 6: front / back, i = 0, 4: sides, odd: diagonals)"""
    ap = 0.3 + 2.2 * ss(CANAL_COLL, 7.0, s_mm)
    side = 2.0 * ss(5.0, 10.0, s_mm)
    diag = 1.8 * ss(7.0, 12.0, s_mm)
    return [side, diag, ap, diag, side, diag, ap, diag]


def star_2d(s_mm):
    """the resting lumen at s: N_SPOKE (x, y) in mm in the canal frame (x along e1, y along e2)"""
    Ls = arm_lengths(s_mm); al = [-math.pi + i * math.pi / 4 for i in range(N_ARM)]
    out = [None] * N_SPOKE
    for i in range(N_ARM):
        ph = np.linspace(al[i] - math.pi / 8, al[i] + math.pi / 8, 1601)
        r = np.full_like(ph, LUMEN_HW)
        for j in range(N_ARM):
            dp = (ph - al[j] + math.pi) % (2 * math.pi) - math.pi
            arm = np.where(np.cos(dp) > 0, np.minimum(Ls[j], LUMEN_HW / np.maximum(np.abs(np.sin(dp)), 1e-9)), 0.0)
            r = np.maximum(r, arm)
        P2 = np.stack([r * np.cos(ph), r * np.sin(ph)], 1)
        cum = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P2, axis=0), axis=1))])
        for jj in range(SEC):
            u = cum[-1] * jj / SEC
            x = float(np.interp(u, cum, P2[:, 0])); y = float(np.interp(u, cum, P2[:, 1]))
            out[(i * SEC - SEC // 2 + jj) % N_SPOKE] = (x, y)
    return out


def canal_ring(s_mm, N):
    """N_SPOKE points at s mm up the canal, then every (N_SPOKE / N)-th: the 2:1 steps keep their angles"""
    t = s_mm / _Lmm; ctr = H.lerp(J, t)
    w_ex = 1 - ss(CANAL_COLL, 3.5, s_mm)          # the skin's own slit, carried up
    w_rd = ss(0.6 * _Lmm, _Lmm, s_mm)             # easing out to the round junction ring
    st = star_2d(s_mm); pts = []
    for kk in range(N_SPOKE):
        th = -math.pi + 2 * math.pi * kk / N_SPOKE
        p_ = ctr + e1 * (st[kk][0] / 1000) + e2 * (st[kk][1] / 1000)
        if w_rd > 0:
            p_ = p_ * (1 - w_rd) + (ctr + (e1 * math.cos(th) + e2 * math.sin(th)) * canal_rc(th, t)) * w_rd
        if w_ex > 0:
            p_ = (hole_loop[_spoke[kk]] + (ctr - H)) * w_ex + p_ * (1 - w_ex)
        pts.append(p_)
    return pts[::N_SPOKE // N]


CANAL_S = [0.15, 0.35, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.25, 5.0, 5.75, 6.5, 7.25, 8.0, 9.0, 10.0, 11.0, 12.0,
           13.5, 15.0, 16.5, 18.0, 19.5, 21.0]
canal_loops = [[hole_loop[_spoke[kk]] for kk in range(N_SPOKE)]]
for s_mm in CANAL_S:
    canal_loops.append(canal_ring(s_mm, N_SPOKE))
for f_, N in [(0.90, N_SPOKE // 2), (0.95, N_SPOKE // 4)]:
    canal_loops.append(canal_ring(f_ * _Lmm, N))
assert N_SPOKE // 8 == NS
J_ring = canal_ring(_Lmm, NS)                    # = canal_r(th, 1.0) around J, as the ampulla expects
canal_loops.append(J_ring)


def build_tube2(name, loops_p):
    """tube through point loops; consecutive loops have equal counts (quads) or halve (2:1 fans)"""
    tb = bmesh.new(); lv = [[tb.verts.new(p_) for p_ in lp] for lp in loops_p]
    for A_, B_ in zip(lv[:-1], lv[1:]):
        if len(A_) == len(B_):
            S = len(A_)
            for k in range(S):
                tb.faces.new((A_[k], B_[k], B_[(k + 1) % S], A_[(k + 1) % S]))
        else:
            assert len(A_) == 2 * len(B_)
            nA, nB = len(A_), len(B_)
            for k in range(nB):
                tb.faces.new((A_[2 * k], B_[k], A_[2 * k + 1]))
                tb.faces.new((A_[2 * k + 1], B_[k], B_[(k + 1) % nB], A_[(2 * k + 2) % nA]))
    bmesh.ops.recalc_face_normals(tb, faces=tb.faces)
    tb.normal_update()
    f = tb.faces[len(tb.faces) // 2]; cc = f.calc_center_median()
    axis_pts = [sum((v.co for v in l_), Vector()) / len(l_) for l_ in lv]
    near_c = min(axis_pts, key=lambda q: (q - cc).length)
    if f.normal.dot(cc - near_c) < 0:            # outward, as build_tube
        bmesh.ops.reverse_faces(tb, faces=tb.faces)
    for f in tb.faces:
        f.smooth = True
    mesh = bpy.data.meshes.new(name); tb.to_mesh(mesh); tb.free()
    ob = bpy.data.objects.new(name, mesh); col_fit.objects.link(ob)
    mesh.materials.append(rect.data.materials[0])
    return ob


canal = build_tube2("AN_AnalCanal", canal_loops)
# the seam: every canal mouth point on a skin vertex; the lips' gap across the slit
from mathutils import kdtree as _kdt
_kd = _kdt.KDTree(len(me.vertices))
for v in me.vertices:
    _kd.insert(v.co, v.index)
_kd.balance()
report["passage"] = {
    "mouth_points": N_SPOKE, "mouth_len_mm": round((hole_loop[0] - hole_loop[N_SPOKE // 2]).length * 1000, 3),
    "seam_max_gap_mm": round(max(_kd.find(canal.data.vertices[i].co)[2] for i in range(N_SPOKE)) * 1000, 4),
    "lip_gap_max_mm": round(max((hole_loop[k] - hole_loop[N_SPOKE - k]).length for k in range(N_SPOKE // 8, 3 * N_SPOKE // 8))
                            * 1000, 4),
    "mouth_depth_below_A_mm": round((A - H).dot(-n) * -1000, 3), "canal_len_from_mouth_mm": round(_Lmm, 2),
    "canal_order_sign": _sg, "canal_k0": _k0,
    "lumen_arms_mm": {str(s_): [round(x_, 2) for x_ in arm_lengths(s_)] for s_ in (2.0, 5.0, 8.0, 12.0)}}

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
