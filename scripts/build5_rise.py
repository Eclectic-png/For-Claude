"""build5_rise: build5_adapted with the funnel turned into a rise - the surface climbs towards the centre (same
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
SQUEEZE = 0.75
FOLD_W = [1.2, 0.8, 1.0, 0.75, 1.05, 0.85, 1.2, 0.85, 1.05, 0.75, 1.0, 0.8]
SIG = 0.09; R_FOLD = 7.2
# rise: build5's 2 mm funnel profile turned upside down - the surface climbs towards the centre and fades to nothing
# at the outer ring
RISE = float(os.environ.get("ANUS_RISE", "0.2"))   # mm (0.5 and up read as an unnatural dome)
ENTRANCE_DIP = 0.0               # mm the pole drops below the rise (build5: 0.6 - on a rise it read as an abrupt sink)


def crease(theta):
    s = 0.0
    for i, w in enumerate(FOLD_W):
        dt = (theta - 2 * math.pi * i / 12 + math.pi) % (2 * math.pi) - math.pi
        s += w * math.exp(-dt * dt / (2 * SIG * SIG))
    return s


def fold_amp(r):
    t = min(max((r - 0.9) / (R_FOLD - 0.9), 0.0), 1.0)
    return 0.42 * math.sin(math.pi * t) ** 0.8


def base_point(r, theta):
    # adapted: same a and (squeezed) lateral offsets as before, but the lateral one is now arc length over the skin,
    # so the sides climb the cheek walls instead of being projected flat across the cleft
    a = r * math.cos(theta)
    b = r * math.sin(theta) * SQUEEZE * (0.45 + 0.55 * ss(0.9, 5.0, r))
    return surf_point(a, b)


rings = []; new_faces = []; info = {}
for r in RINGS:
    S = 72 if r < 7.0 else 36
    ring = []
    for k in range(S):
        th = 2 * math.pi * k / S
        v = bm.verts.new(base_point(r, th)); info[v] = (r, th); ring.append(v)
    rings.append(ring)
pole = bm.verts.new(base_point(0.0, 0.0)); info[pole] = (0.0, 0.0)
for k in range(72):
    new_faces.append(bm.faces.new((rings[0][k], rings[0][(k + 1) % 72], pole)))
for ri in range(len(rings) - 1):
    A_, B_ = rings[ri], rings[ri + 1]
    if len(A_) == len(B_):
        S = len(A_)
        for k in range(S):
            new_faces.append(bm.faces.new((A_[k], B_[k], B_[(k + 1) % S], A_[(k + 1) % S])))
    else:
        for k in range(36):
            new_faces.append(bm.faces.new((A_[2 * k], B_[k], A_[2 * k + 1])))
            new_faces.append(bm.faces.new((A_[2 * k + 1], B_[k], B_[(k + 1) % 36], A_[(2 * k + 2) % 72])))
bm.normal_update()
if sum(f.normal.dot(n) for f in new_faces) < 0:
    for f in new_faces:
        f.normal_flip()
bm.normal_update()
base_pos = {v: v.co.copy() for v in info}
# relief along the (smooth) base normal: radial creases; the closed centre keeps rising with the surface (no sink at
# the end). The rise goes straight out of the body (along n): along the surface normal it would push the steep
# cheek walls sideways into the cleft.
for v, (r, th) in info.items():
    relief = fold_amp(r) * (0.15 - crease(th))
    if v is pole:
        relief = -ENTRANCE_DIP
    rise = RISE * (1 - ss(0, R_FOLD, r))
    v.co = base_pos[v] + v.normal * (relief / 1000) + n * (rise / 1000)
    v[pig_layer] = min(1.0, (1 - ss(3.0, 8.5, r)) * (1 + 0.1 * crease(th) * fold_amp(r) / 0.42))

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
report["pucker_surface_span_mm"] = {
    "AP": round((rings[-1][0].co - rings[-1][18].co).length * 1000, 1),
    "lateral_over_surface": round(sum((rings[-1][k].co - rings[-1][k + 1].co).length for k in range(9, 27)) * 1000, 1)}
bm.to_mesh(me); me.update(); bm.free(); base.free()

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
