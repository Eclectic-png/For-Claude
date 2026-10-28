"""The vagina (female pelvis, reproductive step): replaces the AN_Vagina_Space placeholder with a real organ.
Runs on build_pelvis' output: blender -b out/blend/Hips_urinary.blend --python scripts/build_vagina.py -- <out.blend> <report dir>

  - AN_Vagina: the lumen (mucosal surface), collapsed at rest - a sagittal slit at the introitus, an H in the lower
    and middle vagina, a transverse slit higher up - opening into the vault round the cervix (the fornices). Its
    cross-sections are fitted ring by ring to the placeholder (walls ~3 mm) and to the cervix. Joined:
      * to the skin: a sagittal slit cut into the vestibule's crest (urinary_lib.cut_meatus), the same points;
      * to the cervix: the cervix pierces the vault (the os opens into the lumen), shared loop.
  - AN_Vagina_Wall: the outer surface (the placeholder, its bottom cap removed), rendered.
The placeholder stays (hidden) - build_pelvis' neighbours were fitted round it and do not move into it in any state,
so the vagina needs no shape keys of its own."""
import bpy, bmesh, sys, os, json, math
import numpy as np
from mathutils import Vector, Matrix
from mathutils.bvhtree import BVHTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import urinary_lib as U, tract_lib as T, vagina_lib as V

SAVE_AS, OUT = sys.argv[-2], sys.argv[-1]
O = bpy.data.objects; report = {}
E = O["Internal_Fit_Xform"]; col_fit = bpy.data.collections["Internal_Fit"]
Mw = E.matrix_world.copy(); Mi = Mw.inverted()
hips = O["Hips"]; vs = O["AN_Vagina_Space"]; uterus = O["AN_Uterus"]
STOP = os.environ.get("VAG_STOP", "")

# ---- landmarks (as build_pelvis) ----
X_ = Vector((1, 0, 0))
NECK = Vector(json.loads(os.environ.get("BLADDER_NECK", "[0.0, -0.080, 0.776]")))
_va = math.radians(float(os.environ.get("VAG_DEG", "35")))
VAG_A = Vector((0, math.cos(_va), math.sin(_va))); VAG_M = Vector((0, -math.sin(_va), math.cos(_va)))
VAG_P = NECK - VAG_M * 0.005; VAG_HT = 0.005; VAG_UP = float(os.environ.get("VAG_UP_MM", "40")) / 1000
P1v = VAG_P - VAG_M * (VAG_HT + 0.004); P2v = P1v + VAG_A * VAG_UP
_ib = Vector(json.loads(os.environ.get("VAG_INTROITUS", "[0.0, 0.0275, 0.7797]")))
_skt = U.tree_in(hips, None, cap=False)
_h = _skt.ray_cast(_ib - Vector((0, 0, 0.02)), Vector((0, 0, 1)), 0.04)
INTRO = Mi @ (_h[0] if _h[0] is not None and (_h[0] - _ib).length < 0.005 else _ib)
VC0 = INTRO + (P1v - INTRO).normalized() * 0.010     # the placeholder's first centre point

# ---- 1. the introitus: a sagittal slit in the vestibule's crest (body coords) ----
# the vaginal orifice fills the back of the vestibule, between the meatus (y 13.5 mm, body) and the fourchette
SLIT_Y = float(os.environ.get("VAG_SLIT_Y", "0.0220")); SLIT_HALF = float(os.environ.get("VAG_SLIT_HALF", "0.0045"))
N = 4 * sum(V.SEGS); LIP = 0.00002; G = 0.00001
_h = _skt.ray_cast(Vector((0, SLIT_Y, 0.70)), Vector((0, 0, 1)), 0.2)
M_body = Vector((0.0, SLIT_Y, _h[0].z if _h[0] is not None else 0.78))
REG_U, REG_V = SLIT_HALF + 0.0016, 0.0016
skin_snap = None                                   # cut_meatus carries the skin's shape keys itself
slit_idx, chart, skin_new, rep_s = U.cut_meatus(hips, M_body, Vector((0, 1, 0)), Vector((1, 0, 0)), N, SLIT_HALF, LIP, REG_U, REG_V)
hips.data.update()
slit_body = [hips.data.vertices[i].co.copy() for i in slit_idx]
slit = [Mi @ p for p in slit_body]                 # atlas
Cs = sum(slit, Vector()) / len(slit)
report["introitus"] = {"centre_body": [round(x, 5) for x in M_body], "half_mm": SLIT_HALF * 1000,
                       "length_real_mm": round((slit[0] - slit[N // 2]).length * 1000, 2), **rep_s}
print("INTROITUS", json.dumps(report["introitus"]))

# ---- 2. the centre line: up from the slit into the placeholder's axis, then along it ----
def catmull(pts, n):
    P = [pts[0]] + pts + [pts[-1]]; out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=(i == len(P) - 3)):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    return out
UP = (Mi.to_3x3() @ Vector((0, 0, 1))).normalized()          # into the body at the vestibule
dense = np.array([tuple(p) for p in catmull([Cs, Cs + UP * 0.004, VC0, P1v, P2v], 60)])
seg = np.linalg.norm(np.diff(dense, axis=0), axis=1); cum = np.concatenate([[0.0], np.cumsum(seg)])
S = []; s_ = 0.0
while s_ < cum[-1]:
    S.append(s_); s_ += 0.0010 if s_ < 0.020 else 0.0015
S = np.array(S)
C = np.stack([np.interp(S, cum, dense[:, k]) for k in range(3)], 1)
FR = V.frames_along(C, X_)

# sections of the placeholder (walls ~3 mm inside it) and of the cervix, ring by ring
def tri_arrays(ob):
    me = ob.data; X = np.array([tuple(v.co) for v in me.vertices]); F = []
    for p in me.polygons:
        vv = list(p.vertices)
        F += [(vv[0], vv[k], vv[k + 1]) for k in range(1, len(vv) - 1)]
    return X, np.array(F)
assert vs.matrix_local == Matrix.Identity(4) and uterus.matrix_local == Matrix.Identity(4)
Xv, Fv = tri_arrays(vs); Xu, Fu = tri_arrays(uterus)
env, cer = [], []
for c, (ex, ey), j in zip(C, FR, range(len(C))):
    t = C[min(j + 1, len(C) - 1)] - C[max(j - 1, 0)]; t /= np.linalg.norm(t)
    sv = V.plane_section(Xv, Fv, c, t, ex, ey)
    # only the part of the placeholder's section round the centre (it is one loop where the plane cuts the tube)
    env.append(V.extents(V.near(sv, 0.025)) if len(sv) else None)
    su = V.plane_section(Xu, Fu, c, t, ex, ey)
    cer.append(V.near(su, 0.022) if len(su) else None)
first_env = next(j for j, e in enumerate(env) if e)
last_env = max(j for j, e in enumerate(env) if e)
for j in range(first_env):                         # below the placeholder: its first section
    env[j] = env[first_env]
s_c0 = next(S[j] for j, u in enumerate(cer) if u is not None)          # the first ring that meets the cervix
WALL_LAT, WALL_AP = 0.0035, 0.003
SLIT_HALF_A = (slit[0] - slit[N // 2]).length / 2

# ring 0 = the skin's slit, put into the arm blocks (vertex 0 of the slit is its back tip)
ex0, ey0 = FR[0]
s2d = np.array([((p - Cs).dot(Vector(ex0)), (p - Cs).dot(Vector(ey0))) for p in slit])
sg = 1 if s2d[N // 4][0] > 0 else -1
assert s2d[0][1] < 0, "slit vertex 0 should be the back tip"
ring0 = np.zeros((N, 2)); idx0 = [0] * N
for k in range(N):
    jj = (N - N // 8 + sg * k) % N                  # back tip -> the back arm's tip (index 3N/4 + N/8)
    ring0[jj] = s2d[k]; idx0[jj] = k

rings = []; params = []
for j, s in enumerate(S):
    e = env[j]
    if j > last_env or e is None:
        break
    Lr, Ll = max(e[0] - WALL_LAT, 3 * G), max(e[1] - WALL_LAT, 3 * G)
    Vt = max(min(e[2], e[3]) - WALL_AP, 2 * G)
    f_h = U.ss(0.004, 0.016, s)
    fb = SLIT_HALF_A * (1 - f_h) + 0.0003 * f_h
    lr = (0.0003 * (1 - f_h) + Lr * f_h, 0.0003 * (1 - f_h) + Ll * f_h)
    s_u1 = s_c0 - 0.012; s_u0 = s_u1 - 0.010
    Vh = Vt * U.ss(0.010, 0.022, s) * (1 - U.ss(s_u0, s_u1, s))
    H = V.lumen_ring([(lr[0], Vh), (fb, 0.0), (lr[1], Vh), (fb, 0.0)], G)
    wv = U.ss(s_c0 - 0.010, s_c0 - 0.002, s)
    if wv > 0:
        u = cer[j]
        rx = min(max((np.abs(u[:, 0]).max() + 0.0015) if u is not None else 0.0, min(Lr, Ll)), min(e[0], e[1]) - 0.0015)
        ryb = min(max((-u[:, 1].min() + 0.0015) if u is not None else 0.0, e[3] - WALL_AP), e[3] - 0.0012)
        ryf = e[2] - WALL_AP
        P2 = H * (1 - wv) + V.ellipse_ring(max(rx, 0.001), max(ryf, 0.001), max(ryb, 0.001)) * wv
        params.append((round(s * 1000, 1), round(rx * 1000, 2), round(ryf * 1000, 2), round(ryb * 1000, 2)))
    else:
        P2 = H
    w_ex = 1 - U.ss(0.002, 0.006, s)
    rings.append(ring0 * w_ex + P2 * (1 - w_ex))
NR = len(rings)
report["lumen"] = {"rings": NR, "points_per_ring": N, "length_mm": round(float(S[NR - 1]) * 1000, 1),
                   "envelope_from_mm": round(float(S[first_env]) * 1000, 1), "cervix_from_mm": round(float(s_c0) * 1000, 1),
                   "vault_rings_s_rx_ryf_ryb_mm": params[::2]}
print("LUMEN", json.dumps(report["lumen"]))

# ---- 3. the lumen mesh: rings joined index by index, a dome on top ----
P3 = V.sweep(C[:NR], FR[:NR], rings)
P3[:N] = np.array([tuple(slit[idx0[jj]]) for jj in range(N)])   # ring 0 exactly the skin's slit
tip = C[NR - 1] + (C[NR - 1] - C[NR - 2]) / np.linalg.norm(C[NR - 1] - C[NR - 2]) * 0.0015
P3 = np.vstack([P3, tip[None]])
faces = []
for j in range(NR - 1):
    a, b = j * N, (j + 1) * N
    faces += [(a + k, a + (k + 1) % N, b + (k + 1) % N, b + k) for k in range(N)]
pole = NR * N; a = (NR - 1) * N
faces += [(a + k, a + (k + 1) % N, pole) for k in range(N)]
me = bpy.data.meshes.new("AN_Vagina"); me.from_pydata([tuple(p) for p in P3], [], faces); me.update()
old = O.get("AN_Vagina")
if old is not None:
    O.remove(old)
lumen = O.new("AN_Vagina", me); col_fit.objects.link(lumen); lumen.parent = E; lumen.matrix_parent_inverse = Matrix.Identity(4)
for p in me.polygons:
    p.use_smooth = True
bm = bmesh.new(); bm.from_mesh(me); bm.normal_update()
# normals out of the lumen (into the wall), like the other tubes: away from the centre line
_c = Vector(tuple(C[NR // 2])); f_ = bm.faces[(NR // 2) * N]
if f_.normal.dot(f_.calc_center_median() - _c) < 0:
    bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
bm.to_mesh(me); bm.free(); me.update()
bpy.context.view_layer.update()
if STOP == "lumen":
    U.save(SAVE_AS, OUT, report); sys.exit(0)

# ---- 4. the cervix pierces the vault: the os opens into the lumen (the uterus keeps its tilt keys) ----
snap_u = V.keyed_snapshot(uterus)
_cu = Vector(tuple(C[min(NR - 1, int(np.searchsorted(S, s_c0)) + 6)]))
J = T.open_junction(lumen, uterus, centre=Mw @ _cu, sides=(0,))
moved = V.keys_restore(uterus, snap_u)
report["cervix_join"] = {k: J[k] for k in ("curve_mm", "patch_radius_mm", "killed_faces", "inflate_mm")}
report["cervix_join"]["uterus_new_vertices"] = moved
report["cervix_join"]["opening"] = T.opening_loops(lumen, uterus)
print("CERVIX", json.dumps(report["cervix_join"]))

# ---- 5. the wall: the placeholder's surface without its bottom cap (the lumen leaves through it to the skin) ----
wme = vs.data.copy(); wme.name = "AN_Vagina_Wall"
bm = bmesh.new(); bm.from_mesh(wme)
t0 = (P1v - INTRO).normalized()
cap = [f for f in bm.faces if all((v.co - VC0).dot(t0) < 0.0006 for v in f.verts)]
bmesh.ops.delete(bm, geom=cap, context='FACES')
bmesh.ops.delete(bm, geom=[v for v in bm.verts if not v.link_faces], context='VERTS')
bm.to_mesh(wme); bm.free(); wme.update()
if wme.shape_keys:
    pass
old = O.get("AN_Vagina_Wall")
if old is not None:
    O.remove(old)
wall = O.new("AN_Vagina_Wall", wme); col_fit.objects.link(wall); wall.parent = E; wall.matrix_parent_inverse = Matrix.Identity(4)
wall.display_type = 'TEXTURED'; wall.hide_render = False
for p in wme.polygons:
    p.use_smooth = True
report["wall"] = {"bottom_cap_faces_removed": len(cap)}

# ---- 6. materials, lighting ----
import repro_lib as R
wme.materials.clear(); wme.materials.append(R.tinted(uterus.data.materials[0], "AN_Vagina_Wall", (0.86, 0.50, 0.52)))
lin = bpy.data.materials["AN_AnalCanal_Lining"].copy(); lin.name = "AN_Vagina_Lining"
me.materials.clear(); me.materials.append(lin)
# lining: 0 (skin) at the introitus -> 0.3 (the ramp's pale pink) by 8 mm in
at = me.attributes.get("lining") or me.attributes.new("lining", 'FLOAT', 'POINT')
_ring_of = lambda vi: min(vi // N, NR - 1)
for vi in range(len(me.vertices)):
    at.data[vi].value = 0.3 * float(U.ss(0.0, 0.008, S[_ring_of(vi)])) if vi < NR * N else 0.3
for m_ in (wme.materials[0], lin):
    if m_ and m_.node_tree:
        T.light_insides(m_)
_tr = T.scene_bvh()
report["ao_mean"] = {"AN_Vagina.AO": T.bake_ao(lumen, "AO", False, _tr), "AN_Vagina.AO_in": T.bake_ao(lumen, "AO_in", True, _tr),
                     "AN_Vagina_Wall.AO": T.bake_ao(wall, "AO", False, _tr), "AN_Vagina_Wall.AO_in": T.bake_ao(wall, "AO_in", True, _tr)}

# ---- 7. checks ----
def loops_of(ob):
    bm_ = U.bm_in(ob, None, weld=True); L = U.boundary_loops(bm_); out = [len(l) for l in L]; bm_.free(); return sorted(out)
skin_body = [hips.data.vertices[i].co for i in slit_idx]
lum_w = [Mw @ lumen.data.vertices[jj].co for jj in range(N)]
_wt = U.tree_in(wall, None); _wt_open = U.tree_in(wall, None, cap=False)
_ut = U.tree_in(O["AN_Urethra_Wall"], None) if "AN_Urethra_Wall" in O else None
LW = [Mw @ v.co for v in lumen.data.vertices]
first_in = first_env * N
inside_wall = [U.inside(_wt, p) for p in LW[first_in + 2 * N::7]]
report["checks"] = {
    "skin_seam_mm": round(max((a - b).length for a, b in zip(sorted(skin_body, key=lambda p: tuple(p)), sorted(lum_w, key=lambda p: tuple(p)))) * 1000, 5),
    "lumen_open_loops": loops_of(lumen), "wall_open_loops": loops_of(wall),
    "lumen_inside_wall_share": round(sum(inside_wall) / max(len(inside_wall), 1), 3),
    "lumen_wall_clearance_min_mm": round(min(_wt_open.find_nearest(p)[3] for p in LW[first_in + 2 * N::3]) * 1000, 2),
    "urethra_wall_clearance_min_mm": round(min(_ut.find_nearest(p)[3] for p in LW[::3]) * 1000, 2) if _ut else None}
print("CHECKS", json.dumps(report["checks"]))
U.save(SAVE_AS, OUT, report)
