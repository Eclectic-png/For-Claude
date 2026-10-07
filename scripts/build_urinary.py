"""build_urinary: the female urinary system in the working file (runs on build5_passage's output).

Kidneys, adrenals and ureters come from the atlas (anatomy_ref.blend, read only) into Internal_Fit, under
Internal_Fit_Xform like the other organs. Everything is built in that frame ("atlas space", real-size metres): the
fit squashes it to the body (0.62 / 0.72 / 0.72), so volumes in the body are 0.32x the real ones in the report.
  1. Kidneys: the atlas ureter already carries the renal pelvis and calyces inside the kidney; the kidney's wall is
     opened where the pelvis leaves it at the hilum (tract_lib.open_junction, sides=(0,): pierce, not join).
  2. Bladder (new; the atlas one is a 226-face blob pressed against the rectum): a hollow shell grown by a pressure
     simulation against the pubic bones, the pelvic walls, the inside of the skin (abdominal wall) and the anterior
     vaginal wall plane (room kept for the vagina / uterus). Empty (~70 ml, flattened, lid pressed down by the
     bowel) to full (500 ml, dome rising behind the abdominal wall), in stages so every in-between is a real state.
     The trigone (neck + ureteric orifices) does not move.
  3. Ureters: their last few cm bend to the trigone's corners, entering obliquely; the lumens join the bladder
     (open_junction).
  4. Urethra (new): from the bladder neck (a shared loop) down behind the pubic arch to the external meatus, a
     sagittal slit cut into the skin of the vestibule (a shared loop). Collapsed at rest: a transverse crescent
     turning into the meatus' sagittal slit. Muscular wall sleeve (thicker at the external sphincter) around it.
     The meatus position is URETHRA_MEATUS="x,y,z" (body coordinates, snapped to the skin), else the
     "Urethra_Meatus_Target" empty of URETHRA_MEATUS_FROM=<a previous output .blend> (move the empty there, re-run),
     else the default below; the urethra re-routes to it.
  5. Controls: empty "Urinary_Controls" with Fill (0 empty .. 1 full) and Void (0 closed .. 1 voiding). Fill drives
     the bladder's stage keys and the corrective keys of every neighbour it presses on (the bowel gives way and is
     compressed, bones and the abdominal wall do not); Void opens bladder neck, urethra, sphincter and meatus.
  6. Insides lit + AO / AO_in bakes for every organ and bone in the file.
Run: ANATOMY_REF=anatomy_ref.blend blender -b out/blend/Hips_build5_passage.blend --python scripts/build_urinary.py \
     -- <out.blend> <report folder>"""
import bpy, bmesh, sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tract_lib as T
import urinary_lib as U
from mathutils import Vector, Matrix, kdtree
from mathutils.bvhtree import BVHTree

SAVE_AS, OUT = sys.argv[-2], sys.argv[-1]
REF = os.environ.get("ANATOMY_REF", r"C:\Users\Parker\Anatomy Project\anatomy_ref.blend")
STOP = os.environ.get("URINARY_STOP", "")         # development: stop (and save) after this step
O = bpy.data.objects; report = {}
E = O["Internal_Fit_Xform"]; col_fit = bpy.data.collections["Internal_Fit"]
Mw = E.matrix_world.copy(); Mi = Mw.inverted()   # atlas space <-> body (world)
hips = O["Hips"]


def done(step):
    if STOP == step:
        U.save(SAVE_AS, OUT, report); sys.exit(0)


# ======================= 0. load =======================
LOAD = ["AN_Kidney_L", "AN_Kidney_R", "AN_Adrenal_L", "AN_Adrenal_R", "AN_Ureter_L", "AN_Ureter_R", "AN_Bladder"]
for nm in LOAD:
    ob = O.get(nm)
    if ob:
        me_ = ob.data; O.remove(ob)
        if me_ and me_.users == 0:
            bpy.data.meshes.remove(me_)
with bpy.data.libraries.load(REF) as (src, dst):
    dst.objects = LOAD
for ob in dst.objects:
    col_fit.objects.link(ob); ob.parent = E; ob.matrix_parent_inverse = Matrix.Identity(4)
bpy.context.view_layer.update()
# work in atlas space: each loaded object's mesh in metres, object = identity under E
for ob in dst.objects:
    ob.data.transform(ob.matrix_local); ob.matrix_local = Matrix.Identity(4)
bpy.context.view_layer.update()
old_bladder = O["AN_Bladder"]; old_bladder.name = "_atlas_bladder"

# ======================= 1. kidneys: open the wall where the renal pelvis leaves at the hilum =======================
report["hilum"] = {}
for s in "LR":
    K, Ur = O["AN_Kidney_" + s], O["AN_Ureter_" + s]
    trees = {K.name: T.closed_tree(K), Ur.name: T.closed_tree(Ur)}
    bk = T.world_bm(K); n_in = len(T.overlap_components(bk, trees[Ur.name])); bk.free()
    r = T.open_junction(K, Ur, trees=trees, sides=(0,))
    r["kidney_pieces_inside_pelvis"] = n_in
    r["opening"] = T.opening_loops(K, Ur)
    report["hilum"][s] = r; print("HILUM", s, json.dumps(r))
done("kidneys")

# ======================= 2. bladder: grown against its surroundings, empty .. full =======================
# landmarks (atlas space, m). The neck sits ~19 mm behind the back of the pubic symphysis at the level of its lower
# border (lower than the atlas' male bladder, which rested on the prostate's level); the urethra runs from it to the
# meatus.
X_ = Vector((1, 0, 0))
NECK = Vector(json.loads(os.environ.get("BLADDER_NECK", "[0.0, -0.080, 0.776]")))
# anterior vaginal wall: the bladder base rests on it. A plane through the vagina's axis (VAG_DEG up and back from
# horizontal), VAG_GAP behind the neck; the bladder stays in front of / above it for VAG_LEN along the axis, so the
# vagina (reproductive step) has its room between urethra / bladder and rectum
VAG_DEG = float(os.environ.get("VAG_DEG", "35")); _va = math.radians(VAG_DEG)
VAG_A = Vector((0, math.cos(_va), math.sin(_va))); VAG_M = Vector((0, -math.sin(_va), math.cos(_va)))
VAG_GAP = 0.005; VAG_LEN = 0.075
VAG_P = NECK - VAG_M * VAG_GAP
BONE_GAP = 0.003; WALL = 0.010                  # clear of bone; inside the skin by the abdominal wall
V_STAGES = [float(v) for v in os.environ.get("BLADDER_ML", "70,150,260,380,500").split(",")]
EMPTY_LID = 0.032                                # empty: the bowel presses the lid down to NECK.z + this

bones = [O[n_] for n_ in ("AN_HipBone_L", "AN_HipBone_R", "AN_Sacrum") if n_ in O]
bone_trees = [(U.tree_in(b, Mi), BONE_GAP) for b in bones]
skin_tree = U.tree_in(hips, Mi, cap=False)

# ---- female pelvis: room for the vagina. The atlas is a male body: its rectum / ampulla sit where the vagina goes
# (between urethra / bladder base and rectum). A placeholder for the collapsed vagina (walls included) is laid from the
# introitus up behind the urethra, then back along the bladder base, and the bowel is moved out of it (back toward the
# sacrum; where it cannot go further it flattens); the anal canal stays (it is welded to the anus) and every shape key
# of a moved organ moves with it. The bladder and ureters keep clear of it; the reproductive step builds the vagina in it.
VAG_HW, VAG_HT, VAG_UP = 0.013, 0.005, 0.030     # half width, half thickness (walls, collapsed), upper part length
_ib = Vector(json.loads(os.environ.get("VAG_INTROITUS", "[0.0, 0.0275, 0.7797]")))     # body coords (the dimple)
_sk0 = U.bm_in(hips); _skt0 = BVHTree.FromBMesh(_sk0); _sk0.free()
_h = _skt0.ray_cast(_ib - Vector((0, 0, 0.02)), Vector((0, 0, 1)), 0.04)
INTRO = Mi @ (_h[0] if _h[0] is not None and (_h[0] - _ib).length < 0.005 else _ib)
# behind the neck: on the bladder-base plane at the neck's height, VAG_HT + 1 mm behind it; then up its slope
P1v = VAG_P - VAG_M * (VAG_HT + 0.004)
P2v = P1v + VAG_A * VAG_UP
def catmull(pts, n):
    P = [pts[0]] + pts + [pts[-1]]; out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=(i == len(P) - 3)):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    return out
_t0 = (P1v - INTRO).normalized()
vc = catmull([INTRO + _t0 * 0.010, P1v, P2v], 16)     # starts 10 mm inside the body (clear of the vestibule)
bm = bmesh.new(); rings = []
for j, c in enumerate(vc):
    s_ = j / (len(vc) - 1)
    t_ = (vc[min(j + 1, len(vc) - 1)] - vc[max(j - 1, 0)]).normalized(); e2 = t_.cross(X_).normalized()
    # narrower at the introitus, rounded off at the top (the fornix)
    rnd = math.sqrt(max(0.0, 1 - (max(0.0, s_ - 0.82) / 0.18) ** 2)) * 0.8 + 0.2 if s_ > 0.82 else 1.0
    w_ = VAG_HW * (0.7 + 0.3 * float(U.ss(0.0, 0.2, s_))) * rnd; h_ = VAG_HT * (0.6 + 0.4 * float(U.ss(0.0, 0.1, s_))) * rnd
    rings.append([bm.verts.new(c + X_ * (w_ * math.cos(a)) + e2 * (h_ * math.sin(a))) for a in np.linspace(0, 2 * math.pi, 24, endpoint=False)])
for a_, b_ in zip(rings, rings[1:]):
    for k in range(24):
        bm.faces.new((a_[k], a_[(k + 1) % 24], b_[(k + 1) % 24], b_[k]))
for rr, flip in ((rings[0], True), (rings[-1], False)):
    f_ = bm.faces.new(rr[::-1] if flip else rr)
bm.normal_update()
_cv = sum(vc, Vector()) / len(vc)
if sum((f.calc_center_median() - _cv).dot(f.normal) for f in bm.faces) < 0:
    bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
bmesh.ops.triangulate(bm, faces=bm.faces[:])
vs_me = bpy.data.meshes.new("AN_Vagina_Space"); bm.to_mesh(vs_me); bm.free()
vag_space = O.get("AN_Vagina_Space")
if vag_space is not None:
    O.remove(vag_space)
vag_space = O.new("AN_Vagina_Space", vs_me); col_fit.objects.link(vag_space); vag_space.parent = E
vag_space.matrix_parent_inverse = Matrix.Identity(4); vag_space.display_type = 'WIRE'; vag_space.hide_render = True
bpy.context.view_layer.update()
vag_tree = U.tree_in(vag_space, Mi)
vs_me.update()
_vpts = np.array([tuple(v.co) for v in vs_me.vertices]); _vnrm = np.array([tuple(v.normal) for v in vs_me.vertices])
BOWEL = [n_ for n_ in ("AN_Colon_Descending", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_AnalCanal") if n_ in O]
_before = {n_: np.array([tuple(v.co) for v in O[n_].data.vertices]) for n_ in BOWEL}
_vl = []
# the anal canal is held only near the anus (its sphincter, its weld to the skin): its upper part flexes back with
# the ampulla (in a woman the anorectal junction sits behind the middle of the vagina)
_cb = U.bm_in(O["AN_AnalCanal"], Mi); _loops = U.boundary_loops(_cb)
_anus = np.array([tuple(v.co) for v in min(_loops, key=lambda lp: sum(v.co.z for v in lp) / len(lp))]); _cb.free()
def _pin_v(n_, P):
    if n_ == "AN_AnalCanal":
        return np.array([np.linalg.norm(_anus - p, axis=1).min() < 0.012 for p in P])
    return np.zeros(len(P), bool)
room_v, res_v = U.make_room([O[n_] for n_ in BOWEL], Mi, Mw, [vag_tree], _pin_v,
                            U.Obstacles([(t_, 0.0008) for t_, _ in bone_trees], [(skin_tree, 0.003)], []),
                            gap=0.002, decay=0.996, iters=400, rounds=10, log=_vl, expanders=[(_vpts, _vnrm)])
_sk_w = U.tree_in(hips, None, cap=False)
report["vagina_space"] = {"skin_clearance_mm": round(min(_sk_w.find_nearest(Mw @ Vector(p))[3] for p in _vpts) * 1000, 2),
                          "introitus_atlas": list(INTRO), "top_atlas": list(P2v), "length_mm": round(sum((a - b).length for a, b in zip(vc, vc[1:])) * 1000, 1),
                          "bowel_moved_mm": {}, "residual_mm": res_v}
for n_ in BOWEL:
    ob = O[n_]; D = room_v[n_][0] - _before[n_]
    mv = float(np.linalg.norm(D, axis=1).max())
    report["vagina_space"]["bowel_moved_mm"][n_] = round(mv * (1.0 if ob.parent else 1000.0), 2)
    if mv < 1e-9:
        continue
    keys = ob.data.shape_keys.key_blocks if ob.data.shape_keys else []
    for kb_ in keys:                             # every key moves with the organ (the anus' Open key stays valid)
        K = np.array([tuple(d.co) for d in kb_.data]) + D
        kb_.data.foreach_set("co", K.astype(np.float32).ravel())
    ob.data.vertices.foreach_set("co", room_v[n_][0].astype(np.float32).ravel()); ob.data.update()
bpy.context.view_layer.update()
print("VAGINA", json.dumps(report["vagina_space"]))
done("vagina")
bone_trees_v = bone_trees + [(vag_tree, 0.001)]

# the bowel and anal canal: solid for the EMPTY bladder (it fits round them at rest); filling, they give way
SOFT = [n_ for n_ in ("AN_Colon_Descending", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_AnalCanal", "AN_AnalSphincter")
        if n_ in O]
SOFT_GAP = 0.0015
soft_trees = [(U.tree_in(O[n_], Mi), SOFT_GAP) for n_ in SOFT]
vag_mask = lambda P: ((P - np.array(VAG_P)) @ np.array(VAG_A)) < VAG_LEN
planes = [(tuple(VAG_P), tuple(VAG_M), vag_mask)]

# seed: a 10 mm sphere sitting on the neck, straight above it along the urethra's line (clear of every limit)
MEATUS_DEFAULT = Vector((0.0, 0.0135, 0.7835))   # body coordinates (vestibule crest), see step 4
_mt = O.get("Urethra_Meatus_Target")
if os.environ.get("URETHRA_MEATUS_FROM"):        # a previous output whose target empty was moved by hand
    with bpy.data.libraries.load(os.environ["URETHRA_MEATUS_FROM"]) as (src, dst):
        dst.objects = ["Urethra_Meatus_Target"]
    _mt = dst.objects[0]
if os.environ.get("URETHRA_MEATUS"):
    M_body = Vector(json.loads("[" + os.environ["URETHRA_MEATUS"] + "]"))
elif _mt is not None:
    M_body = _mt.matrix_world.translation.copy()
else:
    M_body = MEATUS_DEFAULT.copy()
# snap to the skin: straight up from below (in the vestibule's slot that reaches the crest), else the nearest point
_sk = U.bm_in(hips); _skt = BVHTree.FromBMesh(_sk); _sk.free()
_hit = _skt.ray_cast(M_body - Vector((0, 0, 0.02)), Vector((0, 0, 1)), 0.04)
M_body = _hit[0] if _hit[0] is not None and (_hit[0] - M_body).length < 0.005 else _skt.find_nearest(M_body)[0]
M_atl = Mi @ M_body
U0 = (NECK - M_atl).normalized()                 # into the bladder, continuing the urethra's line
# retropubic floor: in front of the neck the bladder's anterior wall rises toward the pubis (pubovesical ligaments,
# retropubic fat), FLOOR_DEG up from horizontal, so the neck is the bottom of a V between it and the vaginal wall.
# Both planes pass a little under the neck so its pinned collar sits clear of them.
FLOOR_DEG = float(os.environ.get("FLOOR_DEG", "35")); _fa = math.radians(FLOOR_DEG)
FLOOR_N = Vector((0, math.sin(_fa), math.cos(_fa)))
planes.append((tuple(NECK - U0 * 0.0025), tuple(FLOOR_N), None))
obst = U.Obstacles(bone_trees_v, [(skin_tree, WALL)], planes)
Xs, Fb = U.icosphere(5); r0 = 0.010
Xb = np.array(NECK + U0 * r0) + Xs * r0
dN = np.linalg.norm(Xb - np.array(NECK), axis=1)
mob = U.ss(0.003, 0.008, dN)                     # the neck stays put
log = []
lid = (tuple(NECK + Vector((0, 0, EMPTY_LID))), (0, 0, -1), None)
obst_empty = U.Obstacles(bone_trees_v + soft_trees, [(skin_tree, WALL)], planes + [lid])
_dbg = []
_v0, _ = obst_empty.project(Xb, np.ones(len(Xb), bool), debug=_dbg); print("SEED", _dbg, [b.name for b in bones] + SOFT)
report["seed_violations_mm"] = round(float(np.abs(_v0 - Xb).max()) * 1000, 2)
Xe, Ve = U.grow(Xb, Fb, mob, V_STAGES[0] * 1e-6, obst_empty, log=log)
report["bladder"] = {"neck": list(NECK), "meatus_body": list(M_body), "empty_ml": round(Ve * 1e6, 1), "log_empty": log}
print("BLADDER empty", Ve * 1e6, log[-3:])

# the trigone: the neck and the two ureteric orifices on the base, ~26 mm apart, ~25 mm up the base from the neck
ORI_UP = 0.022; ORI_HALF = 0.013
nrm_e = U.vnormals(Xe, Fb)
def on_shell(p):
    i = int(np.argmin(np.linalg.norm(Xe - np.array(p), axis=1))); return Vector(Xe[i]), Vector(nrm_e[i])
ORI = {s: on_shell(NECK + VAG_A * ORI_UP + X_ * (ORI_HALF if s == "L" else -ORI_HALF)) for s in "LR"}
# ---- the bladder object: the empty stage, opened at the neck onto the urethra's first ring
K_U = 32                                          # vertices round the urethra (and the neck, and the meatus)
NECK_HALF, NECK_BEND, LIP = 0.0025, 0.0006, 0.00002   # resting neck / urethra: a transverse crescent, lips 0.02 mm
R_HOLE = 0.005
T0 = -U0                                         # the urethra leaves the neck straight down its line
E_LAT = X_.copy(); E_AP = T0.cross(E_LAT).normalized()
if E_AP.y < 0:
    E_AP = -E_AP                                 # crescent bulges backward (urethral crest)
neck_rest = U.slit_ring(K_U, NECK, E_LAT, E_AP, NECK_HALF, LIP, NECK_BEND)

bm = bmesh.new()
bv = [bm.verts.new(Vector(p)) for p in Xe]
for f in Fb:
    bm.faces.new([bv[i] for i in f])
bm.verts.ensure_lookup_table()
bmesh.ops.delete(bm, geom=[v for v in bm.verts if (v.co - NECK).length < R_HOLE], context='VERTS')
loops = U.boundary_loops(bm); assert len(loops) == 1, len(loops)
_v = (-U0).cross(E_LAT)
neck_v, n_ann = U.fill_annulus(bm, loops[0], neck_rest, NECK, E_LAT, _v)
bm.normal_update()
# outward normals (the icosphere's faces point out; the annulus follows the same winding test)
_c = Vector(np.array(Xe).mean(0))
if sum((f.calc_center_median() - _c).dot(f.normal) for f in bm.faces) < 0:
    bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
bl_me = bpy.data.meshes.new("AN_Bladder"); bm.to_mesh(bl_me); bm.free()
bladder = O.new("AN_Bladder", bl_me); col_fit.objects.link(bladder); bladder.parent = E
bladder.matrix_parent_inverse = Matrix.Identity(4)
for m_ in old_bladder.data.materials:
    bl_me.materials.append(m_)
bpy.context.view_layer.update()
report["bladder"]["neck_annulus_faces"] = n_ann

# ======================= 3. ureters: their ends bend to the trigone's corners, the lumens join the bladder ===========
URE_IN = 0.003; URE_BLEND = 0.060
report["ureters"] = {}
for s in "LR":
    ob = O["AN_Ureter_" + s]; me_ = ob.data
    P = np.array([tuple(v.co) for v in me_.vertices])
    end = P[P[:, 2] < P[:, 2].min() + 0.003].mean(0)        # the lowest end of the tube
    o_, n_ = ORI[s]
    tgt = np.array(o_ - n_ * URE_IN)
    w = 1 - U.ss(0.0, URE_BLEND, np.linalg.norm(P - end, axis=1))
    P2 = P + np.outer(w, tgt - end)
    me_.vertices.foreach_set("co", P2.ravel()); me_.update()
    report["ureters"][s] = {"end_moved_mm": round(float(np.linalg.norm(tgt - end)) * 1000, 1)}
bpy.context.view_layer.update()
# the atlas' left ureter runs through the sigmoid colon (~44 vertices): it lies behind it, so it is eased out of every
# bowel shell at rest (the colon stays; the ureter's displacement fades along its length, still inside the kidney)
_kt = {s: U.tree_in(O["AN_Kidney_" + s], Mi) for s in "LR"}
_bowel = [U.tree_in(O[n_], Mi) for n_ in ("AN_Colon_Descending", "AN_Rectum") if n_ in O] + [vag_tree]
for s in "LR":
    ob = O["AN_Ureter_" + s]
    for t_ in _bowel:
        res, rr = U.make_room([ob], Mi, Mw, [t_], lambda n_, P: np.array([_kt[n_[-1]].find_nearest(Vector(p))[3] < 0.006 for p in P]),
                              U.Obstacles([(t__, 0.0008) for t__, _ in bone_trees], [], []), gap=0.0015)
        moved = float(np.abs(res[ob.name][0] - np.array([tuple(v.co) for v in ob.data.vertices])).max())
        ob.data.vertices.foreach_set("co", res[ob.name][0].astype(np.float32).ravel()); ob.data.update()
        report["ureters"][s].setdefault("eased_out_of_bowel_mm", []).append(round(moved * 1000, 2))
bpy.context.view_layer.update()
trees = {o_.name: T.closed_tree(o_) for o_ in (bladder, O["AN_Ureter_L"], O["AN_Ureter_R"])}
for s in "LR":
    r = T.open_junction(O["AN_Ureter_" + s], bladder, trees=trees, centre=Mw @ ORI[s][0])
    r["opening"] = T.opening_loops(O["AN_Ureter_" + s], bladder)
    report["ureters"][s]["junction"] = r; print("URETER", s, json.dumps(r)[:400])
# the ureteric openings on the bladder, in atlas space (for the still trigone)
ure_loop = {}
_bb = U.bm_in(bladder, Mi)
for s in "LR":
    _tu = U.tree_in(O["AN_Ureter_" + s], Mi, cap=False)
    ure_loop[s] = np.array([tuple(v.co) for v in _bb.verts if v.is_boundary and _tu.find_nearest(v.co)[3] < 2e-6])
    report["ureters"][s]["opening_from_planned_orifice_mm"] = round(float(np.linalg.norm(ure_loop[s].mean(0) - np.array(ORI[s][0]))) * 1000, 1)
_bb.free()
done("ureters")
# the still trigone round the ACTUAL openings (the ureters cross the wall where their tubes meet it, which can be a
# couple of cm from the planned orifice): the neck and every vertex of both ureteric opening loops
TRIG_FIX, TRIG_FREE = 0.016, 0.050
trig_pts = np.array([tuple(NECK)] + [p for s in "LR" for p in ure_loop[s]])
def dist_trigone(P):
    P = np.atleast_2d(P); out_ = np.empty(len(P))
    kd_ = kdtree.KDTree(len(trig_pts))
    for i, p in enumerate(trig_pts):
        kd_.insert(Vector(p), i)
    kd_.balance()
    for i, p in enumerate(P):
        out_[i] = kd_.find(Vector(p))[2]
    return out_
mob_f = U.ss(TRIG_FIX, TRIG_FREE, dist_trigone(Xe))
stages = [Xe]
for V_ in V_STAGES[1:]:
    log = []
    Xn, Vn = U.grow(stages[-1], Fb, mob_f, V_ * 1e-6, obst, log=log)
    stages.append(Xn); print("BLADDER stage", V_, Vn * 1e6, log[-2:])
    report["bladder"].setdefault("stages_ml", []).append(round(Vn * 1e6, 1))
report["bladder"]["orifices"] = {s: list(ORI[s][0]) for s in "LR"}
report["bladder"]["opening_centres"] = {s: [round(x, 5) for x in np.mean(ure_loop[s], axis=0)] for s in "LR"}

if STOP == "bladder_raw":                         # development: every stage as its own object
    for k, Xk in enumerate(stages):
        me_ = bpy.data.meshes.new(f"_bl{k}"); me_.from_pydata([tuple(p) for p in Xk], [], [tuple(f) for f in Fb])
        ob_ = O.new(f"_bladder_stage{k}", me_); col_fit.objects.link(ob_); ob_.parent = E
    U.save(SAVE_AS, OUT, report); sys.exit(0)

# ======================= 4. urethra: neck -> external meatus (a slit cut into the vestibule's crest) =================
# 4a. the meatus: a sagittal slit, MEATUS_HALF from its middle to each tip (body scale; ~5 mm long real), lips 0.02 mm
MEATUS_HALF = 0.0018; REG_U, REG_V = 0.0032, 0.0015
T_AP = Vector((0, 1, 0)); E_LAT_B = Vector((1, 0, 0))
slit_idx, chart, skin_new, rep_m = U.cut_meatus(hips, M_body, T_AP, E_LAT_B, K_U, MEATUS_HALF, LIP, REG_U, REG_V)
report["meatus"] = {"body": list(M_body), "atlas": list(M_atl), **rep_m}
meatus_w = [hips.data.vertices[i].co.copy() for i in slit_idx]          # Hips has the identity transform
meatus_a = np.array([tuple(Mi @ p) for p in meatus_w])
# Void on the skin: the lips part into a lens (VOID_SKIN to each side at the middle), the rebuilt patch follows,
# fading to nothing at its border; tips stay
VOID_SKIN = 0.0010
skin_void = {}
for i, (u, v) in skin_new.items():
    E_ = math.sqrt(max(0.0, 1 - (u / (MEATUS_HALF * 1.05)) ** 2))
    f_ = 1 - float(U.ss(0.0, REG_V * 0.95, abs(v)))
    sgn = 0.0 if abs(v) < 1e-9 else math.copysign(1.0, v)
    skin_void[i] = hips.data.vertices[i].co + E_LAT_B * (sgn * VOID_SKIN * E_ * f_)
meatus_void_a = np.array([tuple(Mi @ skin_void[i]) for i in slit_idx])

# 4b. centreline: a cubic from the neck (straight down its line) to the meatus (arriving straight up into the body)
def bezier(P, t):
    a, b, c, d = P; m = 1 - t
    return a * m ** 3 + b * 3 * m * m * t + c * 3 * m * t * t + d * t ** 3
T_END = (Mi.to_3x3() @ Vector((0, 0, -1))).normalized()
L0 = (M_atl - NECK).length
CP = [NECK, NECK + T0 * L0 * 0.35, M_atl - T_END * L0 * 0.35, M_atl]
_dense = [bezier(CP, i / 2000) for i in range(2001)]
_cum = [0.0]
for a_, b_ in zip(_dense, _dense[1:]):
    _cum.append(_cum[-1] + (b_ - a_).length)
L_U = _cum[-1]; NR = 56
centre = []; tang = []
for j in range(NR):
    tgt = L_U * j / (NR - 1); k = min(range(len(_cum)), key=lambda i: abs(_cum[i] - tgt))
    centre.append(_dense[k]); tang.append((_dense[min(k + 1, 2000)] - _dense[max(k - 1, 0)]).normalized())
frames = [(E_LAT.copy(), E_AP.copy())]
for j in range(1, NR):                            # parallel transport
    e1, _ = frames[-1]; e1 = (e1 - tang[j] * e1.dot(tang[j])).normalized(); frames.append((e1, tang[j].cross(e1)))
    if frames[-1][1].dot(frames[0][1]) < 0:
        frames[-1] = (e1, -frames[-1][1])
report["urethra"] = {"length_mm": round(L_U * 1000, 1), "length_body_mm": round(sum(((Mw @ a_) - (Mw @ b_)).length for a_, b_ in zip(centre, centre[1:])) * 1000, 1)}

# 4c. rings: rest = transverse crescent turning into the meatus' sagittal slit; open = round, 3.5 -> 2.6 mm
MEATUS_HALF_A = MEATUS_HALF / Mw.to_3x3().col[1].length
R_OPEN_NECK, R_OPEN_MEATUS = 0.0035, 0.0026
rest_r, open_r = [], []
neck_open = U.round_ring(K_U, NECK, E_LAT, E_AP, R_OPEN_NECK)
for j in range(NR):
    s_ = j / (NR - 1); e1, e2 = frames[j]; th = math.pi / 2 * float(U.ss(0.55, 0.95, s_))
    a1 = e1 * math.cos(th) + e2 * math.sin(th); a2 = -e1 * math.sin(th) + e2 * math.cos(th)
    half_ = NECK_HALF + (MEATUS_HALF_A - NECK_HALF) * float(U.ss(0.55, 0.95, s_))
    bend_ = NECK_BEND * (1 - float(U.ss(0.45, 0.85, s_)))
    g = U.slit_ring(K_U, centre[j], a1, a2, half_, LIP, bend_)
    o = U.round_ring(K_U, centre[j], a1, a2, R_OPEN_NECK + (R_OPEN_MEATUS - R_OPEN_NECK) * s_)
    beta = float(U.ss(0.86, 1.0, s_))
    c_ = np.array(centre[j])
    g = g * (1 - beta) + (meatus_a - meatus_a.mean(0) + c_) * beta
    o = o * (1 - beta) + (meatus_void_a - meatus_a.mean(0) + c_) * beta
    rest_r.append(g); open_r.append(o)
rest_r[0] = neck_rest; open_r[0] = neck_open
rest_r[-1] = meatus_a; open_r[-1] = meatus_void_a

bm = bmesh.new()
rv = [[bm.verts.new(Vector(p)) for p in ring] for ring in rest_r]
for j in range(NR - 1):
    for k in range(K_U):
        bm.faces.new((rv[j][k], rv[j][(k + 1) % K_U], rv[j + 1][(k + 1) % K_U], rv[j + 1][k]))
bm.normal_update()
_out = sum((f.calc_center_median() - Vector(centre[min(int(i / K_U), NR - 1)])).dot(f.normal) for i, f in enumerate(bm.faces))
if _out < 0:
    bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
ur_me = bpy.data.meshes.new("AN_Urethra"); bm.to_mesh(ur_me); bm.free()
urethra = O.new("AN_Urethra", ur_me); col_fit.objects.link(urethra); urethra.parent = E
urethra.matrix_parent_inverse = Matrix.Identity(4)
lin_src = bpy.data.materials.get("AN_AnalCanal_Lining")
lin = lin_src.copy(); lin.name = "AN_Urethra_Lining"; ur_me.materials.append(lin)
at_ = ur_me.attributes.new("lining", 'FLOAT', 'POINT')
_dm = [L_U * (1 - j / (NR - 1)) for j in range(NR)]          # distance from the meatus along the urethra
for j in range(NR):
    for k in range(K_U):
        at_.data[j * K_U + k].value = 0.62 * float(U.ss(0.0, 0.008, _dm[j]))   # skin, pale, then mucosa
urethra_open = np.concatenate(open_r)

# 4d. the muscular wall (smooth muscle + the external sphincter, thicker over the middle third): a sleeve round the
# urethra, clear of the open lumen inside, of the bladder and of the skin at its ends
def sleeve_r(s_):
    ri = R_OPEN_NECK + (R_OPEN_MEATUS - R_OPEN_NECK) * s_ + 0.0004
    return ri, ri + 0.0012 + 0.0010 * float(U.ss(0.25, 0.4, s_) * (1 - U.ss(0.6, 0.75, s_)))
bl_tree = U.tree_in(bladder, Mi); skin_tree_a = skin_tree
def ring_clear(j, tree, cap=True):
    ri, ro = sleeve_r(j / (NR - 1)); e1, e2 = frames[j]
    pts = U.round_ring(K_U, centre[j], e1, e2, ro)
    return min(tree.find_nearest(Vector(p))[3] for p in pts)
j_a = next(j for j in range(1, NR // 2) if ring_clear(j, bl_tree) > 0.0012)
j_b = next(j for j in range(NR - 2, NR // 2, -1) if ring_clear(j, skin_tree_a) > 0.0006)
bm = bmesh.new(); S_ = 32
outer_v, inner_v = [], []
for j in range(j_a, j_b + 1):
    ri, ro = sleeve_r(j / (NR - 1)); e1, e2 = frames[j]
    outer_v.append([bm.verts.new(Vector(p)) for p in U.round_ring(S_, centre[j], e1, e2, ro)])
    inner_v.append([bm.verts.new(Vector(p)) for p in U.round_ring(S_, centre[j], e1, e2, ri)])
for a in range(len(outer_v) - 1):
    for k in range(S_):
        k2 = (k + 1) % S_
        bm.faces.new((outer_v[a][k], outer_v[a][k2], outer_v[a + 1][k2], outer_v[a + 1][k]))
        bm.faces.new((inner_v[a][k], inner_v[a + 1][k], inner_v[a + 1][k2], inner_v[a][k2]))
for a, flip in ((0, True), (-1, False)):
    for k in range(S_):
        k2 = (k + 1) % S_; q = (outer_v[a][k], outer_v[a][k2], inner_v[a][k2], inner_v[a][k])
        bm.faces.new(q[::-1] if flip else q)
bm.normal_update()
if sum((f.calc_center_median() - Vector(np.array(centre).mean(0))).dot(f.normal) for f in bm.faces if f.verts[0] in outer_v[len(outer_v) // 2]) < 0:
    bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
w_me = bpy.data.meshes.new("AN_Urethra_Wall"); bm.to_mesh(w_me); bm.free()
uwall = O.new("AN_Urethra_Wall", w_me); col_fit.objects.link(uwall); uwall.parent = E
uwall.matrix_parent_inverse = Matrix.Identity(4)
sph = O.get("AN_AnalSphincter")
if sph and sph.data.materials:
    w_me.materials.append(sph.data.materials[0])
# the urethra lies in the vagina's anterior wall: the placeholder takes a groove round it (its wall, the urethra open)
_um = bpy.data.meshes.new("_uo"); _um.from_pydata([tuple(p) for p in urethra_open], [], [tuple(p.vertices) for p in urethra.data.polygons])
_uob = O.new("_uo", _um); _uot = U.tree_in(_uob, None); O.remove(_uob); bpy.data.meshes.remove(_um)
_wt = U.tree_in(uwall, Mi); w_me.update()
_wp = np.array([tuple(v.co) for v in w_me.vertices]); _wn = np.array([tuple(v.normal) for v in w_me.vertices])
_Xv = np.array([tuple(v.co) for v in vag_space.data.vertices])
_gr, _ = U.make_room([vag_space], Mi, Mw, [_wt], lambda n_, P: np.zeros(len(P), bool), U.Obstacles([(_uot, 0.0005)], [], []),
                     gap=0.0008, decay=0.99, iters=200, rounds=8, expanders=[(_wp, _wn)])
vag_space.data.vertices.foreach_set("co", _gr["AN_Vagina_Space"][0].astype(np.float32).ravel()); vag_space.data.update()
report["vagina_space"]["urethral_groove_mm"] = round(float(np.linalg.norm(_gr["AN_Vagina_Space"][0] - _Xv, axis=1).max()) * 1000, 2)
report["urethra"].update({"rings": NR, "wall_rings": [j_a, j_b], "wall_from_neck_mm": round(_cum[-1] * j_a / (NR - 1) * 1000, 1),
                          "wall_to_meatus_mm": round(L_U * (1 - j_b / (NR - 1)) * 1000, 1)})
bone_min = min(t_.find_nearest(Vector(p))[3] for t_, _ in bone_trees for ring in rest_r for p in ring[::4])
report["urethra"]["bone_clearance_mm"] = round(bone_min * 1000, 2)
print("URETHRA", json.dumps(report["urethra"]), json.dumps(report["meatus"]))
done("urethra")

# ======================= 5. controls, shape keys, neighbours that give way =======================
ctl = O.get("Urinary_Controls") or O.new("Urinary_Controls", None)
if ctl.name not in col_fit.objects:
    col_fit.objects.link(ctl)
ctl.empty_display_type = 'SPHERE'; ctl.empty_display_size = 0.006
ctl.location = Mw @ (NECK + U0 * 0.03)
for prop, desc in (("Fill", "bladder: 0 empty (~70 ml) .. 1 full (500 ml); the bowel gives way"),
                   ("Void", "0 closed .. 1 voiding: bladder neck, urethra and meatus open")):
    ctl[prop] = 0.0
    ui = ctl.id_properties_ui(prop); ui.update(min=0.0, max=1.0, soft_min=0.0, soft_max=1.0, description=desc)
NS = len(stages) - 1                             # Fill_1 .. Fill_NS, each relative to the one before

# 5a. bladder: the stages on the opened, junction-cut mesh (new vertices sit in the still trigone)
Pb = np.array([tuple(v.co) for v in bladder.data.vertices])
kd = kdtree.KDTree(len(Xe))
for i, p in enumerate(Xe):
    kd.insert(Vector(p), i)
kd.balance()
src_i = []; newp = []; near_i = []
for p in Pb:
    _, i, dd = kd.find(Vector(p))
    if dd < 2e-6:                                # (open_junction round-trips through body coordinates)
        src_i.append(i)
    else:
        src_i.append(-1); newp.append(p); near_i.append(i)
n_new = len(newp)
worst = float(U.ss(TRIG_FIX, TRIG_FREE, dist_trigone(np.array(newp))).max()) if newp else 0.0
report["bladder"]["new_vertex_worst_mobility"] = round(worst, 4)
_near = iter(near_i)
nearest_src = [i if i >= 0 else next(_near) for i in src_i]     # new vertices follow their nearest old one
names_fill = [f"Fill_{k}" for k in range(1, NS + 1)]
prev = "Basis"
for k in range(1, NS + 1):
    pos = np.array([stages[k][i] if i >= 0 else Pb[j] + stages[k][nearest_src[j]] - Xe[nearest_src[j]] for j, i in enumerate(src_i)])
    kb = U.add_key(bladder, names_fill[k - 1], pos, relative=prev); prev = kb.name
    U.drive(kb, ctl, "Fill", f"f*{NS}-{k - 1}" if k > 1 else f"f*{NS}")
report["bladder"]["new_vertices_in_trigone"] = n_new

# 5b. Void: the neck funnels open (round, R_OPEN_NECK) and sinks a little; urethra round; meatus lips part
neck_idx = []
_kdb = kdtree.KDTree(len(Pb))
for i, p in enumerate(Pb):
    _kdb.insert(Vector(p), i)
_kdb.balance()
neck_idx = [_kdb.find(Vector(p))[1] for p in neck_rest]
FUNNEL_R = 0.012
void_b = {}
for i, p in enumerate(Pb):
    q = Vector(p) - NECK; h = q.dot(U0); rad = q - U0 * h; rr = rad.length
    if rr > FUNNEL_R or h > 0.012:
        continue
    w = 1 - float(U.ss(R_HOLE * 0.8, FUNNEL_R, rr))
    push = (R_OPEN_NECK - min(rr, R_OPEN_NECK)) * w
    void_b[i] = Vector(p) + (rad.normalized() * push if rr > 1e-9 else Vector())
for k, i in enumerate(neck_idx):
    void_b[i] = Vector(neck_open[k])
kb = U.add_key(bladder, "Void", None); kb.relative_key = bladder.data.shape_keys.key_blocks["Basis"]
for i, p in void_b.items():
    kb.data[i].co = p
U.drive(kb, ctl, "Void", "f")
kb = U.add_key(urethra, "Void", urethra_open); U.drive(kb, ctl, "Void", "f")
kb = U.add_key(hips, "Void", None)
for i, p in skin_void.items():
    kb.data[i].co = p
U.drive(kb, ctl, "Void", "f")
report["void"] = {"neck_open_mm": R_OPEN_NECK * 2000, "meatus_lips_apart_mm": VOID_SKIN * 2000,
                  "neck_seam_mm": round(max((Vector(neck_open[k]) - urethra.data.shape_keys.key_blocks["Void"].data[k].co).length for k in range(K_U)) * 1000, 5)}

# 5c. neighbours: the bowel (and the ureters away from the bladder) give way as the bladder fills
SOFT_N = [n_ for n_ in ("AN_Colon_Descending", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_AnalCanal", "AN_Ureter_L", "AN_Ureter_R") if n_ in O]
kid_trees = [U.tree_in(O["AN_Kidney_" + s], Mi) for s in "LR"]
loop_kd = {}
for s in "LR":
    loop_kd[s] = kdtree.KDTree(len(ure_loop[s]))
    for i, p in enumerate(ure_loop[s]):
        loop_kd[s].insert(Vector(p), i)
    loop_kd[s].balance()
def pinned(name, P):
    if name == "AN_AnalCanal":
        return np.ones(len(P), bool)                 # sits in the skin: stays (it has its own Open key)
    if name.startswith("AN_Ureter"):                 # still where it joins the bladder (and inside the kidney)
        s = name[-1]; near_b = np.array([loop_kd[s].find(Vector(p))[2] < 0.015 for p in P])
        near_k = np.array([kid_trees["LR".index(s)].find_nearest(Vector(p))[3] < 0.006 for p in P])
        return near_b | near_k
    return np.zeros(len(P), bool)
# rigid for them: bones, the inside of the skin, and the organs that stay (kidneys, adrenals)
rigid_soft = U.Obstacles([(t_, 0.0008) for t_, _ in bone_trees] +
                         [(U.tree_in(O[n_], Mi), 0.0005) for n_ in ("AN_Kidney_L", "AN_Kidney_R", "AN_Adrenal_L", "AN_Adrenal_R")] +
                         [(vag_tree, 0.0015)],
                         [(skin_tree, 0.003)], [])
stage_trees = []; stage_pts = []
for k in range(1, NS + 1):
    me_ = bpy.data.meshes.new("_st"); me_.from_pydata([tuple(bladder.data.shape_keys.key_blocks[names_fill[k - 1]].data[j].co) for j in range(len(Pb))], [], [tuple(p.vertices) for p in bladder.data.polygons])
    tmp = O.new("_st", me_); tmp.matrix_world = Matrix.Identity(4)
    stage_trees.append(U.tree_in(tmp, None))
    me_.update(); stage_pts.append((np.array([tuple(v.co) for v in me_.vertices]), np.array([tuple(v.normal) for v in me_.vertices])))
    O.remove(tmp); bpy.data.meshes.remove(me_)
mr_log = []
room, resid = U.make_room([O[n_] for n_ in SOFT_N], Mi, Mw, stage_trees, pinned, rigid_soft, gap=0.0025, log=mr_log, rounds=10, decay=0.996, iters=400,
                         carry_from=np.array([tuple(v.co) for v in bladder.data.vertices]),
                         expanders=stage_pts)
report["make_room"] = {"organs": SOFT_N, "residual_inside_mm": resid, "log": mr_log}
print("MAKE_ROOM", resid, mr_log[-6:])
for n_ in SOFT_N:
    ob = O[n_]; base = np.array([tuple(v.co) for v in ob.data.vertices]); prev = "Basis"
    moved = max(float(np.linalg.norm(room[n_][-1] - base, axis=1).max()), 0.0)
    report["make_room"].setdefault("max_move_mm", {})[n_] = round(moved * 1000 * ob.matrix_world.to_3x3().median_scale if False else moved * (1000 if ob.parent is None else 1.0), 2)
    if moved < 1e-7:
        continue
    for k in range(1, NS + 1):
        kb = U.add_key(ob, f"Bladder_{names_fill[k - 1]}", room[n_][k - 1], relative=prev); prev = kb.name
        U.drive(kb, ctl, "Fill", f"f*{NS}-{k - 1}" if k > 1 else f"f*{NS}")
# the bladder settles against what could not move further (bowel trapped against the sacrum, the pelvic bones): a
# shallow dent in its own wall where it still presses in, the trigone stays
bl_faces = [tuple(p.vertices) for p in bladder.data.polygons]
Pb_now = np.array([tuple(v.co) for v in bladder.data.vertices])
fixed_b = U.ss(TRIG_FIX, TRIG_FREE, dist_trigone(Pb_now)) < 1e-6
report["settle"] = {}
for k in range(1, NS + 1):
    kb = bladder.data.shape_keys.key_blocks[names_fill[k - 1]]
    Xk = np.array([tuple(d.co) for d in kb.data])
    trees_k = [t_ for t_, _ in bone_trees]
    for n_ in SOFT_N:
        ob = O[n_]; Pk = room[n_][k - 1]; Mk = Mi @ ob.matrix_world
        me_ = bpy.data.meshes.new("_sk"); me_.from_pydata([tuple(Mk @ Vector(p)) for p in Pk], [], [tuple(p.vertices) for p in ob.data.polygons])
        tmp = O.new("_sk", me_); trees_k.append(U.tree_in(tmp, None)); O.remove(tmp); bpy.data.meshes.remove(me_)
    Xs_, w0, w1 = U.settle(Xk, bl_faces, fixed_b, trees_k)
    kb.data.foreach_set("co", Xs_.astype(np.float32).ravel())
    report["settle"][names_fill[k - 1]] = {"pressing_in_mm": round(w0 * 1000, 2), "after_mm": round(w1 * 1000, 2)}
print("SETTLE", json.dumps(report["settle"]))
done("keys")

# ======================= 6. lighting: insides lit everywhere, smooth-lighting bakes =======================
# new organs get the outer "AO" bake (at rest); every organ / bone material lights its inside and every object without
# an inner bake gets "AO_in" (the urethra's with it open, like the anal canal's); the bowel keeps build5_passage's
AO_DIST, AO_RAYS = 0.008, 32
def setc(fill, void):
    ctl["Fill"] = fill; ctl["Void"] = void; ctl.update_tag(); bpy.context.view_layer.update()
setc(0.0, 0.0)
_tr = T.scene_bvh(); report["ao"] = {}
for ob in (bladder, urethra, uwall):
    report["ao"][ob.name + ".AO"] = T.bake_ao(ob, "AO", False, _tr, AO_DIST, AO_RAYS)
for o_ in (urethra,):
    if o_.data.color_attributes.get("AO_in"):
        o_.data.color_attributes.remove(o_.data.color_attributes["AO_in"])
meshes = [o_ for o_ in O if o_.type == 'MESH' and o_.name in bpy.context.view_layer.objects and o_ is not hips
          and o_ is not urethra and not o_.name.startswith("_")]
report["ao"].update({k + ".AO_in": v for k, v in T.light_all(meshes, _tr, skip_baked=True).items()})
for m_ in urethra.data.materials:
    T.light_insides(m_)
setc(0.0, 1.0); _tr = T.scene_bvh()
report["ao"]["AN_Urethra.AO_in(open)"] = T.bake_ao(urethra, "AO_in", True, _tr, AO_DIST, AO_RAYS)
setc(0.0, 0.0)
O.remove(old_bladder)

# ======================= 7. checks =======================
def wpos(ob):
    dg = bpy.context.evaluated_depsgraph_get(); ev = ob.evaluated_get(dg); me_ = ev.to_mesh()
    P = [ob.matrix_world @ v.co for v in me_.vertices]; ev.to_mesh_clear(); return P
def ev_tree(ob):
    dg = bpy.context.evaluated_depsgraph_get(); bm_ = bmesh.new(); bm_.from_object(ob, dg); bm_.transform(ob.matrix_world)
    bmesh.ops.remove_doubles(bm_, verts=bm_.verts, dist=1e-7)
    bnd = [e for e in bm_.edges if e.is_boundary]
    if bnd:
        bmesh.ops.holes_fill(bm_, edges=bnd, sides=0)
    t_ = BVHTree.FromBMesh(bm_); bm_.free(); return t_
chk = {}
# drivers
setc(0.5, 0.0)
chk["drivers_at_fill_0.5"] = {kb.name: round(kb.value, 3) for kb in bladder.data.shape_keys.key_blocks[1:]}
# seams in every state
_kidx = {s: None for s in "LR"}
states = [(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
chk["seams_mm"] = {}; chk["volume_ml"] = {}; chk["clipping"] = {}
pairs_ob = [bladder, urethra, uwall] + [O[n_] for n_ in SOFT_N if n_ != "AN_AnalCanal"] + [O[n_] for n_ in ("AN_HipBone_L", "AN_HipBone_R", "AN_Sacrum")] + \
           [O["AN_Kidney_L"], O["AN_Kidney_R"], vag_space]
for fill, void in states:
    setc(fill, void); tag = f"fill{fill}_void{void}"
    Pb_ = wpos(bladder); Pu_ = wpos(urethra); Ph_ = wpos(hips)
    neck = max((Pb_[neck_idx[k]] - Pu_[k]).length for k in range(K_U))
    meat = max((Ph_[slit_idx[k]] - Pu_[(NR - 1) * K_U + k]).length for k in range(K_U))
    ure = {}
    for s in "LR":
        Pr = wpos(O["AN_Ureter_" + s]); kdb = kdtree.KDTree(len(Pb_))
        for i, p in enumerate(Pb_):
            kdb.insert(p, i)
        kdb.balance()
        bm_ = bmesh.new(); bm_.from_mesh(O["AN_Ureter_" + s].data)
        loopv = {v.index for e in bm_.edges if e.is_boundary for v in e.verts}; bm_.free()
        near = [kdb.find(Pr[i])[2] for i in loopv]
        ure[s] = round(max(d_ for d_ in near if d_ < 0.004) * 1000, 4) if near else None
    chk["seams_mm"][tag] = {"neck": round(neck * 1000, 4), "meatus": round(meat * 1000, 4), "ureter_orifices": ure}
    bmv = bmesh.new(); bmv.from_object(bladder, bpy.context.evaluated_depsgraph_get()); bmv.transform(Mi @ bladder.matrix_world)
    bnd = [e for e in bmv.edges if e.is_boundary]
    bmesh.ops.holes_fill(bmv, edges=bnd, sides=0); bmesh.ops.triangulate(bmv, faces=bmv.faces[:])
    chk["volume_ml"][tag] = round(bmv.calc_volume(signed=True) * 1e6, 1); bmv.free()
    if void == 0.0 or fill == 1.0:
        trees_ = {o_.name: ev_tree(o_) for o_ in pairs_ob}
        clip = {}
        for a in pairs_ob:
            Pa = wpos(a)
            for b in pairs_ob:
                if a is b or {a.name, b.name} in ({"AN_Urethra", "AN_Urethra_Wall"},):
                    continue
                tb = trees_[b.name]; deep = []
                for p in Pa:
                    dd = tb.find_nearest(p)[3]
                    if 2e-5 < dd < 0.004 and U.inside(tb, p):
                        deep.append(dd)
                if deep:                         # vertices inside, deepest (mm, body scale)
                    clip[f"{a.name} in {b.name}"] = [len(deep), round(max(deep) * 1000, 2)]
        chk["clipping"][tag] = clip
    print("CHECK", tag, json.dumps(chk["seams_mm"][tag]), chk["volume_ml"][tag], json.dumps(chk["clipping"].get(tag, {}))[:600])
setc(0.0, 0.0)
report["checks"] = chk
# the meatus target, for moving it later (re-run with this file's target: see HANDOFF)
mt = O.get("Urethra_Meatus_Target") or O.new("Urethra_Meatus_Target", None)
if mt.name not in bpy.context.scene.collection.objects:
    bpy.context.scene.collection.objects.link(mt)
mt.empty_display_type = 'SINGLE_ARROW'; mt.empty_display_size = 0.004; mt.location = M_body
mt.rotation_euler = (math.pi, 0, 0)
U.save(SAVE_AS, OUT, report)
print("DONE")
