"""build_pelvis: the female pelvis in the working file - urinary system, bony pelvis, internal reproductive organs
(runs on build5_passage's output; formerly build_urinary.py).

Kidneys, adrenals and ureters come from the atlas (anatomy_ref.blend, read only) into Internal_Fit, under
Internal_Fit_Xform like the other organs. Everything is built in that frame ("atlas space", real-size metres): the
fit squashes it to the body (0.62 / 0.72 / 0.72), so volumes in the body are 0.32x the real ones in the report.
  0b. The atlas' male bony pelvis reshaped female (one smooth field for both hip bones and the sacrum).
  1. Kidneys: the atlas ureter already carries the renal pelvis and calyces inside the kidney; the kidney's wall is
     opened where the pelvis leaves it at the hilum (tract_lib.open_junction, sides=(0,): pierce, not join).
  2. Room for the vagina (AN_Vagina_Space, a placeholder with a vault); the bowel moved out of it.
  2b. Internal reproductive organs (repro_lib): uterus (outer surface + lumen: cervical canal and cavity), uterine
     tubes opening into the cavity and onto the ovaries, ovaries, ovarian / round / suspensory ligaments; the cervix
     opens into the vagina placeholder's vault.
  2c. Bladder (new; the atlas one is a 226-face blob pressed against the rectum): a hollow shell grown by a pressure
     simulation against the pubic bones, the pelvic walls, the inside of the skin (abdominal wall), the anterior
     vaginal wall plane and (empty) the uterus. Empty (~70 ml) to full (500 ml), in stages so every in-between is a
     real state. The trigone (neck + ureteric orifices) does not move.
  3. Ureters: their last few cm bend to the trigone's corners, entering obliquely; the lumens join the bladder.
  4. Urethra: from the bladder neck (a shared loop) to the external meatus, a sagittal slit cut into the skin of the
     vestibule (a shared loop); muscular wall sleeve. The meatus position is URETHRA_MEATUS="x,y,z" (body
     coordinates, snapped to the skin), else the "Urethra_Meatus_Target" empty of URETHRA_MEATUS_FROM=<a previous
     output .blend>, else the default; the urethra re-routes to it.
  5. Controls: empty "Pelvic_Controls" with Bladder_Fill, Void and Rectum_Fill (0..1). Filling drives the organ's
     stage keys and the corrective keys of every neighbour it presses on (the bowel, ureters and reproductive organs
     give way - the uterus tilts on its cervix - bones and the abdominal wall do not); both full has its own
     corrective (Both_Full = Bladder_Fill x Rectum_Fill); a polish pass clears the last small contacts.
  6. Insides lit + AO / AO_in bakes for every organ and bone in the file.
  7. Checks over 14 states (seams, volumes, contacts new against rest).
Run: ANATOMY_REF=anatomy_ref.blend blender -b out/blend/Hips_build5_passage.blend --python scripts/build_pelvis.py \
     -- <out.blend> <report folder>      (URINARY_STOP=<step> saves after a step: pelvis kidneys vagina repro ...)"""
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

# ======================= 0b. the pelvis made female =======================
# the atlas' bones are a male pelvis (narrow subpubic arch, ischial tuberosities ~85 mm apart, heart-shaped inlet,
# curved sacrum). One smooth displacement field (urinary_lib.female_pelvis_field) for both hip bones and the sacrum -
# so the symphysis and the sacroiliac joints stay matched - widens the arch and the outlet, moves the ischial spines
# out, widens the brim and the sacrum, straightens the lower sacrum / coccyx back and lowers the iliac wings a little.
# FEMALE_PELVIS (default 1) scales it (0 = keep the atlas pelvis); it is backed off if a bone would come within
# PELVIS_SKIN mm (body) of the skin.
FEMALE_PELVIS = float(os.environ.get("FEMALE_PELVIS", "1")); PELVIS_SKIN = 0.0025
PELVIS = [n_ for n_ in ("AN_HipBone_L", "AN_HipBone_R", "AN_Sacrum") if n_ in O]
def _bone_pts(n_):
    M_ = Mi @ O[n_].matrix_world; return np.array([tuple(M_ @ v.co) for v in O[n_].data.vertices])
_P0 = {n_: _bone_pts(n_) for n_ in PELVIS}
report["pelvis"] = {"before": U.pelvis_metrics(_P0)}
_skw = U.tree_in(hips, None, cap=False)
amt = FEMALE_PELVIS
while FEMALE_PELVIS > 0:
    _P1 = {n_: _P0[n_] + U.female_pelvis_field(_P0[n_], amt) for n_ in PELVIS}
    clear = min(_skw.find_nearest(Mw @ Vector(p))[3] for n_ in PELVIS for p in _P1[n_][::2])
    if clear >= PELVIS_SKIN or amt < 0.1:
        break
    amt *= 0.8
if FEMALE_PELVIS > 0:
    for n_ in PELVIS:
        M_ = (Mi @ O[n_].matrix_world).inverted()
        O[n_].data.vertices.foreach_set("co", np.array([tuple(M_ @ Vector(p)) for p in _P1[n_]], np.float32).ravel())
        O[n_].data.update()
    bpy.context.view_layer.update()
    report["pelvis"].update({"after": U.pelvis_metrics(_P1), "amount": round(amt, 3), "skin_clearance_mm": round(clear * 1000, 2)})
print("PELVIS", json.dumps(report["pelvis"]))
done("pelvis")

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
VAG_HW, VAG_HT = 0.013, 0.005
VAG_UP = float(os.environ.get("VAG_UP_MM", "40")) / 1000   # past the neck up the bladder base to the vault (~75 mm in all)
VAULT_HW, VAULT_HT = 0.016, 0.011                # the vault round the cervix     # half width, half thickness (walls, collapsed), upper part length
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
    # the vault: the upper part widens into the fornices round the cervix (VAULT_HT / VAULT_HW at the top)
    vlt = float(U.ss(0.55, 0.85, s_))
    w_ = (VAG_HW * (0.7 + 0.3 * float(U.ss(0.0, 0.2, s_))) * (1 - vlt) + VAULT_HW * vlt) * rnd
    h_ = (VAG_HT * (0.6 + 0.4 * float(U.ss(0.0, 0.1, s_))) * (1 - vlt) + VAULT_HT * vlt) * rnd
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
bone_trees_v = bone_trees + [(vag_tree, 0.001)]
done("vagina")

# ======================= 2b. internal reproductive organs =======================
# Uterus (anteverted, anteflexed): the cervix runs from the external os, inside the vagina's vault, square to the
# vagina (its axis C_AX), the body bends forward on it by UT_FLEX_DEG and rests on the bladder. Two surfaces, like the
# urethra and its wall: AN_Uterus (the outer surface) and AN_Uterus_Lumen (cervical canal + the flat triangular
# cavity); they meet at the external os. The uterine tubes open into the cavity's two upper corners (cornua), pass
# through the wall, run out over the ovaries and end in fimbriated funnels, open to the pelvis. Ovaries on the side
# walls; ovarian, round and suspensory ligaments as cords. The cervix pierces the vagina placeholder's vault: the
# passage runs vagina -> os -> canal -> cavity -> tubes -> fimbriae. Real-size measures (atlas space).
import repro_lib as R
UT_FLEX = math.radians(float(os.environ.get("UT_FLEX_DEG", "35")))
_nv = len(vc); VC = vc[int(0.86 * (_nv - 1))]   # the vault's centre
C_AX = VAG_M.copy()                               # cervix: square to the vagina (anteversion)
B_AX = (Matrix.Rotation(UT_FLEX, 3, X_) @ C_AX).normalized()   # body: bent forward (anteflexion)
N_AX = B_AX.cross(X_).normalized()                # the body's thickness direction (front - back)
P_os = VC + C_AX * 0.002; P_io = P_os + C_AX * 0.025
ut_el = [R.capsule_el(P_os, P_io + C_AX * 0.004, 0.0105),
         R.ellipsoid_el(P_io + B_AX * 0.020, B_AX, X_, (0.024, 0.021, 0.015)),
         R.ellipsoid_el(P_io + B_AX * 0.033, B_AX, X_, (0.012, 0.024, 0.014))]
ut_me = R.decimate(R.metaball_mesh("AN_Uterus", ut_el, resolution_mm=0.8), 0.3)
# the lumen: the cervical canal (a stub through the os, trimmed at the join) and the flat triangular cavity
APX = P_io + B_AX * 0.004; CORN = {s_: P_io + B_AX * 0.031 + X_ * (0.0145 if s_ == "L" else -0.0145) for s_ in "LR"}
lu_el = [R.capsule_el(P_os - C_AX * 0.004, APX, 0.0022)]
_tri = [APX, CORN["L"], CORN["R"]]
for i_ in range(9):
    for j_ in range(9 - i_):
        a_, b_ = i_ / 8, j_ / 8; p_ = _tri[0] + (_tri[1] - _tri[0]) * a_ + (_tri[2] - _tri[0]) * b_
        lu_el.append(R.ellipsoid_el(p_, B_AX, X_, (0.0026, 0.0026, 0.0011)))
lu_me = R.decimate(R.metaball_mesh("AN_Uterus_Lumen", lu_el, resolution_mm=0.35), 0.25)
REPRO = {}
def _new(name, me_, mat_name, rgb):
    ob_ = O.get(name)
    if ob_ is not None:
        O.remove(ob_)
    ob_ = O.new(name, me_); col_fit.objects.link(ob_); ob_.parent = E; ob_.matrix_parent_inverse = Matrix.Identity(4)
    me_.materials.append(R.tinted(old_bladder.data.materials[0], mat_name, rgb)); REPRO[name] = ob_; return ob_
uterus = _new("AN_Uterus", ut_me, "AN_Uterus", (0.80, 0.36, 0.40))
ulumen = _new("AN_Uterus_Lumen", lu_me, "AN_Uterus_Lumen", (0.86, 0.45, 0.50))
# ovaries: on the side walls (the ovarian fossa), level with the top of the uterus, behind it
_ht = U.tree_in(O["AN_HipBone_L"], Mi), U.tree_in(O["AN_HipBone_R"], Mi)
OVC = {}; OV_LONG = Vector((0, 0.35, 1)).normalized(); OV_AP = OV_LONG.cross(X_).normalized()
for s_, t_ in zip("LR", _ht):
    sg = 1 if s_ == "L" else -1
    o_ = Vector((0, P_io.y + 0.014, P_io.z + 0.022)); hit = t_.ray_cast(o_, X_ * sg, 0.2)
    wall = hit[0] if hit[0] is not None else o_ + X_ * sg * 0.065
    OVC[s_] = Vector((wall.x - sg * 0.0095, o_.y, o_.z))
    _new(f"AN_Ovary_{s_}", R.ellipsoid_mesh(f"AN_Ovary_{s_}", OVC[s_], (OV_LONG, OV_AP, X_), (0.017, 0.010, 0.0065), bumps=0.04, seed=3 if s_ == "L" else 7),
         "AN_Ovary", (0.93, 0.80, 0.74))
# uterine tubes: in from the cavity's corner, out through the wall, over the ovary, a fimbriated funnel facing it
TUBE = {}
for s_ in "LR":
    sg = 1 if s_ == "L" else -1; Cc = CORN[s_]; ov = OVC[s_]
    pts = [Cc - X_ * sg * 0.002, Cc + X_ * sg * 0.014, Cc + X_ * sg * 0.036 + B_AX * 0.004,
           ov + Vector((sg * 0.004, -0.004, 0.019)), ov + Vector((sg * 0.009, 0.012, 0.010)), ov + Vector((sg * 0.006, 0.016, -0.002))]
    cl, L_t = R.resample(R.catmull([Vector(p) for p in pts], 24), 0.0012)
    s_isth = ((Cc + X_ * sg * 0.036) - Cc).length / L_t
    rad = lambda s, s_isth=s_isth: 0.0016 + 0.0016 * float(U.ss(s_isth, s_isth + 0.3, s)) + 0.0040 * float(U.ss(0.86, 1.0, s))
    fr = lambda th, s: 1.0 + float(U.ss(0.95, 1.0, s)) * 0.55 * max(0.0, math.cos(7 * th)) ** 3
    TUBE[s_] = _new(f"AN_Uterine_Tube_{s_}", R.sweep(f"AN_Uterine_Tube_{s_}", cl, rad, segs=20, cap_start=True, cap_end=False, rim_fn=fr),
                    "AN_Uterine_Tube", (0.88, 0.52, 0.52))
    report.setdefault("reproductive", {}).setdefault("tube_length_mm", {})[s_] = round(L_t * 1000, 1)
# ligaments (cords, ends sunk 1 mm into what they hold): ovarian (ovary -> uterus behind the tube), round (uterus in
# front of the tube -> the deep inguinal ring), suspensory (ovary -> the pelvic brim)
CORDS = {}
for s_, t_ in zip("LR", _ht):
    sg = 1 if s_ == "L" else -1; ov = OVC[s_]
    u_back = P_io + B_AX * 0.028 + X_ * sg * 0.020 - N_AX * 0.006
    u_front = P_io + B_AX * 0.028 + X_ * sg * 0.020 + N_AX * 0.006
    ov_low = ov - OV_LONG * 0.016 - X_ * sg * 0.002; ov_up = ov + OV_LONG * 0.016
    # deep inguinal ring: behind the abdominal wall, above the pubis, out to the side
    o_ = Vector((sg * 0.050, -0.085, P_io.z + 0.030)); h_ = skin_tree.ray_cast(o_, Vector((0, -1, 0)), 0.2)
    ring = (h_[0] + Vector((0, 1, 0)) * (WALL + 0.002)) if h_[0] is not None else o_ + Vector((0, -0.04, 0))
    brim = t_.find_nearest(ov_up + Vector((sg * 0.012, 0.004, 0.022)))[0]
    brim = brim + (ov_up - brim).normalized() * 0.0015
    for nm_, pts, r_, a_org, b_org in (
            (f"AN_Ovarian_Ligament_{s_}", [ov_low, (ov_low + u_back) / 2 - Vector((0, 0, 0.003)), u_back], 0.0016, f"AN_Ovary_{s_}", "AN_Uterus"),
            (f"AN_Round_Ligament_{s_}", [u_front, u_front + Vector((sg * 0.022, -0.010, 0.010)), (u_front + ring) / 2 + Vector((sg * 0.010, 0, 0.012)), ring], 0.0019, "AN_Uterus", None),
            (f"AN_Suspensory_Ligament_{s_}", [ov_up, (ov_up + brim) / 2 + Vector((sg * 0.003, 0, 0.002)), brim], 0.0022, f"AN_Ovary_{s_}", None)):
        cl, L_c = R.resample(R.catmull([Vector(p) for p in pts], 16), 0.0015)
        _new(nm_, R.sweep(nm_, cl, lambda s, r_=r_: r_, segs=12, cap_start=True, cap_end=True), "AN_Ligament", (0.92, 0.86, 0.80))
        CORDS[nm_] = (a_org, Vector(pts[0]), b_org, Vector(pts[-1]), len(cl), 12)
bpy.context.view_layer.update()
# joins: the external os (outer surface <-> lumen), each tube through the wall and into the cavity, the cervix into
# the vagina's vault
J = {}
J["os"] = T.open_junction(uterus, ulumen, centre=Mw @ P_os, sides=(0, 1), outside=(1,))
for s_ in "LR":
    sg = 1 if s_ == "L" else -1
    J["wall_" + s_] = T.open_junction(uterus, TUBE[s_], centre=Mw @ (CORN[s_] + X_ * sg * 0.008), sides=(0,))
    J["ostium_" + s_] = T.open_junction(TUBE[s_], ulumen, centre=Mw @ CORN[s_])
J["vagina"] = T.open_junction(vag_space, uterus, centre=Mw @ (P_os + C_AX * 0.008), sides=(0,))
report["reproductive"]["joins"] = {k: {kk: v[kk] for kk in ("curve_mm", "patch_radius_mm", "killed_faces")} for k, v in J.items()}
report["reproductive"]["uterus_dims_mm"] = [round(float(x) * 1000, 1) for x in (np.array([tuple(v.co) for v in ut_me.vertices]).max(0) - np.array([tuple(v.co) for v in ut_me.vertices]).min(0))]
# the tubes' fimbriated ends, the ovaries and the suspensory ligaments lie against the pelvic side wall: they give
# way to the bone (2 mm clear), the tubes held where they leave the uterus, the ligaments at their wall ends
report["reproductive"]["off_bone_mm"] = {}
for n_, o_ in REPRO.items():
    if n_ in ("AN_Uterus", "AN_Uterus_Lumen"):
        continue
    P_ = np.array([tuple(v.co) for v in o_.data.vertices])
    if n_.startswith("AN_Uterine_Tube"):
        fx = np.linalg.norm(P_ - np.array(CORN[n_[-1]]), axis=1) < 0.015
    elif n_ in CORDS and CORDS[n_][2] is None:
        fx = np.linalg.norm(P_ - np.array(CORDS[n_][3]), axis=1) < 0.003
    else:
        fx = np.zeros(len(P_), bool)
    Xs_, w0_, w1_ = U.settle(P_, [tuple(p.vertices) for p in o_.data.polygons], fx, [t_ for t_, _ in bone_trees], gap=0.002, inward=False)
    if w0_ > 0:
        o_.data.vertices.foreach_set("co", Xs_.astype(np.float32).ravel()); o_.data.update()
        report["reproductive"]["off_bone_mm"][n_] = [round(w0_ * 1000, 2), round(w1_ * 1000, 2)]
bpy.context.view_layer.update()
vag_tree = U.tree_in(vag_space, Mi)
REPRO_SOLID = [REPRO[n_] for n_ in REPRO if n_ != "AN_Uterus_Lumen"]
# the bowel (and the ureters) make room for them, as for the vagina (again after the uterus has settled on the
# empty bladder, 2c)
_kidt = {s_: U.tree_in(O["AN_Kidney_" + s_], Mi) for s_ in "LR"}
def _pin_r(n_, P):
    if n_.startswith("AN_Ureter"):
        return np.array([_kidt[n_[-1]].find_nearest(Vector(p))[3] < 0.006 for p in P])
    return _pin_v(n_, P)
def bowel_room(tag, extra=()):
    global repro_tree
    repro_tree = U.union_tree(REPRO_SOLID, Mi)
    _before = {n_: np.array([tuple(v.co) for v in O[n_].data.vertices]) for n_ in BOWEL + ["AN_Ureter_L", "AN_Ureter_R"]}
    _rp_pts = np.concatenate([np.array([tuple((Mi @ o_.matrix_world) @ v.co) for v in o_.data.vertices]) for o_ in REPRO_SOLID])
    _rp_nrm = np.concatenate([np.array([tuple(((Mi @ o_.matrix_world).to_3x3() @ v.normal).normalized()) for v in o_.data.vertices]) for o_ in REPRO_SOLID])
    room_g, res_g = U.make_room([O[n_] for n_ in _before], Mi, Mw, [repro_tree], _pin_r,
                                U.Obstacles([(t_, 0.0008) for t_, _ in bone_trees] + [(vag_tree, 0.001)] + list(extra), [(skin_tree, 0.003)], []),
                                gap=0.002, decay=0.996, iters=400, rounds=10, expanders=[(_rp_pts, _rp_nrm)])
    rep_ = report["reproductive"].setdefault(tag, {"bowel_moved_mm": {}})
    for n_ in _before:
        ob = O[n_]; D = room_g[n_][0] - _before[n_]; mv = float(np.linalg.norm(D, axis=1).max())
        rep_["bowel_moved_mm"][n_] = round(mv * (1.0 if ob.parent else 1000.0), 2)
        if mv < 1e-9:
            continue
        for kb_ in (ob.data.shape_keys.key_blocks if ob.data.shape_keys else []):
            K = np.array([tuple(d.co) for d in kb_.data]) + D; kb_.data.foreach_set("co", K.astype(np.float32).ravel())
        ob.data.vertices.foreach_set("co", room_g[n_][0].astype(np.float32).ravel()); ob.data.update()
    rep_["bowel_residual_mm"] = res_g
    bpy.context.view_layer.update()
bowel_room("bowel_room")
bone_trees_v = bone_trees + [(vag_tree, 0.001), (repro_tree, 0.0015)]
# ---- how they move when the bladder / rectum fill: the cervix is held (in the vault, by its ligaments); the body
# tilts on it about a transverse axis at the isthmus (no turn below 10 mm up the cervix, all of it above 28 mm), the
# tubes and ovaries follow part of the way (less the further they are from the uterus), each ligament between its two
# ends (a pelvic-wall end stays). The turn at each stage is the least that clears the filling organ (repro_tilt),
# stopped by bone; whatever is left is soft give (make_room).
_ut_P = np.array([tuple(v.co) for v in uterus.data.vertices])
_ut_tree = U.tree_in(uterus, Mi)
_Pos = np.array(P_os); _Cax = np.array(C_AX); _piv = np.array(P_os + C_AX * 0.010)
def tilt_disp(P, th):
    P = np.atleast_2d(P); a = th * U.ss(0.010, 0.028, (P - _Pos) @ _Cax); rel = P - _piv
    c, s = np.cos(a), np.sin(a); q = rel.copy()
    q[:, 1] = rel[:, 1] * c - rel[:, 2] * s; q[:, 2] = rel[:, 1] * s + rel[:, 2] * c
    return q - rel
BACK = 1.0 if tilt_disp(np.array(P_io + B_AX * 0.045), 0.1)[0, 1] > 0 else -1.0      # +y: toward the sacrum
REST_R = {n_: np.array([tuple(v.co) for v in o_.data.vertices]) for n_, o_ in REPRO.items()}
_att = {}
for n_, P_ in REST_R.items():
    if n_ in ("AN_Uterus", "AN_Uterus_Lumen"):
        _att[n_] = np.ones(len(P_))
    elif n_ not in CORDS:
        _att[n_] = 1 - 0.6 * U.ss(0.0, 0.04, np.array([_ut_tree.find_nearest(Vector(p))[3] for p in P_]))
_cord_s = {}
for n_, (a_org, pa, b_org, pb, nr, sg_) in CORDS.items():
    s_ = np.zeros(len(REST_R[n_])); s_[:nr * sg_] = np.repeat(np.arange(nr) / (nr - 1), sg_); s_[nr * sg_ + 1:] = 1.0
    _cord_s[n_] = s_
def repro_disp(th):
    """{name: displacement of each rest vertex} for a tilt th (rad)"""
    D = {n_: _att[n_][:, None] * tilt_disp(REST_R[n_], th) for n_ in _att}
    def end(org, p):
        return np.zeros(3) if org is None else D[org][int(np.argmin(np.linalg.norm(REST_R[org] - np.array(p), axis=1)))]
    for n_, (a_org, pa, b_org, pb, nr, sg_) in CORDS.items():
        s_ = _cord_s[n_]; D[n_] = np.outer(1 - s_, end(a_org, pa)) + np.outer(s_, end(b_org, pb))
    return D
def repro_tilt(exp_tree, signs, blockers, gap=0.0025, max_deg=60):
    """the least turn (rad) that takes the uterus `gap` clear of exp_tree, within the blockers [(tree, min distance)];
    if none clears, the one that leaves the least. Returns (angle, depth left)"""
    Ps = _ut_P[::3]; best = (1e9, 0.0, 0.0)
    for sg in signs:
        for deg in range(0, max_deg + 1):
            th = sg * math.radians(deg); Q = Ps + tilt_disp(Ps, th)
            if deg and any(float(U.signed_dist(t_, Q, band=0.012)[0].min()) < g_ for t_, g_ in blockers):
                break
            pen = max(0.0, gap - float(U.signed_dist(exp_tree, Q, band=0.08)[0].min()))
            best = min(best, (round(pen, 4), deg, th))
            if pen < 1e-4:
                break
    return best[2], best[0]
def repro_pre(thetas):
    """make_room pre_stage: at stage k the tilt goes from thetas[k-1] to thetas[k]"""
    Ds = [repro_disp(t_) for t_ in [0.0] + list(thetas)]
    def pre(k, X, own_n, own_i):
        X = X.copy(); h = np.zeros(len(X), bool)
        for n_ in REPRO:
            m = own_n == n_
            if m.any():
                X[m] += Ds[k + 1][n_][own_i[m]] - Ds[k][n_][own_i[m]]; h |= m
        return X, h
    return pre
def pin_repro(name, P):
    if name in ("AN_Uterus", "AN_Uterus_Lumen"):        # firm: it only tilts (pre_stage), it is never pushed out of
        return np.ones(len(P), bool)                    # shape - its neighbours give way to it, or dent
    if name in CORDS and CORDS[name][2] is None:        # a ligament's end on the pelvic / abdominal wall
        return np.linalg.norm(P - np.array(CORDS[name][3]), axis=1) < 0.003
    return np.zeros(len(P), bool)
REPRO_N = list(REPRO)
report["reproductive"]["vertices"] = {n_: len(o_.data.vertices) for n_, o_ in REPRO.items()}
print("REPRO", json.dumps(report["reproductive"]))
done("repro")

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
obst = U.Obstacles(bone_trees + [(vag_tree, 0.001)], [(skin_tree, WALL)], planes)   # filling: the uterus gives way (5c)
Xs, Fb = U.icosphere(5); r0 = 0.010
Xb = np.array(NECK + U0 * r0) + Xs * r0
dN = np.linalg.norm(Xb - np.array(NECK), axis=1)
mob = U.ss(0.003, 0.008, dN)                     # the neck stays put
log = []
lid = (tuple(NECK + Vector((0, 0, EMPTY_LID))), (0, 0, -1), None)
# the empty bladder fits under the uterus (its top concave where the uterus lies on it), EMPTY_ON_UTERUS=0: grown
# without it, the uterus then tilts back onto it
EMPTY_ON_UTERUS = os.environ.get("EMPTY_ON_UTERUS", "1") == "1"
obst_empty = U.Obstacles(bone_trees + [(vag_tree, 0.001)] + ([(repro_tree, 0.0015)] if EMPTY_ON_UTERUS else []) + soft_trees,
                         [(skin_tree, WALL)], planes + [lid])
_dbg = []
_v0, _ = obst_empty.project(Xb, np.ones(len(Xb), bool), debug=_dbg); print("SEED", _dbg, [b.name for b in bones] + SOFT)
report["seed_violations_mm"] = round(float(np.abs(_v0 - Xb).max()) * 1000, 2)
Xe, Ve = U.grow(Xb, Fb, mob, V_STAGES[0] * 1e-6, obst_empty, log=log)
report["bladder"] = {"neck": list(NECK), "meatus_body": list(M_body), "empty_ml": round(Ve * 1e6, 1), "log_empty": log}
print("BLADDER empty", Ve * 1e6, log[-3:])
# the uterus comes to rest on the empty bladder: it tilts back on its cervix until it clears the bladder's lid (its
# anteflexion is what is left), the bowel gives way to it again; whatever the tilt cannot clear, the bladder settles
_blk0 = [(t_, 0.002) for t_, _ in bone_trees] + [(U.union_tree([O["AN_Rectum"], O["AN_Rectum_LowerAmpulla"]], Mi), -0.004)]
th0, left0 = repro_tilt(U.tree_np(Xe, Fb), (BACK,), _blk0, gap=0.0015)
D0 = repro_disp(th0)
for n_, o_ in REPRO.items():
    P_ = REST_R[n_] + D0[n_]; o_.data.vertices.foreach_set("co", P_.astype(np.float32).ravel()); o_.data.update()
    REST_R[n_] = P_
for n_, (a_org, pa, b_org, pb, nr, sg_) in list(CORDS.items()):
    CORDS[n_] = (a_org, pa + Vector(D0[n_][nr * sg_]), b_org, pb + Vector(D0[n_][nr * sg_ + 1]), nr, sg_)
bpy.context.view_layer.update()
_ut_P = REST_R["AN_Uterus"].copy(); _ut_tree = U.tree_in(uterus, Mi)
bowel_room("bowel_room_after_rest_tilt", [(U.tree_np(Xe, Fb), 0.0015)])
Xe_s, w0_, w1_ = U.settle(Xe, Fb, mob < 1e-6, [repro_tree])
Ve2 = U.volume(Xe_s, Fb); Xe = Xe_s
report["reproductive"]["rest_tilt"] = {"deg": round(math.degrees(abs(th0)), 1), "left_mm": round(left0 * 1000, 2),
                                       "bladder_settled_mm": [round(w0_ * 1000, 2), round(w1_ * 1000, 2)], "empty_ml_after": round(Ve2 * 1e6, 1)}
print("REST_TILT", json.dumps(report["reproductive"]["rest_tilt"]))

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
loops = U.boundary_loops(bm)
if len(loops) != 1:                              # (development: keep the failed empty shell to look at)
    me_ = bpy.data.meshes.new("_bad_empty"); me_.from_pydata([tuple(p) for p in Xe], [], [tuple(f) for f in Fb])
    ob_ = O.new("_bad_empty", me_); col_fit.objects.link(ob_); ob_.parent = E; U.save(SAVE_AS, OUT, report)
    raise AssertionError(f"neck cut: {len(loops)} loops")
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
_bowel = [U.tree_in(O[n_], Mi) for n_ in ("AN_Colon_Descending", "AN_Rectum") if n_ in O] + [vag_tree, repro_tree]
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
# each stage: grown freely, the uterus tilts back on its cervix to clear it (the least turn, stopped by bone or 8 mm
# into the rectum, at most 40 deg), then the stage is grown again with the tilted uterus solid - the bladder takes
# its volume elsewhere (the dome rising along the abdominal wall) instead of being dented by a third
rect_rest_tree = U.union_tree([O["AN_Rectum"], O["AN_Rectum_LowerAmpulla"]], Mi)
_blk = [(t_, 0.002) for t_, _ in bone_trees] + [(rect_rest_tree, -0.008)]   # the rectum gives way up to 8 mm
stages = [Xe]; th_grow = []; obst_st = []
for V_ in V_STAGES[1:]:
    log = []
    Xn, Vn = U.grow(stages[-1], Fb, mob_f, V_ * 1e-6, obst, log=log)
    th_, left_ = repro_tilt(U.tree_np(Xn, Fb), (BACK,), _blk, max_deg=40)
    if th_grow and abs(th_grow[-1]) > abs(th_):
        th_ = th_grow[-1]
    th_grow.append(th_)
    # (the tubes, ovaries and ligaments with it, where the tilt takes them: tied to the uterus and the pelvic wall,
    # they cannot get out of the bladder's way, so it grows round them)
    # (two shells - the uterus, and the rest: one shell of overlapping organs - the tubes run through the uterine
    # wall - would fool the ray-parity inside test)
    _Dt = repro_disp(th_); _pos = {o_.name: REST_R[o_.name] + _Dt[o_.name] for o_ in REPRO_SOLID}
    _rp_t = [(U.union_tree([uterus], Mi, positions=_pos), 0.0015),
             (U.union_tree([o_ for o_ in REPRO_SOLID if o_ is not uterus], Mi, positions=_pos), 0.0015)]
    obst_st.append(U.Obstacles(obst.solid + _rp_t, obst.container, obst.planes))
    Xn, Vn = U.grow(stages[-1], Fb, mob_f, V_ * 1e-6, obst_st[-1], log=log)
    stages.append(Xn); print("BLADDER stage", V_, Vn * 1e6, round(math.degrees(abs(th_)), 1), left_, log[-2:])
    report["bladder"].setdefault("stages_ml", []).append(round(Vn * 1e6, 1))
    report["bladder"].setdefault("uterus_tilt_deg", []).append(round(math.degrees(abs(th_)), 1))
# the crease round the still trigone (the still patch against the expanding wall): Taubin passes (smoothing that does
# not shrink) on each stage, only in the band where the wall's mobility ramps up; a stage keeps them only if its
# volume changes < 2 %
CREASE_IT = int(os.environ.get("CREASE_IT", "20"))
_dT = dist_trigone(Xe); _wc = U.ss(TRIG_FIX, TRIG_FIX + 0.006, _dT) * (1 - U.ss(0.034, 0.046, _dT))
_Lc = U.laplacian_matrix(U.neighbours(len(Xe), Fb)); report["bladder"]["crease_smoothing_ml"] = []
for k in range(1, len(stages)):
    X_ = stages[k].copy(); v0 = U.volume(X_, Fb)
    for it in range(CREASE_IT):
        for lam in (0.55, -0.58):
            X_ += U.umbrella(X_, _Lc) * (lam * _wc)[:, None]
    X_, _ = obst_st[k - 1].project(X_, mob_f > 1e-6); v1 = U.volume(X_, Fb)
    if abs(v1 - v0) / v0 > 0.005:                # a short regrow puts back what the smoothing took
        X_, v1 = U.grow(X_, Fb, mob_f, v0, obst_st[k - 1], iters=400)
    keep = abs(v1 - v0) / v0 < 0.02
    if keep:
        stages[k] = X_
    report["bladder"]["crease_smoothing_ml"].append([round(v0 * 1e6, 1), round(v1 * 1e6, 1), keep])
print("CREASE", report["bladder"]["crease_smoothing_ml"])
report["bladder"]["orifices"] = {s: list(ORI[s][0]) for s in "LR"}
report["bladder"]["opening_centres"] = {s: [round(x, 5) for x in np.mean(ure_loop[s], axis=0)] for s in "LR"}

if STOP == "bladder_raw":                         # development: every stage as its own object
    for k, Xk in enumerate(stages):
        me_ = bpy.data.meshes.new(f"_bl{k}"); me_.from_pydata([tuple(p) for p in Xk], [], [tuple(f) for f in Fb])
        ob_ = O.new(f"_bladder_stage{k}", me_); col_fit.objects.link(ob_); ob_.parent = E
    for k, th_ in enumerate(th_grow):                # and the uterus as tilted for each stage
        P_ = REST_R["AN_Uterus"] + tilt_disp(REST_R["AN_Uterus"], th_)
        me_ = bpy.data.meshes.new(f"_ut{k + 1}"); me_.from_pydata([tuple(p) for p in P_], [], [tuple(p.vertices) for p in uterus.data.polygons])
        ob_ = O.new(f"_uterus_stage{k + 1}", me_); col_fit.objects.link(ob_); ob_.parent = E
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
for _old in ("Urinary_Controls",):
    if O.get(_old):
        O.remove(O[_old])
ctl = O.get("Pelvic_Controls") or O.new("Pelvic_Controls", None)
if ctl.name not in col_fit.objects:
    col_fit.objects.link(ctl)
ctl.empty_display_type = 'SPHERE'; ctl.empty_display_size = 0.006
ctl.location = Mw @ (NECK + U0 * 0.03)
for prop, desc in (("Bladder_Fill", "bladder: 0 empty (~70 ml) .. 1 full (500 ml); the bowel gives way"),
                   ("Void", "0 closed .. 1 voiding: bladder neck, urethra and meatus open"),
                   ("Rectum_Fill", "rectum: 0 resting .. 1 full (+200 ml); the bowel and the bladder give way")):
    ctl[prop] = 0.0
    ui = ctl.id_properties_ui(prop); ui.update(min=0.0, max=1.0, soft_min=0.0, soft_max=1.0, description=desc)
NS = len(stages) - 1                             # Fill_1 .. Fill_NS, each relative to the one before
def setc(fill, void, rect=0.0):
    ctl["Bladder_Fill"] = fill; ctl["Void"] = void; ctl["Rectum_Fill"] = rect; ctl.update_tag(); bpy.context.view_layer.update()

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
    U.drive(kb, ctl, "Bladder_Fill", f"f*{NS}-{k - 1}" if k > 1 else f"f*{NS}")
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
SOFT_N = [n_ for n_ in ("AN_Colon_Descending", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_AnalCanal", "AN_Ureter_L", "AN_Ureter_R") if n_ in O] + REPRO_N
kid_trees = [U.tree_in(O["AN_Kidney_" + s], Mi) for s in "LR"]
loop_kd = {}
for s in "LR":
    loop_kd[s] = kdtree.KDTree(len(ure_loop[s]))
    for i, p in enumerate(ure_loop[s]):
        loop_kd[s].insert(Vector(p), i)
    loop_kd[s].balance()
def pinned(name, P):
    if name in REPRO:
        return pin_repro(name, P)
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
# the uterus tilts back on its cervix as the dome rises under it, as found while growing the stages (2c)
th_f = list(th_grow)
for k in range(NS):
    report.setdefault("uterus_tilt", {})[names_fill[k]] = {"deg": round(math.degrees(abs(th_f[k])), 1)}
print("TILT fill", json.dumps(report["uterus_tilt"]))
mr_log = []
room, resid = U.make_room([O[n_] for n_ in SOFT_N], Mi, Mw, stage_trees, pinned, rigid_soft, gap=0.0025, log=mr_log, rounds=10, decay=0.996, iters=400,
                         carry_from=np.array([tuple(v.co) for v in bladder.data.vertices]),
                         expanders=stage_pts, pre_stage=repro_pre(th_f))
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
        U.drive(kb, ctl, "Bladder_Fill", f"f*{NS}-{k - 1}" if k > 1 else f"f*{NS}")
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
        if n_ == "AN_Uterus_Lumen":
            continue
        ob = O[n_]; Pk = room[n_][k - 1]; Mk = Mi @ ob.matrix_world
        me_ = bpy.data.meshes.new("_sk"); me_.from_pydata([tuple(Mk @ Vector(p)) for p in Pk], [], [tuple(p.vertices) for p in ob.data.polygons])
        tmp = O.new("_sk", me_); trees_k.append(U.tree_in(tmp, None)); O.remove(tmp); bpy.data.meshes.remove(me_)
    Xs_, w0, w1 = U.settle(Xk, bl_faces, fixed_b, trees_k)
    kb.data.foreach_set("co", Xs_.astype(np.float32).ravel())
    report["settle"][names_fill[k - 1]] = {"pressing_in_mm": round(w0 * 1000, 2), "after_mm": round(w1 * 1000, 2)}
print("SETTLE", json.dumps(report["settle"]))

# ======================= 5d. the rectum fills (Rectum_Fill) =======================
# the rectum and the lower ampulla are one reservoir: welded, its two openings (onto the sigmoid, onto the anal canal)
# capped and held still, pressure-grown in RS stages to +RECTUM_ADD_ML (real ml) against the bones (sacrum, coccyx),
# the abdominal wall and the vagina's room. The anal canal stays closed. Then the sigmoid, the bladder (trigone still)
# and the ureters give way to it, as they do to the bladder.
RECTUM_ADD_ML = float(os.environ.get("RECTUM_ADD_ML", "200")); RS = 3
REC = [O["AN_Rectum"], O["AN_Rectum_LowerAmpulla"]]
Xr_, Fr_, idx_r = U.weld_union(REC, Mi)
Xr, Tr, loops_r, cen_r = U.cap_holes(Xr_, Fr_)
if U.volume(Xr, Tr) < 0:
    Tr = [t[::-1] for t in Tr]
V_r0 = U.volume(Xr, Tr)
_rim = [Xr[i] for lp in loops_r for i in lp]
_kr = kdtree.KDTree(len(_rim))
for i, p in enumerate(_rim):
    _kr.insert(Vector(p), i)
_kr.balance()
d_rim = np.array([_kr.find(Vector(p))[2] for p in Xr])
mob_r = U.ss(0.008, 0.025, d_rim); mob_r[cen_r] = 0.0
# on its own the full rectum does not move the bladder or the ureters (in a woman the vagina / uterus take that
# push): they are solid for it; the bladder gives way only when both are full (5e)
obst_r = U.Obstacles([(t_, 0.0035) for t_, _ in bone_trees] + [(vag_tree, 0.0015), (U.tree_in(bladder, Mi), 0.002)] +
                     [(U.tree_in(O["AN_Ureter_" + s_], Mi, cap=False), 0.0012) for s_ in "LR"], [(skin_tree, 0.005)], [])
stages_r = [Xr]; report["rectum_fill"] = {"rest_ml": round(V_r0 * 1e6, 1), "stages_ml": []}
for k in range(1, RS + 1):
    Xn, Vn = U.grow(stages_r[-1], Tr, mob_r, V_r0 + RECTUM_ADD_ML * 1e-6 * k / RS, obst_r)
    stages_r.append(Xn); report["rectum_fill"]["stages_ml"].append(round(Vn * 1e6, 1))
print("RECTUM", json.dumps(report["rectum_fill"]))
names_rf = [f"Rectum_Fill_{k}" for k in range(1, RS + 1)]
# neighbours: the sigmoid, the ureters and the reproductive organs (the uterus tips forward onto the bladder)
SOFT_R = [n_ for n_ in ("AN_Colon_Descending", "AN_Ureter_L", "AN_Ureter_R") if n_ in O] + REPRO_N
_rim_np = np.array(_rim)
def pinned_r(name, P):
    if name == "AN_Bladder":
        return dist_trigone(P) < TRIG_FIX
    if name == "AN_Colon_Descending":                # its seam with the rectum stays with the rectum's still rim
        return np.array([_kr.find(Vector(p))[2] < 0.008 for p in P])
    if name.startswith("AN_Ureter"):                 # held only right at its bladder opening and in the kidney
        s_ = name[-1]; near_b = np.array([loop_kd[s_].find(Vector(p))[2] < 0.006 for p in P])
        near_k = np.array([kid_trees["LR".index(s_)].find_nearest(Vector(p))[3] < 0.006 for p in P])
        return near_b | near_k
    return pinned(name, P)
rigid_r = U.Obstacles([(t_, 0.0008) for t_, _ in bone_trees] +
                      [(U.tree_in(O[n_], Mi), 0.0005) for n_ in ("AN_Kidney_L", "AN_Kidney_R", "AN_Adrenal_L", "AN_Adrenal_R")] +
                      [(vag_tree, 0.0008)], [(skin_tree, 0.003)], [])
bl_rest_tree = U.tree_in(bladder, Mi)
rigid_r1 = U.Obstacles(rigid_r.solid + [(bl_rest_tree, 0.001)], rigid_r.container, [])   # the bladder stays (alone)
st_trees_r = [U.tree_np(Xk, Tr) for Xk in stages_r[1:]]
st_pts_r = [(Xk, U.vnormals(Xk, Tr)) for Xk in stages_r[1:]]
# the uterus stays (the empty bladder in front of it leaves it no room to tip forward): the full rectum settles
# against it, the tubes and ovaries give way
th_r = [0.0] * RS
rlog = []
room_r, resid_r = U.make_room([O[n_] for n_ in SOFT_R], Mi, Mw, st_trees_r, pinned_r, rigid_r1, gap=0.0025, log=rlog,
                              rounds=10, decay=0.996, iters=400, carry_from=Xr, expanders=st_pts_r,
                              pre_stage=repro_pre(th_r))
report["rectum_fill"]["neighbours_residual_mm"] = resid_r
# the rectum settles where the uterus (held at its cervix, against the bladder) could not go further
fixed_r = mob_r < 1e-6; report["rectum_fill"]["settle"] = {}
for k in range(1, RS + 1):
    trees_k = []
    for n_ in REPRO_N:
        if n_ == "AN_Uterus_Lumen":
            continue
        trees_k.append(U.union_tree([O[n_]], Mi, positions={n_: room_r[n_][k - 1]}))
    Xs_, w0, w1 = U.settle(stages_r[k], Tr, fixed_r, trees_k)
    stages_r[k] = Xs_; report["rectum_fill"]["settle"][names_rf[k - 1]] = {"pressing_in_mm": round(w0 * 1000, 2), "after_mm": round(w1 * 1000, 2)}
    report["rectum_fill"].setdefault("settled_ml", []).append(round(U.volume(Xs_, Tr) * 1e6, 1))
for ob in REC:
    Mloc = (Mi @ ob.matrix_world).inverted(); prev = "Basis"
    for k in range(1, RS + 1):
        pos = np.array([tuple(Mloc @ Vector(p)) for p in stages_r[k][idx_r[ob.name]]])
        kb = U.add_key(ob, names_rf[k - 1], pos, relative=prev); prev = kb.name
        U.drive(kb, ctl, "Rectum_Fill", f"f*{RS}-{k - 1}" if k > 1 else f"f*{RS}")
for n_ in SOFT_R:
    ob = O[n_]; base = np.array([tuple(v.co) for v in ob.data.vertices]); prev = "Basis"
    if float(np.abs(room_r[n_][-1] - base).max()) < 1e-9:
        continue
    for k in range(1, RS + 1):
        kb = U.add_key(ob, names_rf[k - 1], room_r[n_][k - 1], relative=prev); prev = kb.name
        U.drive(kb, ctl, "Rectum_Fill", f"f*{RS}-{k - 1}" if k > 1 else f"f*{RS}")
print("RECTUM_ROOM", resid_r, rlog[-4:])

# ======================= 5e. both full (Bladder_Fill x Rectum_Fill) =======================
# the two sets of keys add, so both at 1 would collide: from that sum, the bladder gives way where the full rectum
# presses (it rises and holds less - as it does: a loaded rectum lowers the bladder's capacity), then the bowel and
# ureters are cleared of both. The difference is one more key per organ, driven by Bladder_Fill * Rectum_Fill.
setc(1.0, 0.0, 1.0)
_dg = bpy.context.evaluated_depsgraph_get()
def ev_local(ob):
    e_ = ob.evaluated_get(_dg); m_ = e_.to_mesh(); P_ = np.array([tuple(v.co) for v in m_.vertices]); e_.to_mesh_clear(); return P_
ADD = {n_: ev_local(O[n_]) for n_ in ["AN_Bladder", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_Colon_Descending", "AN_Ureter_L", "AN_Ureter_R"] + REPRO_N}
def shell_tree(names_, extra=None):
    bm_ = bmesh.new()
    for n_ in names_:
        M_ = Mi @ O[n_].matrix_world; vv = [bm_.verts.new(M_ @ Vector(p)) for p in ADD[n_]]
        for poly in O[n_].data.polygons:
            try:
                bm_.faces.new([vv[i] for i in poly.vertices])
            except ValueError:
                pass
    bmesh.ops.remove_doubles(bm_, verts=bm_.verts, dist=1e-6)
    bnd = [e for e in bm_.edges if e.is_boundary]
    if bnd:
        bmesh.ops.holes_fill(bm_, edges=bnd, sides=0)
    t_ = BVHTree.FromBMesh(bm_); bm_.free(); return t_
# the full rectum settles against the uterus as the full bladder has tilted it (Both_Full keys on the rectum too)
Xr_add = stages_r[-1].copy()
def _to_atlas(n_, P_):
    M_ = Mi @ O[n_].matrix_world; return np.array(P_) @ np.array(M_.to_3x3()).T + np.array(M_.translation)
def _to_local(n_, P_):
    M_ = (Mi @ O[n_].matrix_world).inverted(); return np.array(P_) @ np.array(M_.to_3x3()).T + np.array(M_.translation)
for n_ in ("AN_Rectum", "AN_Rectum_LowerAmpulla"):
    Xr_add[idx_r[n_]] = _to_atlas(n_, ADD[n_])      # (the rectum keeps its own object transform; caps stay)
_ut_add = [U.union_tree([O[n_]], Mi, positions={n_: ADD[n_]}) for n_ in ("AN_Uterus",)]
Xr_both, wr0, wr1 = U.settle(Xr_add, Tr, fixed_r, _ut_add)
for n_ in ("AN_Rectum", "AN_Rectum_LowerAmpulla"):
    ADD_R = _to_local(n_, Xr_both[idx_r[n_]]); ob = O[n_]; base = np.array([tuple(v.co) for v in ob.data.vertices])
    if float(np.abs(ADD_R - ADD[n_]).max()) > 1e-9:
        kb = U.add_key(ob, "Both_Full", base + (ADD_R - ADD[n_]), relative="Basis")
        U.drive2(kb, ctl, ["Bladder_Fill", "Rectum_Fill"], "a*b")
    ADD[n_] = ADD_R
report["rectum_fill"]["both_on_uterus_mm"] = [round(wr0 * 1000, 2), round(wr1 * 1000, 2)]
rect_full = shell_tree(["AN_Rectum", "AN_Rectum_LowerAmpulla"])
Xb_add = np.array([tuple(Vector(p)) for p in ADD["AN_Bladder"]])     # bladder local = atlas
Xb_both, w0, w1 = U.settle(Xb_add, bl_faces, fixed_b, [rect_full, shell_tree(["AN_Colon_Descending"]), vag_tree] + [t_ for t_, _ in bone_trees])
ADD_B = ADD["AN_Bladder"]; ADD["AN_Bladder"] = Xb_both
both_expander = shell_tree(["AN_Bladder", "AN_Rectum", "AN_Rectum_LowerAmpulla"])
ADD["AN_Bladder"] = ADD_B
# everything starts from the sum of the two single fillings, so Both_Full is the least correction from it: half-way
# states (one full, the other half) blend two near-identical arrangements instead of two separately solved ones
SOFT_B = [n_ for n_ in ("AN_Colon_Descending", "AN_Ureter_L", "AN_Ureter_R") if n_ in O] + REPRO_N
room_b, resid_b = U.make_room([O[n_] for n_ in SOFT_B], Mi, Mw, [both_expander], pinned_r, rigid_r, gap=0.0025,
                              rounds=10, decay=0.996, iters=400, start={n_: ADD[n_] for n_ in SOFT_B})
# the bladder settles again where the uterus could not go further
_w2 = U.settle(Xb_both, bl_faces, fixed_b, [U.union_tree([O[n_]], Mi, positions={n_: room_b[n_][0]}) for n_ in REPRO_N if n_ != "AN_Uterus_Lumen"])
Xb_both = _w2[0]
setc(0.0, 0.0, 0.0)
report["both_full"] = {"bladder_pressing_mm": round(w0 * 1000, 2), "after_mm": round(w1 * 1000, 2), "neighbours_residual_mm": resid_b,
                       "bladder_on_uterus_mm": [round(_w2[1] * 1000, 2), round(_w2[2] * 1000, 2)]}
for n_, final in [("AN_Bladder", Xb_both)] + [(n_, room_b[n_][0]) for n_ in SOFT_B]:
    ob = O[n_]; base = np.array([tuple(v.co) for v in ob.data.vertices])
    delta = final - ADD[n_]
    if float(np.abs(delta).max()) < 1e-9:
        continue
    kb = U.add_key(ob, "Both_Full", base + delta, relative="Basis")
    U.drive2(kb, ctl, ["Bladder_Fill", "Rectum_Fill"], "a*b")
_bmv = U.volume(np.array(Xb_both), [t for f in bl_faces for t in ([f[0], f[i], f[i + 1]] for i in range(1, len(f) - 1))])
report["both_full"]["bladder_ml"] = round(_bmv * 1e6, 1)
print("BOTH", json.dumps(report["both_full"]))

# ======================= 5f. polish: the last small contacts =======================
# where two organs' corrections meet, a few vertices can still sit a millimetre or two inside a neighbour in a full
# state. For each full state, every organ whose last key is active there steps out of whatever it is NEWLY inside (not
# inside at rest: seams, joins, the ureters in the kidneys), the correction fading over a few mm; written into that key.
def ev_np(ob):
    dg = bpy.context.evaluated_depsgraph_get(); ev = ob.evaluated_get(dg); me_ = ev.to_mesh()
    P = np.array([tuple(v.co) for v in me_.vertices]); ev.to_mesh_clear()
    return P @ np.array(ob.matrix_world.to_3x3()).T + np.array(ob.matrix_world.translation)
def ev_tree(ob):
    dg = bpy.context.evaluated_depsgraph_get(); bm_ = bmesh.new(); bm_.from_object(ob, dg); bm_.transform(ob.matrix_world)
    bmesh.ops.remove_doubles(bm_, verts=bm_.verts, dist=1e-7)
    bnd = [e for e in bm_.edges if e.is_boundary]
    if bnd:
        bmesh.ops.holes_fill(bm_, edges=bnd, sides=0)
    t_ = BVHTree.FromBMesh(bm_); bm_.free(); return t_
def boxes_meet(A, B, pad=0.005):
    return bool(np.all(A.min(0) - pad < B.max(0)) and np.all(B.min(0) - pad < A.max(0)))
CHK_OB = [bladder, urethra, uwall] + [O[n_] for n_ in SOFT_N if n_ != "AN_AnalCanal"] + [O["AN_AnalCanal"], O["AN_AnalSphincter"]] + \
         [O[n_] for n_ in ("AN_HipBone_L", "AN_HipBone_R", "AN_Sacrum")] + [O["AN_Kidney_L"], O["AN_Kidney_R"], vag_space]
# overlaps by design, left out of the polish and the checks: urethra in its wall, the lumen inside the uterus (and
# with it the cervix inside the vagina's vault), the tubes through the uterine wall into the cavity, each ligament's
# ends sunk in the organs it holds
NO_PAIR = [{"AN_Urethra", "AN_Urethra_Wall"}, {"AN_Uterus", "AN_Uterus_Lumen"}, {"AN_Uterus", "AN_Vagina_Space"},
           {"AN_Uterus_Lumen", "AN_Vagina_Space"}]
for s_ in "LR":
    NO_PAIR += [{f"AN_Uterine_Tube_{s_}", "AN_Uterus"}, {f"AN_Uterine_Tube_{s_}", "AN_Uterus_Lumen"}]
for n_, (a_org, _, b_org, _, _, _) in CORDS.items():
    NO_PAIR += [{n_, o_} for o_ in (a_org, b_org) if o_]
setc(0.0, 0.0, 0.0)
REST_W = {o_.name: ev_np(o_) for o_ in CHK_OB}; REST_T = {o_.name: ev_tree(o_) for o_ in CHK_OB}
_rest_in = {}
def rest_in(a, b, gap=0.0001):
    """vertices of a inside b at rest (or on it: seams)"""
    if (a, b) not in _rest_in:
        P_ = REST_W[a]; m_ = np.all((P_ > REST_W[b].min(0) - 0.006) & (P_ < REST_W[b].max(0) + 0.006), axis=1)
        out_ = np.zeros(len(P_), bool); ci = np.nonzero(m_)[0]
        if len(ci):
            out_[ci] = U.signed_dist(REST_T[b], P_[ci], band=0.006)[0] < gap
        _rest_in[(a, b)] = out_
    return _rest_in[(a, b)]
def fixed_of(ob):
    n_ = ob.name; P_ = _to_atlas(n_, [tuple(v.co) for v in ob.data.vertices])     # (the pin tests work in atlas space)
    if ob is bladder:
        return fixed_b
    if n_ in idx_r:
        return mob_r[idx_r[n_]] < 1e-6
    if n_ in ("AN_AnalCanal", "AN_AnalSphincter", "AN_Urethra", "AN_Urethra_Wall"):
        return np.ones(len(P_), bool)
    return pinned_r(n_, P_)
# every stage state, in order, then both full: key k (chained, relative to k-1) is fully on only at its own stage
# and the sum telescopes, so a correction written into it touches only the states between its neighbours
POLISH = [((k / NS, 0.0, 0.0), [names_fill[k - 1], "Bladder_" + names_fill[k - 1]]) for k in range(1, NS + 1)] + \
         [((0.0, 0.0, k / RS), [names_rf[k - 1]]) for k in range(1, RS + 1)] + [((1.0, 0.0, 1.0), ["Both_Full"])]
report["polish"] = {}
for st_, knames in POLISH:
    setc(*st_); tag = f"bladder{round(st_[0], 3)}_void{st_[1]}_rectum{round(st_[2], 3)}"; rep_ = {}
    W = {o_.name: ev_np(o_) for o_ in CHK_OB}; TR = {o_.name: ev_tree(o_) for o_ in CHK_OB}
    for a in CHK_OB:
        keys = a.data.shape_keys.key_blocks if a.data.shape_keys else {}
        kn = next((k_ for k_ in knames if k_ in keys), None)
        if kn is None:
            continue
        others = [b for b in CHK_OB if b is not a and {a.name, b.name} not in NO_PAIR and boxes_meet(W[a.name], W[b.name])]
        if not others:
            continue
        fx_ = fixed_of(a)
        if fx_.all():                                # (held: the uterus - the others step off it)
            continue
        Xn, w0, w1 = U.polish(W[a.name], [tuple(p.vertices) for p in a.data.polygons], fx_,
                              [TR[b.name] for b in others], [rest_in(a.name, b.name) for b in others],
                              boxes=[(W[b.name].min(0), W[b.name].max(0)) for b in others],
                              others=[(W[b.name], rest_in(b.name, a.name)) for b in others])
        if w0 <= 0.0:
            continue
        dl = (Xn - W[a.name]) @ np.array(a.matrix_world.inverted().to_3x3()).T
        kb = keys[kn]; K = np.array([tuple(d.co) for d in kb.data]) + dl
        kb.data.foreach_set("co", K.astype(np.float32).ravel()); a.data.update(); bpy.context.view_layer.update()
        W[a.name] = ev_np(a); TR[a.name] = ev_tree(a)
        rep_[a.name] = {"key": kn, "pressing_mm": round(w0 * 1000, 2), "after_mm": round(w1 * 1000, 2)}
        print("  polish", tag, a.name, rep_[a.name], flush=True)
    report["polish"][tag] = rep_; print("POLISH", tag, json.dumps(rep_))
setc(0.0, 0.0, 0.0)
done("keys")

# ======================= 6. lighting: insides lit everywhere, smooth-lighting bakes =======================
# new organs get the outer "AO" bake (at rest); every organ / bone material lights its inside and every object without
# an inner bake gets "AO_in" (the urethra's with it open, like the anal canal's); the bowel keeps build5_passage's
AO_DIST, AO_RAYS = 0.008, 32
setc(0.0, 0.0)
_tr = T.scene_bvh(); report["ao"] = {}
for ob in [bladder, urethra, uwall] + [REPRO[n_] for n_ in REPRO_N]:
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
    return [Vector(p) for p in ev_np(ob)]
chk = {}
# drivers
setc(0.5, 0.0)
chk["drivers_at_fill_0.5"] = {kb.name: round(kb.value, 3) for kb in bladder.data.shape_keys.key_blocks[1:]}
# seams in every state
_kidx = {s: None for s in "LR"}
states = [(0.0, 0.0, 0.0), (0.25, 0.0, 0.0), (0.5, 0.0, 0.0), (0.75, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
          (1.0, 1.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.0, 1.0), (0.5, 0.0, 0.5), (1.0, 0.0, 0.5), (0.5, 0.0, 1.0),
          (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)]
chk["seams_mm"] = {}; chk["volume_ml"] = {}; chk["clipping"] = {}
pairs_ob = CHK_OB
chk["clipping_at_rest"] = {}; _base = {}
for fill, void, rect in states:
    setc(fill, void, rect); tag = f"bladder{fill}_void{void}_rectum{rect}"
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
    # vertices of one organ inside another (count, deepest mm, body scale): at rest in full (seams, joins, the atlas'
    # own contacts); in every other state only those that are new or more than 0.3 mm deeper than at rest
    trees_ = {o_.name: ev_tree(o_) for o_ in pairs_ob}; W_ = {o_.name: ev_np(o_) for o_ in pairs_ob}
    clip = {}
    for a in pairs_ob:
        for b in pairs_ob:
            if a is b or {a.name, b.name} in NO_PAIR or not boxes_meet(W_[a.name], W_[b.name], 0.004):
                continue
            tb = trees_[b.name]; deep = {}; Pa = W_[a.name]
            lo, hi = W_[b.name].min(0) - 0.004, W_[b.name].max(0) + 0.004
            for i in np.nonzero(np.all((Pa > lo) & (Pa < hi), axis=1))[0]:
                v_ = Vector(Pa[i]); dd = tb.find_nearest(v_)[3]
                if dd is not None and 2e-5 < dd < 0.004 and U.inside(tb, v_):
                    deep[i] = dd
            key_ = f"{a.name} in {b.name}"
            if tag == "bladder0.0_void0.0_rectum0.0":
                _base[key_] = deep
                if deep:
                    chk["clipping_at_rest"][key_] = [len(deep), round(max(deep.values()) * 1000, 2)]
                continue
            b0 = _base.get(key_, {}); new_ = [d_ for i, d_ in deep.items() if i not in b0 or d_ > b0[i] + 0.0003]
            if new_:
                clip[key_] = [len(new_), round(max(new_) * 1000, 2)]
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
