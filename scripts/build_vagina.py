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
# (its re-triangulated patch must stay clear of the meatus' own: the meatus patch reaches y 16.7 mm)
SLIT_Y = float(os.environ.get("VAG_SLIT_Y", "0.0230")); SLIT_HALF = float(os.environ.get("VAG_SLIT_HALF", "0.0042"))
N = 4 * sum(V.SEGS); LIP = float(os.environ.get("VAG_LIP", "0.00002")); G = 0.00001   # lips 0.02 mm apart at the opening (as the meatus)
_h = _skt.ray_cast(Vector((0, SLIT_Y, 0.70)), Vector((0, 0, 1)), 0.2)
M_body = Vector((0.0, SLIT_Y, _h[0].z if _h[0] is not None else 0.78))
# the patch reaches y 17.25 mm: it takes in the crest's odd back-to-back quads at y 17.4 - 19.4 (already in the
# skin; left at the patch's edge they shaded dark) and stays 0.55 mm clear of the meatus' patch
REG_U, REG_V = SLIT_HALF + 0.00155, 0.0016
def _skin_loops():
    bm_ = U.bm_in(hips, None); L = U.boundary_loops(bm_)
    out = sorted(len(lp) for lp in L if (sum((v.co for v in lp), Vector()) / len(lp) - Vector((0, 0.0135, 0.7835))).length < 0.004)
    bm_.free(); return out
_meatus_before = _skin_loops()
# the slit is cut with exactly the lumen's first ring (ring model: long front / back arms, stub side arms), so the
# tube needs no blend out of a differently spaced slit (that folded it inside the closed slit). Chart: u along +y
# (body, backwards), v across; body = atlas scaled by Internal_Fit_Xform
SX_, SY_ = Mw.to_3x3().col[0].length, Mw.to_3x3().col[1].length
G0 = LIP / 2 / SX_                                 # atlas
STUB0 = 3 * G0                                     # side stubs at the opening: next to nothing (in the skin they
                                                   # point down the labial walls, one ring up sideways - a twist)
RING0 = V.lumen_ring([(STUB0, 0.0), (SLIT_HALF / SY_, 0.0), (STUB0, 0.0), (SLIT_HALF / SY_, 0.0)], G0)
slit_uv = [(-y * SY_, x * SX_) for x, y in RING0]   # front (+y ring) = -y body
slit_idx, chart, skin_new, rep_s = U.cut_meatus(hips, M_body, Vector((0, 1, 0)), Vector((1, 0, 0)), N, SLIT_HALF, LIP,
                                                REG_U, REG_V, slit_uv=slit_uv, ref_radius=REG_U + 0.0009)
hips.data.update()
_meatus_after = _skin_loops()
if _meatus_after != _meatus_before:                 # the urethra's opening must not be touched
    raise RuntimeError(f"introitus cut changed the meatus: {_meatus_before} -> {_meatus_after}")
slit_body = [hips.data.vertices[i].co.copy() for i in slit_idx]
slit = [Mi @ p for p in slit_body]                 # atlas
Cs = sum(slit, Vector()) / len(slit)
report["introitus"] = {"centre_body": [round(x, 5) for x in M_body], "half_mm": SLIT_HALF * 1000,
                       "length_real_mm": round((slit[N // 4 + N // 8] - slit[3 * N // 4 + N // 8]).length * 1000, 2), **rep_s}
print("INTROITUS", json.dumps(report["introitus"]))

# the slit's own plane (the crest it is cut into slopes along its length). A closed slit is a line: its plane holds
# the tip-to-tip line and the side-to-side direction; facing into the body
_tip = np.array(tuple(slit[N // 4 + N // 8] - slit[3 * N // 4 + N // 8]))
n_slit = np.cross(np.array(tuple(X_)), _tip); n_slit /= np.linalg.norm(n_slit)
n_slit *= 1 if n_slit @ np.array(tuple((Mi.to_3x3() @ Vector((0, 0, 1))).normalized())) > 0 else -1

# ---- 2. the centre line: up from the slit into the placeholder's axis, then along it ----
def catmull(pts, n):
    P = [pts[0]] + pts + [pts[-1]]; out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=(i == len(P) - 3)):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    return out
UP = (Mi.to_3x3() @ Vector((0, 0, 1))).normalized()          # into the body at the vestibule
# out of the slit square to its plane for the first mm (rings stay parallel to the slit there), then curving into
# the placeholder's axis while the slit's front / back arms shrink
# one cubic from the slit (leaving square to its plane) to P1v (arriving along the placeholder's axis), then the
# axis on to P2v: no kinks (a Catmull-Rom through the placeholder's start point turned up to 30 deg per ring)
_d1 = (P1v - Cs).length; _ta = (P2v - P1v).normalized(); _n = Vector(tuple(n_slit))
_B = [Cs, Cs + _n * (_d1 * 0.45), P1v - _ta * (_d1 * 0.35), P1v]
_bez = [(1 - t) ** 3 * _B[0] + 3 * (1 - t) ** 2 * t * _B[1] + 3 * (1 - t) * t * t * _B[2] + t ** 3 * _B[3] for t in np.linspace(0, 1, 200)]
dense = np.array([tuple(p) for p in _bez] + [tuple(P1v.lerp(P2v, t)) for t in np.linspace(0, 1, 120)[1:]])
seg = np.linalg.norm(np.diff(dense, axis=0), axis=1); cum = np.concatenate([[0.0], np.cumsum(seg)])
S = []; s_ = 0.0
while s_ < cum[-1]:                                # fine where the rugae run, coarser in the vault
    S.append(s_); s_ += 0.0010 if s_ < 0.004 else float(os.environ.get("VAG_STEP_MM", "0.7")) / 1000 if s_ < 0.050 else 0.0012
                                                   # (1 mm where the rings ease out of the curved slit: closer, they cross)
S = np.array(S)
C = np.stack([np.interp(S, cum, dense[:, k]) for k in range(3)], 1)
# ring planes: parallel to the slit's plane at first (stacked, they cannot cross while the slit's long front / back
# arms shrink), turning square to the centre line between 4 and 12 mm in
TT = V.tangents(C)
report["introitus"]["slit_plane_tilt_deg"] = round(math.degrees(math.acos(min(1.0, float(n_slit @ np.array(tuple(UP)))))), 1)
NJ = []
for s_, t_ in zip(S, TT):
    w_ = float(U.ss(0.004, 0.012, s_)); n_ = n_slit * (1 - w_) + t_ * w_; NJ.append(n_ / np.linalg.norm(n_))
NJ = np.array(NJ)
FR = V.frames_along(C, X_, NJ)

# sections of the placeholder (walls ~3 mm inside it) and of the cervix, ring by ring
def tri_arrays(ob):
    """(vertices in atlas space, triangles) of an object parented under Internal_Fit_Xform"""
    me = ob.data; ML = ob.matrix_local; X = np.array([tuple(ML @ v.co) for v in me.vertices]); F = []
    for p in me.polygons:
        vv = list(p.vertices)
        F += [(vv[0], vv[k], vv[k + 1]) for k in range(1, len(vv) - 1)]
    return X, np.array(F)
Xv, Fv = tri_arrays(vs); Xu, Fu = tri_arrays(uterus)
env, cer, shift, envP = [], [], [], []
for c, (ex, ey), j in zip(C, FR, range(len(C))):
    t = NJ[j]
    sv = V.plane_section(Xv, Fv, c, t, ex, ey)
    # only the part of the placeholder's section round the centre (it is one loop where the plane cuts the tube)
    pn = V.near(sv, 0.025) if len(sv) else None
    # the placeholder counts only once the ring's centre is inside its section (below that the plane meets it
    # obliquely); its section's centre (x, y) is where the lumen ring is centred there
    if pn is not None:
        lo_, hi_ = pn.min(0), pn.max(0); cx_, cy_ = (lo_ + hi_) / 2
        ok_ = lo_[0] + 0.001 < 0 < hi_[0] - 0.001 and lo_[1] + 0.001 < 0 < hi_[1] - 0.001
    env.append(V.extents(pn, (cx_, cy_)) if pn is not None and ok_ else None)
    envP.append((pn, (cx_, cy_)) if pn is not None and ok_ else None)
    shift.append((cx_, cy_) if pn is not None and ok_ else None)
    su = V.plane_section(Xu, Fu, c, t, ex, ey)
    cer.append(V.near(su, 0.022) if len(su) else None)
first_env = next(j for j, e in enumerate(env) if e)
last_env = max(j for j, e in enumerate(env) if e)
# where the cervix pierces the placeholder its section has a hole: a direction with no wall there takes the reach
# interpolated between the nearest rings that have one (a smooth front wall for the cervix to pass through)
for d_ in range(4):
    have = [j for j in range(first_env, last_env + 1) if env[j] and env[j][d_] > 0.002]
    for j in range(first_env, last_env + 1):
        if env[j] and env[j][d_] <= 0.002 and have:
            lo = max([h for h in have if h < j], default=None); hi = min([h for h in have if h > j], default=None)
            if lo is not None and hi is not None:
                env[j] = list(env[j]); env[j][d_] = env[lo][d_] + (env[hi][d_] - env[lo][d_]) * (j - lo) / (hi - lo)
            elif lo is not None or hi is not None:
                env[j] = list(env[j]); env[j][d_] = env[lo if lo is not None else hi][d_]
for j in range(first_env):                         # below the placeholder: its first section
    env[j] = env[first_env]
_gaps = [j for j in range(first_env, last_env + 1) if env[j] is None]
for j in _gaps:                                    # a ring whose centre slipped out of the section: the nearest one's
    near_ = min((k for k in range(first_env, last_env + 1) if env[k] is not None), key=lambda k: abs(k - j))
    env[j] = env[near_]
s_c0 = next(S[j] for j, u in enumerate(cer) if u is not None)          # the first ring that meets the cervix
# ring centres follow the placeholder's: blended in over 6 mm once inside it; the front / back offset is held from
# 15 mm below the cervix on (the pierce hole skews the section's centre there)
SH = np.zeros((len(S), 2)); _last = np.zeros(2)
for j in range(len(S)):
    if j >= first_env and shift[j] is not None:
        _cur = np.array(shift[j], float)
        if S[j] > s_c0 - 0.015:
            _cur[1] = _last[1]
        _last = _cur
    SH[j] = _last * float(U.ss(S[first_env], S[first_env] + 0.006, S[j])) if j >= first_env else 0.0
SH[:, 0] = 0.0                                     # symmetric: no side shift (the partial sections near the
for _ in range(4):                                 # placeholder's bottom made it jump); front / back smoothed
    SH[1:-1, 1] = 0.25 * SH[:-2, 1] + 0.5 * SH[1:-1, 1] + 0.25 * SH[2:, 1]
WALL_LAT, WALL_AP = 0.0035, 0.003
RUG_A = float(os.environ.get("VAG_RUGAE_MM", "0.6")) / 1000; RUG_L = 0.0035; RUG_SIDE = 0.4   # fold height, spacing
SLIT_HALF_A = SLIT_HALF / SY_                      # atlas

# ring 0 = the skin's slit, vertex for vertex (cut from RING0). The chart's v may run either way across: then the
# lumen's rings are mirrored in x (MIRX) so that ring point k always sits over slit vertex k
ex0, ey0 = FR[0]
s2d = np.array([((p - Cs).dot(Vector(ex0)), (p - Cs).dot(Vector(ey0))) for p in slit])
_pv = chart.to3d(0.0, 0.0008) - chart.to3d(0.0, 0.0)          # which side of the midline the chart's +v wall is on
MIRX = 1.0 if _pv.x > 0 else -1.0
assert s2d[N // 4 + N // 8][1] > 0, "slit vertex at the front arm's tip should lie in front"
ring0 = s2d; idx0 = list(range(N))

# the T-bars sit at the arms' ends: the placeholder's thickness there (minus a 2.5 mm wall), smoothed along the rings
VT_S = []
for j in range(len(S)):
    e = env[j]
    if e is None or j > last_env:
        VT_S.append(2 * G); continue
    Lr_, Ll_ = max(e[0] - WALL_LAT, 3 * G), max(e[1] - WALL_LAT, 3 * G)
    vt_ = max(min(e[2], e[3]) - WALL_AP, 2 * G)
    if envP[j] is not None:
        _rr = [V.reach_y_at(envP[j][0], sgn_ * (L_ - 0.0005), envP[j][1]) for sgn_, L_ in ((1, Lr_), (-1, Ll_))]
        _rr = [r_ for r_ in _rr if r_ is not None]
        if _rr:
            vt_ = max(min(min(r_) for r_ in _rr) - 0.0025, 2 * G)
    VT_S.append(vt_)
VT_S = np.array(VT_S)
for _ in range(6):
    VT_S[1:-1] = np.minimum(VT_S[1:-1], 0.25 * VT_S[:-2] + 0.5 * VT_S[1:-1] + 0.25 * VT_S[2:])

rings = []; params = []; RUG_Y = []; RUG_ON = []
for j, s in enumerate(S):
    e = env[j]
    if j > last_env or e is None:
        break
    Lr, Ll = max(e[0] - WALL_LAT, 3 * G), max(e[1] - WALL_LAT, 3 * G)
    Vt = max(min(e[2], e[3]) - WALL_AP, 2 * G)
    Vt = VT_S[j]
    f_h = U.ss(0.001, 0.006, s)                # the placeholder's bottom is flat (3 mm half-thick): the slit's
                                               # front / back arms are short by the time the lumen enters it
    fb = SLIT_HALF_A * (1 - f_h) + 0.0003 * f_h
    f_b = U.ss(0.001, 0.005, s)                # the back arm (inside of the curve into the axis) goes sooner
    fbb = SLIT_HALF_A * (1 - f_b) + 0.0003 * f_b
    lr = (STUB0 * (1 - f_h) + Lr * f_h, STUB0 * (1 - f_h) + Ll * f_h)
    # H (the lateral recesses' T-bars) over the lower third, narrowing to the transverse slit by 34 mm
    Vh = Vt * U.ss(0.010, 0.017, s) * (1 - U.ss(0.026, 0.034, s))   # after the side arms have stopped lengthening
    g_ = G0 * (1 - float(U.ss(0.004, 0.012, s))) + G * float(U.ss(0.004, 0.012, s))   # walls: as the lips, then 0.02 mm
    H = V.lumen_ring([(lr[0], Vh), (fb, 0.0), (lr[1], Vh), (fbb, 0.0)], g_)
    wv = U.ss(s_c0 - 0.006, s_c0 - 0.001, s)    # only the vault round the cervix is open
    if wv > 0:
        u = cer[j] - SH[j] if cer[j] is not None else None
        rx = min(max((np.abs(u[:, 0]).max() + 0.0015) if u is not None else 0.0, min(Lr, Ll)), min(e[0], e[1]) - 0.0015)
        ryb = min(max((-u[:, 1].min() + 0.0015) if u is not None else 0.0, e[3] - WALL_AP), e[3] - 0.0012)
        ryf = e[2] - (0.0015 if u is not None else WALL_AP)   # thinner where the cervix passes: it protrudes
        P2 = H * (1 - wv) + V.ellipse_ring(max(rx, 0.001), max(ryf, 0.001), max(ryb, 0.001)) * wv
        params.append((round(s * 1000, 1), round(rx * 1000, 2), round(ryf * 1000, 2), round(ryb * 1000, 2)))
    else:
        P2 = H
    w_ex = 1 - U.ss(0.0, 0.002, s)             # the cut slit (on the curved skin) eases into the planar rings
    P2 = P2 * np.array([MIRX, 1.0])
    # rugae: transverse folds; the touching walls fold together (a corrugated slit), most in the midline columns
    # (fading to RUG_SIDE at the side recesses), over the lower part, fading in after the opening and before the vault
    a_r = RUG_A * float(U.ss(0.004, 0.009, s)) * (1 - float(U.ss(s_c0 - 0.018, s_c0 - 0.010, s))) * (1 - 0.35 * float(U.ss(0.02, 0.04, s)))
    RUG_Y.append(a_r * math.sin(2 * math.pi * s / RUG_L) * (1 - w_ex)); RUG_ON.append(a_r > 0)   # (the ring's centre
                                                                                                    # moves with the fold)
    if a_r > 0:
        xr = np.abs(P2[:, 0]) / max(np.abs(P2[:, 0]).max(), 1e-6)
        P2 = P2 + np.stack([np.zeros(len(P2)), a_r * (1 - (1 - RUG_SIDE) * xr ** 2) * math.sin(2 * math.pi * s / RUG_L)], 1)
    rings.append(ring0 * w_ex + (P2 + SH[j]) * (1 - w_ex))
    if os.environ.get("VAG_DEBUG") and j < 10:
        print("DBG ring", j, "s %.1f" % (s * 1000), "env", [round(x * 1000, 2) for x in e], "Lr %.2f Ll %.2f fb %.2f fbb %.2f Vh %.3f g %.3f" % (lr[0] * 1000, lr[1] * 1000, fb * 1000, fbb * 1000, Vh * 1000, g_ * 1000),
              "turn_deg %.1f" % math.degrees(math.acos(min(1, float(NJ[j] @ NJ[max(j - 1, 0)])))), "SH", [round(x * 1000, 2) for x in SH[j]])
NR = len(rings)
# guard: an arm with a point outside the placeholder (or within CLEAR of it) is shortened toward the ring's centre,
# a little at a time, until it clears (the section's reach is read at the axis; walls pinch in places)
CLEAR = 0.0004; PER = sum(V.SEGS)
_vt = U.tree_in(vs, Mi); _vo = U.tree_in(vs, Mi, cap=False); n_guard = 0
for j in range(first_env + 2, NR):
    ex, ey = (Vector(tuple(a_)) for a_ in FR[j]); c = Vector(tuple(C[j]))
    for blk in range(4):
        sl = slice(blk * PER, (blk + 1) * PER)
        for _ in range(40):
            pts = [c + ex * x + ey * y for x, y in rings[j][sl]]
            def _out(p_):                       # by the nearest face's side (the placeholder's normals face out;
                q_, n_, _, d_ = _vo.find_nearest(p_)   # ray parity is unreliable next to its capped cervix hole)
                return d_ < CLEAR or (p_ - q_).dot(n_) > 0
            bad = any(_out(p_) for p_ in pts)
            if not bad:
                break
            ctr = (SH[j] if j < len(SH) else np.zeros(2)) + np.array([0.0, RUG_Y[j]])   # (about the fold-shifted centre:
                                                                                         # the arms' roots stay joined)
            rings[j][sl] = ctr + (rings[j][sl] - ctr) * 0.96; n_guard += 1
report["lumen_guard_shrinks"] = n_guard
report["lumen"] = {"rings": NR, "points_per_ring": N, "length_mm": round(float(S[NR - 1]) * 1000, 1),
                   "envelope_from_mm": round(float(S[first_env]) * 1000, 1), "envelope_gaps_filled": len(_gaps), "cervix_from_mm": round(float(s_c0) * 1000, 1),
                   "vault_rings_s_rx_ryf_ryb_mm": params[::2]}
print("LUMEN", json.dumps(report["lumen"]))

# ---- 3. the lumen mesh: rings joined index by index, a dome on top ----
P3 = V.sweep(C[:NR], FR[:NR], rings)
P3[:N] = np.array([tuple(slit[idx0[jj]]) for jj in range(N)])   # ring 0 exactly the skin's slit
_td = (C[NR - 1] - C[NR - 2]) / np.linalg.norm(C[NR - 1] - C[NR - 2])
tip = C[NR - 1] + np.append(SH[NR - 1], 0.0)[:2] @ np.array([FR[NR - 1][0], FR[NR - 1][1]])
for _d in (0.0015, 0.001, 0.0005, 0.0):            # the dome's top: as far up as the placeholder allows (CLEAR)
    _q, _n, _, _dd = _vo.find_nearest(Vector(tuple(tip + _td * _d)))
    if _dd >= CLEAR and (Vector(tuple(tip + _td * _d)) - _q).dot(_n) < 0:
        break
tip = tip + _td * _d
report["lumen_dome_mm"] = round(_d * 1000, 2)
P3 = np.vstack([P3, tip[None]])
# where the rugae run, triangles: each quad split from its outer corner on ring j to the inner one on ring j + 1
# (outer: farther from the ring's centre) - the collapsed walls lie 0.04 mm apart and the folds twist their quads by
# more than that; split by index, facing walls would take opposite diagonals and cross. Quads elsewhere
faces = []
for j in range(NR - 1):
    a, b = j * N, (j + 1) * N
    if not (RUG_ON[j] or RUG_ON[j + 1]):
        faces += [(a + k, a + (k + 1) % N, b + (k + 1) % N, b + k) for k in range(N)]
        continue
    rad = np.linalg.norm(np.asarray(rings[j]) - ((SH[j] if j < len(SH) else 0.0) + np.array([0.0, RUG_Y[j]])), axis=1)
    for k in range(N):
        k1 = (k + 1) % N
        if rad[k] >= rad[k1]:
            faces += [(a + k, a + k1, b + k1), (a + k, b + k1, b + k)]
        else:
            faces += [(a + k, a + k1, b + k), (a + k1, b + k1, b + k)]
pole = NR * N; a = (NR - 1) * N
faces += [(a + k, a + (k + 1) % N, pole) for k in range(N)]
me = bpy.data.meshes.new("AN_Vagina"); me.from_pydata([tuple(p) for p in P3], [], faces); me.update()
old = O.get("AN_Vagina")
if old is not None:
    O.remove(old)
lumen = O.new("AN_Vagina", me); col_fit.objects.link(lumen); lumen.parent = E; lumen.matrix_parent_inverse = Matrix.Identity(4)
# the ring planes (atlas centre, normal, arc length mm) for tools: the cervix join renumbers the vertices
lumen["ring_planes"] = [list(map(float, list(C[j]) + list(NJ[j]) + [S[j] * 1000 + (SH[j][1] * 0)])) for j in range(NR)]
lumen["ring_shift_front"] = [float(SH[j][1]) for j in range(NR)]
for p in me.polygons:
    p.use_smooth = True
bm = bmesh.new(); bm.from_mesh(me); bm.normal_update(); bm.faces.ensure_lookup_table()
# normals out of the lumen (into the wall), like the other tubes: away from the centre line
_c = Vector(tuple(C[NR // 2]))                    # (the test face: the one farthest out round that ring - unambiguous)
f_ = max((f for f in bm.faces if all(v.index // N in (NR // 2, NR // 2 + 1) for v in f.verts)), key=lambda f: (f.calc_center_median() - _c).length)
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
_ckd = __import__("mathutils").kdtree.KDTree(NR)
for j in range(NR):
    _ckd.insert(Vector(tuple(C[j])), j)
_ckd.balance()
for vi, v in enumerate(me.vertices):          # by the nearest centre point (vertex numbers shift at the cervix join)
    at.data[vi].value = 0.3 * float(U.ss(0.0, 0.008, S[_ckd.find(lumen.matrix_local @ v.co)[1]]))
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
_beyond = [p for p in LW if _ckd.find(Mi @ p)[1] >= first_env + 2]
def _side(p):
    q_, n_, _, d_ = _wt_open.find_nearest(p); return (p - q_).dot(n_) <= 0
inside_wall = [_side(p) for p in _beyond]
report["checks"] = {
    "skin_seam_mm": round(max((a - b).length for a, b in zip(sorted(skin_body, key=lambda p: tuple(p)), sorted(lum_w, key=lambda p: tuple(p)))) * 1000, 5),
    "lumen_open_loops": loops_of(lumen), "wall_open_loops": loops_of(wall),
    "lumen_outside_wall_points": len(inside_wall) - sum(inside_wall),
    "outside_at_mm": sorted({round(float(S[_ckd.find(Mi @ p)[1]]) * 1000, 1) for p, ok in zip(_beyond, inside_wall) if not ok}),
    "outside_dist_mm": round(max([_wt_open.find_nearest(p)[3] for p, ok in zip(_beyond, inside_wall) if not ok] or [0]) * 1000, 3),
    "wall_bottom_at_mm": round(float(S[_ckd.find(VC0)[1]]) * 1000, 1),
    "lumen_wall_clearance_min_mm": round(min(_wt_open.find_nearest(p)[3] for p in _beyond) * 1000, 2),
    "urethra_wall_clearance_min_mm": round(min(_ut.find_nearest(p)[3] for p in LW[::3]) * 1000, 2) if _ut else None}
print("CHECKS", json.dumps(report["checks"]))
U.save(SAVE_AS, OUT, report)
