"""Connected digestive passage from the atlas (anatomy_ref.blend -> a new file; the atlas itself is not changed).

The organs are closed shells meeting by overlap; tract_lib.open_junction() opens each meeting so the lumens join,
changing only a small patch round the junction. Before that:
  - the ileum carried a duplicate of the whole cecum (a separate closed copy, all 288 faces): removed;
  - the appendix stood 4.6 mm clear of the cecum: moved onto it (sunk APPX_PENETRATE into the wall, at least
    APPX_CLEAR from the ileum's entrance) but not opened - it is a dead end, not part of the passage;
  - junctions that only touch get a smooth local bulge of both walls (BULGE) so they overlap by a lumen-sized lens;
  - jejunum / ileum share coincident walls (split from one source tube): open_junction deflates the jejunum's
    patch by 1 mm (fading to nothing at the patch border) so the two walls cross cleanly.
Afterwards every material (organs, bones, cartilage, lungs) lights its inside (tract_lib.light_insides) and every
object gets an "AO_in" bake.
Run: blender -b anatomy_ref.blend --python scripts/build_tract_passage.py -- <out.blend> <report folder>"""
import bpy, bmesh, sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tract_lib as T
from mathutils import Vector, kdtree
from mathutils.bvhtree import BVHTree

SAVE_AS, OUT = sys.argv[-2], sys.argv[-1]
O = bpy.data.objects; report = {"junctions": []}
JUNCTIONS = [("AN_Esophagus", "AN_Stomach"), ("AN_Stomach", "AN_Duodenum"), ("AN_Duodenum", "AN_Jejunum"),
             ("AN_Jejunum", "AN_Ileum"), ("AN_Ileum", "AN_Cecum"),
             ("AN_Cecum", "AN_Colon_Ascending"), ("AN_Colon_Ascending", "AN_Colon_Transverse"),
             ("AN_Colon_Transverse", "AN_Colon_Descending"), ("AN_Colon_Descending", "AN_Rectum")]
# "A|B": (h_mm, R_mm). The duodenum and jejunum only touch (one face, 0.5 mm deep): a 2.5 mm bulge of both walls
# over 15 mm gives an 8 mm-radius opening
BULGE = {k: tuple(v) for k, v in json.loads(os.environ.get("TRACT_BULGE", '{"AN_Duodenum|AN_Jejunum": [2.5, 15]}')).items()}
APPX_PENETRATE = 0.001                           # m the appendix base sinks into the cecum wall

# 1. the duplicate cecum inside the ileum
I, C = O["AN_Ileum"], O["AN_Cecum"]
kd = kdtree.KDTree(len(C.data.vertices))
for i, v in enumerate(C.data.vertices):
    kd.insert(C.matrix_world @ v.co, i)
kd.balance()
bm = bmesh.new(); bm.from_mesh(I.data); Mi = I.matrix_world
dup = [f for f in bm.faces if all(kd.find(Mi @ v.co)[2] < 1e-7 for v in f.verts)]
bmesh.ops.delete(bm, geom=dup, context='FACES')
bmesh.ops.delete(bm, geom=[v for v in bm.verts if not v.link_faces], context='VERTS')
bm.to_mesh(I.data); bm.free(); I.data.update()
report["ileum_duplicate_cecum_faces_removed"] = len(dup)


# 2. the appendix onto the cecum
def nearest_pair(A, B):
    tB = T.closed_tree(B); best = None
    for v in A.data.vertices:
        p = A.matrix_world @ v.co; q, _, _, dd = tB.find_nearest(p)
        if best is None or dd < best[0]:
            best = (dd, p, q)
    return best


# the appendix is a dead end (nothing passes through it), so it is attached - sunk APPX_PENETRATE into the cecum
# wall - but not opened. It goes to the cecum point nearest its base that is at least APPX_CLEAR from the
# ileum's entrance (the appendix arises a little below the ileocecal valve); the nearest point overall sat right
# on the ileal opening.
APPX_CLEAR = 0.015
Ap = O["AN_Appendix"]; tC = T.closed_tree(C)
base = min((Ap.matrix_world @ v.co for v in Ap.data.vertices), key=lambda p_: tC.find_nearest(p_)[3])
gap0 = tC.find_nearest(base)[3]
_bI = T.world_bm(I); _cs = T.overlap_components(_bI, tC); icj = T.comp_sphere(_bI, _cs[0])[0]; _bI.free()
_bC = T.world_bm(C, weld=True)
cands = [v for v in _bC.verts if (v.co - icj).length >= APPX_CLEAR]
tgt = min(cands, key=lambda v: (v.co - base).length)
mv = (tgt.co - tgt.normal * APPX_PENETRATE) - base; tgt_co = tgt.co.copy(); _bC.free()
Ap.location += mv if Ap.parent is None else Ap.parent.matrix_world.inverted().to_3x3() @ mv
bpy.context.view_layer.update()
report["appendix"] = {"gap_before_mm": round(gap0 * 1000, 2), "moved_mm": round(mv.length * 1000, 2),
                      "base_to_ileal_entrance_mm": round((tgt_co - icj).length * 1000, 1), "opened": False}


# 3. local bulges where neighbours only touch
def bulge(A, B, h, R):
    dd, p, q = nearest_pair(A, B); P = (p + q) / 2; moved = {}
    for ob in (A, B):
        bw = T.world_bm(ob, weld=True)               # one normal per position (seam-split copies move together)
        kw = kdtree.KDTree(len(bw.verts)); nrm = []
        for i, v in enumerate(bw.verts):
            kw.insert(v.co, i); nrm.append(v.normal.copy())
        kw.balance(); bw.free()
        bm = T.world_bm(ob); n_ = 0
        for v in bm.verts:
            d = (v.co - P).length
            if d < R:
                w = 0.5 * (1 + math.cos(math.pi * d / R))
                v.co += nrm[kw.find(v.co)[1]] * (h * w); n_ += 1
        bm.transform(ob.matrix_world.inverted()); bm.to_mesh(ob.data); bm.free(); ob.data.update(); moved[ob.name] = n_
    return {"contact": [round(x, 4) for x in P], "gap_before_mm": round(dd * 1000, 2), "h_mm": h * 1000,
            "R_mm": R * 1000, "verts_moved": moved}


report["bulges"] = []
for (a, b) in JUNCTIONS:
    key = f"{a}|{b}"
    if key in BULGE:
        h, R = BULGE[key]
        report["bulges"].append({"pair": [a, b], **bulge(O[a], O[b], h / 1000, R / 1000)})

def n_loops(ob):
    bm = T.world_bm(ob, weld=True); E = {e for e in bm.edges if e.is_boundary}; seen = set(); n = 0
    for e in E:
        if e in seen:
            continue
        n += 1; st = [e]
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x)
            for v in x.verts:
                st += [y for y in v.link_edges if y in E and y not in seen]
    bm.free(); return n


TREES = {nm: T.closed_tree(O[nm]) for nm in {x for p_ in JUNCTIONS for x in p_}}   # closed shells, before opening
print("LOOPS start", {nm: n_loops(O[nm]) for nm in {x for p_ in JUNCTIONS for x in p_}})
# 4. open every junction
for a, b in JUNCTIONS:
    try:
        r = T.open_junction(O[a], O[b], trees=TREES)
        r["opening"] = T.opening_loops(O[a], O[b])
    except Exception as ex:
        r = {"pair": [a, b], "error": str(ex)}
    report["junctions"].append(r); print("JUNCTION", json.dumps(r))
    print("LOOPS after", a, b, n_loops(O[a]), n_loops(O[b]))

# 5. light the insides of every organ, bone and cartilage (the shaders painted backfaces one flat colour) and bake
# "AO_in", the inner smooth-lighting, against the whole atlas; the outer "AO" of the rebuilt patches came from the
# nearest old vertex in open_junction
_meshes = [o_ for o_ in O if o_.type == 'MESH' and o_.name in bpy.context.view_layer.objects]
report["ao_in_mean"] = T.light_all(_meshes)

os.makedirs(OUT, exist_ok=True)
json.dump(report, open(os.path.join(OUT, "report.json"), "w"), indent=1)
print("REPORT", json.dumps({k: v for k, v in report.items() if k != "junctions"}))
bpy.ops.wm.save_as_mainfile(filepath=SAVE_AS, relative_remap=False)
