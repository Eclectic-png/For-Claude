"""For consecutive tract organs (atlas), the zones where one organ's surface lies INSIDE the other (welded copies,
ray-parity inside test): connected face components with size, centre, depth. Args: -- <out.json>"""
import bpy, bmesh, sys, json
from mathutils import Vector
from mathutils.bvhtree import BVHTree
PAIRS = [("AN_Esophagus", "AN_Stomach"), ("AN_Stomach", "AN_Duodenum"), ("AN_Duodenum", "AN_Jejunum"),
         ("AN_Jejunum", "AN_Ileum"), ("AN_Ileum", "AN_Cecum"), ("AN_Cecum", "AN_Colon_Ascending"),
         ("AN_Cecum", "AN_Appendix"), ("AN_Ileum", "AN_Colon_Ascending"),
         ("AN_Colon_Ascending", "AN_Colon_Transverse"), ("AN_Colon_Transverse", "AN_Colon_Descending"),
         ("AN_Colon_Descending", "AN_Rectum")]
DIRS = [Vector(d).normalized() for d in ((0.31, 0.52, 0.79), (-0.6, 0.2, -0.77), (0.1, -0.95, 0.3))]


def welded(nm):
    o = bpy.data.objects[nm]; bm = bmesh.new(); bm.from_mesh(o.data); bm.transform(o.matrix_world)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7); bm.faces.ensure_lookup_table(); bm.normal_update()
    return bm, BVHTree.FromBMesh(bm)


def inside(tree, p):
    votes = 0
    for d in DIRS:
        c = 0; q = p.copy()
        while True:
            h = tree.ray_cast(q, d, 5.0)
            if h[0] is None: break
            c += 1; q = h[0] + d * 1e-7
        votes += c % 2
    return votes >= 2


cache = {}
def get(nm):
    if nm not in cache: cache[nm] = welded(nm)
    return cache[nm]


R = []
for a, b in PAIRS:
    for x, y in ((a, b), (b, a)):
        bx, tx = get(x); by, ty = get(y)
        ins = {f.index for f in bx.faces if inside(ty, f.calc_center_median())}
        comps = []; seen = set()
        for fi in ins:
            if fi in seen: continue
            st = [fi]; comp = []
            while st:
                k = st.pop()
                if k in seen: continue
                seen.add(k); comp.append(k)
                for e in bx.faces[k].edges:
                    st += [g.index for g in e.link_faces if g.index in ins and g.index not in seen]
            vs = {v for k in comp for v in bx.faces[k].verts}
            c = sum((v.co for v in vs), Vector()) / len(vs)
            depth = max(ty.find_nearest(v.co)[3] for v in vs)
            area = sum(bx.faces[k].calc_area() for k in comp)
            comps.append({"faces": len(comp), "area_mm2": round(area * 1e6, 1), "depth_mm": round(depth * 1000, 2),
                          "centre": [round(t, 4) for t in c]})
        comps.sort(key=lambda d: -d["area_mm2"])
        R.append({"organ": x, "inside": y, "components": comps})
        print("OV", x, "in", y, [(d["faces"], d["area_mm2"], d["depth_mm"]) for d in comps][:6])
json.dump(R, open(sys.argv[-1], "w"), indent=1)
