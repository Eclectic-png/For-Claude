"""Digestive tract continuity (atlas): for each organ, its open-end loops (boundary loops); for each consecutive pair,
the closest-matching loops, their centroid distance, mean point-to-loop distance, radii, and whether the meshes
overlap. Args: -- <out.json>"""
import bpy, bmesh, sys, json
from mathutils import Vector
from mathutils.bvhtree import BVHTree
CHAIN = ["AN_Esophagus", "AN_Stomach", "AN_Duodenum", "AN_Jejunum", "AN_Ileum", "AN_Cecum", "AN_Colon_Ascending",
         "AN_Colon_Transverse", "AN_Colon_Descending", "AN_Rectum", "AN_AnalCanal"]
EXTRA = [("AN_Cecum", "AN_Appendix"), ("AN_Ileum", "AN_Colon_Ascending")]
import os
if os.environ.get("TRACT_CHAIN"):                # e.g. the pelvis subset in Hips.blend
    CHAIN = os.environ["TRACT_CHAIN"].split(","); EXTRA = []


def loops_of(ob):
    bm = bmesh.new(); bm.from_mesh(ob.data); bm.transform(ob.matrix_world)
    seen = set(); out = []
    for e in bm.edges:
        if not e.is_boundary or e in seen:
            continue
        st = [e]; comp = []
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x); comp.append(x)
            for v in x.verts:
                st += [y for y in v.link_edges if y.is_boundary and y not in seen]
        pts = list({v.co.copy().freeze() for x in comp for v in x.verts})
        c = sum(pts, Vector()) / len(pts)
        out.append({"n": len(comp), "c": c, "r": sum((p - c).length for p in pts) / len(pts), "pts": pts})
    tree = BVHTree.FromBMesh(bm); bm.free()
    return sorted(out, key=lambda l: -l["n"]), tree


info = {}; R = {"organs": {}, "pairs": []}
for nm in set(CHAIN) | {b for _, b in EXTRA}:
    ob = bpy.data.objects.get(nm)
    if ob is None:
        R["organs"][nm] = "MISSING"; continue
    L, T = loops_of(ob); info[nm] = (L, T)
    R["organs"][nm] = {"loops": len(L), "big_loops": [(l["n"], round(l["r"] * 1000, 1)) for l in L if l["n"] >= 12][:6]}
pairs = list(zip(CHAIN[:-1], CHAIN[1:])) + EXTRA
for a, b in pairs:
    if a not in info or b not in info:
        R["pairs"].append({"pair": [a, b], "status": "missing"}); continue
    (La, Ta), (Lb, Tb) = info[a], info[b]
    best = None
    for la in [l for l in La if l["n"] >= 8]:
        for lb in [l for l in Lb if l["n"] >= 8]:
            dc = (la["c"] - lb["c"]).length
            if best is None or dc < best[0]:
                best = (dc, la, lb)
    ov = len(Ta.overlap(Tb))
    # closest approach of the two surfaces (sampled from a's vertices)
    ob = bpy.data.objects[a]; M = ob.matrix_world
    dmin = min(Tb.find_nearest(M @ v.co)[3] for v in ob.data.vertices)
    rec = {"pair": [a, b], "surface_min_gap_mm": round(dmin * 1000, 2), "overlapping_face_pairs": ov}
    if best:
        dc, la, lb = best
        mp = sum(min((p - q).length for q in lb["pts"]) for p in la["pts"]) / len(la["pts"])
        rec.update({"end_loops": [la["n"], lb["n"]], "end_radii_mm": [round(la["r"] * 1000, 1), round(lb["r"] * 1000, 1)],
                    "end_centroid_gap_mm": round(dc * 1000, 2), "end_mean_point_gap_mm": round(mp * 1000, 2),
                    "end_centroids": [[round(x, 4) for x in la["c"]], [round(x, 4) for x in lb["c"]]]})
    else:
        rec["end_loops"] = "no open ends"
    R["pairs"].append(rec)
json.dump(R, open(sys.argv[-1], "w"), indent=1)
for p in R["pairs"]:
    print("PAIR", json.dumps({k: v for k, v in p.items() if k != "end_centroids"}))
for k, v in R["organs"].items():
    print("ORG", k, v)
