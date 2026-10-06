"""After build_tract_passage: every organ's open loops (welded at 1e-7 m), and for each whether a neighbour has the
same loop (max vertex distance). An unmatched loop is a real hole. Args: -- <out.json>"""
import bpy, bmesh, sys, json
from mathutils import Vector, kdtree
NAMES = ["AN_Esophagus", "AN_Stomach", "AN_Duodenum", "AN_Jejunum", "AN_Ileum", "AN_Cecum", "AN_Appendix",
         "AN_Colon_Ascending", "AN_Colon_Transverse", "AN_Colon_Descending", "AN_Rectum", "AN_Rectum_LowerAmpulla",
         "AN_AnalCanal", "Hips"]
loops = {}
for nm in NAMES:
    o = bpy.data.objects.get(nm)
    if o is None:
        continue
    bm = bmesh.new(); bm.from_mesh(o.data); bm.transform(o.matrix_world)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7); seen = set(); L = []
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
        pts = [v.co.copy() for v in {v for x in comp for v in x.verts}]
        c = sum(pts, Vector()) / len(pts)
        L.append({"pts": pts, "edges": len(comp), "r_mm": round(sum((p - c).length for p in pts) / len(pts) * 1000, 2),
                  "c": c})
    loops[nm] = L; bm.free()
res = []
for nm, L in loops.items():
    for l in L:
        best = None
        for nm2, L2 in loops.items():
            if nm2 == nm:
                continue
            allp = [p for l2 in L2 for p in l2["pts"]]
            if not allp:
                continue
            kd = kdtree.KDTree(len(allp))
            for i, p in enumerate(allp):
                kd.insert(p, i)
            kd.balance(); d = max(kd.find(p)[2] for p in l["pts"])
            if best is None or d < best[1]:
                best = (nm2, d)
        rec = {"organ": nm, "edges": l["edges"], "r_mm": l["r_mm"], "centre": [round(x, 4) for x in l["c"]],
               "matched_with": best[0] if best else None, "max_gap_mm": round(best[1] * 1000, 4) if best else None}
        res.append(rec)
        print("SEAMCHK", json.dumps(rec))
json.dump(res, open(sys.argv[-1], "w"), indent=1)
