import bpy, bmesh
from mathutils import Vector
A = Vector((0.0, 0.0552, 0.792)); u = Vector((0.0, 0.8776, 0.4793))
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data)
for e in bm.edges:
    if len(e.link_faces) > 2:
        c = (e.verts[0].co + e.verts[1].co) / 2
        fs = e.link_faces
        sets = [frozenset(v.index for v in f.verts) for f in fs]
        print("NM edge at a=%+.2f mm lat=%+.3f  faces=%d dupfaces=%d areas=%s" % ((c - A).dot(u) * 1000, c.x * 1000, len(fs),
              len(sets) - len(set(sets)), [round(f.calc_area() * 1e6, 4) for f in fs]))
bnd = [e for e in bm.edges if e.is_boundary and all((v.co - A).length < 0.006 for v in e.verts)]
seen = set(); comps = []
for e0 in bnd:
    if e0 in seen: continue
    comp = []; st = [e0]
    while st:
        e = st.pop()
        if e in seen: continue
        seen.add(e); comp.append(e); st += [g for v in e.verts for g in v.link_edges if g in bnd and g not in seen]
    comps.append(comp)
for c in comps:
    pts = [(v.co - A).dot(u) * 1000 for e in c for v in e.verts]
    print("LOOP edges=%d  a from %+.2f to %+.2f mm" % (len(c), min(pts), max(pts)))
