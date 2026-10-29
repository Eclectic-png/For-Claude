import bpy, bmesh
from mathutils import Vector
A = Vector((0.0, 0.0552, 0.792)); u = Vector((0.0, 0.8776, 0.4793))
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data)
bv = {v for e in bm.edges if e.is_boundary for v in e.verts if (v.co - A).length < 0.006}
pts = [((v.co - A).dot(u) * 1000, v.co.x * 1000) for v in bv]
print("HOLE boundary verts near anus:", len(pts), " AP range %.2f..%.2f mm" % (min(p[0] for p in pts), max(p[0] for p in pts)))
for a0 in [x * 0.25 for x in range(-9, 10)]:
    xs = [x for a, x in pts if abs(a - a0) < 0.13]
    if xs:
        print("HOLE a=%+5.2f  width %.3f mm  (n=%d)" % (a0, max(xs) - min(xs), len(xs)))
