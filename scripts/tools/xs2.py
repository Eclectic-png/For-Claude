import bpy, bmesh
from mathutils import Vector
from mathutils.bvhtree import BVHTree
O = Vector((0.0, 0.0552, 0.792)); n = Vector((0.0, 0.4793, -0.8776)); u = Vector((0.0, 0.8776, 0.4793))
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data); t = BVHTree.FromBMesh(bm)
row = []
for b in [x * 0.5 for x in range(-12, 13)]:
    P = O + Vector((b / 1000, 0, 0)); h = t.ray_cast(P - n * 0.03, n, 0.1)[0]
    row.append("%+.1f" % ((h - O).dot(n) * 1000) if h else "--")
print("XS2", " ".join(row))
