import bpy, bmesh
from mathutils import Vector
from mathutils.bvhtree import BVHTree
O = Vector((0.0, 0.0565, 0.790)); n = Vector((0.0, 0.4793, -0.8776)); u = Vector((0.0, 0.8776, 0.4793))
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data); t = BVHTree.FromBMesh(bm)
row=[]
for a in range(-20, 21, 2):
  for b in (0.0, 0.003):
    P = O + u*(a/1000) + Vector((b,0,0)); h = t.ray_cast(P - n*0.03, n, 0.1)[0]
    row.append("%+.1f" % ((h-O).dot(n)*1000) if h else " -- ")
print("MID", " ".join(row[0::2])); print("B3 ", " ".join(row[1::2]))
