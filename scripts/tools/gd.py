import bpy, bmesh
from mathutils import Vector
from mathutils.bvhtree import BVHTree
A = Vector((0.0, 0.0552, 0.792)); n = Vector((0.0, 0.4793, -0.8776)); u = Vector((0.0, 0.8776, 0.4793))
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data); t = BVHTree.FromBMesh(bm)
out = []
for s in (0.25, 0.5, 0.8, 1.2, 1.8, 2.6):
    for side in (1, -1):
        hs = []
        for i in range(-50, 51):
            P = A + u * (i * 0.05 / 1000) + Vector((side * s / 1000, 0, 0))
            h = t.ray_cast(P - n * 0.02, n, 0.05)[0]
            if h: hs.append((h - A).dot(n) * 1000)
        # local relief: height minus a 1 mm moving average, then range
        import statistics
        rel = [hs[k] - statistics.mean(hs[max(0, k-10):k+11]) for k in range(len(hs))]
        out.append("%s%.2f:%.2f" % ("+" if side > 0 else "-", s, max(rel) - min(rel)))
print("GD", "  ".join(out))
