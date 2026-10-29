import bpy, bmesh
from mathutils import Vector, kdtree
from mathutils.bvhtree import BVHTree
hips = bpy.data.objects["Hips"]; bm = bmesh.new(); bm.from_mesh(hips.data); tree = BVHTree.FromBMesh(bm)
dirs = [Vector(d).normalized() for d in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,1),(-1,1,-1),(1,-1,1))]
def votes(p):
    t = 0
    for d in dirs:
        c = 0; q = p.copy()
        while True:
            h = tree.ray_cast(q, d, 2.0)
            if h[0] is None: break
            c += 1; q = h[0] + d * 1e-6
        t += c % 2
    return t
o = bpy.data.objects["AN_AnalCanal"]
for i, v in enumerate(o.data.vertices):
    p = o.matrix_world @ v.co; vt = votes(p)
    if vt < 5:
        print("CO vert", i, "ring", i // 36 if i >= 144 else "bottom(on skin)", "votes", vt, "dist_to_skin_mm %.3f" % (tree.find_nearest(p)[3] * 1000))
O = Vector((0.0, 0.0565, 0.790))
me = hips.data; kd = kdtree.KDTree(len(me.vertices))
for v in me.vertices: kd.insert(v.co, v.index)
kd.balance()
near = [v for v in me.vertices if (v.co - O).length < 0.03]
bad = [v for v in near if kd.find(Vector((-v.co.x, v.co.y, v.co.z)))[2] > 2e-5]
print("CO cleft(<30mm) unmatched mirror verts", len(bad), "of", len(near), "max err mm %.3f" % (max([kd.find(Vector((-v.co.x, v.co.y, v.co.z)))[2] for v in near]) * 1000))
