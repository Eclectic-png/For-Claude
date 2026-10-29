import bpy, bmesh
from mathutils import Vector
from mathutils.bvhtree import BVHTree
hips = bpy.data.objects["Hips"]; bm = bmesh.new(); bm.from_mesh(hips.data); bm.normal_update()
bad = 0
for e in bm.edges:
    if len(e.link_faces) == 2:
        f, g = e.link_faces
        def dirn(face):
            vs = [l.vert for l in face.loops]; i = vs.index(e.verts[0]); return vs[(i + 1) % len(vs)] is e.verts[1]
        if dirn(f) == dirn(g): bad += 1
print("ORI inconsistent_edges", bad)
tree = BVHTree.FromBMesh(bm)
dirs = [Vector(d).normalized() for d in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,1,1),(-1,1,-1),(1,-1,1))]
def inside(p):
    votes = 0
    for d in dirs:
        cnt = 0; q = p.copy()
        while True:
            hit = tree.ray_cast(q, d, 2.0)
            if hit[0] is None: break
            cnt += 1; q = hit[0] + d * 1e-6
        votes += cnt % 2
    return votes
for nm in ("AN_AnalSphincter", "AN_AnalCanal", "AN_Rectum_LowerAmpulla", "AN_Rectum", "AN_Colon_Descending"):
    o = bpy.data.objects[nm]; pts = [o.matrix_world @ v.co for v in o.data.vertices][::3]
    outs = sum(1 for p in pts if inside(p) < 5)
    print("ORI %-24s outside_by_ray_parity=%d/%d" % (nm, outs, len(pts)))
