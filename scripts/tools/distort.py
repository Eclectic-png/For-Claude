"""How far each built version moved the skin, per ORIGINAL Hips vertex (distance to the version's final skin, mm).
Run on the untouched original: blender -b Hips.blend --python distort.py -- <version.blend> [...] <out.json>"""
import bpy, bmesh, sys, json, numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree
OUTF = sys.argv[-1]; paths = sys.argv[sys.argv.index("--") + 1:-1]
O = Vector((0.0, 0.0565, 0.790))
orig = bpy.data.objects["Hips"]; assert orig.matrix_world.is_identity if hasattr(orig.matrix_world, "is_identity") else True
P = [orig.matrix_world @ v.co for v in orig.data.vertices]
res = {"xyz": [list(p) for p in P], "d_anus_mm": [(p - O).length * 1000 for p in P], "versions": {}}
for path in paths:
    with bpy.data.libraries.load(path) as (src, dst):
        dst.objects = ["Hips"]
    ob = dst.objects[0]
    bm = bmesh.new(); bm.from_mesh(ob.data); bm.transform(ob.matrix_world); t = BVHTree.FromBMesh(bm)
    res["versions"][path] = [t.find_nearest(p)[3] * 1000 for p in P]
    bm.free(); bpy.data.objects.remove(ob)
json.dump(res, open(OUTF, "w"))
print("DONE", len(P))
