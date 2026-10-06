"""Self-intersections of a mesh object (default AN_AnalCanal): face pairs that overlap and share no vertex.
ANUS_OPEN env = Open key value. Args: -- [object name]"""
import bpy, bmesh, sys
from mathutils.bvhtree import BVHTree
nm = sys.argv[-1] if len(sys.argv) > 1 and "--" in sys.argv and sys.argv[-1] != "--" else "AN_AnalCanal"
import os
ob = bpy.data.objects[nm]
if ob.data.shape_keys and "Open" in ob.data.shape_keys.key_blocks:     # ANUS_OPEN env: test that open state
    ob.data.shape_keys.key_blocks["Open"].value = float(os.environ.get("ANUS_OPEN", "0"))
bm = bmesh.new(); bm.from_object(ob, bpy.context.evaluated_depsgraph_get()); bm.faces.ensure_lookup_table()
t = BVHTree.FromBMesh(bm)
fv = [set(v.index for v in f.verts) for f in bm.faces]
bad = [(a, b) for a, b in t.overlap(t) if a < b and not (fv[a] & fv[b])]
print("SELFX", nm, "faces", len(bm.faces), "intersecting_pairs", len(bad), bad[:10])
from mathutils import Vector
O = Vector((0.0, 0.0565, 0.790))
near = [p for p in bad if min((bm.faces[i].calc_center_median() - O).length for i in p) < 0.03]
print("SELFX near_anus(<30mm)", len(near), near[:6])
