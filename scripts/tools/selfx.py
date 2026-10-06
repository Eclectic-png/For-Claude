"""Self-intersections of a mesh object (default AN_AnalCanal): face pairs that overlap and share no vertex.
Args: -- [object name]"""
import bpy, bmesh, sys
from mathutils.bvhtree import BVHTree
nm = sys.argv[-1] if len(sys.argv) > 1 and "--" in sys.argv and sys.argv[-1] != "--" else "AN_AnalCanal"
ob = bpy.data.objects[nm]; bm = bmesh.new(); bm.from_mesh(ob.data); bm.faces.ensure_lookup_table()
t = BVHTree.FromBMesh(bm)
fv = [set(v.index for v in f.verts) for f in bm.faces]
bad = [(a, b) for a, b in t.overlap(t) if a < b and not (fv[a] & fv[b])]
print("SELFX", nm, "faces", len(bm.faces), "intersecting_pairs", len(bad), bad[:10])
