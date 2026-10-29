import bpy, bmesh
from mathutils.bvhtree import BVHTree
hips = bpy.data.objects["Hips"]; bm = bmesh.new(); bm.from_mesh(hips.data); bm.normal_update(); tree = BVHTree.FromBMesh(bm)
# outward check: the skin normal at the navel-side midline should point forward (-y)
for nm in ("AN_AnalSphincter", "AN_AnalCanal", "AN_Rectum_LowerAmpulla", "AN_Rectum", "AN_Colon_Descending"):
    o = bpy.data.objects[nm]; out = 0; mind = 9; tot = 0
    for v in o.data.vertices:
        p = o.matrix_world @ v.co; loc, nrm, idx, d = tree.find_nearest(p)
        if d < 0.0002: continue            # the canal's bottom ring sits exactly on the pucker's inner ring
        tot += 1; mind = min(mind, d)
        if (p - loc).dot(nrm) > 0: out += 1
    print("INS %-24s outside_skin=%d/%d  min_gap_to_skin=%.1fmm" % (nm, out, tot, mind * 1000))
