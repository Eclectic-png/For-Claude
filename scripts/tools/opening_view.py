"""EEVEE view of one organ alone, looking into its opening towards a neighbour (anatomy_tract_passage.blend).
Args: -- <organ> <neighbour> <distance_m> <out.png>"""
import bpy, bmesh, sys
from mathutils import Vector
from mathutils.bvhtree import BVHTree
nm, nb, dist, OUT = sys.argv[-4:]; dist = float(dist); sc = bpy.context.scene
for o in sc.objects:
    if o.type == 'MESH':
        o.hide_render = o.name != nm
ob = bpy.data.objects[nm]; M = ob.matrix_world
bm = bmesh.new(); bm.from_mesh(ob.data); bm.transform(M); bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7)
nbm = bmesh.new(); nbm.from_mesh(bpy.data.objects[nb].data); nbm.transform(bpy.data.objects[nb].matrix_world)
tN = BVHTree.FromBMesh(nbm)
lv = [v.co.copy() for v in bm.verts if v.is_boundary and tN.find_nearest(v.co)[3] < 2e-6]
c = sum(lv, Vector()) / len(lv); oc = sum((v.co for v in bm.verts), Vector()) / len(bm.verts)
d = (c - oc).normalized()
cam = bpy.data.objects.new("OV", bpy.data.cameras.new("OV")); cam.data.lens = 50; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam; cam.location = c + d * dist
cam.rotation_euler = (-d).to_track_quat('-Z', 'Z').to_euler()
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = 700
sc.render.filepath = OUT; bpy.ops.render.render(write_still=True)
