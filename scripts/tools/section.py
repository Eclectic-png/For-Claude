"""Workbench cut-away on the midline (x = 0): Hips + anal canal + lower ampulla with their x > 0 half removed, seen
from the +x side, so the passage from skin into canal shows in section. Args: -- <outdir> <prefix> <half_width_m>"""
import bpy, bmesh, sys
from mathutils import Vector
OUT, PRE, HW = sys.argv[-3], sys.argv[-2], float(sys.argv[-1])
sc = bpy.context.scene
O = Vector((0.0, 0.0565, 0.790)); d = Vector((0.0, -0.6031, 0.7977))
keep = {"Hips": (0.85, 0.75, 0.72), "AN_AnalCanal": (0.75, 0.35, 0.35), "AN_Rectum_LowerAmpulla": (0.6, 0.3, 0.3)}
for ob in list(bpy.data.objects):
    if ob.type == 'MESH' and ob.name not in keep:
        ob.hide_render = True
for nm, col in keep.items():
    ob = bpy.data.objects[nm]; bm = bmesh.new(); bm.from_mesh(ob.data); M = ob.matrix_world
    bmesh.ops.delete(bm, geom=[v for v in bm.verts if (M @ v.co).x > 1e-6], context='VERTS')
    bm.to_mesh(ob.data); bm.free(); ob.color = (*col, 1.0)
c = O + d * float(__import__("os").environ.get("SEC_ALONG", "0.006"))
cam = bpy.data.objects.new("SC", bpy.data.cameras.new("SC")); cam.data.type = 'ORTHO'; cam.data.ortho_scale = 2 * HW
cam.data.clip_start = 0.001; sc.collection.objects.link(cam); sc.camera = cam
cam.location = c + Vector((0.1, 0, 0)); cam.rotation_euler = Vector((-1, 0, 0)).to_track_quat('-Z', 'Z').to_euler()
sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading
sh.light = 'STUDIO'; sh.color_type = 'OBJECT'; sh.show_backface_culling = False; sh.show_cavity = False
sc.render.resolution_x = sc.render.resolution_y = 1200; sc.render.resolution_percentage = 100
sc.render.filepath = f"{OUT}/{PRE}.png"; bpy.ops.render.render(write_still=True)
