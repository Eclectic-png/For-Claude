"""View along the cleft at a low angle, to judge rises and dips. Args: -- <outdir> <prefix> [shaded]
(default: skin-only Workbench; "shaded": EEVEE with the real material, all objects)"""
import bpy, sys
from mathutils import Vector
SHADED = sys.argv[-1] == "shaded"
OUT, PRE = (sys.argv[-3], sys.argv[-2]) if SHADED else (sys.argv[-2], sys.argv[-1])
sc = bpy.context.scene
O = Vector((0.0, 0.0565, 0.790)); n = Vector((0.0, 0.6035, -0.7973)); u = Vector((0.0, 0.7973, 0.6035))
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 70; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
cam.location = O - u * 0.035 + n * 0.011                # from the perineum side, low over the surface
cam.rotation_euler = (O - cam.location).to_track_quat('-Z', 'Y').to_euler()
if SHADED:
    li = bpy.data.objects.new("GL", bpy.data.lights.new("GL", 'SUN')); li.data.energy = 3; sc.collection.objects.link(li)
    li.rotation_euler = (-(n + Vector((0.4, 0, 0.5)))).to_track_quat('-Z', 'Y').to_euler()
    sc.render.engine = "BLENDER_EEVEE_NEXT"
else:
    for o in sc.objects:
        if o.type == 'MESH' and o.name != "Hips": o.hide_render = True
    sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading; sh.light = 'STUDIO'
    sh.color_type = 'SINGLE'; sh.single_color = (0.8, 0.75, 0.72)
sc.render.resolution_x = sc.render.resolution_y = 900
sc.render.filepath = f"{OUT}/{PRE}.png"; bpy.ops.render.render(write_still=True)
