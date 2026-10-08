import bpy, sys
from mathutils import Vector
OUT, PRE = sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
A = Vector((0.0, 0.0552, 0.792)); u = Vector((0.0, 0.8776, 0.4793)); n = Vector((0.0, 0.4793, -0.8776))
tgt = A + u * 0.0025
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 85; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
d = (n + Vector((0, -0.5, -0.2))).normalized(); cam.location = tgt + d * 0.03
cam.rotation_euler = (tgt - cam.location).to_track_quat('-Z', 'Y').to_euler()
for o in sc.objects:
    if o.type == 'MESH' and o.name != "Hips": o.hide_render = True
sc.render.engine = "BLENDER_WORKBENCH"; sc.display.shading.light = 'STUDIO'; sc.display.shading.color_type = 'SINGLE'
sc.render.resolution_x = sc.render.resolution_y = 700
sc.render.filepath = f"{OUT}/{PRE}_tip.png"; bpy.ops.render.render(write_still=True)
