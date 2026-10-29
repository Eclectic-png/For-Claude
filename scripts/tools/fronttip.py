import bpy, sys
from mathutils import Vector
OUT, PRE, MODE = sys.argv[-3], sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
A = Vector((0.0, 0.0552, 0.792)); n = Vector((0.0, 0.6035, -0.7973)); u = Vector((0.0, 0.8776, 0.4793))
tgt = A - u * 0.0017                       # front (bottom-in-image) tip of the slit
cam = bpy.data.objects.new("RC", bpy.data.cameras.new("RC")); cam.data.lens = 85; cam.data.clip_start = 0.0005
sc.collection.objects.link(cam); sc.camera = cam
li = bpy.data.objects.new("RL", bpy.data.lights.new("RL", 'SUN')); li.data.energy = 3; sc.collection.objects.link(li)
li.rotation_euler = (-(n + Vector((0.4, 0, 0.5)))).to_track_quat('-Z', 'Y').to_euler()
cam.location = tgt + n * 0.014; cam.rotation_euler = (tgt - cam.location).to_track_quat('-Z', 'Y').to_euler()
if "nocanal" in MODE:
    bpy.data.objects["AN_AnalCanal"].hide_render = True
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = 700
sc.render.filepath = f"{OUT}/{PRE}_{MODE}.png"; bpy.ops.render.render(write_still=True); print("OK", MODE)
