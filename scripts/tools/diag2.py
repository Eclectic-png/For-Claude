import bpy, sys
from mathutils import Vector
OUT, PRE, MODE = sys.argv[-3], sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
A = Vector((0.0, 0.0568, 0.7898)); n = Vector((0.0, 0.6035, -0.7973))
hips = bpy.data.objects["Hips"]
cam = bpy.data.objects.new("RC", bpy.data.cameras.new("RC")); cam.data.lens = 85; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
li = bpy.data.objects.new("RL", bpy.data.lights.new("RL", 'SUN')); li.data.energy = 3; sc.collection.objects.link(li)
li.rotation_euler = (-(n + Vector((0.4, 0, 0.5)))).to_track_quat('-Z', 'Y').to_euler()
tgt, dist = A, 0.11
if MODE.startswith("zoom"):
    tgt, dist = A + Vector((0, 0.0045, 0.0035)) * 1.0, 0.035
cam.location = tgt + n * dist; cam.rotation_euler = (tgt - cam.location).to_track_quat('-Z', 'Y').to_euler()
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = 1200
if "wire" in MODE:
    w = hips.copy(); w.data = hips.data.copy(); sc.collection.objects.link(w)
    m = w.modifiers.new("W", 'WIREFRAME'); m.thickness = 0.00008 if MODE.startswith("zoom") else 0.00015; m.use_replace = True
    mat = bpy.data.materials.new("wire"); mat.use_nodes = True
    bs = mat.node_tree.nodes["Principled BSDF"]; bs.inputs["Base Color"].default_value = (0, 0.2, 1, 1)
    bs.inputs["Emission Color"].default_value = (0, 0.3, 1, 1); bs.inputs["Emission Strength"].default_value = 2
    w.data.materials.clear(); w.data.materials.append(mat)
if "subsurf" in MODE:
    m = hips.modifiers.new("S", 'SUBSURF'); m.levels = 2; m.render_levels = 2
sc.render.filepath = f"{OUT}/{PRE}_{MODE}.png"; bpy.ops.render.render(write_still=True); print("OK", MODE)
