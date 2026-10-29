"""Skin-only Workbench close-up of a point on the slit, looking along the outward normal.
Args: -- <outdir> <prefix> <mm along the slit, + = back tip> <distance_m> <wire 0/1>. Image top is the back end."""
import bpy, sys
from mathutils import Vector
OUT, PRE, OFF, DIST, WIRE = sys.argv[-5], sys.argv[-4], float(sys.argv[-3]), float(sys.argv[-2]), sys.argv[-1] == "1"
sc = bpy.context.scene
O = Vector((0.0, 0.0565, 0.790)); n = Vector((0.0, 0.6035, -0.7973)); u = Vector((0.0, 0.7973, 0.6035))
tgt = O + u * (OFF / 1000)
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 85; cam.data.clip_start = 0.0005
sc.collection.objects.link(cam); sc.camera = cam
cam.location = tgt + (n + Vector((0.15, 0, 0))).normalized() * DIST
cam.rotation_euler = (tgt - cam.location).to_track_quat('-Z', 'Y').to_euler()
for o in sc.objects:
    if o.type == 'MESH' and o.name != "Hips": o.hide_render = True
sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading; sh.light = 'STUDIO'
sh.color_type = 'SINGLE'; sh.single_color = (0.8, 0.75, 0.72)
sc.render.resolution_x = sc.render.resolution_y = 900
if WIRE:
    h = bpy.data.objects["Hips"]; w = h.copy(); w.data = h.data.copy(); sc.collection.objects.link(w)
    m = w.modifiers.new("W", "WIREFRAME"); m.thickness = 0.000004; m.use_replace = True
    mt = bpy.data.materials.new("k"); mt.diffuse_color = (0.1, 0.2, 1, 1); w.data.materials.clear(); w.data.materials.append(mt)
    sh.color_type = "MATERIAL"; h.data.materials.clear()
sc.render.filepath = f"{OUT}/{PRE}.png"; bpy.ops.render.render(write_still=True)
