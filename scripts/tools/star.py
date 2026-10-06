"""Geometry-only (Workbench) view of the whole crease star, skin only. Args: -- <outdir> <prefix> <dist_m> <wire 0/1>"""
import bpy, sys
from mathutils import Vector
OUT, PRE, DIST, WIRE = sys.argv[-4], sys.argv[-3], float(sys.argv[-2]), sys.argv[-1] == "1"
sc = bpy.context.scene
A = Vector((0.0, 0.0559, 0.7909)); n = Vector((0.0, 0.6035, -0.7973))
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.lens = 85; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
d = (n + Vector((0.15, 0, 0))).normalized(); cam.location = A + d * DIST
cam.rotation_euler = (A - cam.location).to_track_quat('-Z', 'Y').to_euler()
for o in sc.objects:
    if o.type == 'MESH' and o.name != "Hips": o.hide_render = True
sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading; sh.light = 'STUDIO'
sh.color_type = 'SINGLE'; sh.single_color = (0.8, 0.75, 0.72)
sc.render.resolution_x = sc.render.resolution_y = 1200
if WIRE:
    h = bpy.data.objects["Hips"]; w = h.copy(); w.data = h.data.copy(); sc.collection.objects.link(w)
    m = w.modifiers.new("W", "WIREFRAME"); m.thickness = 0.00001; m.use_replace = True
    mt = bpy.data.materials.new("k"); mt.diffuse_color = (0.1, 0.2, 1, 1); w.data.materials.clear(); w.data.materials.append(mt)
    sh.color_type = "MATERIAL"; h.data.materials.clear()
sc.render.filepath = f"{OUT}/{PRE}.png"; bpy.ops.render.render(write_still=True)
