"""Diagnostic: geometry + wireframe close-ups of the anus from several angles. Args: -- <outdir> <prefix>"""
import bpy, sys
from mathutils import Vector
OUT, PRE = sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
A = Vector((0.0, 0.0568, 0.7898)); n = Vector((0.0, 0.6035, -0.7973))
hips = bpy.data.objects["Hips"]
wire = hips.copy(); wire.data = hips.data.copy(); sc.collection.objects.link(wire)
wm = wire.modifiers.new("W", 'WIREFRAME'); wm.thickness = 0.00012; wm.use_replace = True
mw = bpy.data.materials.new("wire"); mw.diffuse_color = (0.05, 0.05, 0.08, 1)
mb = bpy.data.materials.new("skin"); mb.diffuse_color = (0.8, 0.75, 0.72, 1)
wire.data.materials.clear(); wire.data.materials.append(mw)
hips.data.materials.clear(); hips.data.materials.append(mb)
for o in sc.objects:
    if o.type == 'MESH' and o not in (hips, wire):
        o.hide_render = True
cam = bpy.data.objects.new("RC", bpy.data.cameras.new("RC")); cam.data.lens = 85; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading; sh.light = 'STUDIO'; sh.color_type = 'MATERIAL'
sc.render.resolution_x = sc.render.resolution_y = 800
views = {"n": (n, 0.065), "oblique": ((n + Vector((0.7, 0, 0))).normalized(), 0.065),
         "front": ((n + Vector((0, -0.6, -0.2))).normalized(), 0.07)}
for vname, (dirv, dist) in views.items():
    cam.location = A + dirv * dist
    cam.rotation_euler = (A - cam.location).to_track_quat('-Z', 'Y').to_euler()
    for tag, show in (("wire", True), ("geo", False)):
        wire.hide_render = not show
        sc.render.filepath = f"{OUT}/{PRE}_{vname}_{tag}.png"; bpy.ops.render.render(write_still=True)
