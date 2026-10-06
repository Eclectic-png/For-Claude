"""Full-resolution EEVEE close-up centred on the anus (no downscaling). ANUS_OPEN env = Open key value. Args: -- <outdir> <prefix> <distance_m> <px>"""
import bpy, sys
from mathutils import Vector
OUT, PRE, DIST, PX = sys.argv[-4], sys.argv[-3], float(sys.argv[-2]), int(sys.argv[-1])
sc = bpy.context.scene
_open = float(__import__("os").environ.get("ANUS_OPEN", "0"))   # value for every "Open" shape key (build5_passage)
for _o in bpy.data.objects:
    if _o.type == 'MESH' and _o.data.shape_keys and "Open" in _o.data.shape_keys.key_blocks:
        _o.data.shape_keys.key_blocks["Open"].value = _open
O = Vector((0.0, 0.0565, 0.790)); n = Vector((0.0, 0.6035, -0.7973))
cam = bpy.data.objects.new("HQ", bpy.data.cameras.new("HQ")); cam.data.lens = 85; cam.data.clip_start = 0.001
sc.collection.objects.link(cam); sc.camera = cam
li = bpy.data.objects.new("HQL", bpy.data.lights.new("HQL", 'SUN')); li.data.energy = 3; sc.collection.objects.link(li)
li.rotation_euler = (-(n + Vector((0.4, 0, 0.5)))).to_track_quat('-Z', 'Y').to_euler()
cam.location = O + n * DIST; cam.rotation_euler = (O - cam.location).to_track_quat('-Z', 'Y').to_euler()
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = PX
sc.render.resolution_percentage = 100
sc.render.filepath = f"{OUT}/{PRE}.png"; bpy.ops.render.render(write_still=True)
