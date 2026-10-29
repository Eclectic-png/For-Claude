"""Close-up renders of the anus: geometry-only (Workbench) and shaded (EEVEE). Args: -- <outdir> <prefix>"""
import bpy, sys, math
from mathutils import Vector
OUT, PRE = sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
A = Vector((0.0, 0.0568, 0.7898))
n = Vector((0.0, 0.6035, -0.7973))           # outward, opposite the canal axis
cam_d = bpy.data.cameras.new("RC"); cam_d.lens = 85; cam_d.clip_start = 0.001
cam = bpy.data.objects.new("RC", cam_d); sc.collection.objects.link(cam); sc.camera = cam
li_d = bpy.data.lights.new("RL", 'SUN'); li_d.energy = 3
li = bpy.data.objects.new("RL", li_d); sc.collection.objects.link(li)
sc.render.resolution_x = sc.render.resolution_y = 1200
views = {"close": (n, 0.11), "below": ((n + Vector((0, 0, -0.9))).normalized(), 0.30)}
for vname, (dirv, dist) in views.items():
    cam.location = A + dirv * dist
    cam.rotation_euler = (A - cam.location).to_track_quat('-Z', 'Y').to_euler()
    li.rotation_euler = (-(dirv + Vector((0.4, 0, 0.5)))).to_track_quat('-Z', 'Y').to_euler()
    for eng in ("BLENDER_WORKBENCH", "BLENDER_EEVEE_NEXT"):
        sc.render.engine = eng
        if eng == "BLENDER_WORKBENCH":
            sh = sc.display.shading; sh.light = 'STUDIO'; sh.color_type = 'SINGLE'; sh.single_color = (0.8, 0.75, 0.72)
            tag = "wb"
        else:
            tag = "std"
        sc.render.filepath = f"{OUT}/{PRE}_{vname}_{tag}.png"
        try:
            bpy.ops.render.render(write_still=True); print("RENDERED", sc.render.filepath)
        except Exception as ex:
            print("FAILED", eng, ex)
