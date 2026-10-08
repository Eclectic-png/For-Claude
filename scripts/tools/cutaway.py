"""EEVEE cut-away of the digestive tract (anatomy_tract_passage.blend): the meshes with their front half (y below the
plane through the point) removed, seen from the front, everything else hidden. Args: -- <x> <y> <z> <half_width_m> <out.png>"""
import bpy, bmesh, sys
from mathutils import Vector
x, y, z, hw, OUT = sys.argv[-5:]
P = Vector((float(x), float(y), float(z))); HW = float(hw); sc = bpy.context.scene
TR = {"AN_Esophagus", "AN_Stomach", "AN_Duodenum", "AN_Jejunum", "AN_Ileum", "AN_Cecum", "AN_Appendix",
      "AN_Colon_Ascending", "AN_Colon_Transverse", "AN_Colon_Descending", "AN_Rectum"}
for o in list(sc.objects):
    if o.type != 'MESH':
        continue
    if o.name not in TR:
        o.hide_render = True; continue
    bm = bmesh.new(); bm.from_mesh(o.data); M = o.matrix_world
    bmesh.ops.delete(bm, geom=[f for f in bm.faces if (M @ f.calc_center_median()).y < P.y], context='FACES')
    bm.to_mesh(o.data); bm.free()
cam = bpy.data.objects.new("CW", bpy.data.cameras.new("CW")); cam.data.type = 'ORTHO'; cam.data.ortho_scale = 2 * HW
sc.collection.objects.link(cam); sc.camera = cam; cam.location = P + Vector((0, -0.5, 0))
cam.rotation_euler = Vector((0, 1, 0)).to_track_quat('-Z', 'Z').to_euler(); cam.data.clip_end = 2
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = 900
sc.render.filepath = OUT; bpy.ops.render.render(write_still=True)
