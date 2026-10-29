"""Workbench view of one atlas organ with its boundary (hole) edges drawn as red tubes. Args: -- <obj> <out.png> [view]"""
import bpy, bmesh, sys
from mathutils import Vector
NM, OUT = sys.argv[-2], sys.argv[-1]
sc = bpy.context.scene
for o in sc.objects:
    if o.type == 'MESH' and o.name != NM: o.hide_render = True
ob = bpy.data.objects[NM]; M = ob.matrix_world
bm = bmesh.new(); bm.from_mesh(ob.data)
cm = bpy.data.meshes.new("holes"); hb = bmesh.new()
for e in bm.edges:
    if e.is_boundary:
        a, b = [hb.verts.new(M @ v.co) for v in e.verts]; hb.edges.new((a, b))
hb.to_mesh(cm); ho = bpy.data.objects.new("holes", cm); sc.collection.objects.link(ho)
sk = ho.modifiers.new("s", "SKIN"); 
for v in cm.skin_vertices[0].data: v.radius = (0.0012, 0.0012)
mt = bpy.data.materials.new("r"); mt.diffuse_color = (1, 0, 0, 1); cm.materials.append(mt); ho.color = (1, 0, 0, 1)
P = [M @ v.co for v in ob.data.vertices]; c = sum(P, Vector()) / len(P)
ext = max((p - c).length for p in P)
cam = bpy.data.objects.new("C", bpy.data.cameras.new("C")); cam.data.type = 'ORTHO'; cam.data.ortho_scale = ext * 2.2
sc.collection.objects.link(cam); sc.camera = cam
cam.location = c + Vector((0, -1, 0)); cam.rotation_euler = Vector((0, 1, 0)).to_track_quat('-Z', 'Z').to_euler()
cam.data.clip_end = 10
sc.render.engine = "BLENDER_WORKBENCH"; sh = sc.display.shading; sh.light = 'STUDIO'; sh.color_type = 'OBJECT'
ob.color = (0.8, 0.75, 0.7, 1); sh.show_backface_culling = False
sc.render.resolution_x = sc.render.resolution_y = 1000; sc.render.filepath = OUT; bpy.ops.render.render(write_still=True)
