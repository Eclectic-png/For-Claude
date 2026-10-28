"""EEVEE views of the pelvis with real materials. Modes:
  below <x> <y> <z> <dist>   camera under the point looking up (the vulva), everything as is;
  cut <x> <y> <z> <half_w>   every mesh's x > 0 half removed, orthographic from +x (midline cut-away).
Args: -- <mode> x y z size <out.png>   (SAG_KEYS-like PELVIS_KEYS='Bladder_Fill=1' sets Pelvic_Controls first)"""
import bpy, bmesh, sys, os
from mathutils import Vector
argv = sys.argv[sys.argv.index("--") + 1:]
mode, P, size, OUT = argv[0], Vector(tuple(float(a) for a in argv[1:4])), float(argv[4]), argv[5]
sc = bpy.context.scene; ctl = bpy.data.objects.get("Pelvic_Controls")
for kv in filter(None, os.environ.get("PELVIS_KEYS", "").split(",")):
    k, v = kv.split("="); ctl[k] = float(v); ctl.update_tag()
bpy.context.view_layer.update()
cam = bpy.data.objects.new("PC", bpy.data.cameras.new("PC")); sc.collection.objects.link(cam); sc.camera = cam
cam.data.clip_start = 0.0005; cam.data.clip_end = 3
for n in [n for n in os.environ.get("PCUT_HIDE", "").split(",") if n]:
    if n in bpy.data.objects:
        bpy.data.objects[n].hide_render = True
if mode == "below":
    d = Vector((0, 0.25, -1)).normalized(); cam.data.lens = 60
    cam.location = P + d * size; cam.rotation_euler = (-d).to_track_quat('-Z', 'Y').to_euler()
else:
    only = [n for n in os.environ.get("PCUT_ONLY", "").split(",") if n]
    if only:
        for o in sc.objects:
            if o.type == 'MESH':
                o.hide_render = o.name not in only
    dg = bpy.context.evaluated_depsgraph_get()
    for o in list(sc.objects):
        if o.type != 'MESH' or o.hide_render:
            continue
        bm = bmesh.new(); bm.from_object(o, dg); M = o.matrix_world
        bmesh.ops.delete(bm, geom=[f for f in bm.faces if (M @ f.calc_center_median()).x > 0], context='FACES')
        me = bpy.data.meshes.new(o.name + "_cut"); bm.to_mesh(me); bm.free()
        for m in o.data.materials:
            me.materials.append(m)
        if o.data.shape_keys:
            o.shape_key_clear()
        o.data = me
    cam.data.type = 'ORTHO'; cam.data.ortho_scale = 2 * size
    cam.location = P + Vector((0.3, 0, 0)); cam.rotation_euler = Vector((-1, 0, 0)).to_track_quat('-Z', 'Y').to_euler()
sc.render.engine = "BLENDER_EEVEE_NEXT"; sc.render.resolution_x = sc.render.resolution_y = 1000
sc.render.filepath = OUT; bpy.ops.render.render(write_still=True); print("PCUT", OUT)
