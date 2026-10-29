"""2D section of many objects through one plane, drawn as coloured outlines on a 5 mm grid (labels in metres).
Shape keys / modifiers are evaluated, so set them first (SAG_KEYS='Bladder_Fill=1,Rectum_Fill=1' sets custom props on
Pelvic_Controls before drawing).
Args: -- <axis x|y|z> <offset_m> <centre_u> <centre_v> <half_width_m> <out.png> [obj ...]   (no objs = all meshes)"""
import bpy, bmesh, sys, os, colorsys
from mathutils import Vector
from PIL import Image, ImageDraw

argv = sys.argv[sys.argv.index("--") + 1:]
ax, off, cu, cv, hw, OUT = argv[0], float(argv[1]), float(argv[2]), float(argv[3]), float(argv[4]), argv[5]
names = argv[6:] or [o.name for o in bpy.data.objects if o.type == 'MESH' and not o.hide_render]
ctl = bpy.data.objects.get("Pelvic_Controls") or bpy.data.objects.get("Urinary_Controls")
for kv in filter(None, os.environ.get("SAG_KEYS", "").split(",")):
    k, v = kv.split("="); ctl[k] = float(v); ctl.update_tag()
bpy.context.view_layer.update()
k = "xyz".index(ax); i, j = [(1, 2), (0, 2), (0, 1)][k]
S = 1600; sc = S / (2 * hw)
img = Image.new("RGB", (S, S + 40), "white"); dr = ImageDraw.Draw(img)
to_px = lambda u, v: ((u - cu) * sc + S / 2, S / 2 - (v - cv) * sc)
g = 0.005
u0 = (cu - hw) // g * g
while u0 < cu + hw:
    x_, _ = to_px(u0, cv); dr.line([(x_, 0), (x_, S)], fill=(235, 235, 235)); u0 += g
v0 = (cv - hw) // g * g
while v0 < cv + hw:
    _, y_ = to_px(cu, v0); dr.line([(0, y_), (S, y_)], fill=(235, 235, 235))
    dr.text((2, y_ + 1), "%.3f" % v0, fill=(150, 150, 150)); v0 += g
dg = bpy.context.evaluated_depsgraph_get()
for n_, nm in enumerate(names):
    ob = bpy.data.objects[nm]; bm = bmesh.new(); bm.from_object(ob, dg); bm.transform(ob.matrix_world)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    col = tuple(int(255 * c) for c in colorsys.hsv_to_rgb((n_ * 0.618) % 1, 0.8, 0.75))
    for f in bm.faces:
        ds = [v.co[k] - off for v in f.verts]; pts = []
        for m in range(3):
            d0, d1 = ds[m], ds[(m + 1) % 3]
            if (d0 < 0) != (d1 < 0):
                q = f.verts[m].co.lerp(f.verts[(m + 1) % 3].co, d0 / (d0 - d1)); pts.append(to_px(q[i], q[j]))
        if len(pts) == 2:
            dr.line(pts, fill=col, width=2)
    bm.free()
    dr.text((8 + (n_ % 6) * 265, S + 4 + (n_ // 6) * 12), nm, fill=col)
img.save(OUT); print("SAGITTAL", OUT)
