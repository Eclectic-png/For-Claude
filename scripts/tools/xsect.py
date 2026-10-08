"""2D cross-sections of two organs through a point (planes normal to x, y, z), drawn as outlines.
Args: -- <objA> <objB> <x> <y> <z> <half_width_m> <out.png>"""
import bpy, bmesh, sys
from mathutils import Vector
from PIL import Image, ImageDraw
a, b, x, y, z, hw, OUT = sys.argv[-7:]
P0 = Vector((float(x), float(y), float(z))); HW = float(hw)
cols = {a: (40, 110, 220), b: (220, 80, 60)}
img = Image.new("RGB", (1830, 620), "white"); dr = ImageDraw.Draw(img)
axes = [(Vector((1, 0, 0)), 1, 2), (Vector((0, 1, 0)), 0, 2), (Vector((0, 0, 1)), 0, 1)]
for k, (nrm, i, j) in enumerate(axes):
    ox = k * 615 + 305; oy = 310; sc = 300 / HW
    dr.text((k * 615 + 8, 6), "plane normal %s" % "xyz"[k], fill="black")
    for nm in (a, b):
        ob = bpy.data.objects[nm]; bm = bmesh.new(); bm.from_mesh(ob.data); bm.transform(ob.matrix_world)
        bmesh.ops.triangulate(bm, faces=bm.faces[:])
        for f in bm.faces:
            ds = [(v.co - P0).dot(nrm) for v in f.verts]
            pts = []
            for m in range(3):
                d0, d1 = ds[m], ds[(m + 1) % 3]
                if (d0 < 0) != (d1 < 0):
                    t = d0 / (d0 - d1); q = f.verts[m].co.lerp(f.verts[(m + 1) % 3].co, t) - P0
                    pts.append((ox + q[i] * sc, oy - q[j] * sc))
            if len(pts) == 2 and all(abs(p[0] - ox) < 305 and abs(p[1] - oy) < 305 for p in pts):
                dr.line(pts, fill=cols[nm], width=2)
        bm.free()
    dr.ellipse((ox - 3, oy - 3, ox + 3, oy + 3), fill="black")
img.save(OUT); print("XSECT", OUT)
