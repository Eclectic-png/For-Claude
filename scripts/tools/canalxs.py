"""Cross-sections of AN_AnalCanal (build5_passage layout: 288-point rings from the mouth, then 144 / 72 / 36) drawn
in the canal frame (x lateral, y front/back), each scaled to fit its own panel. Args: -- <out.png>"""
import bpy, sys, math
from mathutils import Vector
from PIL import Image, ImageDraw
OUT = sys.argv[-1]
me = bpy.data.objects["AN_AnalCanal"].data; V = [v.co.copy() for v in me.vertices]
sizes = [288] * 26 + [144, 72, 36]; S_MM = [0, 0.15, 0.35, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.25, 5.0, 5.75, 6.5, 7.25,
                                         8.0, 9.0, 10.0, 11.0, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.2, 23.4, 24.7]
rings = []; i = 0
for n in sizes:
    rings.append(V[i:i + n]); i += n
assert i == len(V)
ctr = [sum(r, Vector()) / len(r) for r in rings]
d = (ctr[-1] - ctr[0]).normalized(); X = Vector((1, 0, 0)); e1 = (X - d * X.dot(d)).normalized(); e2 = d.cross(e1)
pick = [0, 6, 9, 11, 13, 15, 17, 19, 21, 23, 25, 28]
P = 300; img = Image.new("RGB", (P * 6, P * 2 + 40), "white"); dr = ImageDraw.Draw(img)
for j, ri in enumerate(pick):
    r = rings[ri]; c = ctr[ri]
    pts = [((p - c).dot(e1) * 1000, (p - c).dot(e2) * 1000) for p in r]
    ext = max(max(abs(x), abs(y)) for x, y in pts) * 1.15
    ox, oy = (j % 6) * P + P / 2, (j // 6) * (P + 20) + P / 2 + 20
    sc = (P / 2) / ext
    dr.line([(ox + x * sc, oy - y * sc) for x, y in pts] + [(ox + pts[0][0] * sc, oy - pts[0][1] * sc)], fill="black", width=2)
    dr.text((ox - P / 2 + 6, oy - P / 2 - 16), "s=%.1f mm  (panel %.1f mm wide)" % (S_MM[ri], 2 * ext), fill="black")
img.save(OUT); print("XS saved", OUT)
