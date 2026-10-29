"""Cross-sections of the vagina (build_vagina output): lumen ring (black), wall (pink) and cervix (blue) cut by the
ring's plane, every panel at the same scale (20 mm grid ticks). Args: -- <out.png> [ring numbers comma-separated]"""
import bpy, sys, os
import numpy as np
from PIL import Image, ImageDraw
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import vagina_lib as V
argv = sys.argv[sys.argv.index("--") + 1:]
OUT = argv[0]
O = bpy.data.objects; lum = O["AN_Vagina"]; N = 4 * sum(V.SEGS)
RP = [list(r) for r in lum["ring_planes"]]; SF = list(lum.get("ring_shift_front", [0.0] * len(RP))); NR = len(RP)
picks = [int(x) for x in argv[1].split(",")] if len(argv) > 1 else list(np.linspace(0, NR - 1, 12).astype(int))
def tri(ob):
    me = ob.data; ML = ob.matrix_local; P = np.array([tuple(ML @ v.co) for v in me.vertices]); F = []
    for p in me.polygons:
        vv = list(p.vertices); F += [(vv[0], vv[k], vv[k + 1]) for k in range(1, len(vv) - 1)]
    return P, np.array(F)
L_ = tri(lum); W = tri(O.get("AN_Vagina_Wall") or O["AN_Vagina_Space"]); U_ = tri(O["AN_Uterus"])
P = 300; HW = 0.020; sc = P / (2 * HW); cols = 6; rows = (len(picks) + cols - 1) // cols
img = Image.new("RGB", (cols * P, rows * (P + 18)), "white"); dr = ImageDraw.Draw(img)
for n, j in enumerate(picks):
    c = np.array(RP[j][:3]); t = np.array(RP[j][3:6]); t /= np.linalg.norm(t)
    ex = np.array([1.0, 0, 0]) - t * t[0]; ex /= np.linalg.norm(ex); ey = np.cross(ex, t)
    c = c + ey * SF[j]                              # the ring's centre (shifted to the placeholder's)
    ox, oy = (n % cols) * P + P / 2, (n // cols) * (P + 18) + 18 + P / 2
    to = lambda q: (ox + q[0] * sc, oy - q[1] * sc)
    for (Pm, Fm), col, wd in ((W, (230, 120, 140), 2), (U_, (60, 110, 220), 2), (L_, (0, 0, 0), 1)):
        for a, b in V.plane_section(Pm, Fm, c, t, ex, ey):
            if max(abs(a).max(), abs(b).max()) < HW:
                dr.line([to(a), to(b)], fill=col, width=wd)
    dr.text(((n % cols) * P + 4, (n // cols) * (P + 18) + 2), "ring %d, %.0f mm in (front up)" % (j, RP[j][6]), fill="black")
    dr.line([(ox - 0.01 * sc, oy + P / 2 - 6), (ox + 0.01 * sc, oy + P / 2 - 6)], fill=(150, 150, 150), width=2)
img.save(OUT); print("VAGXS", OUT, NR)
