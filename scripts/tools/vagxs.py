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
X = np.array([tuple(v.co) for v in lum.data.vertices]); NR = (len(X) - 1) // N
picks = [int(x) for x in argv[1].split(",")] if len(argv) > 1 else list(np.linspace(0, NR - 1, 12).astype(int))
cent = np.array([X[j * N:(j + 1) * N].mean(0) for j in range(NR)])
def tri(ob):
    me = ob.data; P = np.array([tuple(v.co) for v in me.vertices]); F = []
    for p in me.polygons:
        vv = list(p.vertices); F += [(vv[0], vv[k], vv[k + 1]) for k in range(1, len(vv) - 1)]
    return P, np.array(F)
W = tri(O["AN_Vagina_Wall"]); U_ = tri(O["AN_Uterus"])
cum = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(cent, axis=0), axis=1))])
P = 300; HW = 0.020; sc = P / (2 * HW); cols = 6; rows = (len(picks) + cols - 1) // cols
img = Image.new("RGB", (cols * P, rows * (P + 18)), "white"); dr = ImageDraw.Draw(img)
for n, j in enumerate(picks):
    t = cent[min(j + 1, NR - 1)] - cent[max(j - 1, 0)]; t /= np.linalg.norm(t)
    ex = np.array([1.0, 0, 0]) - t * t[0]; ex /= np.linalg.norm(ex); ey = np.cross(ex, t)
    ox, oy = (n % cols) * P + P / 2, (n // cols) * (P + 18) + 18 + P / 2
    to = lambda q: (ox + q[0] * sc, oy - q[1] * sc)
    for (Pm, Fm), col in ((W, (230, 120, 140)), (U_, (60, 110, 220))):
        for a, b in V.plane_section(Pm, Fm, cent[j], t, ex, ey):
            if max(abs(a).max(), abs(b).max()) < HW:
                dr.line([to(a), to(b)], fill=col, width=2)
    R = [((p - cent[j]) @ ex, (p - cent[j]) @ ey) for p in X[j * N:(j + 1) * N]]
    dr.line([to(q) for q in R] + [to(R[0])], fill="black", width=1)
    dr.text(((n % cols) * P + 4, (n // cols) * (P + 18) + 2), "ring %d, %.0f mm in (front up)" % (j, cum[j] * 1000), fill="black")
    dr.line([(ox - 0.01 * sc, oy + P / 2 - 6), (ox + 0.01 * sc, oy + P / 2 - 6)], fill=(150, 150, 150), width=2)
img.save(OUT); print("VAGXS", OUT, NR)
