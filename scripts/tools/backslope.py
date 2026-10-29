"""Metrics for the slope behind the anus (the "cliff" on the coccyx side). Prints one JSON line "BACKSLOPE {...}".
Run: blender -b <file.blend> --python backslope.py -- <out.json>

Profiles are sampled along the cleft (u, + = toward the coccyx / back) at fixed lateral offsets b = -2, 0, +2 mm, by
casting along the outward normal n. b = 0 runs along the up/down creases; b = +-2 mm runs across the pads.
Reported (all in degrees unless noted; a in mm):
  tilt_front / tilt_back : fitted slope of the anus's front half (a -7..-1.5) and back half (a 1.5..7), b = +-2 mean
  max_slope_back         : steepest 1 mm-smoothed slope for a in 3..20 (b = +-2 mean)
  max_kink_back          : largest change in slope between neighbouring 2 mm windows for a in 3..20 (b = +-2 mean)
  profiles               : {b: [[a, h_mm], ...]} for plotting
"""
import bpy, bmesh, sys, json, math
import numpy as np
from mathutils import Vector
from mathutils.bvhtree import BVHTree

OUT = sys.argv[-1]
O = Vector((0.0, 0.0565, 0.790)); n = Vector((0.0, 0.6035, -0.7973)); u = Vector((0.0, 0.7973, 0.6035))
bv = n.cross(u).normalized()
bm = bmesh.new(); bm.from_mesh(bpy.data.objects["Hips"].data); bm.transform(bpy.data.objects["Hips"].matrix_world)
tree = BVHTree.FromBMesh(bm)
prof = {}
for b in (-2.0, 0.0, 2.0):
    rows = []
    for i in range(-300, 501):
        a = i * 0.05
        P = O + u * (a / 1000) + bv * (b / 1000)
        hit = tree.ray_cast(P - n * 0.025, n, 0.06)[0]
        if hit is not None:
            rows.append([a, (hit - O).dot(n) * 1000])
    prof[b] = np.array(rows)


def slope_deg(arr, a0, a1):
    m = (arr[:, 0] >= a0) & (arr[:, 0] <= a1)
    return math.degrees(math.atan(np.polyfit(arr[m, 0], arr[m, 1], 1)[0]))


def local_slopes(arr, a0, a1, win):
    out = []
    a = a0
    while a + win <= a1 + 1e-9:
        out.append((a + win / 2, slope_deg(arr, a, a + win)))
        a += 0.25
    return out


res = {}
side = [prof[-2.0], prof[2.0]]
res["tilt_front"] = round(float(np.mean([slope_deg(p, -7, -1.5) for p in side])), 2)
res["tilt_back"] = round(float(np.mean([slope_deg(p, 1.5, 7) for p in side])), 2)
ms = [max(s for _, s in local_slopes(p, 3, 20, 1.0)) for p in side]
res["max_slope_back"] = round(float(np.mean(ms)), 2)
kinks = []
for p in side:
    ls = local_slopes(p, 3, 20, 2.0)
    by_a = {round(a, 2): s for a, s in ls}
    k = 0.0
    for a, s in ls:
        a2 = round(a + 2.0, 2)
        if a2 in by_a:
            k = max(k, abs(by_a[a2] - s))
    kinks.append(k)
res["max_kink_back"] = round(float(np.mean(kinks)), 2)
res["profiles"] = {str(b): prof[b].round(3).tolist() for b in prof}
json.dump(res, open(OUT, "w"))
print("BACKSLOPE", json.dumps({k: v for k, v in res.items() if k != "profiles"}))
