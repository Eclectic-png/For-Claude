"""Internal female reproductive organs for build_pelvis: shapes in atlas space (real-size metres).

Metaballs make the smooth organic bodies (uterus, its lumen, ovaries); Blender clamps a metaball's resolution to
>= 5 mm, so they are built in millimetre units and scaled down. Tubes and ligaments are swept along Catmull-Rom
centrelines like the urethra."""
import bpy, bmesh, math
import numpy as np
from mathutils import Vector, Matrix, Quaternion


def _quat_from_axes(ax_x, ax_y):
    """rotation taking local x -> ax_x, local y -> ax_y (z = x cross y)"""
    x = Vector(ax_x).normalized(); y = (Vector(ax_y) - x * x.dot(Vector(ax_y))).normalized(); z = x.cross(y)
    return Matrix((x, y, z)).transposed().to_quaternion()


def metaball_mesh(name, elements, resolution_mm=0.7, threshold=0.6):
    """elements: dicts {type, co (m), radius (m), size (x, y, z: m, for ELLIPSOID / half-length x for CAPSULE),
    ax_x, ax_y (orientation, optional)}. Returns a new mesh (metres) of their blended surface."""
    mb = bpy.data.metaballs.new("_mb_" + name); mb.resolution = resolution_mm; mb.render_resolution = resolution_mm
    mb.threshold = threshold
    for e in elements:
        el = mb.elements.new(type=e["type"]); el.co = tuple(Vector(e["co"]) * 1000); el.radius = e["radius"] * 1000
        if e["type"] == 'CAPSULE':
            el.size_x = e["half_len"] * 1000          # a length (object units: mm)
        elif "size" in e:
            el.size_x, el.size_y, el.size_z = e["size"]   # unitless factors
        if "ax_x" in e:
            el.rotation = _quat_from_axes(e["ax_x"], e.get("ax_y", (0, 0, 1)))
    ob = bpy.data.objects.new("_mb_" + name, mb); bpy.context.scene.collection.objects.link(ob)
    bpy.context.view_layer.update()
    dg = bpy.context.evaluated_depsgraph_get(); ev = ob.evaluated_get(dg); m = ev.to_mesh()
    bm = bmesh.new(); bm.from_mesh(m); ev.to_mesh_clear()
    bm.transform(Matrix.Scale(0.001, 4))
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    bm.normal_update()
    c = sum((v.co for v in bm.verts), Vector()) / len(bm.verts)
    if sum((f.calc_center_median() - c).dot(f.normal) for f in bm.faces) < 0:
        bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
    me = bpy.data.meshes.new(name); bm.to_mesh(me); bm.free()
    bpy.data.objects.remove(ob); bpy.data.metaballs.remove(mb)
    return me


def decimate(me, ratio):
    """the mesh collapsed to `ratio` of its faces (the metaball surface is far denser than the shape needs; every
    later step - room making, checks, bakes - walks its vertices one by one)"""
    name = me.name; ob = bpy.data.objects.new("_dec", me); bpy.context.scene.collection.objects.link(ob)
    md = ob.modifiers.new("dec", 'DECIMATE'); md.decimate_type = 'COLLAPSE'; md.ratio = ratio
    md.use_collapse_triangulate = True
    bpy.context.view_layer.update()
    m = bpy.data.meshes.new_from_object(ob.evaluated_get(bpy.context.evaluated_depsgraph_get()))
    bpy.data.objects.remove(ob); bpy.data.meshes.remove(me); m.name = name
    return m


def catmull(pts, n):
    P = [pts[0]] + list(pts) + [pts[-1]]; out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=(i == len(P) - 3)):
            out.append(0.5 * ((2 * p1) + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    return out


def resample(pts, step):
    """equal arc-length points along a polyline"""
    cum = [0.0]
    for a, b in zip(pts, pts[1:]):
        cum.append(cum[-1] + (b - a).length)
    L = cum[-1]; n = max(2, int(round(L / step)) + 1); out = []
    for k in range(n):
        t = L * k / (n - 1); i = max(0, min(len(cum) - 2, int(np.searchsorted(cum, t)) - 1))
        f = 0.0 if cum[i + 1] == cum[i] else (t - cum[i]) / (cum[i + 1] - cum[i])
        out.append(pts[i].lerp(pts[i + 1], f))
    return out, L


def sweep(name, centre, radius_fn, segs=20, cap_start=True, cap_end=True, rim_fn=None, up=Vector((0, 0, 1))):
    """tube along `centre` (Vectors), radius_fn(s) for s in 0..1 (arc length); rim_fn(theta, s) -> extra radial
    factor (fringes); caps close the ends (a fan), an uncapped end stays open (a free boundary)."""
    n = len(centre); bm = bmesh.new(); rings = []
    tang = [(centre[min(i + 1, n - 1)] - centre[max(i - 1, 0)]).normalized() for i in range(n)]
    e1 = up.cross(tang[0]).normalized() if up.cross(tang[0]).length > 1e-6 else Vector((1, 0, 0)).cross(tang[0]).normalized()
    for i in range(n):
        e1 = (e1 - tang[i] * e1.dot(tang[i])).normalized(); e2 = tang[i].cross(e1)
        s = i / (n - 1); r = radius_fn(s); ring = []
        for k in range(segs):
            th = 2 * math.pi * k / segs; f = rim_fn(th, s) if rim_fn else 1.0
            ring.append(bm.verts.new(centre[i] + (e1 * math.cos(th) + e2 * math.sin(th)) * r * f))
        rings.append(ring)
    for a, b in zip(rings, rings[1:]):
        for k in range(segs):
            bm.faces.new((a[k], a[(k + 1) % segs], b[(k + 1) % segs], b[k]))
    if cap_start:
        c = bm.verts.new(centre[0] - tang[0] * radius_fn(0.0) * 0.3)
        for k in range(segs):
            bm.faces.new((rings[0][(k + 1) % segs], rings[0][k], c))
    if cap_end:
        c = bm.verts.new(centre[-1] + tang[-1] * radius_fn(1.0) * 0.3)
        for k in range(segs):
            bm.faces.new((rings[-1][k], rings[-1][(k + 1) % segs], c))
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    bm.normal_update()
    # outward: faces point away from their ring's centre
    votes = 0
    for f in bm.faces[:200]:
        p = f.calc_center_median(); i = min(range(n), key=lambda j: (centre[j] - p).length)
        votes += 1 if (p - centre[i]).dot(f.normal) > 0 else -1
    if votes < 0:
        bmesh.ops.reverse_faces(bm, faces=bm.faces[:])
    me = bpy.data.meshes.new(name); bm.to_mesh(me); bm.free()
    return me


def ellipsoid_mesh(name, centre, axes, radii, subdiv=4, bumps=0.0, seed=1):
    """an ellipsoid (axes: three orthonormal Vectors, radii in m), optional small surface bumps (follicles)"""
    bm = bmesh.new(); bmesh.ops.create_icosphere(bm, subdivisions=subdiv, radius=1.0)
    rng = np.random.default_rng(seed)
    centres = [Vector(rng.normal(size=3)).normalized() for _ in range(9)]
    for v in bm.verts:
        d = v.co.normalized(); f = 1.0
        for c in centres:
            f += bumps * math.exp(-((d - c).length / 0.25) ** 2)
        v.co = centre + (axes[0] * (d.x * radii[0]) + axes[1] * (d.y * radii[1]) + axes[2] * (d.z * radii[2])) * f
    bm.normal_update()
    me = bpy.data.meshes.new(name); bm.to_mesh(me); bm.free()
    return me


def tinted(src, name, rgb):
    """copy of an atlas toon material with its colours moved to the hue / saturation of rgb (values kept)"""
    import colorsys
    m = src.copy(); m.name = name
    h, s, _ = colorsys.rgb_to_hsv(*rgb)
    for n in m.node_tree.nodes:
        for sock in n.inputs:
            if sock.type == 'RGBA' and not sock.is_linked:
                r, g, b, a = sock.default_value
                _, s0, v0 = colorsys.rgb_to_hsv(r, g, b)
                sock.default_value = (*colorsys.hsv_to_rgb(h, s if s0 > 0.05 else s0, v0), a)
    m["insides_lit"] = src.get("insides_lit", False)
    return m


def ellipsoid_el(c, ax_x, ax_y, half):
    """a metaball ellipsoid with the given half extents (m) along ax_x, ax_y, ax_x x ax_y (calibrated: an isolated
    element's surface sits at 0.574 x radius x size at threshold 0.6)"""
    m = max(half); return {"type": 'ELLIPSOID', "co": c, "radius": m / 0.574, "size": tuple(h / m for h in half),
                           "ax_x": ax_x, "ax_y": ax_y}


def capsule_el(a, b, r):
    """a metaball capsule whose surface runs from a to b (end caps included), radius r (m)"""
    a = Vector(a); b = Vector(b); L = (b - a).length
    return {"type": 'CAPSULE', "co": (a + b) / 2, "radius": r / 0.574, "half_len": max(L / 2 - r, 1e-5),
            "ax_x": (b - a).normalized(), "ax_y": (b - a).normalized().orthogonal()}
