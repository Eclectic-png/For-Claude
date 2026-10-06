"""Helpers for turning the atlas' separate digestive organ shells into one connected passage.

The atlas organs are closed shells (their open-looking edges are only unwelded shading seams) that overlap where
neighbours meet. open_junction() opens one such meeting: both shells are cut along their intersection curve in a
small patch round the junction (so side contacts elsewhere stay closed), the bit of each shell inside the other is
deleted, and the two cut edges share their positions - the lumens join. Nothing outside the patch moves."""
import bpy, bmesh, math
from mathutils import Vector, kdtree
from mathutils.bvhtree import BVHTree

DIRS = [Vector(d).normalized() for d in ((0.31, 0.52, 0.79), (-0.6, 0.2, -0.77), (0.1, -0.95, 0.3))]


def world_bm(ob, weld=False):
    bm = bmesh.new(); bm.from_mesh(ob.data); bm.transform(ob.matrix_world)
    if weld:
        bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-7)
    bm.verts.ensure_lookup_table(); bm.faces.ensure_lookup_table(); bm.normal_update()
    return bm


def closed_tree(ob):
    """BVH of the organ's shell made watertight (seams welded, any open rim - e.g. a cut-off rectum - capped), for
    inside / outside tests"""
    bm = world_bm(ob, weld=True)
    bnd = [e for e in bm.edges if e.is_boundary]
    if bnd:
        bmesh.ops.holes_fill(bm, edges=bnd, sides=0)
    t = BVHTree.FromBMesh(bm); bm.free(); return t


def inside(tree, p):
    votes = 0
    for d in DIRS:
        c = 0; q = p.copy()
        while True:
            h = tree.ray_cast(q, d, 5.0)
            if h[0] is None:
                break
            c += 1; q = h[0] + d * 1e-7
        votes += c % 2
    return votes >= 2


def overlap_components(bm, tree):
    ins = {f.index for f in bm.faces if inside(tree, f.calc_center_median())}
    comps = []; seen = set()
    for fi in ins:
        if fi in seen:
            continue
        st = [fi]; comp = []
        while st:
            k = st.pop()
            if k in seen:
                continue
            seen.add(k); comp.append(k)
            for e in bm.faces[k].edges:
                st += [g.index for g in e.link_faces if g.index in ins and g.index not in seen]
        comps.append(comp)
    return sorted(comps, key=lambda c: -sum(bm.faces[k].calc_area() for k in c))


def comp_sphere(bm, comp):
    vs = {v for k in comp for v in bm.faces[k].verts}
    c = sum((v.co for v in vs), Vector()) / len(vs)
    return c, max((v.co - c).length for v in vs)


def _ctx_edit(ob):
    for o in bpy.context.view_layer.objects:
        o.select_set(False)
    bpy.context.view_layer.objects.active = ob; ob.select_set(True)
    bpy.ops.object.mode_set(mode='EDIT')


def _cut_patches(bmA, bmB, tA, tB, cA, R, eps=0.0, eps_side=1, sides=(0, 1)):
    """intersect the patches of A and B within R of cA; returns (tb, src layer, kill, closed curve?, curve length).
    sides: which sides lose the piece inside the other organ ((0,) = only A: B passes through A's wall)"""
    patchA = [f for f in bmA.faces if (f.calc_center_median() - cA).length < R]
    patchB = [f for f in bmB.faces if (f.calc_center_median() - cA).length < R]
    # temporary object: both patches, tagged by source; cut along the intersection with shared vertices
    tb = bmesh.new(); src = tb.faces.layers.int.new("src"); vmap = {}
    for tag, bm_, patch in ((0, bmA, patchA), (1, bmB, patchB)):
        for f in patch:
            vs = []
            for v in f.verts:
                key = (tag, v.index)
                if key not in vmap:
                    vmap[key] = tb.verts.new(v.co)
                vs.append(vmap[key])
            if len(set(vs)) < 3:
                continue
            try:
                nf = tb.faces.new(vs); nf[src] = tag
            except ValueError:
                pass
    for tag in (0, 1):                             # weld each side's shading seams inside the patch
        vs_ = {v for f in tb.faces if f[src] == tag for v in f.verts}
        bmesh.ops.remove_doubles(tb, verts=list(vs_), dist=1e-7)
    if eps != 0:
        # neighbours split from one source mesh share exactly coincident edges or walls along the split, which the
        # cut cannot resolve (gaps in the curve): one side's patch is inflated (eps > 0) or deflated (eps < 0) by
        # |eps|, fading to nothing at the patch border
        tb.normal_update()
        sv = {v for f in tb.faces if f[src] == eps_side for v in f.verts}
        bv = [v for v in sv if any(len([g for g in e.link_faces if g[src] == eps_side]) == 1 for e in v.link_edges)]
        kb = kdtree.KDTree(max(len(bv), 1))
        for i, v in enumerate(bv):
            kb.insert(v.co, i)
        kb.balance()
        moves = {}
        for v in sv:                               # fades out towards the patch border, which must not move
            db = kb.find(v.co)[2] if bv else 1.0
            w = min(db / (0.3 * R), 1.0)
            moves[v] = v.normal * (eps * w * w * (3 - 2 * w))
        for v, mv in moves.items():
            v.co += mv
    me = bpy.data.meshes.new("_junction"); tb.to_mesh(me); tb.free()
    T = bpy.data.objects.new("_junction", me); bpy.context.scene.collection.objects.link(T)
    _ctx_edit(T)
    # A's faces selected, B's not: cut A against B only (each organ's own coils touching are left alone)
    bpy.ops.object.mode_set(mode='OBJECT')
    _s = [0] * len(T.data.polygons); T.data.attributes["src"].data.foreach_get("value", _s)
    for poly, tag in zip(T.data.polygons, _s):
        poly.select = tag == 0
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.intersect(mode='SELECT_UNSELECT', separate_mode='NONE', solver='EXACT')
    bpy.ops.object.mode_set(mode='OBJECT')
    tb = bmesh.new(); tb.from_mesh(T.data); src = tb.faces.layers.int["src"]; tb.faces.ensure_lookup_table()
    bpy.data.objects.remove(T); bpy.data.meshes.remove(me)
    # the junction = the longest intersection curve C (edges shared by both sides). On each side, the piece of the
    # patch that C encloses and that does not reach the patch's outer border is what lies inside the other organ
    # (the cap pushed into it / the disc of wall it pierces): that piece goes. Same curve on both sides, so the
    # openings match exactly; smaller contacts caught in the patch stay closed.
    tb.edges.ensure_lookup_table()
    cut = {e for e in tb.edges if len({f[src] for f in e.link_faces}) == 2}
    comps = []; seen = set()
    for e0 in cut:
        if e0 in seen:
            continue
        st = [e0]; comp = []
        while st:
            e = st.pop()
            if e in seen:
                continue
            seen.add(e); comp.append(e)
            for v in e.verts:
                st += [g for g in v.link_edges if g in cut and g not in seen]
        comps.append(comp)
    if not comps:
        return tb, src, [], False, 0.0, (patchA, patchB), []
    # the cut can miss short stretches (0.1-3 mm) where the two shells graze each other or share an edge from the
    # source mesh they were split from: loose ends within GAP_MAX are paired, and on each side the cut follows the
    # shortest edge path between them; the sliver between the two sides' paths is filled afterwards (stitch)
    GAP_MAX = 0.003
    deg = {}
    for e in cut:
        for v in e.verts:
            deg[v] = deg.get(v, 0) + 1
    main = max(comps, key=lambda c: sum(e.calc_length() for e in c))
    comp_of = {}
    for i, c in enumerate(comps):
        for e in c:
            for v in e.verts:
                comp_of[v] = i
    ends = [v for v, d_ in deg.items() if d_ == 1]
    pairs = []; used = set()
    for u, v in sorted(((u, v) for i, u in enumerate(ends) for v in ends[i + 1:]), key=lambda uv: (uv[0].co - uv[1].co).length):
        if u in used or v in used or (u.co - v.co).length > GAP_MAX:
            continue
        pairs.append((u, v)); used |= {u, v}
    adj = {}
    for i, c in enumerate(comps):
        adj.setdefault(i, set())
    for u, v in pairs:
        adj[comp_of[u]].add(comp_of[v]); adj[comp_of[v]].add(comp_of[u])
    mi = comps.index(main); grp = {mi}; st = [mi]
    while st:
        k = st.pop()
        for j in adj[k]:
            if j not in grp:
                grp.add(j); st.append(j)
    C = {e for i in grp for e in comps[i]}
    pairs = [(u, v) for u, v in pairs if comp_of[u] in grp]
    pd = {}
    for u, v in pairs:
        pd[u] = pd.get(u, 0) + 1; pd[v] = pd.get(v, 0) + 1
    closed = all(deg[v] + pd.get(v, 0) == 2 for e in C for v in e.verts)

    def side_path(u, v, tag):
        import heapq
        dist = {u: 0.0}; prev = {}; h = [(0.0, id(u), u)]
        while h:
            d_, _, x = heapq.heappop(h)
            if x is v:
                break
            if d_ > dist.get(x, 1e9) or d_ > 4 * GAP_MAX:
                continue
            for e in x.link_edges:
                if not any(f[src] == tag for f in e.link_faces):
                    continue
                y = e.other_vert(x); nd = d_ + e.calc_length()
                if nd < dist.get(y, 1e9):
                    dist[y] = nd; prev[y] = (x, e); heapq.heappush(h, (nd, id(y), y))
        if v not in prev:
            return None, None
        vs = [v]; es = []
        while vs[-1] is not u:
            x, e = prev[vs[-1]]; vs.append(x); es.append(e)
        return vs[::-1], set(es)
    paths = {0: [], 1: []}; barrier = {0: set(C), 1: set(C)}
    for u, v in pairs:
        for tag in (0, 1):
            pv, pe = side_path(u, v, tag)
            if pv is None:
                closed = False; continue
            paths[tag].append(pv); barrier[tag] |= pe
    kill = []
    for tag in sides:
        side = [f for f in tb.faces if f[src] == tag]; seen = set(); pieces = []
        for f0 in side:
            if f0 in seen:
                continue
            st = [f0]; comp = []; border = False; onC = False
            while st:
                f = st.pop()
                if f in seen:
                    continue
                seen.add(f); comp.append(f)
                for e in f.edges:
                    same = [g for g in e.link_faces if g[src] == tag]
                    if e in barrier[tag]:
                        onC = True; continue
                    if len(same) == 1 and e not in cut:
                        border = True
                    st += [g for g in same if g not in seen]
            pieces.append((comp, border, onC))
        # among the pieces bounded by C, the one lying inside the other organ (area-weighted centre test; a piece
        # reaching the patch border is outside by construction unless the whole organ fits in the patch)
        tr = tB if tag == 0 else tA; best = None
        for comp, border, onC in pieces:
            if not onC:
                continue
            a_tot = sum(f.calc_area() for f in comp)
            frac = sum(f.calc_area() for f in comp if inside(tr, f.calc_center_median())) / max(a_tot, 1e-18)
            score = frac - (0.5 if border else 0.0)
            if frac > 0.5 and (best is None or score > best[0]):
                best = (score, comp)
        if best:
            kill.append(best[1])
    stitch = list(zip(paths[0], paths[1]))
    ok = closed and len(kill) == len(sides)
    return tb, src, [f for k in kill for f in k], ok, sum(e.calc_length() for e in C), (patchA, patchB), stitch


def open_junction(A, B, centre=None, margins=(0.012, 0.024, 0.04), report=None, trees=None, sides=(0, 1)):
    """Open the A/B junction. centre: world point inside the junction zone (default: the biggest zone of A inside B).
    The patch round it grows through `margins` until the junction's intersection curve closes inside it.
    trees: {name: BVHTree} of the organs' CLOSED shells taken before any junction was opened (an opened shell
    would break the inside tests); built here if not given.
    sides=(0,) pierces instead: only A's wall is opened where B passes through it, B stays whole (a ureter's renal
    pelvis leaving the kidney at the hilum)"""
    tA = trees[A.name] if trees and A.name in trees else closed_tree(A)
    tB = trees[B.name] if trees and B.name in trees else closed_tree(B)
    bmA, bmB = world_bm(A), world_bm(B)
    compsA = overlap_components(bmA, tB)
    if not compsA:
        raise RuntimeError(f"{A.name} and {B.name} do not overlap")
    if centre is None:
        best = compsA[0]
    else:
        best = min(compsA, key=lambda c: (comp_sphere(bmA, c)[0] - centre).length)
    cA, rA = comp_sphere(bmA, best)
    compsB = overlap_components(bmB, tA) if 1 in sides else []
    if compsB:                                     # B's part inside A, near the same place, sizes the patch too
        cb = min(compsB, key=lambda c: (comp_sphere(bmB, c)[0] - cA).length)
        rA = max(rA, max((v.co - cA).length for k in cb for v in bmB.faces[k].verts))
    tries = [(m_, e_, s_) for e_, s_ in ((0.0, 1), (5e-5, 1), (5e-5, 0), (-5e-5, 1), (-5e-5, 0), (2e-4, 1), (2e-4, 0),
                                          (5e-4, 1), (-5e-4, 0), (1e-3, 1), (-1e-3, 0)) for m_ in margins]
    for margin, eps, eps_side in tries:
        R = rA + margin
        tb, src, kill, ok, clen, (patchA, patchB), stitch = _cut_patches(bmA, bmB, tA, tB, cA, R, eps, eps_side, sides)
        if ok:
            break
        tb.free()
    else:
        raise RuntimeError(f"{A.name}/{B.name}: no closed junction curve up to {R * 1000:.0f} mm")
    bmesh.ops.delete(tb, geom=kill, context='FACES')
    # close the slivers where the two sides' cuts took different paths across a gap: A's opening there follows
    # B's path, so both openings are the same loop
    n_st = 0
    for p0, p1 in stitch:
        poly = p0 + p1[::-1][1:-1]
        if len(set(poly)) < 3 or any(not v.is_valid for v in poly):
            continue
        for k in range(1, len(poly) - 1):
            try:
                f = tb.faces.new((poly[0], poly[k], poly[k + 1])); f[src] = 0; n_st += 1
            except ValueError:
                pass

    # put each side back into its organ, welded to the untouched rest at the patch border
    out = {}
    for tag, ob, bm_, patch in ((0, A, bmA, patchA), (1, B, bmB, patchB)):
        ao = bm_.verts.layers.float_color.get("AO")
        if ao is not None:                         # the baked AO of the new vertices: from the nearest old one
            kd = kdtree.KDTree(len(bm_.verts))
            for v in bm_.verts:
                kd.insert(v.co, v.index)
            kd.balance(); ao_old = [tuple(v[ao]) for v in bm_.verts]
        bmesh.ops.delete(bm_, geom=patch, context='FACES_ONLY')
        loose = [v for v in bm_.verts if not v.link_faces]
        bmesh.ops.delete(bm_, geom=loose, context='VERTS')
        border = [v for v in bm_.verts if v.is_boundary]
        nv = {}
        for f in tb.faces:
            if f[src] != tag:
                continue
            vs = []
            for v in f.verts:
                if v not in nv:
                    nv[v] = bm_.verts.new(v.co)
                vs.append(nv[v])
            try:
                bm_.faces.new(vs)
            except ValueError:
                pass
        if ao is not None:
            for v in nv.values():
                v[ao] = ao_old[kd.find(v.co)[1]]
        bmesh.ops.remove_doubles(bm_, verts=border + list(nv.values()), dist=1e-7)
        bm_.normal_update()
        bm_.transform(ob.matrix_world.inverted()); bm_.to_mesh(ob.data); ob.data.update()
        out[ob.name] = {"patch_faces": len(patch)}
    tb.free()
    res = {"pair": [A.name, B.name], "curve_mm": round(clen * 1000, 1), "inflate_mm": eps * 1000, "inflated": [A.name, B.name][eps_side] if eps else None,
           "gaps_stitched": len(stitch), "stitch_faces": n_st, "centre": [round(x, 4) for x in cA], "patch_radius_mm": round(R * 1000, 1),
           "killed_faces": len(kill), **out}
    if report is not None:
        report.append(res)
    return res


def opening_loops(A, B, tol=2e-6):
    """edge loops of A lying on B's surface after open_junction (= the openings), with their size"""
    bm = world_bm(A); tB = BVHTree.FromBMesh(world_bm(B, weld=True))
    E = [e for e in bm.edges if e.is_boundary and all(tB.find_nearest(v.co)[3] < tol for v in e.verts)]
    seen = set(); loops = []
    for e in E:
        if e in seen:
            continue
        st = [e]; comp = []
        while st:
            x = st.pop()
            if x in seen:
                continue
            seen.add(x); comp.append(x)
            for v in x.verts:
                st += [y for y in v.link_edges if y in E and y not in seen]
        vs = {v.co.copy().freeze() for x in comp for v in x.verts}
        c = sum(vs, Vector()) / len(vs)
        loops.append({"edges": len(comp), "perimeter_mm": round(sum(x.calc_length() for x in comp) * 1000, 1),
                      "mean_radius_mm": round(sum((p - c).length for p in vs) / len(vs) * 1000, 1)})
    bm.free()
    return loops


# ---- lighting (shared with build5_passage) ----
def light_insides(mat):
    """atlas toon shader: light backfaces (the inside) with the flipped normal and the "AO_in" bake instead of the
    one flat backface colour. Works on the three atlas shader families (organ, skeleton / cartilage, lung): in each
    the last Mix before the Emission picks the flat colour (B) on backfaces (Factor = Geometry.Backfacing); its lit
    colour (A) is wired straight to the Emission instead. Returns False if the material is not one of them."""
    nt = mat.node_tree if mat else None
    if nt is None or mat.get("insides_lit"):
        return False
    N = nt.nodes; Lk = nt.links
    if not all(n_ in N for n_ in ("Geometry", "Vector Transform", "Color Attribute", "Emission")):
        return False
    em_in = N["Emission"].inputs["Color"]
    bf = em_in.links[0].from_node if em_in.is_linked else None
    if bf is None or bf.bl_idname != "ShaderNodeMix" or not bf.inputs[0].is_linked or \
            bf.inputs[0].links[0].from_socket.name != "Backfacing" or not bf.inputs[6].is_linked:
        return False
    lit_src = bf.inputs[6].links[0].from_socket
    geo = N["Geometry"]; vt = N["Vector Transform"]
    sg = N.new("ShaderNodeMath"); sg.operation = "MULTIPLY_ADD"; sg.inputs[1].default_value = -2.0
    sg.inputs[2].default_value = 1.0; Lk.new(geo.outputs["Backfacing"], sg.inputs[0])
    fl = N.new("ShaderNodeVectorMath"); fl.operation = "SCALE"
    Lk.new(geo.outputs["Normal"], fl.inputs[0]); Lk.new(sg.outputs[0], fl.inputs["Scale"])
    Lk.new(fl.outputs["Vector"], vt.inputs["Vector"])
    ao = N["Color Attribute"]; ai = N.new("ShaderNodeVertexColor"); ai.layer_name = "AO_in"
    mx = N.new("ShaderNodeMix"); mx.data_type = "RGBA"
    Lk.new(geo.outputs["Backfacing"], mx.inputs["Factor"]); Lk.new(ao.outputs["Color"], mx.inputs[6])
    Lk.new(ai.outputs["Color"], mx.inputs[7])
    tgt = [l_.to_socket for l_ in ao.outputs["Color"].links if l_.to_node is not mx]
    for t_ in tgt:
        Lk.new(mx.outputs[2], t_)
    Lk.new(lit_src, em_in)                         # the flat backface colour is bypassed
    mat["insides_lit"] = True
    return True


def light_all(objs, tree=None, report=None, skip_baked=False):
    """light_insides on every material of `objs` and an "AO_in" bake for each object whose material was lit (now or
    earlier; skip_baked: not for objects that already carry an AO_in bake). tree: scene BVH for the bake (built if
    None). Returns {name: mean AO_in}."""
    out = {}; lit = []
    for ob in objs:
        if ob.type != 'MESH':
            continue
        mats = [m_ for m_ in ob.data.materials if m_]
        for m_ in mats:
            light_insides(m_)
        if mats and all(m_.get("insides_lit") for m_ in mats):
            if not (skip_baked and ob.data.color_attributes.get("AO_in")):
                lit.append(ob)
    tree = tree or scene_bvh()
    for ob in lit:
        out[ob.name] = bake_ao(ob, "AO_in", True, tree)
    if report is not None:
        report.update(out)
    return out


def fib_hemi(nr):
    out = []
    for i in range(nr):                            # cosine-weighted Fibonacci hemisphere (z up)
        u = (i + 0.5) / nr; r_ = math.sqrt(u); ph = i * math.pi * (3 - math.sqrt(5))
        out.append(Vector((r_ * math.cos(ph), r_ * math.sin(ph), math.sqrt(1 - u))))
    return out


def set_open(v_):
    for o_ in bpy.data.objects:
        if o_.type == 'MESH' and o_.data.shape_keys and "Open" in o_.data.shape_keys.key_blocks:
            o_.data.shape_keys.key_blocks["Open"].value = v_
    bpy.context.view_layer.update()


def scene_bvh():
    dg = bpy.context.evaluated_depsgraph_get(); bb = bmesh.new()
    for o_ in bpy.data.objects:
        if o_.type == 'MESH' and not o_.hide_render:
            t_ = bmesh.new()
            try:
                t_.from_object(o_, dg)
            except ValueError:                     # not in the view layer / no evaluated mesh
                t_.free(); continue
            t_.transform(o_.matrix_world)
            me_ = bpy.data.meshes.new("_tmp"); t_.to_mesh(me_); t_.free(); bb.from_mesh(me_); bpy.data.meshes.remove(me_)
    tr = BVHTree.FromBMesh(bb); bb.free(); return tr


def bake_ao(ob, name, inner, tree, dist=0.008, rays=32):
    """per-vertex ambient occlusion ("smooth lighting") into colour attribute `name`: the share of a cosine-weighted
    hemisphere (around the normal, or the inverted normal for the inside) clear of `tree` within `dist`"""
    hemi = fib_hemi(rays)
    dg = bpy.context.evaluated_depsgraph_get(); t_ = bmesh.new(); t_.from_object(ob, dg)
    t_.transform(ob.matrix_world); t_.normal_update(); vals = []
    for v in t_.verts:
        nrm = -v.normal if inner else v.normal
        q = nrm.to_track_quat('Z', 'Y'); o_ = v.co + nrm * 2e-5; acc = 0
        for dd in hemi:
            if tree.ray_cast(o_, q @ dd, dist)[0] is None:
                acc += 1
        vals.append(acc / rays)
    t_.free()
    ca = ob.data.color_attributes.get(name) or ob.data.color_attributes.new(name, 'FLOAT_COLOR', 'POINT')
    for i, a_ in enumerate(vals):
        ca.data[i].color = (a_, a_, a_, 1.0)
    return round(sum(vals) / max(len(vals), 1), 3)
