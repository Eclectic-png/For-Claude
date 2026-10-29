"""Renders of the urinary build: overview, sagittal cut-aways at Fill 0 / 0.5 / 1, the meatus closed / voiding, and an
EEVEE (toon materials) cut-away of the bladder showing its lit inside. Args: -- <outdir>"""
import bpy, sys, os, math
from mathutils import Vector
OUT = sys.argv[-1]; os.makedirs(OUT, exist_ok=True)
sc = bpy.context.scene; O = bpy.data.objects; ctl = O["Pelvic_Controls"]
COL = {"Hips": (0.85, 0.75, 0.72), "AN_Bladder": (0.95, 0.8, 0.3), "AN_Urethra": (0.95, 0.5, 0.55),
       "AN_Urethra_Wall": (0.75, 0.35, 0.35), "AN_Ureter_L": (0.9, 0.55, 0.1), "AN_Ureter_R": (0.9, 0.55, 0.1),
       "AN_Kidney_L": (0.6, 0.2, 0.2), "AN_Kidney_R": (0.6, 0.2, 0.2), "AN_Adrenal_L": (0.8, 0.6, 0.3),
       "AN_Adrenal_R": (0.8, 0.6, 0.3), "AN_HipBone_L": (0.93, 0.93, 0.9), "AN_HipBone_R": (0.93, 0.93, 0.9),
       "AN_Sacrum": (0.93, 0.93, 0.9), "AN_Rectum": (0.55, 0.65, 0.45), "AN_Colon_Descending": (0.45, 0.6, 0.4),
       "AN_Rectum_LowerAmpulla": (0.55, 0.65, 0.45), "AN_AnalCanal": (0.6, 0.4, 0.4), "AN_AnalSphincter": (0.7, 0.3, 0.3),
       "AN_Uterus": (0.80, 0.36, 0.42), "AN_Uterus_Lumen": (0.98, 0.70, 0.75), "AN_Ovary_L": (0.95, 0.85, 0.70),
       "AN_Ovary_R": (0.95, 0.85, 0.70), "AN_Uterine_Tube_L": (0.90, 0.50, 0.60), "AN_Uterine_Tube_R": (0.90, 0.50, 0.60)}
for s_ in "LR":
    for nm_ in ("Ovarian", "Round", "Suspensory"):
        COL[f"AN_{nm_}_Ligament_{s_}"] = (0.75, 0.72, 0.62)
REP = [n for n in COL if n.startswith(("AN_Uter", "AN_Ov", "AN_Round", "AN_Susp"))]
for o in O:
    if o.type == 'MESH' and o.name in COL:
        o.color = (*COL[o.name], 1)
cam = O.new("_cam", bpy.data.cameras.new("_cam")); cam.data.type = 'ORTHO'; sc.collection.objects.link(cam); sc.camera = cam


def setc(fill, void, rect=0.0):
    ctl["Bladder_Fill"] = fill; ctl["Void"] = void; ctl["Rectum_Fill"] = rect; ctl.update_tag(); bpy.context.view_layer.update()


def show(names):
    for o in O:
        if o.type == 'MESH':
            o.hide_render = o.name not in names


def shot(name, target, d, scale, clip=0.001, res=(1100, 1100), engine="BLENDER_WORKBENCH"):
    sc.render.engine = engine
    if engine == "BLENDER_WORKBENCH":
        sh = sc.display.shading; sh.light = 'STUDIO'; sh.color_type = 'OBJECT'; sh.show_backface_culling = False
        sh.show_cavity = False
    d = Vector(d).normalized(); cam.location = Vector(target) - d * 0.4; cam.data.ortho_scale = scale
    cam.data.clip_start = clip; cam.data.clip_end = 2.0
    cam.rotation_euler = d.to_track_quat('-Z', 'Y').to_euler()
    sc.render.resolution_x, sc.render.resolution_y = res
    sc.render.filepath = f"{OUT}/{name}.png"; bpy.ops.render.render(write_still=True)


org = [n for n in COL if n not in ("Hips",)]
setc(0, 0); show(org)
shot("overview_front", (0, 0.0, 0.88), (0, 1, 0), 0.30)
shot("overview_side", (0, 0.02, 0.86), (-1, 0, 0), 0.24)
# midline cut-aways (camera on +x, clipped at x = 0): bladder filling, bowel giving way
for f in (0.0, 0.5, 1.0):
    setc(f, 0); show(org + ["Hips"])
    shot(f"cut_fill{f}", (0, 0.01, 0.83), (-1, 0, 0), 0.13, clip=0.4)
# rectum full, and both full
for f, r, nm in ((0.0, 1.0, "cut_rectum_full"), (1.0, 1.0, "cut_both_full")):
    setc(f, 0, r); show(org + ["Hips"])
    shot(nm, (0, 0.01, 0.83), (-1, 0, 0), 0.13, clip=0.4)
# meatus from below, closed / voiding (labia cut away by a clip plane just under the crest)
for v in (0.0, 1.0):
    setc(0, v); show(["Hips", "AN_Urethra"])
    shot(f"meatus_void{v}", (0, 0.0135, 0.7835), (0, 0, 1), 0.012, clip=0.4 - 0.0008)
    shot(f"meatus_side_void{v}", (0, 0.0135, 0.792), (-1, 0, 0), 0.026, clip=0.4)
# EEVEE, real materials: bladder (full) cut open from the front, inside lit
setc(1.0, 1.0); show(org)
sc.eevee.taa_render_samples = 16
shot("eevee_cut_inside", (0, 0.0, 0.84), (0, 1, 0.2), 0.12, clip=0.4 - 0.01, engine="BLENDER_EEVEE_NEXT")
setc(0, 0); show(org)
shot("eevee_front", (0, 0.0, 0.88), (0, 1, 0), 0.30, engine="BLENDER_EEVEE_NEXT")
# the reproductive organs: front / top / side with the bladder and the bony pelvis, the vagina placeholder (pink);
# sagittal cut-aways (uterus cut open: canal + cavity) at rest, bladder full, rectum full, both full
if "AN_Vagina_Space" in O:
    O["AN_Vagina_Space"].color = (0.85, 0.55, 0.62, 1)
setc(0, 0)
_bb = [O["AN_Uterus"].matrix_world @ Vector(c) for c in O["AN_Uterus"].bound_box]
cu = sum(_bb, Vector()) / 8
ctx = REP + ["AN_Bladder", "AN_Vagina_Space", "AN_HipBone_L", "AN_HipBone_R", "AN_Sacrum", "AN_Rectum", "AN_Colon_Descending"]
show(ctx)
shot("repro_front", cu + Vector((0, 0, 0.01)), (0, 1, 0), 0.17)
shot("repro_top", cu, (0, 0.25, -1), 0.17)
show(REP + ["AN_Bladder", "AN_Vagina_Space", "AN_Rectum", "AN_Rectum_LowerAmpulla", "AN_Colon_Descending", "AN_Sacrum"])
shot("repro_side", cu, (-1, 0, 0), 0.13)
for f, r, nm in ((0.0, 0.0, "rest"), (1.0, 0.0, "bladder_full"), (0.0, 1.0, "rectum_full"), (1.0, 1.0, "both_full")):
    setc(f, 0, r); show(org + REP + ["Hips", "AN_Vagina_Space"])
    shot(f"cut_repro_{nm}", (0, cu.y + 0.005, cu.z - 0.01), (-1, 0, 0), 0.12, clip=0.4)
setc(0, 0); show(REP + ["AN_Bladder"])
shot("eevee_repro_front", cu + Vector((0, 0, 0.01)), (0, 1, 0.15), 0.15, engine="BLENDER_EEVEE_NEXT")
print("VIEWS DONE")
