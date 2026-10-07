# NEXT SESSION START HERE

**Project goal:** functioning organ systems from the base of the neck down, in Blender, toon shader, female model
first (`Hips.blend`), a male version later. Done: the digestive system (esophagus -> anus, real openable passage).
**Current: the female pelvis - `scripts/build_pelvis.py` (was `build_urinary.py`): urinary system, female bony pelvis,
internal reproductive organs (first pass), built on build5_passage's output.**

**Urinary system, state (details under "Urinary system" below):**
1. Kidneys / ureters / adrenals from the atlas, in `Internal_Fit` under `Internal_Fit_Xform`. The atlas ureter already
   carries the renal pelvis + calyces inside the kidney; the kidney wall is opened where the pelvis leaves it at the
   hilum (`tract_lib.open_junction(..., sides=(0,))` = pierce).
2. New bladder `AN_Bladder` (the atlas one was a 226-face blob pressed on the rectum): pressure-grown shell, empty
   ~70 ml -> full 500 ml in 4 stages, against the pubic bones, sacrum, abdominal wall (skin - 10 mm), a retropubic
   floor and the anterior vaginal wall plane (room kept for the vagina). Trigone (neck + ureteric orifices) still.
3. Ureters bend to the trigone corners and open into the bladder (shared loops).
4. New `AN_Urethra` (38 mm real / 27 mm in the body) from the bladder neck (shared loop) to the external meatus: a
   3.6 mm sagittal slit cut into the crest of the vestibule's slot (shared loop with the skin). Collapsed at rest
   (transverse crescent -> sagittal slit), `AN_Urethra_Wall` sleeve (thicker at the external sphincter).
5. One control object, **`Pelvic_Controls`** (custom properties **Bladder_Fill**, **Void**, **Rectum_Fill**, each 0..1) drives every key:
   bladder `Fill_1..4` (chained stages) + `Void`; urethra `Void`; skin `Void`; and the neighbours' corrective keys
   `Bladder_Fill_1..4` (descending/sigmoid colon, rectum, lower ampulla, ureters) - they give way, nothing clips.
6. Insides lit for every organ / bone material (`tract_lib.light_all`, all three atlas shader families), AO_in
   baked; also in the atlas tract build (all 122 meshes).
7. Internal reproductive organs (first pass, see "Internal reproductive organs" below): uterus with its canal and
   cavity, uterine tubes, ovaries, ovarian / round / suspensory ligaments, the cervix opening into the vagina
   placeholder's vault. The uterus tilts back on its cervix as the bladder fills (22 deg at full); bladder, rectum,
   bowel, tubes and ovaries give way around it. Renders: `reproductive_views.png`, `reproductive_cutaways.png`,
   `pelvis_sections.png` (midline sections: rest / bladder full / rectum full / both full).

**Builds (blends gitignored, rebuild locally; renders + report in `out/urinary/`):**
```
cd drive
ANATOMY_REF=$PWD/anatomy_ref.blend blender -b Hips.blend --python ../scripts/build5_passage.py -- \
    $PWD/../out/blend/Hips_build5_passage.blend <report dir>          # anus / tract working file, ~25 s
ANATOMY_REF=$PWD/anatomy_ref.blend blender -b ../out/blend/Hips_build5_passage.blend \
    --python ../scripts/build_pelvis.py -- $PWD/../out/blend/Hips_urinary.blend ../out/urinary   # ~25-30 min
blender -b anatomy_ref.blend --python ../scripts/build_tract_passage.py -- \
    $PWD/../out/blend/anatomy_tract_passage.blend <report dir>        # full atlas tract, ~70 s
```
Setup: Blender 4.2.3 at /opt (see "Project context"), `apt-get update` before the EEVEE libs, Pillow in Blender's
Python (`/opt/blender-4.2.3-linux-x64/4.2/python/bin/python3.11 -m pip install pillow`). Drive files via gdown into
`drive/` (the network policy must allow download.blender.org and Google Drive).

**Open questions for the user / next steps (urinary):**
- The vulva is rudimentary; the reproductive step rebuilds it. The meatus then needs re-cutting: re-run
  build_pelvis with the new skin (and `URETHRA_MEATUS` / `URETHRA_MEATUS_FROM` if the vestibule moved).
- The organ / bone atlas is a male reference body (the skin is female). Done for the female layout: the bony pelvis
  reshaped (0b), room for the vagina (`AN_Vagina_Space`), the bowel moved out of it.
- Pelvis: done (female reshaping, step 0b). Rectum filling and both full: done (5d, 5e). The bones were reshaped
  only in the Hips working file; the full-atlas tract file (anatomy_tract_passage) still has the male pelvis.
- Male version later: no vaginal plane, prostate round the urethra's first 3 cm, a ~20 cm urethra through the penis.
- Fill is linear between the stages: 70 / 149 / 257 / 361 / 471 ml at Bladder_Fill 0 / .25 / .5 / .75 / 1 (real ml;
  grown to 150 / 260 / 380 / 500 targets, ~5 % lost to the crease smoothing, the polish and the room the tilted
  uterus takes - 470 ml still is a normal full bladder). Rectum 109 -> 176 / 243 / 308 ml. Both full: bladder 470 ml.
  The checks run at those values; in-between values blend two clean states.
- Remaining contact (`checks` -> `clipping` in report.json, 14 states; only contacts new or > 0.3 mm deeper than at
  rest, intended overlaps excluded): rest and Void states clean; bladder-only and rectum-only states <= 0.8 mm except
  one colon / rectum vertex against the uterus at Bladder_Fill 0.5 (2.4-2.7 mm); both full <= 0.2 mm. The half-way
  mixes are the worst: bladder 0.5 + rectum 1 has the uterus 1.8 mm into the rectum (370 vertices) and the same single
  vertices at 2.4-2.9 mm. Tubes / ovaries touch the hip bone by <= 0.6 mm in some states. Seams (neck, meatus,
  ureteric openings) <= 0.0005 mm everywhere.
- Crease round the still trigone: smoothed (Taubin passes in the 16-46 mm band, then a short regrow to put the
  volume back); see `eevee_cut_inside.png`.
- Lighting of the full atlas tract file: every material is now lit inside (122 meshes); `AO_in` of small closed
  shells (bronchial trees, sphincter) is dark (0.04-0.07), as a fully enclosed inside would be.

## Urinary system (female) - `scripts/build_pelvis.py`, `scripts/urinary_lib.py`

Runs on build5_passage's output (`out/blend/Hips_build5_passage.blend`) and writes `out/blend/Hips_urinary.blend` +
`out/urinary/report.json`. Everything is built in "atlas space" (the frame under `Internal_Fit_Xform`, real-size
metres): the fit squashes it to the body (0.62 / 0.72 / 0.72), so body-scale lengths are ~0.7x and volumes 0.32x the
real ones. Objects added / replaced in `Internal_Fit`: AN_Kidney_L/R, AN_Adrenal_L/R, AN_Ureter_L/R, AN_Bladder,
AN_Urethra, AN_Urethra_Wall, AN_Vagina_Space (placeholder); empties `Pelvic_Controls` and `Urethra_Meatus_Target`.

**How to use it (for the user):** select `Pelvic_Controls` -> Object properties -> Custom Properties: drag
**Bladder_Fill** (0 empty ~70 ml .. 1 full 500 ml, real volumes), **Void** (0 closed .. 1 voiding) and
**Rectum_Fill** (0 resting ~124 ml .. 1 full ~325 ml). Everything that has to
move follows through drivers (simple expressions, no Python / auto-run needed). A shape key is a stored version of a
mesh's shape that can be blended in by a slider; each control blends several of them.

**Steps (numbers from the last build):**
1. Kidneys: the atlas ureter mesh already includes the renal pelvis and its calyces inside the kidney (top 4-5 cm,
   up to 9 mm wide). The kidney's wall is opened where the pelvis leaves it at the hilum: `open_junction(K, Ur,
   sides=(0,))` (new `sides` option: pierce - only A loses the piece inside B). Openings r 7.4 / 10.5 mm.
2. Bladder, pressure-grown (`urinary_lib.grow`: membrane pushed out along its normals with PI control on the
   volume, tangential + fairing smoothing, projected out of `Obstacles` every few passes):
   - neck `BLADDER_NECK` [0, -0.080, 0.776] (atlas): ~19 mm behind the back of the symphysis at its lower border;
   - limits: hip bones + sacrum (3 mm), inside the skin by `WALL` 10 mm (abdominal wall), a retropubic floor
     (`FLOOR_DEG` 35, rising forward from the neck), the anterior vaginal wall plane (`VAG_DEG` 35 up and back, 5 mm
     behind the neck, for `VAG_LEN` 75 mm): room for the vagina between urethra / bladder and rectum;
   - empty stage (70 ml): seed 10 mm sphere on the neck, lid pressed down to neck + 32 mm (the bowel), the bowel and
     anal canal solid too; stages 150 / 260 / 380 / 500 ml (`BLADDER_ML`) grow from it with the bowel soft;
   - the trigone stays still in every stage: within `TRIG_FIX` 16 mm of the neck or of any vertex of the two
     ureteric openings, free by `TRIG_FREE` 50 mm.
   - Neck: the icosphere is cut 5 mm round the neck and the hole filled (constrained Delaunay) onto the urethra's
     first ring (32 vertices, a transverse crescent, lips 0.02 mm).
3. Ureters: last 60 mm bent (translation blend) so the end sits 3 mm inside the bladder at the planned orifice
   (22 mm up the base from the neck, 13 mm either side); the atlas' left ureter ran through the sigmoid (44 vertices)
   and is first eased out of the bowel at rest; then `open_junction` with the bladder. The tubes cross the wall
   where they meet it (the right one ~2 cm from the planned spot), so the still trigone is built round the actual
   opening loops, BEFORE the fill stages are grown.
4. Urethra: a cubic from the neck (leaving straight down the bladder's axis) to the meatus (arriving straight up),
   38 mm real / 27 mm in the body, 56 rings of 32. Rest: transverse crescent turning (55-95 % of the way) into the
   meatus' sagittal slit; Void: round 7 mm at the neck -> 5.2 mm. "lining" attribute 0 at the meatus -> 0.62 by 8 mm on
   a copy of the canal lining material (`AN_Urethra_Lining`). `AN_Urethra_Wall`: a closed sleeve round it (0.4 mm
   clear of the open lumen, +1 mm over the external sphincter's third), trimmed clear of bladder and skin.
   - Meatus (`urinary_lib.cut_meatus`): the vestibule here is a 5 mm deep slot whose labial walls meet in a knife
     edge crest (0.1-0.4 mm apart). `CrestChart` unfolds it: u along the crest, v = arc length down either wall
     from the crest, through cross-section polylines; the skin round the target is refined (edges split), the
     crest's unwelded duplicates welded (5 um) and zero-thickness midline flaps removed (they made the crest
     branch), a (u, v) box is re-triangulated (CDT with a 0.22 mm point grid) round a K=32 slit 3.6 mm long. UVs are
     sampled on the same side of the midline seam, every existing shape key (Open) keeps Basis there.
   - Target: `URETHRA_MEATUS="x,y,z"` (body coords, snapped up onto the crest), or `URETHRA_MEATUS_FROM=<previous
     output .blend>` (its `Urethra_Meatus_Target` empty, move it by hand), or the default (0, 0.0135, crest):
     ~31 mm behind the clitoral hood, ~13 mm in front of the introitus dimple. Tested at y 0.0110 / 0.0175 and
     0.5 mm off the midline: all build.
5. Keys + drivers: bladder `Fill_1..4` (each relative to the one before: driver f*4-(k-1), clamped 0..1 by the key
   range) + `Void` (neck funnels open, round); urethra `Void`; Hips `Void` (the lips part into a 2 mm lens, fading out
   over the rebuilt patch). Neighbours (`urinary_lib.make_room`): descending/sigmoid colon, rectum, lower ampulla,
   ureters (anal canal pinned; ureters pinned 15 mm round their bladder openings and inside the kidney) get
   `Bladder_Fill_1..4`: per stage, soft vertices inside (or within 2.5 mm of) the grown bladder are pushed out,
   faces the bladder bulges through are pushed past it, bowel pressed into bowel steps back (both half), the pushes
   spread as a decaying average over each organ (welded seams move together), bones / skin rigid. Then the bladder
   settles (`urinary_lib.settle`): where it still presses into bowel trapped against bone it takes a shallow dent.
   Before the pushes, soft tissue next to the wall is carried along with the wall's own motion (`carry_from`).
   At full the sigmoid moves up to ~46 mm (real scale) up / back, the rectum ~37 mm, the lower ampulla ~6 mm (body).
6. Lighting: new organs get the outer "AO"; `tract_lib.light_all` lights every organ / bone material (all three
   atlas shader families now: organ, skeleton / cartilage, lung) and bakes AO_in where missing (the urethra's with it
   open); build5_passage's bowel bakes are kept.
7. Checks in `report.json` -> `checks`, at Fill 0 / 0.25 / 0.5 / 0.75 / 1 x Void 0 / 1: seams (neck, meatus, both
   ureteric openings), bladder volume, and every pair of organs for vertices inside another (count, deepest mm).

**Gotchas:**
- Setting a custom property from Python does not re-run drivers: `ctl.update_tag()` then `view_layer.update()`
  (the UI does this itself).
- `open_junction` round-trips vertices through body coordinates (~1e-7 m): match old vertices within 2 um.
- `signed_dist` decides inside by ray parity (the nearest-face normal lies near capped rims and collapsed walls).
- The bladder stages and the bladder object have different topology (icosphere vs opened / junction-cut mesh): keys
  are mapped by position; new vertices take their nearest old vertex's motion.
- The atlas is a male body: its rectum sits where the vagina must go. The reproductive step will need the rectum /
  ampulla moved back (or compressed) by roughly the vagina's thickness; the vaginal plane already keeps the bladder
  and urethra in front of it.

**Female bony pelvis (step 0b, first thing after loading):** `urinary_lib.female_pelvis_field` - one smooth
displacement for both hip bones and the sacrum (joints stay matched): below the acetabula the inferior pubic /
ischial rami and the tuberosities spread out (12 mm a side, nothing at the symphysis), the ischial spines move out
6 mm, the brim and the sacral wings widen 5 mm a side, the lower sacrum / coccyx straighten back 6 mm, the iliac
wings sit a little lower. Measured (atlas mm, `report.json` -> pelvis): intertuberous 85 -> 109, brim width 123 ->
133, sacrum width 118 -> 123; skin clearance stays 5.3 mm (body; backed off automatically below 2.5 mm).
`FEMALE_PELVIS=0` keeps the atlas pelvis, other values scale it. The `subpubic_angle_deg` metric is unreliable (it
read 91 deg on the male arch) - judge the arch from `out/urinary/pelvis_male_vs_female.png`.

**Rectum filling (5d):** AN_Rectum + AN_Rectum_LowerAmpulla welded into one reservoir (`weld_union`), its openings
onto the sigmoid and the anal canal capped (`cap_holes`) and held still (8 mm, free by 25 mm), grown like the bladder
in 3 stages, +`RECTUM_ADD_ML` 200 ml: 124 -> 191 / 258 / 325 ml (real). Solid for it: bones (3.5 mm), the abdominal
wall (5 mm), the vagina space, and - on its own - the bladder and the ureters (in a woman the vagina / uterus take that
push). Keys `Rectum_Fill_1..3` on both (driver f*3-(k-1)); the sigmoid and the ureters give way (`make_room`, carried
along, ureters held only 6 mm round their bladder openings) with keys of the same names.

**Both full (5e):** with Bladder_Fill = Rectum_Fill = 1 the keys add and would collide; from that sum the bladder
settles away from the full rectum (it holds less: 495 instead of 500 ml), then the sigmoid and ureters are cleared
of both. The difference is a `Both_Full` key on bladder, colon and ureters, driven by `a*b` (Bladder_Fill x
Rectum_Fill), so it fades in with either.

**Internal reproductive organs (step 2b, `scripts/repro_lib.py`; first pass, connected to the vagina placeholder):**
all real-size in atlas space, parented to `Internal_Fit_Xform`, materials tinted copies of the atlas organ material.
- `AN_Uterus` (outer surface) + `AN_Uterus_Lumen` (cervical canal + flat triangular cavity), metaballs (built in mm -
  Blender clamps metaball resolution to >= 5 mm - then decimated to ~5-6k vertices each). Cervix 25 mm, square to the
  vagina (anteversion), body bent forward on it by `UT_FLEX_DEG` 35 (anteflexion); ~52 x 65 x 50 mm. The two meet at
  the external os (shared loop, `open_junction(..., outside=(1,))`: the lumen's stub outside the uterus is removed);
  the cervix pierces the vagina placeholder's vault (shared loop, 40 mm round).
- `AN_Uterine_Tube_L/R` (~127 mm): from the cavity's corners (shared loops), through the wall, over the ovary, ending
  in an open fimbriated funnel. `AN_Ovary_L/R` (34 x 20 x 13 mm, small surface bumps) on the side walls.
  Ligaments as cords: `AN_Ovarian_Ligament_*` (ovary -> uterus), `AN_Round_Ligament_*` (uterus -> deep inguinal
  ring, inside the abdominal wall), `AN_Suspensory_Ligament_*` (ovary -> pelvic brim); ends sunk 1 mm.
- The bowel and ureters give way to them (`bowel_room`, keys moved too); tubes / ovaries are kept 2 mm off bone.
- The empty bladder is grown under the uterus (solid), then the uterus settles on it (`repro_tilt`: 7 deg back);
  its top is dented where the uterus lies on it, as in life. Each filling stage is grown freely, the uterus tilt that
  clears it is found, and the stage is grown again with the tilted uterus, tubes, ovaries and ligaments solid (two
  shells: the uterus, the rest - one shell of overlapping organs fools the ray-parity inside test), so the bladder
  takes its volume elsewhere (the dome rising along the abdominal wall) instead of being dented.
- Filling: the uterus is firm - it is never pushed out of shape, it only tilts on its cervix (`tilt_disp`: no turn
  below 10 mm up the cervix, all of it above 28 mm; the least turn that clears the filling bladder, stopped by bone
  or 8 mm into the rectum, at most 40 deg). Its neighbours give way to it (make_room: a held organ's face makes the
  other take the whole push) and the bladder settles (dents) where it presses on. Tubes, ovaries, ligaments follow
  part of the tilt (less the further from the uterus; ligaments blend between their ends) and are soft. Rectum full:
  the uterus stays (no room to tip forward with the bladder in front), the rectum settles against it. Both full:
  every soft organ starts from the sum of the two fillings (the least correction), the rectum and bladder settle
  against the uterus (`Both_Full` keys on the rectum and ampulla too). Settling now dents inward (`settle(...,
  inward=True)`: a vertex goes back along its own inward normal to where it leaves the neighbour; the nearest way
  out wrapped the bladder round the uterus).
- Not done yet: the vagina itself (still the placeholder `AN_Vagina_Space`), the vulva, broad ligament / peritoneum,
  blood supply; the male side.
- What was tried: the uterus soft like the bowel -> its fine mesh was torn into shards when squeezed 10-15 mm (the
  vertex pushes are fine on coarse bowel meshes only); the empty bladder grown without the uterus, then the uterus
  tilted onto it -> the cervix sat in the bladder, settling wrapped the bladder round it; Both_Full solved from rest
  -> half-way states (one full, one half) blended two unrelated arrangements and clipped 4 mm; the bladder grown
  without the tilted uterus and dented afterwards -> it held only 330-420 ml at full.
- Gotcha: AN_Rectum / AN_Colon_Descending keep their own object transforms (the new organs are identity under
  `Internal_Fit_Xform`): convert local <-> atlas (`_to_atlas` / `_to_local`) before mixing them with atlas data.
- Speed: `urinary_lib.lap_avg` (np.bincount) replaced np.add.at in every diffusion loop (~11x); the polish and the
  checks test only vertices inside a neighbour's bounding box. Blender ignores PYTHON* env vars: add
  `--python-use-system-env` (with PYTHONUNBUFFERED=1) to see the log live.

**Polish (5f) and checks:** after all keys, every stage state (bladder 0.25 .. 1, rectum 1/3 .. 1, both full) is
polished: every organ with a key active there steps out of anything it is newly inside (not inside at rest), and
pushes off fixed vertices (bone, held organs) poking through its faces, with a short falloff - written into that
stage's key (chained keys telescope, so a correction touches only the neighbouring states). The checks report, per
state, only contacts that are new or > 0.3 mm deeper than at rest (`clipping`), the rest ones once
(`clipping_at_rest`); overlaps by design are left out (lumen in uterus, cervix in the vault, tubes through the wall,
ligament ends).

**Female pelvis (step 2, before the bladder):** `AN_Vagina_Space` (wire display, not rendered) is a placeholder for
the collapsed vagina with its walls: a flattened tube (26 x 10 mm real, narrower at the bottom, rounded fornix at the
top) from 10 mm above the introitus dimple (`VAG_INTROITUS`, body coords) up behind the urethra to the bladder-base
plane at the neck's height, then 40 mm (`VAG_UP_MM`) up that plane (~64 mm long, ~75 mm from the introitus; 30 mm put
the cervix right on the bladder neck). The bowel (sigmoid, rectum, lower
ampulla, anal canal except its lowest 12 mm at the anus) is moved out of it with `make_room` (2 mm septum); every
shape key of a moved organ moves with it, so the anus' Open key still works. Moved: ampulla 2.7 mm, canal top 2.7 mm
(body), rectum 0.8 mm (atlas). The bladder grows clear of it (solid, 1 mm), the ureters are eased out of it, the
bowel treats it as solid while the bladder fills, and it takes a groove round the urethra (urethra in the anterior
vaginal wall, 2.5 mm). The reproductive step should build the vagina in it (and may let it compress).

**What was tried for the neighbours (so it is not repeated):** pushing soft vertices to the nearest bladder point
only -> the sigmoid crumpled (its parts pushed every which way) and the bladder bulged through coarse rectum faces;
adding bulge-through and soft-soft pushes -> better, but the colon slid into the rectum as a whole (a 3 cm falloff
moves both walls of the tube); shorter falloff -> helps; carrying the bowel along with the wall's own motion
(`carry_from`) -> the big improvement (no crumpling); the bladder settling against what is trapped -> clears most of
the rest. 3-ray inside tests gave false answers near grazing edges: 5 rays now.

**Tools:** `scripts/tools/sagittal.py` (any objects, any axis-aligned plane; `SAG_KEYS="Bladder_Fill=1,Rectum_Fill=1"`),
`scripts/tools/urinary_views.py` (overview, midline cut-aways at Fill 0 / 0.5 / 1, rectum full, both full, meatus
closed / voiding, EEVEE cut-away of the full bladder's lit inside, reproductive organs front / top / side and midline
cut-aways at rest / bladder full / rectum full / both full). Renders in `out/urinary/`.

---

# Digestive system (previous task) - state when the urinary work started

**State: all five handoff items are DONE** (details further down, under "Progress, item N"):
1. Real passage: skin slit edge loop = canal mouth (288 points, 0 mm seam). `scripts/build5_passage.py`.
2. "Open" shape key (0..1) on Hips, AN_AnalCanal, AN_Rectum_LowerAmpulla, AN_AnalSphincter (12 mm, `ANUS_OPEN_D`).
3. Collapsed resting canal: slit, then a closed 8-armed star (anal columns).
4. Lining: skin -> anoderm -> mucosa ramp (`AN_AnalCanal_Lining`); organ insides lit via `AO_in` bakes.
5. Whole tract connected: `scripts/build_tract_passage.py` + `scripts/tract_lib.py` open all 9 junctions
   (esophagus -> rectum) with exact shared openings; Hips build also opens descending colon -> rectum.

**Builds (both outputs gitignored, rebuild locally):**
```
cd drive
ANATOMY_REF=$PWD/anatomy_ref.blend blender -b Hips.blend --python ../scripts/build5_passage.py -- \
    $PWD/../out/blend/Hips_build5_passage.blend <report dir>          # working file, ~20 s
blender -b anatomy_ref.blend --python ../scripts/build_tract_passage.py -- \
    $PWD/../out/blend/anatomy_tract_passage.blend <report dir>        # full tract, a few minutes
```
Setup: Blender 4.2.3 at /opt (see "Project context"), `apt-get update` before the EEVEE libs, Pillow in Blender's
Python for `canalxs.py` (`/opt/blender-4.2.3-linux-x64/4.2/python/bin/python3.11 -m pip install pillow`). Drive files
via gdown into `drive/` (the network policy must allow download.blender.org and Google Drive).

**Open questions for the user / possible next steps:**
- Appendix was moved 21.9 mm to attach to the cecum >= 15 mm from the ileal entrance (not opened: dead end). The
  user may want a smaller move (e.g. `APPX_CLEAR` 12 mm) or a different attachment spot.
- Canal: it opens from the star over the top 40 % (from ~15 mm); anatomically it could stay closed over the
  sphincter's full height and open more sharply at the anorectal junction (offered, not done).
- Lining colour at the canal mouth renders a little greyer than the rim skin (~180,141,142 vs 204,162,170); user
  hasn't asked to change it.
- Only tract organs have lit insides; other atlas organs keep the flat backface colour.
- Where the full-tract file should live in the user's real project (they have the full atlas on Windows at
  `C:\Users\Parker\Anatomy Project\`) - the atlas-derived file is a new file, anatomy_ref.blend is untouched.

**Tools added this session (`scripts/tools/`):** `section.py` (midline cut-away), `selfx.py` (self-intersections,
`ANUS_OPEN`), `canalxs.py` (canal cross-sections), `tract.py`, `overlaps.py`, `tract_seams.py` (every open loop must
be matched by a neighbour), `xsect.py`, `cutaway.py`, `opening_view.py`, `holes_view.py`. `hq.py` / `star.py` take
`ANUS_OPEN`.

**Gotchas:** atlas "holes" are unwelded shading seams (weld at 1e-7 m to test); the ampulla loop in build5_passage
reassigns `e1, e2`; `mesh.transform` does not move shape keys (key the sphincter after re-seating it); inside tests
need closed shells taken BEFORE junctions are opened (`trees=`); neighbours split from one source mesh share
coincident walls (open_junction's inflate / deflate retry).

---

# Earlier state (history)

**Current best build:** `scripts/build5_rise.py` (all defaults). The built file is
`out/blend/Hips_build5_rise_corridor.blend` (gitignored; rebuild with the command under "How to build"). It has:
- build5's radial creases: `ANUS_SIG` 0.06 and 2x depth;
- a 0.2 mm rise, then a smooth 1.5 mm plunge into the centre;
- levelled along the cleft, with the corridor fix for the slope behind the anus;
- the normal map faded over the cleft, and the cleft mesh refined.

The user approved the look. Alignment checks pass: the canal bottom ring sits on the skin (0.00 mm gap) and is
centred within 1.3 mm. Canal 37 deg from vertical, anorectal angle 108 deg, clearances unchanged.

**Progress, item 1 DONE (`scripts/build5_passage.py`, file `out/blend/Hips_build5_passage.blend`):**
- The pole and the rings inside `ANUS_HOLE` (0.3 mm) are gone. The 0.3 mm ring (288 points) is the skin's edge loop.
  Inside `ANUS_COLL` (0.75 mm) the rings collapse sideways, so the edge loop is a closed slit 0.6 mm long with the
  lips 0.020 mm apart, at the bottom of the plunge.
- `AN_AnalCanal` starts ON that loop: the same 288 points (seam gap 0.0000 mm, in `report.json` -> `passage`). It
  carries the slit straight up for `ANUS_CANAL_COLL` (2 mm), opens to build5_rise's section by 7.5 mm, and steps 2:1
  down to the 36-point junction ring. Still two objects (skin / canal) sharing the seam, not one joined mesh.
- Gotcha fixed: the cleft levelling (Jacobi, 3000 passes) had relied on the pole; without it the middle of each lip
  lagged ~1 mm below the slit tips. It now runs on build5_rise's closed layout with the removed rings + pole as
  virtual points.
- The corridor no longer holds the mouth loop fixed (it is a boundary now).
- Checks: skin within 0.05 mm of build5_rise beyond 2 mm from the centre, 0 non-manifold, 0 degenerate faces, canal
  angle / anorectal angle / clearances unchanged. Renders: `out/variants/build5_passage/` (`cmp_*` shaded,
  `section_cmp_*` midline cut-away via the new `scripts/tools/section.py`).
- A brownish hairline shows in the slit (canal material) - item 4.

**Progress, item 3 DONE (same script):** the canal is collapsed at rest. The lumen is the skin's slit carried up 2 mm,
then an AP slit lengthening to 2.5 mm half-length by 7 mm, then a closed 8-armed star (sides 2.0 mm from 5-10 mm,
diagonals 1.8 mm from 7-12 mm; walls 0.01 mm off each arm's midline). The tissue between arms = the anal columns.
Over the top 40 % it eases out to the round 36-point junction ring. Each arm owns a fixed block of 36 points (tip in
the middle), and point k sits near angle -pi + 2 pi k / 288, ready for the Open key to spread them round a circle.
Checks: 0 self-intersections (new `scripts/tools/selfx.py`), every canal point beyond the seam inside the skin
(canalout), 0 non-manifold. Plots: `canal_sections_rest.png` (new `scripts/tools/canalxs.py`; needs Pillow in
Blender's Python), `star_section_cmp.png`. Next: item 2, the Open shape key (proposed ~12 mm across at the skin,
subtle columns, unless the user says otherwise).

**Progress, item 2 DONE (same script):** an "Open" shape key (0..1) on Hips, AN_AnalCanal, AN_Rectum_LowerAmpulla and
AN_AnalSphincter. Skin: each ring follows a profile (flat from the 8.6 mm ring, a rounded rim of 1 mm, then a wall
2.7 mm down into the canal), old rings spread along it by arc length, so the old plunge becomes the opening's wall;
the opening is round, `ANUS_OPEN_D` (12 mm) across; creases flatten to 35 %; the rim rolls out 0.4 mm. Canal: a round
tube from the skin's opened mouth (same points: seam 0.0000 mm) to the junction, which widens x1.35 with the ampulla
easing (fades to 1 at the rectum rim; junction gap 0.0001 mm); anal columns as 5 % ridges. Sphincter pushed out
uniformly (0.64 mm) to clear the open canal by 0.4 mm. Checks at open 0 / 0.5 / 1: canal 0 self-intersections, skin
none within 30 mm of the anus (the 157 elsewhere are in the original too), no internal object inside the open canal.
Gotcha: the ampulla loop reassigns `e1, e2`; 9b resets them. `hq.py`, `star.py`, `selfx.py` take `ANUS_OPEN`
(`star.py` also `ANUS_STAR_ALL=1`). Renders: `open_states_0.03.png` (shaded), `geo_open_states.png` (geometry).
Next: item 4 (lining material: the open canal shows as a flat brown disc in the toon shader).

**Progress, item 4 DONE + insides lit (section 9c):** the organ shader (emission toon, AN_Rectum and
AN_Colon_Descending) painted every backface one flat brown (`Mix.003`, colour 0.527/0.159/0.118), and the generated
canal / ampulla had no baked `AO`. `light_insides()` now lights backfaces with the flipped normal and switches the AO
input to a new `AO_in` attribute on backfaces (the flat colour is bypassed). AO is baked by ray casting (32 cosine
rays, 8 mm reach, all meshes): `AO` (outside, rest state) for canal + ampulla, `AO_in` (inside, OPEN state) for
canal, ampulla, rectum, descending colon. The user asked for this ("smooth lighting" like Minecraft / Surgeon
Simulator also on the insides). New material `AN_AnalCanal_Lining` (canal + ampulla) = the lit organ shader with base
and rim colours from ramps on a `lining` attribute (height up the canal 0..1; ampulla 1): skin (calibrated against
the rendered rim skin) to 2 mm, pale anoderm 5-11 mm, red mucosa (the rectum's colour) past ~15 mm (dentate line).
The open junction is round (the 8-fold ripple fades back in up the ampulla). The closed slit's hairline is now
dark instead of brown. Renders: `lining_open_states.png`, `lin_open1_close.png`, `hairline_before_after.png`.
Next: item 5 (mouth-to-anus continuity check).

**Item 5, checked (`scripts/tools/tract.py`, results in `out/tract/`):**
- The atlas has no mouth / pharynx (it starts at the neck, C5) and no sigmoid colon (the descending colon runs down
  to the rectum). Chain: Esophagus, Stomach, Duodenum, Jejunum, Ileum, Cecum (+ Appendix), Colon Asc / Trans / Desc,
  Rectum, then (Hips.blend) LowerAmpulla, AnalCanal, skin.
- Surfaces: every consecutive pair touches (closest gap 0.00-0.07 mm, overlapping faces at each junction), except the
  APPENDIX, 4.1 mm (atlas units) clear of the cecum. Hips.blend: Desc colon - rectum 0.02 mm; rectum - ampulla -
  canal - skin share exact loops (0.00 mm).
- Lumens: the atlas organs are separate shells (many with dozens of small holes: Jejunum 76 boundary loops, Ileum 85),
  overlapping at the junctions; open ends face each other only at ileum -> cecum (shared loop, 0 mm) and roughly at
  duodenum -> jejunum (1.6 mm). The stomach is closed (no inlet / outlet opening), the esophagus does not open into
  it, and jejunum->ileum, cecum->ascending, asc->trans, trans->desc, desc->rectum have no facing openings (7-50 mm off).
  So the tract is visually continuous, but a hollow passage only from the rectum down (and ileum->cecum).
- Internal_Fit_Xform scales the atlas by (0.62, 0.72, 0.72), so atlas mm are ~0.7x in the fitted body.
- User's answers: no mouth / pharynx needed (the model starts at the base of the neck); fix holes if inaccurate;
  look for missing parts; make the full connected passage only if the models are not changed much.

**Item 5 DONE - connected passage (`scripts/tract_lib.py`, `scripts/build_tract_passage.py`):**
- Holes: there were none. Every "hole" was an unwelded shading seam; each organ is a closed shell once welded at
  1e-7 m. Not changed.
- Missing parts: none in the depicted region. The "descending colon" model includes the sigmoid's course (curves
  into the pelvis and overlaps the rectum 11 mm deep, `out/tract/lower_colon.png`); it is just not split out.
- Atlas defect fixed: the ileum carried a complete duplicate of the cecum (all 288 faces, a separate closed copy).
- `open_junction(A, B)`: cuts both shells along their intersection inside a patch round the junction (bpy
  mesh.intersect, EXACT, A against B only), opens on each side the piece the junction curve encloses that lies
  inside the other organ, and puts the patch back welded to the untouched rest. Both openings are the same vertex
  loop. Side contacts between coils stay closed. Patch grows 12 / 24 / 40 mm until the curve closes; loose curve
  ends <= 3 mm apart are bridged and the slivers stitched; if still open, one side is inflated / deflated slightly
  (fading to 0 at the patch border). Rebuilt vertices take "AO" from the nearest old one.
- Atlas result (`out/blend/anatomy_tract_passage.blend`, gitignored; anatomy_ref.blend itself untouched): all 9
  junctions open, radii mm: esoph-stomach 8.0, stomach-duod 10.6, duod-jejunum 8.1 (needed a 2.5 mm / 15 mm bulge:
  they only touched), jejunum-ileum 12.7 (jejunum deflated 1 mm locally: coincident walls), ileum-cecum 8.9,
  cecum-asc 12.4, asc-trans 11.6, trans-desc 10.7, desc-rectum 25.4. `tools/tract_seams.py`: 18 openings, all matched
  (0.0001 mm), no other open edges; 0 non-manifold, 0 degenerate faces.
- Appendix: attached (moved 21.9 mm so its base sinks 1 mm into the cecum >= 15 mm from the ileal entrance; the
  nearest spot was on the ileal opening), NOT opened (dead end). The user may prefer a smaller move.
- Insides: every tract organ's material is lit inside + "AO_in" baked (atlas build step 5).
- Hips.blend (`build5_passage.py` 9c): descending colon -> rectum opened the same way (r 16.9 mm); so the working
  file runs colon -> rectum -> ampulla -> canal -> skin with exact seams (`out/tract/seams_hips_passage.json`).
- Renders: `out/tract/openings.png` (organs alone looking into their openings), `junction_cutaways.png` (top: atlas
  before, bottom: after, insides now lit).
- Tools: `tract_seams.py`, `overlaps.py`, `xsect.py` (2D sections), `cutaway.py`, `opening_view.py`, `holes_view.py`.

**What is missing (the next task, in order):**
1. **Real passage.** The skin is closed at a point (the pole) and the canal (`AN_AnalCanal`) starts 0.08 mm under the
   0.9 mm ring, so they are separate surfaces.
   - Open the skin at the centre (remove the pole and the inner rings below ~0.3-0.9 mm) and weld its edge loop to the
     canal's bottom ring, like the earlier slit versions (build8_v5/v6: "canal opening = skin slit edge loop").
   - Keep it visually closed at rest: lips touching, no visible hole.
2. **Opening shape key** (e.g. "Open") on Hips + AN_AnalCanal:
   - the creases flatten and spread, and the rim widens and rolls outward;
   - the canal dilates (and the lower ampulla eases); a 0..1 slider.
3. **Resting canal as a collapsed slit with anal columns**, not a round tube. Currently it is only lightly flattened.
4. **Lining colour:** the material should go from skin to mucosa inside the opening.
5. **Bigger picture:** check the rest of the tract from mouth to anus is continuous, without gaps between organ meshes
   in `anatomy_ref.blend` / the Internal_Fit collection.

**Don't redo:** the long history below covers the crease and slope experiments and their results. Key tools:
`backslope.py`, `distort.py`, `grazing.py` (TOP = front), `hq.py`, `tipcam.py`, `chk.py`, `canalout.py`.

---

# Session handoff: anus detail on `Hips.blend`

Branch: `claude/optimistic-gates-jorged`. Pull request: https://github.com/Eclectic-png/For-Claude/pull/1

## Project context

The user is building a detailed 3D anatomy model of the upper body's internal organs in Blender. This work
connects the digestive system to the exterior anus and details the anus itself. The style is a toon / anime shader.

- **Source files:** in the user's Google Drive folder
  https://drive.google.com/drive/folders/1oAWugWowv7qyMIMWI6KbRiJFBTXoJlFS
  - `Hips.blend`: the skin plus the pelvis subset of the internals.
  - `anatomy_ref.blend`: the full atlas. It's read-only and the scripts only copy `AN_AnalSphincter` from it.
- **Download:** `pip install gdown && gdown --folder <url> -O drive`. `drive/` and `*.blend` are gitignored.
- **Blender:** 4.2.3 LTS, headless. Install it with:
  `curl -sSLO https://download.blender.org/release/Blender4.2/blender-4.2.3-linux-x64.tar.xz && tar xf ... -C /opt && ln -sf /opt/blender-4.2.3-linux-x64/blender /usr/local/bin/blender`
- **EEVEE renders:** these need
  `apt-get install -y libegl1 libgl1 libgl1-mesa-dri libegl-mesa0 libxkbcommon0 libsm6 libxi6 libxxf86vm1 libxfixes3 libxrender1`
- **Local path:** the user's local chat originally wrote `build7.py`. On Windows it lives in
  `C:\Users\Parker\Anatomy Project\`.

## How to build

Every build script runs on the original `Hips.blend`:

```
cd drive
ANATOMY_REF=$PWD/anatomy_ref.blend blender -b Hips.blend --python ../scripts/build8_vN.py -- <out.blend or INPLACE> <report folder>
```

On Windows `ANATOMY_REF` defaults to the user's path. Each script writes `report.json` into the report folder.

## Versions

| Script | Status | What it is |
|---|---|---|
| `scripts/build5.py` | original (user's local chat) | Pole-closed anus with 12 radial creases of fixed angular width whose depth grows from the centre (`sin(pi t)^0.8`), so they all narrow into one point. The user likes its lines. It stood up like a lump (biharmonic fill bridging the cleft, flat projection) and reshaped the cheeks. Renders: `out/variants/build5/`. |
| `scripts/build5_adapted.py` | **current** | The user's build5 itself, with only three changes (each marked `adapted:`): the rebuilt region is cut narrow in surface coordinates (`surf_coords`, `D_AP` / `D_LAT`), the rings are cast on a Phong-curved base, and the lateral offsets are arc length over the skin (`surf_point`). There is also a portability fix (env var, report path, texture lookup). The cross-section (`xs2.py`) shows the centre at +2.3 mm (lump) in build5 and -0.1 mm here, with the cheek walls ~3 mm further out in build5. Renders: `out/variants/build5_adapted/`. |
| `scripts/build5_rise.py` | **in progress** | build5_adapted with the 2 mm funnel turned into a rise toward the centre: same `1 - ss(0, R_FOLD, r)` profile, applied along n, with build5's 0.6 mm entrance dip kept at the pole. Set the height with the `ANUS_RISE` env var (default 1.0 mm). Measured rises above the surrounding floor: 0.5 gives ~0.7 mm, 1.0 gives ~1 mm, 2.0 gives ~1.8 mm. Renders for all three: `out/variants/build5_rise/`. |
| `scripts/build7.py` | original | The local chat's last build (v7): a near-circular pucker. |
| `scripts/build8.py` | **saved v1** (`out/versions/closed_fissure_v1/`) | Squeezed cleft, closed fissure, lips peaking at the edge, cusped wrinkles. Uses the older confocal ring layout. |
| `scripts/build8_raised.py` | variant | v1 with a ~1 mm raised dome. |
| `scripts/build8_painted_lines.py` | variant | Crisp painted crease lines from a per-pixel mask texture. The user rejected painted lines and wants 3D geometry. |
| `scripts/build8_v2.py` | superseded | Straight wrinkles running into the crease. Beaded and jagged. |
| `scripts/build8_v2_wip.py` | parked | A failed experiment. |
| `scripts/build8_v3.py` | superseded | Straight grid "spokes" with creases on them, a short slit (3.6 mm), a closed slit, clean tip topology. Has a U-shaped artifact and sliver faces at the tips. |
| `scripts/build8_v4.py` | **saved v4** (`out/versions/even_mesh_v4/`) | **Current best.** The inner pucker is an even ~0.07 mm triangle mesh made by constrained Delaunay (`mathutils.geometry.delaunay_2d_cdt`). The lip ring plus four outer rings keep the join to the skin. The surface is a full membrane relax, plus a shaped lift that is smoothed and limited to 2.5 mm of the midline. `PUFF = 0`, creases on 20° spokes. The user approved its overall shape. |
| `scripts/build8_v5.py` | superseded by v6 | Like v4, but the pads between creases are domes meeting at the creases (`crease()` rewritten). The user still sees grooves. See next steps. |
| `scripts/build8_v6.py` | **in progress** | v5, plus a crease straight up the middle from the back tip (requested by the user), and the tip tidy and lip-curl taper shrunk to `TIP_R = 0.35` mm so the creases no longer fade out near the tips. |
| `scripts/build8_v7.py` | **experiment, rejected** | v6 with crease depth capped at `ASPECT` x half the gap to the neighbouring crease, and creases blended into the lip curl. It removed the tip knots and the spikes along the slit, but the creases fade near both tips. The user rejected it. |
| `scripts/build8_v8.py` | abandoned | Rays from the slit centre; crumpled at the centre. Superseded by build9. |
| `scripts/build9.py` | **in progress, current** | The user asked to base the work on build5. It is build5's anus (pole, 288 spokes, rings stepping 0.04 to 7.2 mm, build5's `crease()` and `fold_amp` measured from the centre) on v4's surroundings: narrow hole, arc-length layout riding up the cheek walls (`K_LAT = 0.77`), squeezed cleft, membrane relax plus level lift, and `FUNNEL = 0`. The canal opens under the centre (the 0.9 mm ring, tucked 0.08 mm under). Renders: `out/variants/build9/`. |
| `scripts/build10.py` | **in progress, current** | build9 plus build5's centre: `FUNNEL` 2 mm and a smooth `PIT` of 0.6 mm, both applied along -n (along the surface normal they dug a pocket into the steep left cheek wall). The rings are squeezed sideways near the centre (`PINCH` 0.3, easing out by `R_PINCH` 5.5 mm) for build5's pinched, gathered look. Creases run from the centre. Renders: `out/variants/build10/`. |

The renders for each version are in `out/variants/<name>/` and `out/versions/<name>/`.

## Design decisions the user has settled

### Crease and fissure shape

- **Closed fissure:** the anus is a closed fissure. It's a passage into the canal that can open, but no hole may be visible.
  - The lips touch, about 0.02 mm apart.
  - The canal opening is exactly the skin's final slit edge loop.
  - The canal stays collapsed for its first ~2 mm.
- **Short slit:** about 3.6 mm (`SLIT = 1.8` half-length), so the creases converge like a star. The user's reference
  image is a radial star of thin creases.
- **Creases are 3D geometry, not painted lines.** The painted-line strength node is 0.
- **Clean lines:** each crease is a clean line of constant width that runs all the way into the centre crease without
  flattening. No wobble.
- **Randomness is allowed only in:**
  - which spokes carry a crease
  - how long each runs
  - how deep it is

  Never mirrored left/right.
- **Crease layout:**
  - Front half (toward the perineum): approved as is, a crease every 20°.
  - Back half: the sides are offset by 10° so no crease runs straight up the cleft.
  - `_FRONT` / `_BACK` in v3+.
- **Profile across the creases:** the target is the user's sketch.
  - Each pad is a dome curving across its whole width.
  - Neighbouring domes meet in a sharp cusp one point wide.
  - Not flat pads with a U-groove, and not "ridges with grooves inside".

### Overall surface

- **Level:** no mound and no pit. The area around the slit should be nearly level with the surrounding skin.
  - A rise toward the slit was tried: 0.7 mm was too much, and 0.35 mm still too much.
  - The current state is no rise and half the lift.
- **Squeeze:** the squeezed cleft (`SQ = 0.35`) and the pigment riding slightly up the cheek walls were requested early on.
- **Smoothness:** no abrupt lines, rims or U shapes anywhere. The user zooms in closely.

## What we learned (avoid repeating)

### Toon shader and normal map

- **Hard two-tone cut:** the toon shader's active ramp (`Full range test`, a narrow 0.147–0.192 band) is effectively a
  two-tone cut, so small relief only shows where the shadow edge crosses it.
  - v3+ adds an "Anus crease shading (3D)" node group.
  - It darkens the toon result by the light ratio of the creased surface to the smooth surface underneath.
- **Body normal map:** it's 1024 px and was magnified ~15× at the anus, which drew texel spikes.
  - It's faded out over the cleft through RawShade's `Use Normals?` input.
  - The fade uses the `cleft_nm_off` attribute.

### Mesh and topology

- **Cleft refinement:** the low-poly game mesh made the toon shadow edge zigzag. The cleft is refined with each edge
  split into 4, and the new vertices are placed on a Phong surface.
- **Vertex snapping:** creases cut diagonally across the grid look jagged and beaded. On an unstructured mesh a crease
  falls between vertices and reads as a narrow trough.
  - **Next step:** force edges along every crease line in the CDT so each crease bottom is a single row of vertices,
    and never smooth across a crease.
- **Tips:**
  - Proximity welds and point merges at the tips create non-manifold fans.
  - Instead, open the full slit and weld the lips **pairwise** over the last 0.35 mm (`OPEN_MARGIN`).
- **Lift:** a lift that falls off evenly in every direction left a U-shaped ridge, because the old dimple is lopsided.
  Use the per-AP profile lift, smoothed over neighbouring vertices and applied along `n`.
- **Measure in real millimetres:** the old confocal layout's lateral offset jumps near the slit, so measure lateral
  distances via `ring_xy`.
- **Smoothing passes erase narrow grooves.** Measure the groove depth with `scripts/tools/gd.py` instead of judging by eye.

### Process

- **Frame the actual spot** in renders before claiming something is fixed. One "confirmed" render missed the spot entirely.
- **Variants:** the user wants changes made as new versions or files, so saved versions are never disturbed.

## Diagnostic tools (`scripts/tools/`)

Run each as `blender -b <file.blend> --python scripts/tools/X.py [-- args]`.

| Tool | What it does |
|---|---|
| `tip.py` / `tipw.py` | Skin-only geometry close-up of the slit. `tipw` adds a wireframe. Args: `-- <outdir> <prefix>`. In this view the image top is the front (perineum) end of the slit. |
| `diag2.py` | EEVEE shaded views. Modes: `close`, `zoom`, `wire`, `subsurf`. In the `close` view the image top is the back (coccyx) end. |
| `fronttip.py` | EEVEE close-up of the front tip, with `canal` / `nocanal` modes. Edit `tgt` to look at other points. |
| `chk.py` | Vertex, boundary and non-manifold counts. |
| `nm.py` | Lists non-manifold edges and the slit edge loops. |
| `canalout.py` | Canal vertices outside the skin. Ones that are "on" sit exactly on the skin edge and are fine. |
| `orient.py` | Ray-parity inside/outside test. The colon reads as partly "outside" because the Hips mesh is cut at the waist, which is expected. |
| `xs2.py` | Lateral cross-section through the anus. |
| `mid.py` | AP profile along the midline. |
| `gd.py` | Groove depth at various distances from the slit. |
| `hole.py` | Width of the slit hole along its length. |
| `inside.py` | Nearest-normal inside test. Unreliable at sharp edges. |
| `star.py` | Skin-only Workbench view of the whole crease star, looking along the outward normal. Args: `-- <outdir> <prefix> <distance_m> <wire 0/1>`. 0.022 frames the star, 0.014 with wire shows the crease rows, 0.007 frames the slit edge. In this view the image top is the back (coccyx) end. |
| `grazing.py` | Skin-only Workbench view low along the cleft from the perineum side, for judging rises and dips. Args: `-- <outdir> <prefix>`. |

`scripts/render_close.py` and `scripts/render_wire.py` produce the standard close / below / geometry renders.

## Current state and next step

- **Back-side slope, ADOPTED (latest):** the corridor fix is now ON by default in `build5_rise.py`
  (`ANUS_CORRIDOR=3000`, `ANUS_CORR_ANCHOR=0.03`, `ANUS_CORR_LAT=0.5`; the solver is from the corridor builder). It smooths
  the cleft-floor profile along the cleft from -10 to +28 mm, with the cheek walls held and the anus relief taken off
  and put back.
  - **Bare-surface metrics** (built with `ANUS_CREASE_SCALE=0`): max_slope_back 40.1 -> 21.4 and max_kink_back
    24.6 -> 9.7. Anus tilt is about the same (it sits on the rising floor).
  - **Metric caveat:** with creases on, backslope.py's slope and kink include the crease walls, so compare bare-surface
    builds.
  - **Renders:** `backslope_corridor_vs_current.png`.
  - **File for the user:** `out/blend/Hips_build5_rise_corridor.blend`.
- **Back-side slope, comparison run (stopped early to save credits):** the user confirmed the problem side is the
  back (coccyx, red marker). Both the back half's tilt and the steep floor behind it are involved. Metrics come from
  `scripts/tools/backslope.py`.
  - **Current build:** tilt_front 3.8, tilt_back 15.2, max_slope_back 40, max_kink_back 23.
  - **Region candidate** (`scripts/build5_rise_region.py`, `ANUS_BACK_REACH=20`): tilt_back 11.3, max_slope_back 26.5,
    max_kink_back 27. It adds a hard seam around the anus and a fold behind it
    (`backslope_region_candidate_vs_current.png`), so it is not adopted.
  - **Corridor candidate** (`scripts/build5_rise_corridor.py`): no improvement yet.
  - **Target-curve candidate:** never started.
  - **Default:** unchanged (`build5_rise.py` with levelling on).
  - **To resume:** `Workflow({scriptPath: <session workflows/scripts/back-slope-fix-*.js>, resumeFromRunId: "wf_64543cfe-4a4"})`.
    The run IDs are session-local, so in a new session rerun the approaches instead.
- **Cliff behind the anus:** the user saw a steep climb just behind (toward the back of) the anus.
  - **Measured** (`cliff_profiles.png`): the narrow fill sat low near the original pit, then climbed ~4 mm within
    2.5 mm at the back edge to meet the cleft floor. The smoothing pass also leaves that floor ~1.3 mm above the
    original.
  - **Fix:** `ANUS_LEVEL` (3000 iterations) smooths heights along the cleft only: edge weight (da / length)^2, mobility
    fading up the walls from 3 to 6 mm lateral, and the outer ring fixed.
  - **Tried and dropped:** a straight-line lift (left ridges at the 7.6 mm ring), and turning the smoothing pass off
    (`ANUS_TAUBIN=0`). Without that pass the original's low-poly facets and a hard crease show in front of the anus,
    so it stays on (default 20).
  - **Renders:** `cliff_lowview_*`.
- **Cheek distortion check:** `scripts/tools/distort.py` measures, for each of the original's 1,580 skin
  vertices, its distance to a version's final skin. The plot is `cheek_distortion_map.png`.
  - **Mean / max movement in mm:** at 12-20 mm from the anus, current 0.88 / 2.31, build7 1.13 / 2.73, build5 1.13 /
    2.80, build8 1.22 / 3.60. At 20-30 mm, every version is ~0.25-0.32 / 1.5-2.05.
  - **Beyond 45 mm:** zero for all versions.
  - **Cause:** nearly all the movement outside the rebuilt area is step 2, the Taubin smoothing of the cleft crease
    (in the user's build7 from the start). On its own it moves vertices up to 4 mm out to ~30 mm, mostly along the cleft
    toward the back.
- **build5_rise, glitch check:** the user picked `ANUS_SIG` 0.06 (now the default, no taper) and asked for a
  check of textures and shadows.
  - **Diagnosis** (renders with one change each): turning the body normal map off removed both glitches, the spiky
    shadow edge on the left cheek wall and the crease ends splitting into streaks. Resetting custom normals changed
    nothing.
  - **Fix:** ported build8's three cleft fixes: 2c refinement, custom normals reset over the refined region, and the
    `cleft_nm_off` attribute driving RawShade "Use Normals?".
  - `ANUS_RING_STEP` (extra outer rings) was a wrong guess for the streaks. It is kept as an option, off by default.
- **build5_rise, crease width:** the user felt the grooves took about as much area as the pads. They want the
  creases to keep widening with distance from the centre, just less.
  - **The user's way:** a smaller `ANUS_SIG` (build5: 0.09; tried 0.06 and 0.045).
  - **Alternative:** `ANUS_CREASE_TAPER` K shrinks the angular width as `(R_W / r)^K` past `R_W` = 1.5 mm (tried 0.4).
  - Both need finer spokes, so `N_SPOKE` is 288 (build5: 72), with 2:1 steps to 144 / 72 at 7.6 / 8.6 mm.
  - Only the midline spokes snap to x = 0 now (distance snapping would fold the dense centre).
  - Defaults: `ANUS_SIG` 0.09, taper 0.6, and `ANUS_DIP` now 1.5 (the file the user opened).
  - Comparison: `crease_width_geometry.png`.
- **build5_rise, plunge:** the user's sketch shows shoulders rising toward the centre that roll over smoothly into
  a narrow plunge. That plunge is `DIP * (1 - sqrt(r / R_DIP))^2` with `R_DIP` 2.5 mm, along n; set it with the
  `ANUS_DIP` env var (default 1.2). The build adds rings inside 0.9 mm and at 1.15 / 1.65 / 2.35 mm so it curves
  smoothly, and the canal still opens at the 0.9 mm ring. 0.8 and 1.5 were rendered, with cross-section plots
  (`profiles_plunge_0.8_1.5.png`) and 1600 px shaded renders (`scripts/tools/hq.py`).
  - **Render quality:** the softness came from my downscaled contact sheets and the 1024 px textures. The jagged crease
    edges in the shaded view come from build5's 72 spokes, which make each crease only ~2 vertices wide under the
    toon cut. More spokes would fix it but would slightly change the build5 look, so ask first.
- **build5_rise, crease depth:** `ANUS_CREASE_SCALE` (default 2.0, the user's pick of 1x / 1.5x / 2x) multiplies build5's crease relief. The creases were
  measured as deep as build5's, but they read faint once the anus follows the curved cleft (toon two-tone cut). The
  rise default is now 0.2 mm with no pole sink (`ENTRANCE_DIP = 0`). 1x, 1.5x and 2x were rendered at a 0.2 mm rise
  (`crease_depth_1x_1.5x_2x_rise0.2.png`).
- **build5_rise (earlier):** the user asked for the surface to rise toward the centre. Three heights were rendered side
  by side, and the user is choosing one. Note that "level" in build5_adapted already included the 2 mm funnel dive.

- **build5_adapted (latest):** the user rejected build9 and build10 as rebuilds rather than build5. They asked to go back
  to build5 and only fix its problems: the lump, the edges not riding up the cheeks, and the deformed cheeks.
  - **Changed:** only those fixes. Everything else (72 spokes, pole, `crease()`, `fold_amp`, `FUNNEL`, pigment, canal)
    is build5 verbatim.
  - **Open:** none of the later additions are in it (`SQ` cleft squeeze, cleft refinement, normal-map fade, crease
    shading node). With the narrow cut, the hole loop is only 14 coarse vertices bridged to the 36-vertex outer ring.

- **build10 (latest):** the user wanted build5's dive toward the centre and its "pinched together" look.
  - **Tried and rejected:** build5's numbers literally. The creases stopped at 0.9 mm and a cone pit sat inside, which
    on the dense mesh became a hard-rimmed oval crater.
  - **Now:** creases from the centre, a smooth pit, a stronger pinch, and the dive along the cleft axis.
  - **Open:** a small dark notch remains just left of the centre, and the user hasn't reviewed it yet.

- **build9 (latest):** the user asked to use build5 as the base and fix its lump, its separate-disc edge and its
  cheek deformation.
  - **Done:** the creases meet in one point, the pigment and rings climb the cheek walls, there is no funnel or lump,
    there are no non-manifold edges and no degenerate faces.
  - **Fixed along the way:**
    - Snapping by distance to the midline folded the dense centre rings. Only the two midline spokes snap now.
    - With `K_LAT = 0.85` the outer ring came so close to the hole that the joining strip creased (slivers on the
      cheeks). It is now 0.77.
  - **Open:**
    - The anus is closed at a point (build5 style) rather than the earlier 3.6 mm closed slit. Ask whether the slit
      should come back.
    - The crease layout is build5's 12 creases at 30 deg (with one straight up and one straight down the cleft), not the
      v3+ `_FRONT` / `_BACK` layout.
    - The pads on the steep left cheek wall read as raised lobes in the geometry view.

### Earlier (v5-v8)

- **v4:** saved and committed.
- **v5 (`build8_v5.py`):** has the domed-pad `crease()`. The last render the user saw still read as grooves, which led
  to the question of how to be sure the creases meet in a line a single point wide.
- **Answer, implemented in `build8_v5.py` (commit ad776fa):** crease edges are forced into the triangulation, and the
  build verifies that each crease bottom is one vertex wide:
  1. **`crease_rows`:** a row of vertices sits exactly on each crease line, from a lip vertex outwards. The rows are
     passed as constraint edges:
     `delaunay_2d_cdt(co2, crease_edges, [outer_poly, lip_poly], 2, 1e-7, True)`.
  2. **`side_pts`:** rows at ±0.6 and ±1.2 `MESH_H` run parallel to each crease. Lattice points within 1.6 `MESH_H` of
     a crease are removed.
  3. **Smoothing:** `RELIEF_SMOOTH = 0`, so nothing smooths across a crease.
  4. **Verification:** the build writes `crease_bottom_single_vertex` to `report.json`. The latest report says
     "682 of 682 crease vertices are the strict lowest point across", with `inner_mesh_verts` 23600,
     `slit_edge_loops` 1 and `canal_opening_points` 84.
- **v5 rendered** (`out/variants/v5/`: `star_geo`, `star_wire`, `slit_geo`, `v5_close_*`, `v5_below_*`, `v5_zoom`).
  A rebuild reproduced `report.json` exactly, and `hole.py` reads 0.000 mm along the whole slit.
  - **Good:** in `star_wire` every crease is a straight single row of vertices. In `star_geo` the pads are domes meeting
    in a sharp cusp, as in the user's sketch.
  - **Defect 1, spikes at the slit (`slit_geo`):** each crease ends at the slit in a small spike that juts into the slit.
    The crease row starts on a lip vertex, so that lip vertex gets pulled down while the pads either side stay up.
  - **Defect 2, abrupt outer ends (`star_geo`, left side):** some creases stop at their outer end with a step, because
    the pad on one side stands higher than the other.
  - **Defect 3, crowded tips:** where the creases converge at both slit tips, the triangles are small and irregular, and
    the surface looks lumpy.
- **v6 (`build8_v6.py`, renders in `out/variants/v6/`, `v5_vs_v6.png` side by side):**
  - **Requested:** one more crease at the top middle, i.e. straight up from the back tip at 0 deg. This overrides
    the earlier "no crease straight up the cleft" rule. It has a fixed depth and length, so `_rng` and every other
    crease stay as in v5.
  - **Fading fixed:** the user asked why creases fade out near the centre. The v4 tip tidy (40 normal-smoothing passes
    within 1.2 mm of each tip) flattened them. The lip curl's taper (last 0.7 mm) also left a lump on the right lip
    near the front tip. Both now use `TIP_R = 0.35` (the zipped `OPEN_MARGIN`). The report shows 725 of 725
    single-vertex crease bottoms, `slit_edge_loops` 1, hole 0.000 mm, and no non-manifold edges.
  - **Still open:** a small knot at each tip where the creases converge (`TIP_R = 0.5` was tried and was lumpier).
    The spikes where creases meet the slit (defect 1) and the steps at outer ends (defect 2) are not yet addressed.
- **Tip problem (open):** at each tip, three to five creases meet the slit within ~0.1-0.3 mm of one point. At full
  depth, the pads between them become thinner than the mesh (`MESH_H` ~0.07 mm) and collapse into a knot (v6). With
  the depth reduced, the creases visibly fade (v7). Turning off the tip tidy and raising `ASPECT` to 3 helped but did not
  fix it; the user confirmed that none of these is right.
- **Next:** the user is choosing between three options:
  1. (Recommended.) Spread where the tip creases meet the slit, at least ~0.3 mm apart, so only the 0/180 deg creases
     enter at the very tips.
  2. Join converging creases in Y-junctions before the tip.
  3. Refine the mesh at the tips. This risks reading as a pit.
- **Tool:** `scripts/tools/tipcam.py` renders a skin-only close-up at a point along the slit.
- **Smaller open items:**
  - A thin line runs from the slit's back tip along the cleft; not yet addressed.
  - The saved v1 still has a couple of non-manifold edges.
