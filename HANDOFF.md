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

- **build5_rise, glitch check (latest):** the user picked `ANUS_SIG` 0.06 (now the default, no taper) and asked for a
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
