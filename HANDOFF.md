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
| `scripts/build7.py` | original | The local chat's last build (v7): a near-circular pucker. |
| `scripts/build8.py` | **saved v1** (`out/versions/closed_fissure_v1/`) | Squeezed cleft, closed fissure, lips peaking at the edge, cusped wrinkles. Uses the older confocal ring layout. |
| `scripts/build8_raised.py` | variant | v1 with a ~1 mm raised dome. |
| `scripts/build8_painted_lines.py` | variant | Crisp painted crease lines from a per-pixel mask texture. The user rejected painted lines and wants 3D geometry. |
| `scripts/build8_v2.py` | superseded | Straight wrinkles running into the crease. Beaded and jagged. |
| `scripts/build8_v2_wip.py` | parked | A failed experiment. |
| `scripts/build8_v3.py` | superseded | Straight grid "spokes" with creases on them, a short slit (3.6 mm), a closed slit, clean tip topology. Has a U-shaped artifact and sliver faces at the tips. |
| `scripts/build8_v4.py` | **saved v4** (`out/versions/even_mesh_v4/`) | **Current best.** The inner pucker is an even ~0.07 mm triangle mesh made by constrained Delaunay (`mathutils.geometry.delaunay_2d_cdt`). The lip ring plus four outer rings keep the join to the skin. The surface is a full membrane relax, plus a shaped lift that is smoothed and limited to 2.5 mm of the midline. `PUFF = 0`, creases on 20° spokes. The user approved its overall shape. |
| `scripts/build8_v5.py` | **in progress** | Like v4, but the pads between creases are domes meeting at the creases (`crease()` rewritten). The user still sees grooves. See next steps. |

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

`scripts/render_close.py` and `scripts/render_wire.py` produce the standard close / below / geometry renders.

## Current state and next step

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
- **Next:** render v5 with `render_close.py` and `tipw.py`, framing the creases. Show the user, then save it as a version
  if they approve.
- **Smaller open items:**
  - A thin line runs from the slit's back tip along the cleft; not yet addressed.
  - The saved v1 still has a couple of non-manifold edges.
