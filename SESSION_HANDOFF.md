# Session handoff - female pelvis (urinary + reproductive), 2026-10-08

Read this first, then `HANDOFF.md` for the full technical record (every step, parameter, measurement and what
was tried and dropped).

## Where the work is
- Repo: `Eclectic-png/For-Claude`, branch **`claude/optimistic-gates-jorged`** (latest commit: "Anus tilt default
  3 deg"). Everything is committed and pushed. To continue on another branch, merge or cherry-pick from it;
  the work is all in `scripts/` + `HANDOFF.md` (+ renders / report in `out/urinary/`).
- Blends are gitignored. Rebuild locally (Blender 4.2.3, drive files in `drive/`):
  ```
  cd drive
  ANATOMY_REF=$PWD/anatomy_ref.blend blender -b Hips.blend --python ../scripts/build5_passage.py -- \
      $PWD/../out/blend/Hips_build5_passage.blend <report dir>                  # anus / tract, ~25 s
  ANATOMY_REF=$PWD/anatomy_ref.blend blender -b --python-use-system-env ../out/blend/Hips_build5_passage.blend \
      --python ../scripts/build_pelvis.py -- $PWD/../out/blend/Hips_urinary.blend ../out/urinary   # ~27 min
  blender -b ../out/blend/Hips_urinary.blend --python ../scripts/tools/urinary_views.py -- ../out/urinary
  ```
  (`--python-use-system-env` + `PYTHONUNBUFFERED=1` to see the log live; `URINARY_STOP=<step>` saves early:
  pelvis, anus, kidneys, vagina, repro, ureters, bladder_raw, urethra, keys.)

## What is done (female model)
1. Urinary system: kidneys / ureters / adrenals (atlas), new fillable bladder (70 -> ~470 ml), ureter openings,
   urethra to a movable meatus cut into the vestibule, Void.
2. Female bony pelvis (the atlas pelvis is male; reshaped by one smooth field).
3. Rectum filling and both-full correctives.
4. Internal reproductive organs, first pass: uterus (outer + canal / cavity), tubes, ovaries, ovarian / round /
   suspensory ligaments; cervix opens into the vault of `AN_Vagina_Space` (a hidden placeholder tube, NOT a real
   vagina). The uterus is firm: it tilts on its cervix (22 deg back at full bladder); neighbours give way.
5. Perineum: anus tilted **3 deg** about the back end of its opening (anchored on the skin), the mound + incline
   in front replaced by an even ramp (`urinary_lib.anus_tilt_ramp`, `ANUS_TILT_DEG`).
6. Lighting (lit insides + AO bakes), polish pass, checks over 14 states.
Controls: select **`Pelvic_Controls`** -> Custom Properties: **Bladder_Fill**, **Void**, **Rectum_Fill** (0..1).

## State of the last build
- Seams (bladder neck, meatus, ureter openings) <= 0.0005 mm in every state; rest / Void states clip-free.
- Other states <= ~1.4 mm, a few isolated spots up to ~3.5 mm (rectosigmoid seam pressed into the uterus).
- Tubes / ovaries touch the hip bone by <= ~0.9 mm in some states.

## Update 2026-10-28
- **Vagina, first pass: DONE** as a separate, swappable step: `scripts/build_vagina.py` on build_pelvis' output ->
  `out/blend/Hips_vagina.blend` (~1 min; see HANDOFF.md "Vagina"). `AN_Vagina` (collapsed lumen: introitus slit ->
  H -> transverse slit -> vault round the cervix, os opening into it) + `AN_Vagina_Wall`; one cut in the skin. The
  user has not said whether their full-body model already has a vaginal canal - if it does, delete these two objects
  (and the introitus cut) when merging.
- **Vagina, accuracy pass: DONE** - portio ~10 mm into the canal (was ~2.5), posterior fornix (~5 mm past the
  cervix), squarer lower walls, 2.5 mm rectovaginal septum, rugae. Needs build_pelvis re-run (done; ~13 min) then
  build_vagina. Details and what was tried in HANDOFF.md "Vagina".
- Anal canal: now held shut its whole length, opening only over the last 5 mm into the rectum
  (`ANUS_CANAL_OPENLEN`, build5_passage); build_pelvis was re-run on it (report updated).

## Next steps (agreed order)
1. ~~The vagina~~ (first pass + accuracy pass done; limits in HANDOFF.md: short H recesses).
2. **The vulva** rebuild (rudimentary now); then re-cut the urethral meatus (`URETHRA_MEATUS` /
   `URETHRA_MEATUS_FROM`) and bring the clitoris closer (meatus is ~43 mm behind the hood now, real ~25-30).
3. Smaller: the lower ampulla's folds where it is pushed back for the vagina (160 self-intersecting face pairs;
   relax / remesh), broad ligament / peritoneum. (Septum done: 2.5 mm.)
4. Later: the male version.

## Working notes (user preferences learned)
- Show before/after sections and close renders before committing to a full 27-minute build; the user judges
  shapes from the midline section and the view from below.
- Keep the user's approved anus look (creases, cleft) intact; change only what was asked.
- Commit and push after each step (a stop hook enforces it); GitHub pushes sometimes fail with a transient 500 -
  retry after a minute.
- Body scale ~0.72x real (atlas organs are real size under `Internal_Fit_Xform`); report real-size numbers.
