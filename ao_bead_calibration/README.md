# AO bead calibration — scope2k acquisition

Acquisition code for the holographic adaptive-optics (AO) bead calibration
(NeAT, bioRxiv 2024.10.20.619284). Génesis provides the bead slide, the SLM
target pattern, and the analysis; **this folder is the scope2k side** — Phase 1
(SLM→back-pupil mapping parameters) and Phase 2 (the 20-spot grid z-stack on
beads, plus references and metadata). See the full protocol in
`~/Downloads/README_AO_bead_calibration.md`.

Detection is on the **Basler alignment camera**: the holographic grid is the
excitation, and the camera images its two-photon fluorescence directly, so each
spot *is* a PSF on the camera. Laser power is controlled over **Holochat by the
DAQ computer** (`alignment/alignCodeDAQ2K.m`), exactly like
`align_slm_to_camera_scope2k.m`: the hologram is gated **ON only during each
grab** and OFF in between, which gives a true beam-off background and limits
bleaching.

`alignCodeDAQ2K` reads the wavelength **once** at startup and then services power
until it receives `'end'` (which each script sends when it finishes). So
**restart `alignCodeDAQ2K` on the DAQ before each of the two scripts.**

## Files

| File | Role |
|---|---|
| `slm_pupil_mapping.m` | **Phase 1.** Measures rotation, mirror flip, pupil diameter (px), pupil center (px), and confirms the SLM→image affine. Leaves a `pupil_mapping` struct in the workspace. |
| `acquire_bead_grid_stack.m` | **Phase 2.** Feeds the fixed N-spot grid, sweeps the objective in z, grabs per plane, saves the raw stack + references + metadata. Auto-includes `pupil_mapping` if present. |
| `make_grid_coords.m` | Builds/loads the N×4 `[x y z power]` grid (Génesis's file, or a synthetic fallback). |
| `save_bead_stack.m` | Writes the float32 multipage TIFF + JSON/`.mat` metadata sidecar (no TIFF writer existed in `holography2k`). |

## Run order

1. **DAQ computer:** run `alignCodeDAQ2K`. It waits for the wavelength, then sets
   laser power on request and acks each with `'gotit'`. It exits when the script
   sends `'end'`, so **restart it before each of the two scripts below.**
2. **Holography computer:** run `slm_pupil_mapping.m` cell-by-cell on the bead
   slide. Its setup cell sends the wavelength to the DAQ and previews a central
   spot (gated on) so you can tune `pwr` to ~80% of camera max. Confirm the four
   printed parameters look sane; it saves a `pupil_mapping` `.mat` and sends
   `'end'`.
3. **Restart `alignCodeDAQ2K`**, then run `acquire_bead_grid_stack.m` cell-by-cell:
   - preview → find a bead field with beads well separated vs. the grid spacing;
   - **set-power / no-saturation cell:** raise `pwr` until the brightest bead is
     well under the 8-bit ceiling (255 DN) — the value is recorded in the metadata;
   - verify each grid spot lands on a bead (nudge the grid / coordinate file if
     not — spots off beads measure nothing);
   - acquire the stack, let the bleaching re-check pass; it loads the Phase-1
     summary from the `.mat`, saves, and sends `'end'`.
4. **Handoff** (see below).

Sutter settle times (`settle_first_s` = 3 s after the big initial jump,
`settle_step_s` = 0.3 s per step) are tuned to match
`alignment/alignSLMtoCam/align_slm_to_camera_scope2k.m`. If planes look blurred or
z-misregistered, the stage is being read before it settles — raise these, don't
lower them.

## Two rules that matter for the analysis

- **No saturation.** Both phase retrieval and deconvolution fail on clipped
  pixels. The no-saturation cell warns; do not proceed past a saturated frame.
- **No smoothing / no background subtraction / no clipping in the saved stack.**
  The saved TIFF is the raw frame-average, float32. `get_psf_no_power` smooths and
  background-subtracts *for its own FWHM display only*; we deliberately do not,
  because those steps corrupt the data Génesis deconvolves. Background and
  reference frames are saved separately for her to use as she sees fit.

## What to hand Génesis

Everything lands in `<data_root>/ao_bead_calibration/<yymmdd>/` under one stem:

- `*_grid_stack.tif` — the z-stack (float32, planes = TIFF directories, axis
  order Y,X,Z).
- `*_background.tif`, `*_reference_blank.tif` — beam-off background and blank-SLM
  reference (README item 14a). `*_reference_nolut.tif` if you captured the
  correction-disabled pass (item 14b — see note below).
- `*_pattern.mat` — the exact `gridCoords` projected.
- `*_metadata.json` / `*_metadata.mat` — z-step/range, pixel size (`pxPerMu`),
  wavelength, power, frames averaged, camera exposure, and the four Phase-1
  mapping parameters. **NA, objective magnification, and pupil-fill are not
  knowable from the rig code** — set `na`, `objective_mag`, `pupil_fill` at the
  top of `acquire_bead_grid_stack.m`, or any field left as `NEEDS_INPUT` in the
  JSON must be filled before handoff.

## Open questions to settle with Génesis (do not block acquisition)

- **Temporal focusing (README item 4).** RESOLVED — scope2k **has TF**, so the
  TF-derived theoretical PSF model applies as-is; no forward-model swap needed.
- **Coordinate frame.** Confirm whether her target list is in normalized-SLM
  coords (`gridOpts.frame='SLM'`) or ScanImage/FOV coords (`'SI'`, converted here
  via `function_SItoSLM` + `find_latest_calib`). `make_grid_coords` errors out
  rather than silently dropping any target outside the SLM range.
- **Existing static correction = the Meadowlark hardware LUT.** There is **no
  software base-phase/Zernike layer** in the CGH path (the wrap in
  `function_Make_3D_SHOT_Holos.m` adds only the target phase). The item-14b
  "correction-disabled" reference therefore means reloading the SLM with a
  **linear LUT** (`...\Blink OverDrive Plus\LUT Files\linear.LUT`), grabbing, and
  restoring the real LUT — a manual step noted in the acquisition script. When
  the correction is later *applied* (Phase 4), a new software base-phase hook
  will be needed at that wrap line; that is out of scope for image collection.
- **Zernike convention** on the returned correction: Noll vs OSA/ANSI, normalized
  or not (README item 17). Confirm it matches her analysis input spec.
