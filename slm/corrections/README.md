# SLM correction / pattern files

Phase masks for `load_slm_correction`, tracked here so the holography computer can
`git pull` them instead of receiving them by hand.

## Files

| File | Panel | Format | Wrapped RMS | Source |
|---|---|---|---|---|
| `slm_pattern_14umaFWHM.bmp` | 1024×1024 | 8-bit BMP, 0–255 DN = 2π | 0.202 waves | Génesis, AO bead calibration |
| `slm_pattern_30umaFWHM.bmp` | 1024×1024 | 8-bit BMP, 0–255 DN = 2π | 0.186 waves | Génesis, AO bead calibration |

Both are for the **900 nm** board. 8 bpp with an identity grayscale palette, so
`imread` returns the DN directly — which is what `load_slm_correction` auto-detects.

Measured, over the inscribed unit pupil, after unwrapping:

- peak-to-valley 1.47 and 1.38 waves; RMS 0.199 and 0.186 waves
- 99.4% described by Zernikes up to n = 4 (residual RMS 0.007 rad), so both were
  rendered from a low-order coefficient fit rather than measured pixel-wise
- the two are **mutually uncorrelated** (r = 0.013) and are not a ±φ pair
  (astigmatism and spherical flip sign between them, coma does not)

That last point is worth knowing before picking one: two estimates of the same
optical system's aberration would look nearly identical, and these do not. The names
refer to a spot's axial FWHM. If you need the authoritative Phase 4 correction,
confirm with Génesis which of these is it — AO protocol Phase 3 item 17 says the
deliverable is *Zernike coefficients with the convention stated*, and
`slm/zernike_phase.m` renders those directly, which keeps the convention auditable
instead of frozen into an image.

## Using one

By hand:

```matlab
slm   = get_slm(900);
Setup = function_loadparameters2();
Setup.SLM.correction = load_slm_correction( ...
    fullfile(holography2k_root, 'slm', 'corrections', 'slm_pattern_14umaFWHM.bmp'), ...
    'Size', [slm.Nx slm.Ny]);
```

From the rig file, point `rig.paths.slm_correction_dir` at this folder and name the
file on the 900 nm opto channel — see `docs/slm-wavefront-correction.md`, which also
carries the `holodaq` patch that adds the `slm_correction` field.

## Before trusting a result

- **Orientation** is transposed by default (`imread` gives (y,x); a hologram is
  indexed (x,y)) and **cannot be verified from a square file**. Feed a known tilt and
  check the spot moves as predicted.
- **Sign** is ambiguous by convention. AO Phase 4 item 19: apply +φ and −φ, keep
  whichever improves the spot.
- **Delivered power will shift.** The diffraction efficiencies the holo box returns
  were measured uncorrected, and `Sequence.calculate_power` divides the laser command
  by them. Re-run `calibrate_DE_powermeter`.

## Adding a file here

`.gitattributes` pins `*.bmp` as `binary`. Do not remove that: the repo-wide
`* text=auto` rule lets git guess, and a misdetected mask gets LF-normalised on
checkout — a silently corrupted wavefront rather than an error. Verify a new file
survives the round trip:

```bash
git show HEAD:slm/corrections/<file>.bmp | shasum -a 256
shasum -a 256 <original>.bmp
```
