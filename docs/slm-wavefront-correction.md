# Static SLM wavefront correction

A per-board phase map added to **every** hologram before the phase is wrapped to the
panel's 8-bit frame. This is Phase 4 item 18 of the AO bead protocol
(`ao_bead_calibration/README.md`), and it is independent of the Meadowlark hardware
LUT — the LUT is a per-pixel voltage response curve, not a wavefront. Both apply.

## Where it enters

`cgh/phase_to_frame.m` is the single place phase becomes a frame. It replaced this
line, which had been copy-pasted into four CGH functions:

```matlab
Hologram = uint8(floor((Setup.SLM.pixelmax*(Holo.phase+pi)/(2*pi))));
```

Those four — `function_Make_3D_SHOT_Holos`, `_disks_KCZ`, `_tuneZ`, and
`function_makeAB` — are upstream of every hologram this repo produces. That is
deliberate: a correction applied to stimulation but *not* to
`calibrate_DE_powermeter` would leave the power model calibrated against a different
optical system than the one delivering light.

With `Setup.SLM.correction` empty, the output is **byte-identical** to the old line,
so this can sit on the rig ahead of any mask.

`function_VolumeIntensity` is still fed the raw, uncorrected `Holo.phase`. It
predicts the focus for the phase you *commanded*; correcting it would show a
simulated aberration that the real optics is supposed to be cancelling.

## Using it by hand (scripts, alignment, PSF, AO acquisition)

```matlab
slm   = get_slm(900);
Setup = slm.add_slm(function_loadparameters2());
Setup.SLM.correction = load_slm_correction( ...
    'C:\slm\corrections\ao_900nm.bmp', 'Size', [slm.Nx slm.Ny]);

[Holo, ~, ~] = function_Make_3D_SHOT_Holos(Setup, [0.5 0.5 0 1]);
slm.feed(Holo);
```

`load_slm_correction` accepts a Meadowlark-style image (`.bmp`, `.png`, `.tif`),
a `.mat`, a numeric array, or Zernike coefficients. It prints the units it detected,
the orientation it applied, and the RMS / peak-to-valley in waves.

`slm.blank()` still feeds true zeros. That is what AO item 14b's reference pass
needs, and it is what "SLM off" should mean.

## Enabling it for experiments (the listener path)

`HoloListener` loads the correction once per board at startup, from the rig's opto
table, and injects it per channel at compile time. Two boards therefore get two
different corrections — a shared `Setup` would have silently given both the same one,
decided by whichever channel compiled last.

The holography2k side is done and inert until a rig file names a file. **The rig-file
half lives in `holodaq` and is not applied**, because that checkout had unrelated
in-flight work in `rigs/Scope2KRig.m` when this was written. Apply it there:

### 1. `holodaq/rigs/opto_channel.m`

Add two optional parameters next to `slm_board` / `slm_lut`:

```matlab
    p.addParameter('slm_correction', '');       % '' = run uncorrected
    p.addParameter('slm_correction_opts', struct());
```

and two fields in the returned struct, **after** `slm_lut` and before `label` (field
name *and* order are fixed by that constructor — MATLAB refuses to concatenate
structs whose fields differ in either, so every entry must come out of it):

```matlab
        'slm_correction',      local_char(o.slm_correction, 'slm_correction'), ...
        'slm_correction_opts', o.slm_correction_opts, ...
```

Do **not** add either to `opto_signature`. The signature is `name@wavelength#board`
on both sides and gates `Experiment.confirm_opto_agreement`; the correction is a
holography-computer concern, and adding it would invalidate every existing agreement.

### 2. `holodaq/rigs/Scope2KRig.m`

A base folder for relative paths, next to `rig.paths.slm_lut_dir`:

```matlab
rig.paths.slm_correction_dir = ['C:\Program Files\Meadowlark Optics\' ...
                                'Blink OverDrive Plus\WFC Files'];
```

and the declaration on the 900 nm channel:

```matlab
rig.opto = [ opto_channel('red',  1100, 'fpc_1100', 'slm_1100'), ...
             opto_channel('blue',  900, 'fpc_900',  'slm_900', ...
                          'slm_correction', 'ao_900nm.bmp', ...
                          'slm_correction_opts', struct('Sign', 1)) ];
```

### 3. On the DAQ

Run `publish_rig_config()`, then restart the holo listener. Its startup banner prints
the correction file and its RMS in waves for each board.

## Settling sign and orientation — the part code cannot do

Two conventions are genuinely ambiguous and must be measured, not reasoned about.

**Sign.** Phase 4 item 19: apply `+φ` and `−φ`, re-image, keep whichever improves the
spot. Set `struct('Sign', -1)` in `slm_correction_opts` to flip it. Two-photon signal
is quadratic in intensity, so peak brightness is the most sensitive readout.

**Orientation.** `imread` returns `(row, column) = (y, x)`; a hologram here is indexed
`(x, y)`, because `function_Make_3D_SHOT_Holos` builds `Masks(sxx, syy)` and hands
that array straight to `Write_image`. An image source is therefore **transposed by
default**. On a 1024×1024 panel no size check can catch a wrong transpose or flip, so
verify it:

1. Feed a known tilt *as the correction* — one full turn across the panel:
   `zernike_phase(1024, 1024, 2*pi, 'Modes', [1 1])`.
2. Grab. The spot must move along the axis you expect, by the amount a one-period
   ramp predicts.
3. If it moves along the other axis, set `struct('Transpose', false)`; if it moves the
   wrong way along the right axis, set `struct('FlipX', true)` or `('FlipY', true)`.

Record the answer in the rig file rather than in a re-exported image — a config field
can be diffed, a baked-in orientation cannot.

## Expect delivered power to shift

The diffraction efficiencies the holography computer returns to the DAQ come from a
CoC calibration measured **without** correction, and `Sequence.calculate_power`
divides the laser command by that DE. A correction that tightens the focus puts more
power on target for the same command.

Re-run `calibrate_DE_powermeter` after enabling it, or at minimum measure the change
before trusting a calibrated experiment (`activation_calibrated`,
`bidirectional_calibrated`, `power_matrix_multiple_power`, `power_screen_graded`).

## Units, and the two traps worth knowing

A Meadowlark-style file stores DN with `0..pixelmax` spanning 2π, and
`load_slm_correction` detects that automatically (integer-valued, within range). It
**refuses to guess** between radians and waves for small non-integer data, because a
wrong guess scales the whole correction by 2π — pass `'Units'` explicitly there.

If you ever add a correction outside this hook, note that the obvious approach —
adding to an already-compiled `uint8` stack — is wrong twice over:

- MATLAB `uint8` arithmetic **saturates**: `uint8(200) + uint8(100)` is 255, not 44.
  Every pixel whose sum passes the top clips into a flat maximum-phase region, which
  looks like a diffuse haze rather than an error.
- the DN modulus is `pixelmax` (255), **not 256**, because the forward map puts 2π
  across 255 counts.

Adding in radians before the wrap, which is what `phase_to_frame` does, avoids both.

## Tests

Both run offline, no rig and no SLM:

```matlab
addpath(fullfile(pwd,'cgh'), fullfile(pwd,'slm'), fullfile(pwd,'tests'))
test_phase_to_frame     % byte-identity without a correction; wrap vs saturate
test_slm_correction     % units, orientation, Noll/OSA index tables
```
