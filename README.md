# holography2k

Computer-generated holography for a 2-photon holographic stimulation rig: compile
holograms, drive the SLM, align the SLM to the imaging field, and calibrate
photostim power. 142 MATLAB files.

This code runs on the **holography computer** — one of the four machines in a full
rig (see [`holodaq`](https://github.com/adesnik-lab/holodaq)'s README §2). It is
not the machine that runs experiments: the DAQ computer owns the rig file and the
NI card and drives this box over a digital trigger plus the holochat broker.

| repo | owns |
|---|---|
| [`holodaq`](https://github.com/adesnik-lab/holodaq) | hardware layer: DAQ, module classes, `TrialManager`, holochat, and the `rigs/` channel map |
| [`holoexpt`](https://github.com/adesnik-lab/holoexpt) | the experiment lifecycle, launcher GUI, satellite priming |
| **`holography2k`** (this repo) | SLM, hologram compilation, alignment, power calibration |
| `scope2k-experiments` | one lab's experiment definitions |

Dependency direction is one-way: this repo uses holodaq, never the reverse.

## Getting on the path — read this first

Every hand-run script here starts with `makePaths()` or `holo_paths()`. Both need
**holodaq checked out beside this repo**, and both *error* rather than continue if
it is not:

```
Documents/GitHub/
├── holography2k/     <- this repo
└── holodaq/          <- must be here
```

```bash
git clone https://github.com/adesnik-lab/holodaq
```

Resolution is a sibling directory named `holodaq`, then the only holodaq-shaped
sibling, then an error with instructions.

**There is no environment variable for this, on purpose.** A `HOLODAQ_HOME` used
to come first and was removed: an env var is invisible in the checkout,
machine-global, and stays silently wrong the moment it outlives the path it names
— and pointing at a stale holodaq means a stale `rigs/` layer, which is the wrong
physical channel map. A sibling directory you can see beats a variable you have to
remember you set.

The error is deliberate. These scripts used to open a hardcoded trio of
`addpath(genpath('C:\Users\holos\...'))` naming one machine's username. On any
other computer `genpath` of a missing folder returns `''` and `addpath('')` is a
silent no-op, so the paths simply were not there and the first SLM or camera call
died somewhere unrelated. Failing at the path step is much cheaper to diagnose.

`holo_paths` is the modern entry point — this checkout plus holodaq, appended
rather than prepended so a stale tree cannot shadow the copy you are editing.
`makePaths` wraps it for the ~13 existing scripts that call it, and additionally
resolves msocket and the SLM vendor SDK. New code should call `holo_paths`
directly and add only what it needs.

## Configuration comes from the rig, not from files here

This machine has **no rig file of its own**. It reads the config the DAQ published
over holochat, through `rig_remote_get('<dotted.path>', <fallback>)`. Every literal
in this repo is a last-resort fallback, not a default — so pointing this code at a
different scope means editing that scope's rig file, not editing anything here.

| what this repo reads | used for |
|---|---|
| `paths.calib_dir` | spatial (SLM↔camera) calibrations — read *and* written |
| `paths.slm_sdk` | SLM vendor SDK, added to the path at listener start |
| `paths.slm_lut_dir` | base for a relative `slm_lut` on an opto channel |
| `paths.holo_scratch` | holeburn `.tif` staging for the live alignment flow |
| `paths.power_calib_dir` | where the `AutoLaserPowerCalib_*` scripts write LUTs |
| `paths.holo_request` | the `holoRequest.mat` folder |
| `paths.si_root` | the ScanImage box's drive, for tiff staging during alignment |
| `holo.cgh_method`, `holo.use_gpu`, `holo.slm_timeout_ms` | hologram compilation |
| `opto` | the rig's ordered opto-channel table |
| `serial.sutter.port` | the MP285 manipulator, which hangs off **this** machine |

If the satellites appear to be running on stale values, someone edited the rig
file and did not re-run `publish_rig_config()` on the DAQ.

**Calibrations are resolved by date, not by filename.** `find_latest_calib(nm)`
picks the newest `*_Calib_<nm>*.mat` in `paths.calib_dir`. It replaced a
hand-maintained `switch` over dated filenames that had to be re-pasted every time
anyone recalibrated.

## Layout

```
holorequest/     start_holo_listener (the live listener), HoloListener (its loop)
                 and HoloListenerMonitor (its status window), the holorequest
                 entry points, find_latest_calib, PlaySequence2K
alignment/       SLM-to-camera alignment, PSF measurement, CoC fitting
  alignSLMtoCam/   align_slm_to_camera_scope2k -> process_holeburns_claude
cgh/             hologram algorithms (Gerchberg-Saxton and friends)
slm/             SLM start/stop, LUT calibration, the SLM class
cameras/         ThorCam and Basler acquisition (vendor DLLs vendored in-tree)
roi_class/       Pattern, Sequence
sutter/          MP285 manipulator control
manual-use/      manual_slm, manual_power_control_setup
tests/           bench checks: SLM range, zero-order location, ScanImage bounds,
                 test_holo_listener_monitor (offline: no rig, no broker, no SLM)
```

## Entry points

| run this | for |
|---|---|
| `holorequest/start_holo_listener` | the live listener during an experiment. Prints the config it resolved — read those lines and confirm they are your scope's. Async by default, with a **Holo Listener** status window (lamp + Start/Stop). Closing that window stops the listener; `'nogui'` skips it, `'blocking'` is the old Ctrl-C loop |
| `holorequest/MsocketHolorequest2K` / `SingleHolorequest2K` | hand-run hologram compilation |
| `alignment/alignSLMtoCam/align_slm_to_camera_scope2k` | the current SLM↔camera alignment |
| `AutoLaserPowerCalib_HWP` / `_EOM` / `_MOD` | photostim power calibration. `_EOM` is the one for a rig that sets power with a modulator rather than a rotator |
| `calibrate_DE_powermeter` | diffraction-efficiency calibration |
| `getImagingPSF` | PSF measurement |

`AutoLaserPowerCalib_HWP` needs the ELL14 bus, which `ScopeController` /
`PowerControllerCalibrated` hold open on the DAQ. **Close those GUIs first** — the
bus is exclusive.

## Known rough edges

- **`alignment/alignCodeDAQ2K.m` is broken and still hardcoded.** It names DAQ
  channels and a dated power-calibration path directly, and it passes a 4th
  positional argument that `LaserPowerControl` stopped accepting — deliberately
  left failing loudly, because making it "work" positionally would reinterpret
  `1250` as a modulator resting voltage. To be integrated with the other
  calibration scripts later.
- **`alignment/alignment-refactor/optotune/AlignOptotune.m` does not parse.**
  Pre-existing, untouched since 2023.
- A superseded alignment chain was deleted on 2026-08-03:
  `alignSLMtoCamMultiTarg2K{,2D}` → `process_holeburns{,_607,2d}`, plus
  `dirtyAlignment2D`, `getPSF` and `alignment_tools`. The live chain is
  `align_slm_to_camera_scope2k` → `process_holeburns_claude`. If you were using one
  of the deleted ones, it is in git history and in
  `_git-backups/holography2k-pre-prune-20260803.bundle`.
- `Holomaker.getLaserEOM` calls `get_the_laser_eom_pulse_stuff()`, which does not
  exist. `modules/LaserEOM.m` in holodaq is likewise a stub with no methods. If you
  want an analog modulator, the working piece is holodaq's `LaserModulator`
  component, reached by declaring `kind = 'eom'` on an `fpc_*` module.
- The repo name says `2k`, but the config layer no longer does: everything above
  comes from whichever rig published it.
