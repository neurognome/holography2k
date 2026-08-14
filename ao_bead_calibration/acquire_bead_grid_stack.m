%% acquire_bead_grid_stack.m
% Phase 2 of the holographic AO bead-calibration protocol (see the README in
% this folder, and ~/Downloads/README_AO_bead_calibration.md).
%
% Projects a fixed N-spot hologram grid onto a dense bead slide and acquires a
% z-stack by TRANSLATING THE OBJECTIVE (Sutter), imaging the two-photon
% fluorescence of the grid on the Basler alignment camera. The SLM phase is held
% FIXED for the entire sweep -- the SLM phase is the quantity being measured, so
% defocus must come from the stage, never from the SLM (README item 13).
%
% Modeled on get_psf_no_power.m: laser power is set MANUALLY by the operator (no
% DAQ power server / no msocket gating) and left on for the acquisition.
% Differences from get_psf_no_power: N-spot grid instead of one spot, SLM fed
% once and never touched again, and the RAW averaged stack is saved to TIFF (no
% smoothing / no background subtraction / no clipping) for Genesis's pipeline.
%
% RUN ORDER (see README): set the laser power by hand (shutter/rotator/EOM as
% usual), then run this script cell-by-cell on the holography computer.
%
% Hand-run script (clears the workspace), matching the get_psf_no_power.m style.

clear
close all
clc

tBegin = tic;
disp('Setting up AO bead-grid acquisition...');

makePaths()

%% ---- user parameters --------------------------------------------------------
wavelength      = 1030;              % SLM wavelength (900 / 1030 / 1100)
mouse_or_slide  = 'beadslide';       % label for the output stem
epoch           = 1;

% Grid coordinate source (see make_grid_coords.m). [] -> synthetic 20-pt grid
% until Genesis's coordinate file is in hand. When it arrives:
%   gridSource = '/path/to/genesis_targets.csv';  gridOpts.frame = 'SLM'; % or 'SI'
gridSource      = [];
gridOpts        = struct('frame','SLM', 'n',20);

% z-sweep: holographic spots have long axial tails -> wide range, fine step.
UZ              = linspace(-30, 30, 61);   % um about focus, ~1 um step (README item 13)
nframesCapture  = 10;                       % frames averaged per plane

% Power is set MANUALLY (no DAQ control here). Record the mW you dialed in for
% the metadata handoff -- leave [] and it is flagged NEEDS_INPUT.
pwr_mW          = [];               % <- fill in the power you set by hand

% Objective / analysis metadata that is NOT knowable from the rig code -- fill
% these in for the handoff (README item 3 / 15). Left here so they are obvious.
na              = [];               % objective NA            e.g. 0.8
objective_mag   = [];               % objective magnification e.g. 16
pupil_fill      = [];               % pupil-fill fraction     e.g. 0.9

%% ---- hardware setup (mirrors get_psf_no_power.m:1-33) -----------------------
Setup = function_loadparameters2();
Setup.CGHMethod = 2;                 % Global GS, as in get_psf_no_power
Setup.GSoffset  = 0;
Setup.verbose   = 0;

if Setup.useGPU
    disp('Getting gpu...');
    g = gpuDevice;
end

slm = get_slm(wavelength);
slm.stop();
slm.wait_for_trigger = 0;
slm.start();

sutter = sutterController();

bas = bascam();
bas.start()

disp('Setup complete.')

%% ---- 1p preview: find the objective / bead field ----------------------------
bas.preview()

%% ---- build & feed the grid hologram (ONCE) ----------------------------------
gridCoords = make_grid_coords(gridSource, wavelength, gridOpts);
nSpots = size(gridCoords, 1);

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, gridCoords);
slm.feed(Holo);
fprintf('Fed %d-spot grid hologram to SLM. It will NOT be changed again.\n', nSpots);

% Reconstruction sanity view (single_slm_patch.m:49): confirm N spots + layout.
figure('Name','Grid reconstruction (FFT of fed hologram)');
imagesc(abs(fftshift(fft2(exp(1i*double(Holo)/255*2*pi)))));
axis image; colorbar; title(sprintf('Expect %d spots', nSpots));

%% ---- set power by hand + no-saturation check (README item 11) ---------------
% Dial the laser power up by hand (shutter/rotator/EOM as usual) while previewing
% the grid, until the brightest bead is ~80% of camera max and NO raw pixel
% clips at the 8-bit ceiling (255). This power is used for the whole stack.
bas.preview()   % close the preview window when the power looks right

% Snapshot at the working power and verify no saturation before acquiring.
raw = bas.grab(nframesCapture);
rawmax = max(raw(:));
frac_sat = mean(raw(:) >= bas.camMax);
fprintf('Raw camera max = %d / %d (%.3f%% of pixels saturated).\n', ...
    rawmax, bas.camMax, 100*frac_sat);
if rawmax >= bas.camMax
    warning('acquire_bead_grid:saturated', ...
        ['Camera is SATURATING at the current power (%d DN). Lower the laser ' ...
         'power and re-run this cell before acquiring -- clipped data breaks ' ...
         'phase retrieval and deconvolution.'], rawmax);
end

figure('Name','Working-power frame'); clf
imagesc(mean(double(raw),3)); axis image; colorbar
title(sprintf('max = %d DN  (ceiling %d)', rawmax, bas.camMax));

%% ---- background (beam left on, as in get_psf_no_power) ----------------------
% NOTE: with manual power there is no beam-off background. This is a same-state
% reference grab; it is saved for Genesis but NOT subtracted from the raw stack.
bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

%% ---- acquire the z-stack (RAW; SLM fixed) -----------------------------------
% Store the raw frame-averaged plane. No imgaussfilt, no bgd subtraction, no
% clipping in the SAVED array (a smoothed display copy is fine).
sutter.setRef()
sz = size(bgd);
dataUZ = zeros([sz numel(UZ)]);

figure('Name','Acquiring z-stack'); clf
disp('Collecting bead-grid z-stack.')

for i = 1:numel(UZ)
    fprintf('Plane %d/%d (z = %+.1f um)\n', i, numel(UZ), UZ(i));
    sutter.moveZ(UZ(i))

    if i == 1, pause(1); else, pause(0.1); end

    data = bas.grab(nframesCapture);
    dataUZ(:,:,i) = mean(double(data), 3);   % RAW averaged plane

    % live view only (does not touch the saved array)
    subplot(1,2,1)
    imagesc(imgaussfilt(dataUZ(:,:,i), 2)); axis image; colorbar
    title(sprintf('Plane %d (z=%+.1f)', i, UZ(i)))
    subplot(1,2,2)
    imagesc(max(dataUZ, [], 3)); axis image; colorbar; title('Max projection')
    drawnow
end

sutter.moveToRef()
disp('Done collecting stack.')

%% ---- bleaching check (README item 11) ---------------------------------------
% Re-grab the first plane; if the field is now dimmer at equivalent focus the
% axial profile is corrupted -> retake at lower power. (Beam is on throughout,
% so this measures real bleaching of the beads.)
sutter.moveZ(UZ(1)); pause(0.5);
recheck = mean(double(bas.grab(nframesCapture)), 3);
sutter.moveToRef()

first_signal  = sum(dataUZ(:,:,1) - mean(bgd(:)), 'all');
recheck_signal = sum(recheck - mean(bgd(:)), 'all');
bleach_ratio = recheck_signal / max(first_signal, eps);
fprintf('Bleaching check: end/start signal ratio = %.3f\n', bleach_ratio);
if bleach_ratio < 0.9
    warning('acquire_bead_grid:bleaching', ...
        ['Field dimmed to %.0f%% of its start value over the stack. Retake at ' ...
         'lower power -- the axial profile is bleach-corrupted.'], 100*bleach_ratio);
end

%% ---- pixel size: pxPerMu via a 50 um Sutter move (get_psf_no_power:163-192) --
muUsed = 50;
disp('Determining pxPerMu...')

sutter.moveTo([0 0 0]); pause(1);
p1 = mean(double(bas.grab(nframesCapture)), 3);

sutter.moveTo([0 muUsed 0]); pause(1);
p2 = mean(double(bas.grab(nframesCapture)), 3);
sutter.moveToRef()

[x1, y1] = function_findcenter(max(p1 - bgd, 0));
[x2, y2] = function_findcenter(max(p2 - bgd, 0));
pxPerMu = pdist([x1 y1; x2 y2]) / muUsed;
pixel_size_um = 1 / pxPerMu;
fprintf('pxPerMu = %.3f  ->  pixel size = %.3f um\n', pxPerMu, pixel_size_um);

%% ---- reference images (README item 14) --------------------------------------
refs = struct();

% (a) blank-SLM raster/frame of the same field, no hologram.
slm.blank();   % feeds zeros(Nx,Ny); dimension-safe
refs.reference_blank = mean(double(bas.grab(nframesCapture)), 3);
slm.feed(Holo);   % restore the grid

% (b) "existing static correction disabled" reference.
% On scope2k the only existing correction is the Meadowlark hardware LUT loaded
% by slm.start() (there is NO software base-phase layer in the CGH path). To
% capture an uncorrected reference, reload the SLM with a LINEAR LUT and re-grab.
% This is a MANUAL step -- reloading the LUT mid-session depends on the Blink
% SDK. To do it: stop the SLM, point Setup.SLM.lut_file at the linear LUT, e.g.
%   'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\linear.LUT'
% restart, feed(Holo), grab, then restore the real LUT. See README, item 14b.
refs.reference_nolut = [];   % populate manually if you capture this pass

%% ---- save (RAW stack + refs + metadata) -------------------------------------
stamp = datestr(now, 'yymmdd');                     %#ok<TNOW1,DATST>
stem  = sprintf('%s_%s_e%d_%dnm', stamp, mouse_or_slide, epoch, wavelength);

try
    data_root = rig_remote_get('paths.data_root', pwd);
catch
    data_root = pwd;
end
outdir = fullfile(data_root, 'ao_bead_calibration', stamp);

% Pattern coordinate file, saved next to the stack so the exact grid is recorded.
pattern_file = fullfile(outdir, [stem '_pattern.mat']);
if ~isfolder(outdir), mkdir(outdir); end
save(pattern_file, 'gridCoords', 'wavelength');

meta = struct();
meta.z_step_um        = mean(diff(UZ));
meta.z_range_um       = [min(UZ) max(UZ)];
meta.z_planes_um      = UZ;
meta.pxPerMu          = pxPerMu;
meta.pixel_size_um    = pixel_size_um;
meta.wavelength_nm    = wavelength;
meta.power_mW         = pwr_mW;      % operator-recorded (set by hand)
meta.frames_averaged  = nframesCapture;
meta.camera_exposure  = bas.exposure;
meta.na               = na;
meta.objective_mag    = objective_mag;
meta.pupil_fill_fraction = pupil_fill;
meta.n_spots          = nSpots;
meta.pattern_file     = pattern_file;
meta.bleach_ratio     = bleach_ratio;
meta.cgh_method       = Setup.CGHMethod;
meta.timestamp        = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<TNOW1,DATST>
% If slm_pupil_mapping.m has been run in this session, fold in only its SCALAR
% summary parameters (the full struct holds raw camera frames -- those stay in
% the pupil-mapping .mat, not in this stack's JSON).
if exist('pupil_mapping', 'var')
    pm = struct();
    for f = {'rotation_deg','mirror_flip','pupil_center_px', ...
             'pupil_diameter_px','defocus_scale_um_per_unit'}
        if isfield(pupil_mapping, f{1}), pm.(f{1}) = pupil_mapping.(f{1}); end
    end
    meta.pupil_mapping = pm;
end

outpaths = save_bead_stack(outdir, stem, dataUZ, bgd, refs, meta);

slm.blank()   % park the SLM, as get_psf_no_power does on exit

fprintf('\nAcquisition complete in %.0f s. Output in:\n  %s\n', toc(tBegin), outdir);
disp(outpaths)
