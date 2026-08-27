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
% Laser power is controlled over Holochat by the DAQ computer, exactly like
% align_slm_to_camera_scope2k.m: the hologram is gated ON only during each grab
% and OFF in between. That gives a true beam-off background and limits bleaching.
%
% RUN ORDER (see README):
%   1. On the DAQ computer: run alignCodeDAQ2K (it waits for the wavelength, then
%      services power requests and acks each with 'gotit').
%   2. Run this script cell-by-cell on the holography computer.
%
% Hand-run script (clears the workspace), matching the get_psf/align conventions.

clear
close all
clc

tBegin = tic;
disp('Setting up AO bead-grid acquisition...');

makePaths()

%% ---- user parameters --------------------------------------------------------
wavelength      = 900;              % SLM wavelength (900 / 1030 / 1100)
mouse_or_slide  = 'beadslide';       % label for the output stem
epoch           = 1;

% Grid coordinate source (see make_grid_coords.m). [] -> synthetic 20-pt grid
% until Genesis's coordinate file is in hand. When it arrives:
%   gridSource = '/path/to/genesis_targets.csv';  gridOpts.frame = 'SLM'; % or 'SI'
gridSource      = [];
gridOpts        = struct('frame','SLM', 'n',36);

% Per-hologram power in mW, commanded to the DAQ. Set LOW and tune it in the
% no-saturation cell below so no raw pixel reaches the 8-bit ceiling (255).

% z-sweep: holographic spots have long axial tails -> wide range, fine step.
UZ              = linspace(-50, 50, 101);   % um about focus, ~1 um step (README item 13)
nframesCapture  = 10;                       % frames averaged per plane

% Sutter settle times. Matched to align_slm_to_camera_scope2k.m's camera z-loop:
% the first plane is a large jump from reference and needs a long settle, every
% subsequent 1 um step a short one. Grabbing before the stage settles gives
% blurred / z-misregistered planes -- do not shorten these without checking.
settle_first_s  = 3;                         % s, after the large initial jump
settle_step_s   = 0.3;                        % s, after each small z step

% Objective / analysis metadata that is NOT knowable from the rig code -- fill
% these in for the handoff (README item 3 / 15). Left here so they are obvious.
na              = 0.8;               % objective NA            e.g. 0.8
objective_mag   = 16;               % objective magnification e.g. 16
pupil_fill      = 0.9;               % pupil-fill fraction     e.g. 0.9

%% ---- hardware setup ---------------------------------------------------------
Setup = function_loadparameters2();
Setup.CGHMethod = 2;                 % Global GS
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

% Power control over Holochat (DAQ computer must be running alignCodeDAQ2K).
comm = HolochatInterface('holo');
comm.send(wavelength, 'daq');        % tells the DAQ which laser to calibrate
comm.flush();

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

%% ---- set power + no-saturation check (README item 11) -----------------------
% Preview with the hologram gated ON so you can watch the beads while tuning; then
% confirm NO raw pixel clips at 255. Adjust `pwr` above and re-run until the
% brightest bead is ~80% of the ceiling.

pwr             = 75;                 % mW
power_gate(comm, pwr/1000);   % beam ON
bas.preview()                 % close the preview window when the power looks right
power_gate(comm, 0);          % beam OFF

raw = grab_gated(comm, bas, pwr/1000, nframesCapture, 'uint8');
rawmax = max(raw(:));
frac_sat = mean(raw(:) >= bas.camMax);
fprintf('Raw camera max = %d / %d (%.3f%% of pixels saturated).\n', ...
    rawmax, bas.camMax, 100*frac_sat);
if rawmax >= bas.camMax
    warning('acquire_bead_grid:saturated', ...
        ['Camera is SATURATING at the current power (%d DN). Lower `pwr` and ' ...
         're-run this cell before acquiring the stack -- clipped data breaks ' ...
         'phase retrieval and deconvolution.'], rawmax);
end

figure('Name','Working-power frame'); clf
imagesc(mean(double(raw),3)); axis image; colorbar
title(sprintf('max = %d DN  (ceiling %d)', rawmax, bas.camMax));

%% ---- beam-off background ----------------------------------------------------
% True dark frame: the DAQ keeps the beam off (power_gate is not called), so this
% is real camera/room background. Saved alongside the stack; NOT subtracted from
% the raw stack.
bgd = mean(double(bas.grab(10)), 3);

%% ---- acquire the z-stack (RAW; SLM fixed; gated per plane) ------------------
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

    if i == 1, pause(settle_first_s); else, pause(settle_step_s); end

    % hologram ON only for the grab, OFF in between (limits bleaching)
    dataUZ(:,:,i) = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double');

    % live view only (does not touch the saved array)
    subplot(1,2,1)
    imagesc(imgaussfilt(dataUZ(:,:,i), 2)); axis image; colorbar
    title(sprintf('Plane %d (z=%+.1f)', i, UZ(i)))
    subplot(1,2,2)
    imagesc(max(dataUZ, [], 3)); axis image; colorbar; title('Max projection')
    drawnow
end

sutter.moveToRef()
pause(0.1)
disp('Done collecting stack.')

%% ---- bleaching check (README item 11) ---------------------------------------
% Re-grab the first plane; if the field is now dimmer at equivalent focus the
% axial profile is corrupted -> retake at lower power.
sutter.moveZ(UZ(1)); pause(settle_first_s);   % large jump back to the first plane
recheck = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double');
sutter.moveToRef()
pause(0.1)

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

sutter.moveTo([0 0 0]); pause(settle_first_s);
p1 = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double');

sutter.moveTo([0 muUsed 0]); pause(settle_first_s);   % 50 um lateral jump
p2 = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double');
sutter.moveToRef()
pause(0.1)

[x1, y1] = function_findcenter(max(p1 - bgd, 0));
[x2, y2] = function_findcenter(max(p2 - bgd, 0));
pxPerMu = pdist([x1 y1; x2 y2]) / muUsed;
pixel_size_um = 1 / pxPerMu;
fprintf('pxPerMu = %.3f  ->  pixel size = %.3f um\n', pxPerMu, pixel_size_um);

%% ---- reference images (README item 14) --------------------------------------
refs = struct();

% (a) blank-SLM raster/frame of the same field, no hologram.
slm.blank();   % feeds zeros(Nx,Ny); dimension-safe
refs.reference_blank = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double');
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
meta.power_mW         = pwr;         % commanded to the DAQ over Holochat
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
% Fold in only the SCALAR mapping summary from slm_pupil_mapping.m (the full
% struct holds raw camera frames -- those stay in the pupil-mapping .mat, not in
% this stack's JSON). This script `clear`s at the top, so read it from the .mat
% slm_pupil_mapping.m saved today, not the workspace.
pm_src = [];
if exist('pupil_mapping', 'var')            % (if run without clearing)
    pm_src = pupil_mapping;
else
    pm_file_in = fullfile(outdir, sprintf('%s_pupil_mapping_%dnm.mat', stamp, wavelength));
    if isfile(pm_file_in), pm_src = load(pm_file_in); end
end
if ~isempty(pm_src)
    pm = struct();
    for f = {'rotation_deg','mirror_flip','pupil_center_px', ...
             'pupil_diameter_px','defocus_scale_um_per_unit'}
        if isfield(pm_src, f{1}), pm.(f{1}) = pm_src.(f{1}); end
    end
    meta.pupil_mapping = pm;
end

outpaths = save_bead_stack(outdir, stem, dataUZ, bgd, refs, meta);

comm.send('end', 'daq');   % release the DAQ power loop (it zeros power, acks 'kthx')
slm.blank()                % park the SLM

fprintf('\nAcquisition complete in %.0f s. Output in:\n  %s\n', toc(tBegin), outdir);
disp(outpaths)

%% ---- local helpers ----------------------------------------------------------
function img = grab_gated(comm, bas, pwr_W, nframes, castTo)
% Gate the hologram ON, grab+average, gate OFF. Returns 'double' (default) for
% the saved stack, or 'uint8' raw frames for the saturation check.
if nargin < 5, castTo = 'double'; end
power_gate(comm, pwr_W);
frames = bas.grab(nframes);
power_gate(comm, 0);
if strcmp(castTo, 'uint8')
    img = frames;                       % raw, for max()/saturation test
else
    img = mean(double(frames), 3);      % averaged plane
end
end

function power_gate(comm, val_W)
% Command the DAQ (alignCodeDAQ2K) to set laser power to val_W (0 = off) and wait
% for its 'gotit' ack. Message is [power_W, divisor, multiplier]; the DAQ computes
% PowerRequest = power_W*mult/div, so [val 1 1] requests exactly val watts.
comm.send([val_W, 1, 1], 'daq');
invar = [];
while ~strcmp(invar, 'gotit')
    invar = comm.read(0.01);
end
end
