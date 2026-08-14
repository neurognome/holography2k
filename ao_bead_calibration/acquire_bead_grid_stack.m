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
% Modeled on get_psf_v2.m. Differences: N-spot grid instead of one spot, SLM fed
% once and never touched again, and the RAW averaged stack is saved to TIFF (no
% smoothing / no background subtraction / no clipping) for Genesis's pipeline.
%
% RUN ORDER (see README):
%   1. On the DAQ computer: start the power responder (alignment/alignCodeDAQ2K.m
%      or equivalent) so the 42130 msocket gate below gets 'gotit' acks.
%   2. Run this script cell-by-cell on the holography computer.
%
% Hand-run script (clears the workspace), matching the get_psf_v2.m convention.

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

% Power: set LOW. This is the per-hologram power in mW; it is gated on only
% during each grab. Tune it in the "no-saturation check" cell below BEFORE the
% stack so that no raw pixel reaches the 8-bit camera ceiling (255).
pwr             = 0.1;               % mW  (matches get_psf_v2 default; verify per rig)

% z-sweep: holographic spots have long axial tails -> wide range, fine step.
UZ              = linspace(-30, 30, 61);   % um about focus, ~1 um step (README item 13)
nframesCapture  = 10;                       % frames averaged per plane

% Objective / analysis metadata that is NOT knowable from the rig code -- fill
% these in for the handoff (README item 3 / 15). Left here so they are obvious.
na              = [];                % objective NA           e.g. 0.8
objective_mag   = [];                % objective magnification e.g. 16
pupil_fill      = [];                % pupil-fill fraction     e.g. 0.9

%% ---- hardware setup (mirrors get_psf_v2.m:1-40) -----------------------------
Setup = function_loadparameters2();
Setup.CGHMethod = 2;                 % Global GS, as in get_psf_v2
Setup.GSoffset  = 0;
Setup.verbose   = 0;
Setup.useGPU    = 1;
Setup.SLM.is_onek = 1;

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

%% ---- connect to DAQ computer (power gating, port 42130) ---------------------
% Identical handshake to get_psf_v2.m:45-63. Start the DAQ-side responder first.
fprintf('Waiting for msocket communication From DAQ... ')
srvsock = mslisten(42130);
masterSocket = msaccept(srvsock, 15);
msclose(srvsock);
mssend(masterSocket, 'A');

invar = [];
while ~strcmp(invar, 'B')
    invar = msrecv(masterSocket, .5);
end
mssend(masterSocket, wavelength);
fprintf('done.\r')

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

%% ---- no-saturation check (README item 11) -----------------------------------
% Project the grid at the working power and confirm NO raw pixel clips at 255.
% Adjust `pwr` above and re-run this cell until the brightest spot is well below
% the ceiling (aim ~80% of camMax on the brightest bead).
bgd_frames = bas.grab(10);
bgd = mean(bgd_frames, 3);

mssend(masterSocket, [pwr/1000 1 1]);   % beam ON
wait_gotit(masterSocket);
raw = bas.grab(nframesCapture);
mssend(masterSocket, [0 1 1]);          % beam OFF
wait_gotit(masterSocket);

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

    mssend(masterSocket, [pwr/1000 1 1]);   % beam ON
    wait_gotit(masterSocket);

    data = bas.grab(nframesCapture);

    mssend(masterSocket, [0 1 1]);          % beam OFF (limits bleaching)
    wait_gotit(masterSocket);

    dataUZ(:,:,i) = mean(double(data), 3);  % RAW averaged plane

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
% axial profile is corrupted -> retake at lower power.
sutter.moveZ(UZ(1)); pause(0.5);
mssend(masterSocket, [pwr/1000 1 1]); wait_gotit(masterSocket);
recheck = mean(double(bas.grab(nframesCapture)), 3);
mssend(masterSocket, [0 1 1]); wait_gotit(masterSocket);
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

%% ---- pixel size: pxPerMu via a 50 um Sutter move (get_psf_v2.m:210-263) -----
muUsed = 50;
disp('Determining pxPerMu...')

% peak plane / spot to track (brightest bead in the grid)
mxProj = max(dataUZ, [], 3);
[cx, cy] = function_findcenter(mxProj);

sutter.moveTo([0 0 0]); pause(1);
mssend(masterSocket, [pwr/1000 1 1]); wait_gotit(masterSocket);
p1 = mean(double(bas.grab(nframesCapture)), 3);
mssend(masterSocket, [0 1 1]); wait_gotit(masterSocket);

sutter.moveTo([0 muUsed 0]); pause(1);
mssend(masterSocket, [pwr/1000 1 1]); wait_gotit(masterSocket);
p2 = mean(double(bas.grab(nframesCapture)), 3);
mssend(masterSocket, [0 1 1]); wait_gotit(masterSocket);
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
mssend(masterSocket, [pwr/1000 1 1]); wait_gotit(masterSocket);
refs.reference_blank = mean(double(bas.grab(nframesCapture)), 3);
mssend(masterSocket, [0 1 1]); wait_gotit(masterSocket);
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
meta.power_mW         = pwr;
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

fprintf('\nAcquisition complete in %.0f s. Output in:\n  %s\n', toc(tBegin), outdir);
disp(outpaths)

%% ---- local helper -----------------------------------------------------------
function wait_gotit(sock)
% Block until the DAQ power responder acks the gate command (get_psf_v2 pattern).
invar = [];
while ~strcmp(invar, 'gotit')
    invar = msrecv(sock, 0.01);
end
end
