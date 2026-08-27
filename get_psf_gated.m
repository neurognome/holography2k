%% get_psf_gated.m
% PSF of a SINGLE hologram, with laser power controlled over Holochat by the DAQ
% computer -- the same gating the AO bead scripts use.
%
% This is get_psf_no_power.m plus power control, and minus three things that bias
% a width measurement (see MEASUREMENT NOTES below). The hologram is gated ON only
% for each grab and OFF in between, which gives a true beam-off background and
% limits hole-burning during a long axial sweep.
%
% RUN ORDER
%   1. DAQ computer: run alignCodeDAQ2K FIRST. Wait for it to print
%      'Waiting for Hologram info'.
%   2. Holography computer: run this file cell by cell.
%
%   The order is not a preference. HolochatInterface's constructor calls
%   io.reset(id), which is a DELETE on the whole 'daq' store -- so starting
%   alignCodeDAQ2K DESTROYS any wavelength already sitting in that mailbox. Send it
%   first and the DAQ wipes it, then times out its comm.read(10) and errors with
%   alignCodeDAQ2K:wavelength.
%
%   If that happens, you do NOT have to re-run the setup cell (which re-opens the
%   SLM, camera and stage). Restart alignCodeDAQ2K, then in the command window:
%       comm.send(wavelength, 'daq')
%
% MEASURING MORE THAN ONCE -- do not run the release cell
%   The last cell sends 'end', which exits the DAQ's power loop for good. Leave it
%   alone until you are finished for the session and everything below is repeatable
%   with nothing restarted: alignCodeDAQ2K's own timeout is 100000 s and resets on
%   every power request, so it will sit there all day.
%
%   To repeat a measurement, or to A/B a correction:
%       slm.feed(Holo)                     % the release cell blanks the SLM; this
%                                          % is only needed if you ran it
%       re-run 'build & feed' if you changed correction_file / correction_sign
%       re-run 'acquire the axial stack' onward
%
%   Only run the release cell at the very end. If you did run it, restart
%   alignCodeDAQ2K and re-send the wavelength as above.
%
% MEASUREMENT NOTES -- why this does not just copy get_psf_no_power
%   * Fits run on the RAW (background-subtracted) stack; smoothing is display-only.
%     get_psf_no_power stores imgaussfilt(frame, 2) and then fits that, and a
%     Gaussian smooth adds to the width you are trying to measure, in quadrature:
%
%         FWHM_measured = sqrt( FWHM_true^2 + (2.355*sigma_px/pxPerMu)^2 )
%
%     How big that is depends entirely on pxPerMu, so measure it before deciding
%     the bias is tolerable -- for a true 2 um FWHM read through a sigma = 2 px
%     blur, +16% at 4 px/um and +86% at 1.5 px/um. Smoothing is still used for
%     centroid finding, where it only helps.
%   * Background is subtracted but NOT clipped at zero for the fit. max(x-bgd,0)
%     rectifies the noise floor into a positive pedestal, which pulls a gauss1 fit
%     (it has no offset term) wider and shallower.
%   * The peak plane is checked for saturation AFTER the sweep. A clipped peak is
%     flat-topped, so every FWHM from it is an overestimate -- and it looks like a
%     clean measurement.
%
% Set correction_file to compare a wavefront correction against none: run once
% with it empty, once with it set, and read the peak-intensity change. Two-photon
% signal is quadratic in intensity, so peak brightness moves further than FWHM
% and is the more sensitive readout (AO protocol Phase 4 item 20). Flip
% correction_sign to settle the sign the same way (item 19).
%
% See also: get_psf_no_power, ao_bead_calibration/acquire_bead_grid_stack,
%           load_slm_correction, docs/slm-wavefront-correction.md

clear
close all
clc

tBegin = tic;
disp('Setting up single-hologram PSF measurement...');

makePaths()

%% ---- user parameters --------------------------------------------------------
wavelength      = 900;               % 900 / 1030 / 1100 / 607
label           = 'psf';             % label for the saved stem

% The single target, in normalized SLM coords: [x y z power].
% Put it near where you actually stimulate, not at the zero order.
slmCoords       = [0.4, 0.4, 0, 1];

% Commanded power, mW. Start LOW and raise it in the no-saturation cell until the
% spot peaks around 80% of the camera ceiling.
pwr             = 5;                 % mW

% Axial sweep. A single holographic spot has long axial tails; keep the range wide
% enough to see the baseline on both sides or the axial fit has nothing to sit on.
UZ              = linspace(-50, 50, 51);   % um about focus, 2 um step
nframesCapture  = 10;                       % frames averaged per plane

% Sutter settle times, matched to acquire_bead_grid_stack.m. Grabbing before the
% stage settles gives blurred, z-misregistered planes. Do not shorten these
% without checking against a known bead.
settle_first_s  = 3;                 % s, after the large initial jump
settle_step_s   = 0.3;               % s, after each small z step

% Static wavefront correction. '' = none (and then this behaves exactly like the
% uncorrected case). A path is loaded by load_slm_correction and applied by
% phase_to_frame to the hologram below.
correction_file = '';                % e.g. fullfile('slm','corrections','slm_pattern_14umaFWHM.bmp')
correction_sign = 1;                 % +1 / -1, test both (Phase 4 item 19)

muUsed          = 50;                % um, lateral Sutter move used for pxPerMu

save_result     = true;              % write a .mat next to the other rig data

%% ---- hardware setup ---------------------------------------------------------
Setup = function_loadparameters2();
Setup.CGHMethod = 2;                 % Global GS
Setup.GSoffset  = 0;
Setup.verbose   = 0;

if Setup.useGPU
    disp('Getting gpu...');
    g = gpuDevice; %#ok<NASGU>
end

slm = get_slm(wavelength);
slm.stop();
slm.wait_for_trigger = 0;
slm.start();

sutter = sutterController();

bas = bascam();
bas.start()

% Power control over Holochat. alignCodeDAQ2K must already be running on the DAQ.
comm = HolochatInterface('holo');
comm.send(wavelength, 'daq');        % which laser to drive; it reads this ONCE
comm.flush();

disp('Setup complete.')

%% ---- 1p preview: find the objective ----------------------------------------
bas.preview()

%% ---- build & feed the single-spot hologram ---------------------------------
assert(size(slmCoords, 1) == 1, 'get_psf_gated:notSingle', ...
    ['slmCoords must be ONE row for a single-hologram PSF, got %d. A multi-spot ' ...
     'hologram\nsplits the power and the fits below assume one spot -- use ' ...
     'ao_bead_calibration/acquire_bead_grid_stack.m for a grid.'], ...
    size(slmCoords, 1));

Setup.SLM.correction = [];
if ~isempty(correction_file)
    Setup.SLM.correction = load_slm_correction(correction_file, ...
        'Size', [Setup.SLM.Nx Setup.SLM.Ny], 'Sign', correction_sign);
    fprintf('Applying wavefront correction (sign %+d).\n', correction_sign);
else
    disp('No wavefront correction (uncorrected reference).');
end

[Holo, Reconstruction, Masksg] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
slm.feed(Holo);
fprintf('Fed single-spot hologram at SLM [%.3f %.3f %.1f].\n', slmCoords(1:3));

%% ---- set power + no-saturation check ---------------------------------------
% Preview gated ON so you can watch the spot while tuning pwr, then confirm the
% raw frames do not clip. Re-run this cell after changing pwr.
power_gate(comm, pwr/1000);
bas.preview()                        % close the window when the power looks right
power_gate(comm, 0);

raw    = grab_gated(comm, bas, pwr/1000, nframesCapture, 'raw');
rawmax = max(raw(:));
fprintf('Raw camera max = %d / %d (%.3f%% of pixels at ceiling).\n', ...
    rawmax, bas.camMax, 100*mean(raw(:) >= bas.camMax));
if rawmax >= bas.camMax
    warning('get_psf_gated:saturated', ...
        ['Camera SATURATES at %d DN. Lower pwr and re-run this cell -- a ' ...
         'flat-topped peak\nmakes every FWHM below an overestimate.'], rawmax);
end

figure('Name','Working-power frame'); clf
imagesc(mean(double(raw), 3)); axis image; colorbar
title(sprintf('max = %d DN (ceiling %d)', rawmax, bas.camMax));

%% ---- beam-off background ---------------------------------------------------
% A true dark frame: power_gate is never called, so the DAQ holds the beam off.
bgd = mean(double(bas.grab(nframesCapture)), 3);
fprintf('Background: mean %.2f DN, std %.2f DN\n', mean(bgd(:)), std(bgd(:)));

%% ---- acquire the axial stack (gated per plane) ------------------------------
% The SLM phase is fixed for the whole sweep: defocus comes from the STAGE. Moving
% it with the SLM would change the thing being measured.
sutter.setRef()
sz     = size(bgd);
dataUZ = zeros([sz numel(UZ)]);

figure('Name','Acquiring PSF stack'); clf
disp('Collecting PSF stack.')

for i = 1:numel(UZ)
    fprintf('Plane %d/%d (z = %+.1f um)\n', i, numel(UZ), UZ(i));
    sutter.moveZ(UZ(i))

    if i == 1, pause(settle_first_s); else, pause(settle_step_s); end

    % Raw frame-average minus background, NOT clipped and NOT smoothed.
    dataUZ(:,:,i) = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double') - bgd;

    subplot(1,2,1)
    imagesc(imgaussfilt(dataUZ(:,:,i), 2)); axis image; colorbar
    title(sprintf('Plane %d (z = %+.1f um)', i, UZ(i)))
    subplot(1,2,2)
    imagesc(max(dataUZ, [], 3)); axis image; colorbar; title('Max projection')
    drawnow
end

sutter.moveToRef()
pause(0.1)
disp('Done collecting stack.')

%% ---- locate the spot and its peak plane ------------------------------------
% Centroid on a SMOOTHED max-projection (robust), everything downstream on raw.
mxProj  = max(dataUZ, [], 3);
[cx, cy] = function_findcenter(imgaussfilt(mxProj, 2));

roi     = 2;                          % +-px about the centre for the axial trace
dimx    = max(cx-roi, 1):min(cx+roi, sz(1));
dimy    = max(cy-roi, 1):min(cy+roi, sz(2));

axialTrace = squeeze(mean(mean(dataUZ(dimx, dimy, :), 1), 2));
[~, peakPlane] = max(axialTrace);
peakFrame = dataUZ(:, :, peakPlane);

fprintf('Spot at (%d, %d); peak plane %d of %d (z = %+.1f um).\n', ...
    cx, cy, peakPlane, numel(UZ), UZ(peakPlane));

% A clipped peak invalidates every width below, so say so now.
peak_raw = peakFrame(cx, cy) + bgd(cx, cy);
if peak_raw >= bas.camMax
    warning('get_psf_gated:peakSaturated', ...
        ['The PEAK PLANE is saturated (%.0f DN at the ceiling of %d). Its ' ...
         'flat top makes\nboth FWHMs overestimates. Retake at lower power.'], ...
        peak_raw, bas.camMax);
end

%% ---- pixel size: pxPerMu from a lateral Sutter move ------------------------
% Measured AT the peak plane, so the spot is in focus for both grabs.
disp('Determining pxPerMu...')

sutter.moveTo([0 0 UZ(peakPlane)]);      pause(settle_first_s);
pos1 = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double') - bgd;

sutter.moveTo([0 muUsed UZ(peakPlane)]); pause(settle_first_s);
pos2 = grab_gated(comm, bas, pwr/1000, nframesCapture, 'double') - bgd;

sutter.moveToRef()
pause(0.1)

[x1, y1] = function_findcenter(imgaussfilt(pos1, 2));
[x2, y2] = function_findcenter(imgaussfilt(pos2, 2));
pxPerMu  = pdist([x1 y1; x2 y2]) / muUsed;
pixel_size_um = 1 / pxPerMu;
fprintf('pxPerMu = %.3f  ->  pixel size = %.3f um\n', pxPerMu, pixel_size_um);

%% ---- fit lateral and axial profiles ----------------------------------------
% MATLAB's gauss1 is a1*exp(-((x-b1)/c1)^2), so FWHM = 2*sqrt(log(2))*c1.
% (get_psf_no_power writes this as 2*sqrt(2*log(2))*c1/sqrt(2) -- same number,
% 1.6651*c1, just harder to check.)
fwhm_of = @(c1) 2*sqrt(log(2)) * c1;

% AXIS CONVENTION. function_findcenter returns (x, y) with x indexing DIMENSION 1
% and y indexing dimension 2 -- the same order function_Make_3D_SHOT_Holos uses for
% Masks(x,y). So the x profile runs DOWN a column and the y profile runs ALONG a row.
% get_psf_no_power labels these the other way round (its "xLine" is peakFrame(x,:),
% which varies along dimension 2); the widths it prints are right, the axis names
% are swapped. Named by dimension here so the two cannot be confused.
xLine = double(peakFrame(:, cy))';            % along dim 1 = x
yLine = double(peakFrame(cx, :));             % along dim 2 = y
xUm   = ((1:numel(xLine)) - cx) / pxPerMu;    % um, centred on the spot
yUm   = ((1:numel(yLine)) - cy) / pxPerMu;

fitX = fit(xUm(:),  xLine(:),      'gauss1');
fitY = fit(yUm(:),  yLine(:),      'gauss1');
fitZ = fit(UZ(:),   axialTrace(:), 'gauss1');

xFWHM = fwhm_of(fitX.c1);
yFWHM = fwhm_of(fitY.c1);
zFWHM = fwhm_of(fitZ.c1);
peakIntensity = max(peakFrame(:));

fprintf('\n--- PSF ---\n');
fprintf('lateral FWHM  x = %.2f um   y = %.2f um\n', xFWHM, yFWHM);
fprintf('axial   FWHM  z = %.2f um   (centre %+.2f um)\n', zFWHM, fitZ.b1);
fprintf('peak intensity  = %.1f DN above background\n', peakIntensity);
if ~isempty(correction_file)
    fprintf('correction      = %s (sign %+d)\n', correction_file, correction_sign);
else
    fprintf('correction      = none\n');
end

%% ---- summary figure --------------------------------------------------------
figure('Name','PSF summary'); clf

subplot(2,2,1)
plot(UZ, axialTrace, 'o'); hold on; plot(fitZ);
title(sprintf('Axial FWHM %.2f \\mum', zFWHM))
xlabel('Z position (\mum)'); ylabel('Intensity (DN)')
legend('Measured','Fit','Location','best'); legend boxoff

subplot(2,2,2)
plot(xUm, xLine, 'o'); hold on; plot(fitX);
xlim([fitX.b1-15 fitX.b1+15])
title(sprintf('Lateral FWHM %.2f \\mum (x)', xFWHM))
xlabel('X position (\mum)'); ylabel('Intensity (DN)')
legend('Measured','Fit','Location','best'); legend boxoff

subplot(2,2,3)
% imagesc(xdata, ydata, C): the horizontal axis runs along dim 2 (= y here).
imagesc(yUm, xUm, peakFrame); axis image
xlim([-15 15]); ylim([-15 15]); colorbar
title(sprintf('Peak frame (z = %+.1f \\mum)', UZ(peakPlane)))
xlabel('y (\mum)'); ylabel('x (\mum)')

subplot(2,2,4)
% Axial section through the spot centre: how the width evolves with z. The slice is
% dim1-by-Z, so the vertical axis is x.
imagesc(UZ, xUm, squeeze(dataUZ(:, cy, :))); colorbar
ylim([-15 15])
title('Axial section (through spot centre)')
xlabel('Z (\mum)'); ylabel('x (\mum)')

%% ---- save ------------------------------------------------------------------
if save_result
    stamp = datestr(now, 'yymmdd');                       %#ok<TNOW1,DATST>
    stem  = sprintf('%s_%s_%dnm', stamp, label, wavelength);
    try
        data_root = rig_remote_get('paths.data_root', pwd);
    catch
        data_root = pwd;
    end
    outdir = fullfile(data_root, 'psf', stamp);
    if ~isfolder(outdir), mkdir(outdir); end

    psf = struct();
    psf.slmCoords       = slmCoords;
    psf.wavelength_nm   = wavelength;
    psf.power_mW        = pwr;
    psf.z_planes_um     = UZ;
    psf.axial_trace     = axialTrace;
    psf.pxPerMu         = pxPerMu;
    psf.pixel_size_um   = pixel_size_um;
    psf.peak_plane      = peakPlane;
    psf.spot_px         = [cx cy];
    psf.fwhm_x_um       = xFWHM;
    psf.fwhm_y_um       = yFWHM;
    psf.fwhm_z_um       = zFWHM;
    psf.peak_intensity  = peakIntensity;
    psf.frames_averaged = nframesCapture;
    psf.camera_exposure = bas.exposure;
    psf.cgh_method      = Setup.CGHMethod;
    psf.correction_file = correction_file;
    psf.correction_sign = correction_sign;
    psf.timestamp       = datestr(now, 'yyyy-mm-ddTHH:MM:SS'); %#ok<TNOW1,DATST>

    outfile = fullfile(outdir, [stem '.mat']);
    save(outfile, 'psf', 'dataUZ', 'bgd', 'peakFrame', '-v7.3');
    fprintf('Saved %s\n', outfile);
end

%% ---- release ---------------------------------------------------------------
comm.send('end', 'daq');   % DAQ zeros the power and acks 'kthx'
slm.blank()                % park the SLM

fprintf('\nPSF measurement complete in %.0f s.\n', toc(tBegin));

%% ---- local helpers --------------------------------------------------------
function img = grab_gated(comm, bas, pwr_W, nframes, castTo)
% Gate the hologram ON, grab, gate OFF. 'double' averages the frames; 'raw'
% returns the uint8 stack for max()/saturation checks.
if nargin < 5, castTo = 'double'; end
power_gate(comm, pwr_W);
frames = bas.grab(nframes);
power_gate(comm, 0);
if strcmp(castTo, 'raw')
    img = frames;
else
    img = mean(double(frames), 3);
end
end

function power_gate(comm, val_W)
% Ask the DAQ (alignCodeDAQ2K) for val_W watts (0 = off) and wait for its ack.
% The message is [power_W divisor multiplier] and the DAQ computes
% PowerRequest = power_W*mult/div, so [val 1 1] requests exactly val watts.
% It must be 3 elements: a shorter numeric message is ignored with a warning
% on the DAQ side, and this would then block forever waiting for 'gotit'.
comm.send([val_W, 1, 1], 'daq');
invar = [];
while ~strcmp(invar, 'gotit')
    invar = comm.read(0.01);
end
end
