%% slm_pupil_mapping.m
% Phase 1 of the holographic AO bead-calibration protocol: measure the four
% SLM -> back-pupil mapping parameters Genesis needs (README items 5-9):
%
%     rotation angle, mirror flip (yes/no), pupil diameter in SLM pixels,
%     pupil center in SLM pixels
%
% plus a confirmation of the existing SLM -> image affine.
%
% These are measured empirically on the Basler camera by feeding single-target
% (or half-blocked) holograms and reading where / how bright the spot lands.
% Reuses the same primitives as get_psf_v2 / acquire_bead_grid_stack:
%   function_Make_3D_SHOT_Holos, slm.feed, bas.grab, function_findcenter.
%
% Run cell-by-cell on the holography computer. Laser power is controlled over
% Holochat by the DAQ computer (run alignCodeDAQ2K there first), so the hologram
% is gated ON only during each grab. The final cell assembles a `pupil_mapping`
% struct that acquire_bead_grid_stack.m folds into the handoff metadata
% automatically if it is present in the workspace.
%
% Several parameters here are semi-quantitative read-offs (README explicitly
% frames them as empirical). Where a value cannot be reduced to a single number
% purely in software, the raw sweep is stored in `pupil_mapping` so it can be
% finalized against Genesis's coordinate convention.

clear
close all
clc
makePaths()

%% ---- setup (same as acquire_bead_grid_stack.m) ------------------------------
wavelength = 900;
nframes    = 10;

Setup = function_loadparameters2();
Setup.CGHMethod = 2; Setup.GSoffset = 0; Setup.verbose = 0;
Setup.useGPU = 1; Setup.SLM.is_onek = 1;
if Setup.useGPU, g = gpuDevice; end

slm = get_slm(wavelength);
slm.stop(); slm.wait_for_trigger = 0; slm.start();
sutter = sutterController();
bas = bascam(); bas.start()


% Power control over Holochat (DAQ computer must be running alignCodeDAQ2K).
comm = HolochatInterface('holo');
comm.send(wavelength, 'daq');
comm.flush();
%%
% Preview a central spot (gated ON) and tune `pwr` to ~80% of camera max before
% running the mapping cells.
pwr        = 0.5;      % mW, commanded to the DAQ; gated ON per grab
slm.feed(function_Make_3D_SHOT_Holos(Setup, [0.55 0.55 0 1]));
power_gate(comm, pwr/1000); bas.preview(); power_gate(comm, 0);

sutter.setRef()

% central reference target used throughout
xc = 0.5; yc = 0.5;

pupil_mapping = struct('wavelength_nm', wavelength);

%% ---- (5) ROTATION + HANDEDNESS from a pure x-tilt ---------------------------
% Displace the target along SLM-x (and separately SLM-y) and measure the spot's
% motion on the camera. The 2x2 map [dCam/dSLM] gives the rotation angle; a
% negative determinant means the SLM<->camera mapping includes a mirror flip.
dxy = 0.15;                          % SLM-coord step for the tilt probe
probes = [ xc-dxy, yc     ;          % -x
           xc+dxy, yc     ;          % +x
           xc,     yc-dxy ;          % -y
           xc,     yc+dxy ];         % +y
cen = zeros(4, 2);
for k = 1:4
    cen(k,:) = grab_spot_center(Setup, slm, bas, comm, pwr, [probes(k,:) 0 1], nframes);
end
% camera displacement per unit SLM move, columns = [d/dx , d/dy]
dCam_dx = (cen(2,:) - cen(1,:))' / (2*dxy);   % 2x1 [drow; dcol] per unit SLM-x
dCam_dy = (cen(4,:) - cen(3,:))' / (2*dxy);   % 2x1 per unit SLM-y
J = [dCam_dx, dCam_dy];                        % 2x2 Jacobian cam<-SLM
rotation_deg = atan2d(J(2,1), J(1,1));
mirror_flip  = det(J) < 0;

pupil_mapping.rotation_deg = rotation_deg;
pupil_mapping.mirror_flip  = mirror_flip;
pupil_mapping.tilt_jacobian_cam_from_slm = J;
fprintf('(5) rotation ~ %.1f deg;  mirror flip = %d;  det(J) = %.3g\n', ...
    rotation_deg, mirror_flip, det(J));

%% ---- (6) MIRROR-FLIP CONFIRM via half-pupil block ---------------------------
% Blank the left half of the SLM and image the resulting PSF asymmetry; repeat
% for the right half. Which half of the PSxF dims (and on which camera side)
% independently confirms orientation / handedness from (5).
spot = [xc yc 0 1];
[Holo, ~, ~] = function_Make_3D_SHOT_Holos(Setup, spot);

Holo_L = Holo; Holo_L(:, 1:512) = 0;     % block left half of SLM
Holo_R = Holo; Holo_R(:, 513:end) = 0;   % block right half

img_full = feed_and_grab(slm, bas, comm, pwr, Holo,   nframes);
img_L    = feed_and_grab(slm, bas, comm, pwr, Holo_L, nframes);
img_R    = feed_and_grab(slm, bas, comm, pwr, Holo_R, nframes);
slm.feed(Holo);

figure('Name','(6) Half-pupil block'); clf
subplot(1,3,1); imagesc(img_full); axis image; title('full pupil'); colorbar
subplot(1,3,2); imagesc(img_L);    axis image; title('left half blocked'); colorbar
subplot(1,3,3); imagesc(img_R);    axis image; title('right half blocked'); colorbar

pupil_mapping.halfblock_full  = img_full;
pupil_mapping.halfblock_left  = img_L;
pupil_mapping.halfblock_right = img_R;
disp('(6) Half-pupil frames stored. Read the PSF-shift direction to confirm flip.')

%% ---- (7) DEFOCUS SCALE: um focal shift per unit SLM-z -----------------------
% Feed a central spot at a range of SLM defocus values; for each, find the
% Sutter-z that maximizes on-axis intensity (the physical focus of that
% hologram). The slope d(focal_um)/d(SLM_z) is the axial scale of the SLM's
% defocus term -- the quantity that, with the Fresnel model in
% function_GenerateFresnelPropagationStack, sets the effective pupil diameter.
zSLM   = [-0.03 -0.015 0 0.015 0.03];   % SLM defocus coeffs to probe
zScan  = linspace(-100, 100, 41);         % Sutter search range (um) per defocus
focus_um = zeros(size(zSLM));

% Sutter settle times, matched to align_slm_to_camera_scope2k.m: long settle
% after the large jump back to the start of each scan, short after each step.
settle_first_s = 3;
settle_step_s  = 0.3;

sutter.setRef()
for m = 1:numel(zSLM)
    [Hm, ~, ~] = function_Make_3D_SHOT_Holos(Setup, [xc-dxy yc-dxy zSLM(m) 1]);
    slm.feed(Hm);
    prof = zeros(size(zScan));
    for j = 1:numel(zScan)
        sutter.moveZ(zScan(j));
        if j==1, pause(settle_first_s); else, pause(settle_step_s); end
        fr = grab_gated(comm, bas, pwr, nframes);
        prof(j) = max(fr(:));
    end
    ff = fit(zScan(:), prof(:), 'gauss1');   % peak = physical focus for this defocus
    focus_um(m) = ff.b1;
    fprintf('(7) SLM z=%+.3f -> focus at Sutter z=%+.1f um\n', zSLM(m), focus_um(m));
end
sutter.moveToRef()
pause(0.1)

pfit = polyfit(zSLM(:), focus_um(:), 1);
defocus_scale_um_per_unit = pfit(1);      % um focal shift per unit SLM-z
pupil_mapping.defocus_zSLM      = zSLM;
pupil_mapping.defocus_focus_um  = focus_um;
pupil_mapping.defocus_scale_um_per_unit = defocus_scale_um_per_unit;
fprintf('(7) defocus scale = %.1f um focal shift per unit SLM-z.\n', defocus_scale_um_per_unit);
disp(['    -> pupil diameter in SLM px follows from this scale + the Fresnel ' ...
      'kernel (psSLM=17um, f=' num2str(Setup.focal_SLM) 'm). Convert against ' ...
      'Genesis''s convention; raw scale stored.']);

%% ---- (8) PUPIL EDGE / CENTER via blaze at increasing radius -----------------
% Sweep a single spot outward along +x from the SLM center; the diffracted spot
% intensity falls as the blaze pushes it past the pupil edge. The radius where
% intensity collapses locates the pupil edge; symmetric sweeps (+/-x, +/-y) give
% the center. Reported in normalized SLM units AND SLM pixels (Nx=1024).
radii = linspace(0, 0.48, 25);          % normalized offset from center
dirs  = [1 0; -1 0; 0 1; 0 -1];         % +x -x +y -y
edge_norm = zeros(1, 4);
sweepI = zeros(numel(radii), 4);
for d = 1:4
    for r = 1:numel(radii)
        tgt = [xc yc] + radii(r)*dirs(d,:);
        fr = feed_and_grab(slm, bas, comm, pwr, ...
                function_Make_3D_SHOT_Holos(Setup, [tgt 0 1]), nframes);
        sweepI(r,d) = max(fr(:));
    end
    % edge = radius where intensity drops below half of its near-center value
    I0 = max(sweepI(1:3,d));
    idx = find(sweepI(:,d) < 0.5*I0, 1, 'first');
    if isempty(idx), edge_norm(d) = radii(end); else, edge_norm(d) = radii(idx); end
end
slm.feed(function_Make_3D_SHOT_Holos(Setup, [xc yc 0 1]));

% center = midpoint of opposing edges; diameter = sum of opposing half-extents
cx_norm = xc + (edge_norm(1) - edge_norm(2))/2;
cy_norm = yc + (edge_norm(3) - edge_norm(4))/2;
diam_norm = mean([edge_norm(1)+edge_norm(2), edge_norm(3)+edge_norm(4)]);
Nx = Setup.SLM.Nx;
pupil_center_px   = [cx_norm, cy_norm] * (Nx - 1) + 1;
pupil_diameter_px = diam_norm * (Nx - 1);

pupil_mapping.blaze_radii_norm    = radii;
pupil_mapping.blaze_intensity     = sweepI;
pupil_mapping.pupil_center_px     = pupil_center_px;
pupil_mapping.pupil_diameter_px   = pupil_diameter_px;

figure('Name','(8) Blaze rolloff'); clf
plot(radii, sweepI); xlabel('offset from center (norm SLM)'); ylabel('spot max DN')
legend({'+x','-x','+y','-y'}); title('Pupil edge = intensity collapse')
fprintf('(8) pupil center ~ [%.0f %.0f] px, diameter ~ %.0f px (of %d).\n', ...
    pupil_center_px(1), pupil_center_px(2), pupil_diameter_px, Nx);

%% ---- (9) CONFIRM existing SLM -> image affine -------------------------------
% Feed one spot at 5-9 SLM positions and compare the camera centroid against the
% existing calibration's prediction. Uses function_SLMtoSI + the latest calib.
% (This validates the SLM<->FOV map the acquisition relies on; README item 9.)
try
    calib_dir = rig_remote_get('paths.calib_dir', '');
catch
    calib_dir = '';
end
CoC = find_latest_calib(wavelength, calib_dir);

probe_pts = [0.35 0.35; 0.65 0.35; 0.5 0.5; 0.35 0.65; 0.65 0.65];
meas_cam = zeros(size(probe_pts));
pred_SI  = zeros(size(probe_pts));
for k = 1:size(probe_pts,1)
    meas_cam(k,:) = grab_spot_center(Setup, slm, bas, comm, pwr, [probe_pts(k,:) 0 1], nframes);
    si = function_SLMtoSI([probe_pts(k,:) 0], CoC);
    pred_SI(k,:) = si(1:2);
end
slm.feed(function_Make_3D_SHOT_Holos(Setup, [xc yc 0 1]));

pupil_mapping.affine_probe_slm   = probe_pts;
pupil_mapping.affine_meas_cam    = meas_cam;
pupil_mapping.affine_pred_si     = pred_SI;
disp('(9) SLM->image probe complete; compare affine_meas_cam vs affine_pred_si.')

%% ---- summary ----------------------------------------------------------------
fprintf('\n===== SLM -> back-pupil mapping (hand to Genesis) =====\n');
fprintf('  rotation_deg        : %.2f\n', pupil_mapping.rotation_deg);
fprintf('  mirror_flip         : %d\n',   pupil_mapping.mirror_flip);
fprintf('  pupil_center_px     : [%.0f, %.0f]\n', pupil_mapping.pupil_center_px);
fprintf('  pupil_diameter_px   : %.0f\n', pupil_mapping.pupil_diameter_px);
fprintf('  defocus_scale_um    : %.1f um / unit SLM-z\n', pupil_mapping.defocus_scale_um_per_unit);
% Persist the FULL struct (incl. raw half-block frames + blaze sweeps) to a .mat,
% so the measurements survive the session. acquire_bead_grid_stack.m folds only
% the scalar summary into the stack metadata from the workspace copy.
try
    data_root = rig_remote_get('paths.data_root', pwd);
catch
    data_root = pwd;
end
pm_dir = fullfile(data_root, 'ao_bead_calibration', datestr(now, 'yymmdd')); %#ok<TNOW1,DATST>
if ~isfolder(pm_dir), mkdir(pm_dir); end
pm_file = fullfile(pm_dir, sprintf('%s_pupil_mapping_%dnm.mat', datestr(now,'yymmdd'), wavelength)); %#ok<TNOW1,DATST>
save(pm_file, '-struct', 'pupil_mapping');
fprintf('Saved full pupil_mapping -> %s\n', pm_file);

comm.send('end', 'daq');   % release the DAQ power loop (alignCodeDAQ2K exits)

disp('Saved pupil_mapping .mat; acquire_bead_grid_stack.m loads its summary from');
disp('there. RESTART alignCodeDAQ2K on the DAQ before running the acquisition.');

%% ============================ local helpers =================================
function c = grab_spot_center(Setup, slm, bas, comm, pwr, coord4, nframes)
% Feed a single target, grab+average (gated), return spot center [row col].
fr = feed_and_grab(slm, bas, comm, pwr, function_Make_3D_SHOT_Holos(Setup, coord4), nframes);
[x, y] = function_findcenter(fr);
imagesc(fr)
c = [x, y];
end

function img = feed_and_grab(slm, bas, comm, pwr, Holo, nframes)
% Feed a hologram, then grab with the beam gated ON only for the grab.
slm.feed(Holo);
img = grab_gated(comm, bas, pwr, nframes);
end

function img = grab_gated(comm, bas, pwr, nframes)
% Gate the hologram ON (pwr in mW), grab+average, gate OFF.
power_gate(comm, pwr/1000);
img = mean(double(bas.grab(nframes)), 3);
power_gate(comm, 0);
end

function power_gate(comm, val_W)
% Command the DAQ (alignCodeDAQ2K) to set laser power to val_W (0 = off) and wait
% for its 'gotit' ack. Message [power_W, div, mult]; DAQ computes power_W*mult/div.
comm.send([val_W, 1, 1], 'daq');
invar = [];
while ~strcmp(invar, 'gotit')
    invar = comm.read(0.01);
end
end
