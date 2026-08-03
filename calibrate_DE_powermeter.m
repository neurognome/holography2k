%% set parameters
% This script did no pathing of its own, so it only ran in a session something
% else had already set up. rig_remote_get below lives in holodaq. See holo_paths.
holo_paths();

wavelength = 1030;

% Not used by the calibration itself: a COPY of this file is written at the
% bottom with the DE info appended. The folder now comes from the rig config the
% DAQ published (this machine has no rig file of its own), so another scope needs
% no edit here; the literal is the fallback, so this rig behaves identically.
try
    calib_dir = rig_remote_get('paths.calib_dir', 'C:\Users\holos\Documents\calibs');
catch
    calib_dir = 'C:\Users\holos\Documents\calibs';
end

% Deliberately NOT find_latest_calib, which is the rig-driven calib resolver
% everywhere else in this repo. Its glob is '*_Calib_<nm>*.mat', and that also
% matches THIS script's own output, '<base>_DE_calib.mat'. Reusing it would make
% a second run take the previous run's product as its base and write
% '<base>_DE_calib_DE_calib.mat', compounding a suffix every time. Excluding the
% outputs is the entire difference, so the lookup is inline rather than a flag
% threaded through the shared resolver.
files = dir(fullfile(calib_dir, sprintf('*_Calib_%d*.mat', wavelength)));
if ~isempty(files)
    files = files(~contains({files.name}, '_DE_calib'));
end
assert(~isempty(files), 'calibrate_DE_powermeter:noBaseCalib', ...
    ['No base *_Calib_%d*.mat in %s (excluding this script''s own _DE_calib ' ...
     'outputs).\nRun the SLM/CoC calibration for %dnm first, or point ' ...
     'rig.paths.calib_dir at the right folder.'], wavelength, calib_dir, wavelength);
[~, idx] = max([files.datenum]);           % newest by file modification date
base_calib = fullfile(files(idx).folder, files(idx).name);
fprintf('base calib %dnm: %s\n', wavelength, base_calib);

slmXrange = [0, 1];
slmYrange = [0, 1];
slmZrange = 0*[-0.2 0.2];

% discretization levels for each dimension:
numX = 31;
numY = 31;
numZ = 1;  % make this odd number

% make sure all odd:
assert(mod(numX,2)==1)
assert(mod(numY,2)==1)
assert(mod(numZ,2)==1)
%% setup SLM

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;

slm = get_slm(wavelength);

slm.stop();
slm.wait_for_trigger = 0;
slm.start();

slmCoords = [.4, .4, 0, 1];
[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
slm.feed(Holo);

%% load powermeter reader
% get_power_meter owns the listdevices/connect handshake -- the bare constructor
% only makes a discovery stub, and set_wavelength/read were never real methods.
tpm = get_power_meter(wavelength);

% One reading up front, purely to populate meterPowerUnit and check it. A head left
% in dBm returns plausible-looking numbers, and this sweep is ~15 min of them.
tpm.updateReading(0);
assert(strcmpi(tpm.meterPowerUnit, 'W'), ...
    'Meter is reporting %s, not W -- DE_calib.powers must be in watts.', ...
    tpm.meterPowerUnit);

disp('Turn on FS50, set FS50 parameters to something not crazy (e.g, a few mW)')
input('Press enter when ready: ')
%% create SLM coordinate meshgrid
x = linspace(slmXrange(1), slmXrange(2), numX);
y = linspace(slmYrange(1), slmYrange(2), numY);
z = linspace(slmZrange(1), slmZrange(2), numZ);

% [x,y,z] = meshgrid(x,y,z);
% x = x(:);
% y = y(:);
% z = z(:);

%% loop through SLM positions and record power
powers = zeros(numX, numY, numZ);
count = 1;
tic
for i = 1:numX
    for j = 1:numY
        for k = 1:numZ
            slmCoords = [x(i), y(j), z(k), 1]; % 0.
            [Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
            slm.feed(Holo);
            
            pause(.25)
            tpm.updateReading(0);   % writes meterPowerReading; returns nothing
            pwr = tpm.meterPowerReading;   % watts, as DE_calib.powers expects

            powers(i, j, k) = pwr;
            
            count = count + 1;

            time_per_iter = toc/count;
            disp(['Elapsed: ', num2str(toc), ' sec'])
            disp(['Remaining (est.): ', num2str((numX*numY*numZ-count)*time_per_iter), ' sec'])
            
            if numZ == 1
                figure(123)
                imagesc(powers)
                colorbar
            end
            
        end
    end
    i
end
tpm.disconnect(); % release the head, so a re-run can connect again

%% save this calibration data:

load(base_calib)

DE_calib.slmXgrid = x;
DE_calib.slmYgrid = y;
DE_calib.slmZgrid = z;
DE_calib.powers = powers;

CoC.DE_calib = DE_calib;

save([base_calib(1:end-4), '_DE_calib.mat'], 'CoC')