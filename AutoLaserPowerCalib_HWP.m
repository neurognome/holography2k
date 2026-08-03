%% What is tthe rep rate?

% Put this checkout and its holodaq on the path. Replaces three hardcoded
% addpaths naming one machine's username; see holo_paths.
addpath(fileparts(mfilename('fullpath')), '-end');
holo_paths();

clear
clc
close all;

% fill this in
wavelength = 900; % 900, 1100, 1030
used_khz = 750;
aom = 1;

gate = 'bi';%'uni'; % or none or normal?

% Where the power->angle LUTs go. From the rig so another scope needs no edit
% here; the literal is the fallback, so this rig behaves identically.
try
    save_base = rig_remote_get('paths.power_calib_dir', ...
        'C:\Users\holos\Documents\power-calibrations');
catch
    save_base = 'C:\Users\holos\Documents\power-calibrations';
end
if ~isfolder(save_base)
    mkdir(save_base);
end

%% Params
nsamplesPM = 1000; % counts (at 1000 Hz I think), produces an average
interStepPause = 0.5; % seconds (0.5 for the fast head)

% The meter's averaging window, in seconds -- setAverageTime takes a DURATION, not
% a sample count. Also the minimum we have to wait after moving the rotator before
% a reading means anything; see the read in the sweep loop below.
avg_time_s = nsamplesPM/1000;

%% initialize the thing
% The power meter is connected after the switch, since connecting needs the
% wavelength. get_power_meter owns the listdevices/connect handshake.

dq = daq('ni');

% Load the rig BEFORE anything reads it: a rig_get above this line silently
% returns its fallback.
try, load_rig(); catch, end

% The rotator bus. Resolved the way PowerControllerCalibrated does it
% (PowerControllerCalibrated.m:242-244) -- via the module's 'serial' field --
% rather than naming a bus here. This script used to name rig.serial.hwp directly,
% which had drifted to a stale COM5 while the GUI drove these very same rotators on
% rig.serial.ell14 (COM4): one physical bus, two declarations, and only the GUI's
% was maintained. Going through the module means the two cannot disagree again.
%
% NOTE the bus is exclusive. ScopeController / PowerControllerCalibrated holds it
% open, along with the shutter line below, so close that GUI before running this.
ell14_cfg = struct('port', 'COM4', 'baud', 9600, 'byte_order', 'big-endian', ...
    'parity', 'none', 'stop_bits', 1, 'data_bits', 8, 'terminator', 'CR/LF');
sname = rig_get(sprintf('modules.fpc_%d.serial', wavelength), 'ell14');
s = open_serial(rig_get(['serial.' sname], ell14_cfg));

% Shutter line and rotator address per laser. The literals are the fallback and are
% what this script has always used; the rig wins wherever it declares the module
% (Scope2KRig declares fpc_900 and fpc_1100, but not fpc_1030).
switch wavelength
    case 900
        shutter_line = 'port0/line5'; default_ch = 1;
    case 1100
        shutter_line = 'port0/line4'; default_ch = 2;
    case 1030
        shutter_line = 'port0/line6'; default_ch = 3;
    otherwise
        error('AutoLaserPowerCalib_HWP:badWavelength', ...
            'No shutter/rotator mapping for %dnm (900, 1100, 1030).', wavelength);
end
dq.addoutput('Dev1', ...
    rig_get(sprintf('modules.fpc_%d.shutter', wavelength), shutter_line), 'Digital');
% The 1100 branch never built this at all, so hwp.set(0) below died on the very
% wavelength this script defaults to.
hwp = ELL14(SerialInterface(s), ...
    rig_get(sprintf('modules.fpc_%d.ell14_channel', wavelength), default_ch), 'hwp');
tpm = get_power_meter(wavelength);
%tpm.setAverageTime(avg_time_s);                        % SECONDS
tpm.setTimeout(1000 * (3 + 1.1*nsamplesPM*3/1000));    % MILLISECONDS
disp('Devices connected.')

%% 
hwp.set(0); % start at 0
dq.write(0); % close shutter

%%
% disp('Do this calibration with the z-block in place')
% disp('On the Holography Computer put a hologram near the zero order')
% disp('Put Laser Gate on Bypass') %added 3/18/21

%% initial search for low and high points
initial_search_queries = linspace(0, 120, 60*2); % 0 to 70

initial_search_values = zeros(size(initial_search_queries));

% start
dq.write(1);
for ii = 1:numel(initial_search_queries)
    hwp.set(initial_search_queries(ii)); % move to deg

    % updateReading measures FIRST and pauses after, so it hands back whatever the
    % meter has already accumulated -- it does not trigger a fresh averaged read the
    % way the raw-VISA 'read?' query did. Wait out the rotator settle AND a full
    % averaging window here, or every point is the previous angle's power.
    pause(interStepPause + avg_time_s);
    tpm.updateReading(0);
    assert(strcmpi(tpm.meterPowerUnit, 'W'), ...
        'Meter is reporting %s, not W -- the calibration would be garbage.', ...
        tpm.meterPowerUnit);

    val = tpm.meterPowerReading * 1000; % in mW
    if val < 0
        val = 0;
    end
    initial_search_values(ii) = val;
    disp(['Deg: ' num2str(initial_search_queries(ii)) ' Power (mW):  ' num2str(val)])
end
dq.write(0);
tpm.disconnect(); % release the head, so re-running this script can connect again

%% find hi-low
[max_pwr, max_idx] = max(initial_search_values(3:end-3)); %exclude ends
[min_pwr, min_idx] = min(initial_search_values(3:end-3));

max_idx = max_idx + 2;
min_idx = min_idx + 2;

segment = linspace(min_idx, max_idx, abs(max_idx-min_idx) + 1);
plot(initial_search_queries(segment), initial_search_values(segment));
fprintf('Initial Search Report\n Max: %0.2fmW at %0.2fdeg \n Min: %0.2fmW and %0.2fdeg\n',...
    max_pwr, initial_search_queries(max_idx), min_pwr, initial_search_queries(min_idx));


initial_search_values = initial_search_values/1000; % convert to W


%% interp only
%% information
calib = struct();
calib.max_deg = initial_search_queries(max_idx);
calib.min_deg = initial_search_queries(min_idx);
calib.max_power = initial_search_values(max_idx);
calib.min_power = initial_search_values(min_idx);

calib.degrees = initial_search_queries(segment);
calib.powers = initial_search_values(segment);
calib.khz = used_khz;

%% save
fn = fullfile(save_base, sprintf('%s_%dnm_%dkHz_%dAOM_%s_gate_calibration_temp.mat', datetime('now', 'Format', 'yyMMdd'), wavelength, used_khz, aom*100, gate));
save(fn, 'calib')
save(sprintf('active_%d.mat', wavelength), 'calib')

fprintf('Saved in %s\n', fn)