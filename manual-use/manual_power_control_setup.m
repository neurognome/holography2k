%%Create DAQ Session
% Put this checkout and its holodaq on the path (see holo_paths).
addpath(fileparts(fileparts(mfilename('fullpath'))), '-end');
holo_paths();

% Load the rig ONCE, up front, before anything reads it. rig_get resolves against
% the rig loaded this session, so a rig_get before this line silently returns its
% fallback -- which on this rig happens to equal the real value, making the mistake
% invisible.
try, load_rig(); catch, end

fprintf('Starting daq...\r')

fprintf('Making MATLAB NIDAQ object... ')
% Rate from the rig. This called getDefaults(), which exists in none of the four
% repos -- it belonged to the K:\ scaffolding this rig file replaced -- so line 5
% would have errored before anything else ran.
dq = daq('ni');
dq.Rate = rig_get('daq.rate', 20000);
fprintf('OK.\n')

fprintf('Making SerialPort object... ')
% The rotator bus, from rig.serial.hwp. open_serial applies every field
% (ByteOrder/Parity/StopBits/DataBits) and the terminator, so this is the same
% device config as before with the port no longer hardcoded.
hwp_cfg = struct('port', 'COM5', 'baud', 9600, 'byte_order', 'big-endian', ...
    'parity', 'none', 'stop_bits', 1, 'data_bits', 8, 'terminator', 'CR/LF');
s = open_serial(rig_get('serial.hwp', hwp_cfg));
fprintf('OK.\n')
%%

% Power LUTs from the rig, so this follows whatever the AutoLaserPowerCalib_*
% scripts most recently produced. All three channels previously shared ONE pinned
% 2023 1030nm LUT, under a comment saying to update them as calibrations arrived --
% so 900 and 1100 were being driven by a 1030nm calibration.
%
% fpc_1030 has no module in Scope2KRig (the rig declares fpc_900 and fpc_1100
% only), so it falls back to the legacy literal, which is the only record of that
% path. Declare a modules.fpc_1030 in the rig file to bring it under the same roof.
LEGACY_LUT = ['C:\Users\holos\Documents\power-calibrations\' ...
              '231031_1030nm_25kHz_30AOM_uni_gate_calibration.mat'];
cal_900  = rig_get('modules.fpc_900.calibration',  LEGACY_LUT);
cal_1100 = rig_get('modules.fpc_1100.calibration', LEGACY_LUT);
cal_1030 = rig_get('modules.fpc_1030.calibration', LEGACY_LUT);

% Say which LUT each channel got, and check it BEFORE constructing: an empty or
% missing LUT makes FiberPowerControl fail inside get_pwr_fun (scale =
% current_khz/calib.khz on an empty calib), which is an opaque way to learn that a
% file is not there.
for c = {'900', cal_900; '1100', cal_1100; '1030', cal_1030}'
    fprintf('  LUT %-4s : %s\n', c{1}, c{2});
    if isempty(c{2}) || exist(c{2}, 'file') ~= 2
        warning('manual_power_control_setup:noLut', ...
            ['Power LUT for %s is not a readable file:\n  %s\n' ...
             'FiberPowerControl will throw when it importdata''s this. Run the ' ...
             'relevant\nAutoLaserPowerCalib_* script, or fix ' ...
             'rig.modules.fpc_%s.calibration.'], c{1}, c{2}, c{1});
    end
end

%initalize contact
fpc_900 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(s), 1, '_Power 900'),...
    cal_900);

fpc_1100 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line4'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 2, 'Power 1100'),...
    cal_1100);

fpc_1030 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    ELL14(SerialInterface(s), 3, 'Power 1030'),...
    cal_1030);

fpc_900.initialize();
fpc_1100.initialize();
fpc_1030.initialize();

%%  examples
% change 900 path power

% if you have calibration
fpc_900.power(10); % set to 10mW

% if you dont
fpc_900.hwp.moveto(10); % move hwp to 10 deg

%
fpc_900.open(); % open shuttre for 900
fpc_900.close_all(); % close all shutters