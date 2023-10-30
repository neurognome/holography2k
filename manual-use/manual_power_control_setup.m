%%Create DAQ Session
fprintf('Starting daq...\r')

fprintf('Loading defaults... ')
setup = getDefaults();  
fprintf('OK.\n')

fprintf('Making MATLAB NIDAQ object... ')
dq = daq('ni');
dq.Rate = setup.daqrate;
fprintf('OK.\n')

fprintf('Making SerialPort object... ')
s = serialport("COM5", ...
    9600,...
    'ByteOrder', 'big-endian',...
    'Parity', 'none',...
    'StopBits', 1,...
    'DataBits', 8);
s.configureTerminator('CR/LF');
fprintf('OK.\n')
%%

%initalize contact
fpc_900 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(s), 1, '_Power 900'),...
    'C:\Users\holos\Documents\power-calibrations\231019_1030nm_50kHz_25AOM_none_gate_calibration.mat'); % update these calibration paths as you get them...

fpc_1100 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line4'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 2, 'Power 1100'),...
    'C:\Users\holos\Documents\power-calibrations\231019_1030nm_50kHz_25AOM_none_gate_calibration.mat');

fpc_1030 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    ELL14(SerialInterface(s), 3, 'Power 1030'),...
    'C:\Users\holos\Documents\power-calibrations\231019_1030nm_50kHz_25AOM_none_gate_calibration.mat');

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