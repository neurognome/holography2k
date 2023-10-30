function alignCodeDAQ2K
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

fpc_1100 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 2, 'Power 1100'),...
    "C:\Users\holos\Documents\power-calibrations\231019_1030nm_50kHz25AOM_none_gate_calibration.mat");

fpc_1030 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    ELL14(SerialInterface(s), 3, 'Power 1100'),...
    "C:\Users\holos\Documents\power-calibrations\231019_1030nm_50kHz_25AOM_none_gate_calibration.mat");
%%

[HoloSocket]=msocketPrep();
calibration_wavelength = msrecv(HoloSocket, 5); % long timeout

fprintf('Calibrating %dnm.\n', calibration_wavelength)
switch calibration_wavelength
    case 900
        pwr = fpc_900;
    case 1100
        pwr = fpc_1100;
    case 1030
        pwr = fpc_1030;
end

pwr.initialize(); % add channels to daq

%% get data
disp('Waiting for Hologram info')
timeoutTime = 100000;

tstart = tic;
go =1;
while go;
    invar=msrecv(HoloSocket,.01); %order: power, DE, nTargets
    if ~isempty(invar) && ~ischar(invar)
        fprintf('Update Power: ')
        PowerRequest = (invar(1)*invar(3))/invar(2);
        
        if PowerRequest == 0
            pwr.zero()
        else
            pwr.power(PowerRequest); % set power, check if mW or W
            pwr.open(); % open shutter
        end
        
        mssend(HoloSocket,'gotit');
        
        fprintf(['Time since last run ' num2str(toc(tstart),2) 's\n']);
        tstart=tic;
    end
    
    if toc(tstart)>timeoutTime || strcmp(invar,'end');
        go=0;
        pwr.zero()
        mssend(HoloSocket,'kthx');
    end
end

pwr.zero()
disp('Ended')