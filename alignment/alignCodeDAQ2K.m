function alignCodeDAQ
%%Create DAQ Session

%initalize contact
% IP='128.32.173.87';
% IP='128.32.173.209';
fpc_900 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(s), 0, 'Power 900'),...
    'C:\Users\holos\Documents\power-calibrations\230918_900nm_100kHz_25AOM_fast_gate_calibration.mat');

fpc_1100 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 0, 'Power 1100'),...
    'C:\Users\holos\Documents\power-calibrations\230918_1100nm_100kHz_25AOM_fast_gate_calibration.mat');

fpc_1030 = FiberPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 0, 'Power 1100'),...
    'C:\Users\holos\Documents\power-calibrations\230918_1100nm_100kHz_25AOM_fast_gate_calibration.mat');

[HoloSocket]=msocketPrep();

% choose the appropriate one here?
calibration_wavelength = msrecv(HoloSocket, 5) % long timeout

switch calibration_wavelength
    case 900
        pwr = fpc_900;
    case 1100
        pwr = fpc_1100;
    case 1030
        pwr = fpc_1030;
end

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