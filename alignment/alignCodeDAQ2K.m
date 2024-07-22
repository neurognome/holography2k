function alignCodeDAQ2K
%%Create DAQ Session
fprintf('Starting daq...\r')

fprintf('Loading defaults... ')
%setup = getDefaults();  
fprintf('OK.\n')

fprintf('Making MATLAB NIDAQ object... ')
dq = daq('ni');
dq.Rate = 20000;%setup.daqrate;
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
power_calibrations;
%initalize contact
fpc_900 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(s), 1, 'Power 900'),...
     power_calibration.calibration_900, 250); % update these calibration paths as you get them...

fpc_1100 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line4'), 'Shutter 1100'),...
    ELL14(SerialInterface(s), 2, 'Power 1100'),...
    power_calibration.calibration_1100, 50);

fpc_1030 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    ELL14(SerialInterface(s), 3, 'Power 1030'),...
    power_calibration.calibration_1030H, 500);
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
            % pwr.hwp.moveto(22);
            pwr.power(min(PowerRequest, pwr.max_pwr)); % set power, check if mW or W
            pwr.open(); % open shutter
        end
        
        fprintf('%0.2fmW\t', PowerRequest*1000);
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