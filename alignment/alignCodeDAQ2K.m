function alignCodeDAQ2K
%%Create DAQ Session
fprintf('Starting daq...\r')

fprintf('Loading defaults... ')
default_setup()
%%
power_calibrations;
%initalize contact
laser_900 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(sp), 1, 'Power 900'),...
     power_calibration.calibration_900, 50); % update these calibration paths as you get them...

laser_1100 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line4'), 'Shutter 1100'),...
    ELL14(SerialInterface(sp), 2, 'Power 1100'),...
    power_calibration.calibration_1100, 50);

%calib_607 = 'C:\Users\holos\Documents\power-calibrations\240927_1030nm_100kHz_35AOM_uni_gate_calibration(607nm_first_order_suboptimal_dichroic).mat';
calib_607 = 'C:\Users\holos\Documents\power-calibrations\240927_1030nm_100kHz_35AOM_uni_gate_calibration(607nm_first_order).mat';
laser_607 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 607'),...
    ELL14(SerialInterface(sp), 3, 'Power 607'), calib_607, 100); 

laser_1030 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    LaserModulator(DAQOutput(dq, 'ao1'), 'Power 1030'),...
    power_calibration.calibration_1030A, 1250); % if we want to use this, make sure to do a more careful calibration on power...

%%
comm = HolochatInterface('daq');
calibration_wavelength = comm.read(10);

fprintf('Calibrating %dnm.\n', calibration_wavelength)
switch calibration_wavelength
    case 900
        pwr = laser_900;
    case 607
        pwr = laser_607;
    case 1100
        pwr = laser_1100;
    case 1030
        pwr = laser_1030;
end

pwr.initialize();

%% get data
disp('Waiting for Hologram info')
timeoutTime = 100000;

tstart = tic;
go =1;
while go;
    invar = comm.read(0.01);
    if ~isempty(invar) && ~ischar(invar)
        fprintf('Update Power: ')
        PowerRequest = (invar(1)*invar(3))/invar(2);
        
        if PowerRequest == 0
            if calibration_wavelength == 1030
                dq.write([0, 0]);
            else
                pwr.zero(); 
            end
        else
            if calibration_wavelength == 1030
                dq.write([0, min(pwr.pwr_fun(PowerRequest), pwr.max_pwr)]);
            else
                dq.write(1);
                setpwr = min([(PowerRequest), pwr.max_pwr]);
                pwr.control.set(pwr.pwr_fun(setpwr));
            end
        end
        
        fprintf('%0.8fmW\t', PowerRequest*1000);
        comm.send('gotit', 'holo');
        
        fprintf(['Time since last run ' num2str(toc(tstart),2) 's\n']);
        tstart=tic;
    end
    
    if toc(tstart)>timeoutTime || strcmp(invar,'end');
        go=0;
        if calibration_wavelength == 1030
            dq.write([0, 0]);
        else
            pwr.zero();
        end
        comm.send('kthx', 'holo');
    end
end

pwr.zero()
disp('Ended')