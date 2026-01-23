function alignCodeDAQ2K
%%Create DAQ Session
fprintf('Starting daq...\r')

fprintf('Loading defaults... ')
default_setup()
%%
power_calibrations;
%initalize contact
lpc = dictionary();
lpc('900') = LaserPowerControl(Output(DAQOutput(dq, 'port0/line5'), 'Shutter 900'),...
    ELL14(SerialInterface(sp), 1, 'Power 900'),...
     power_calibration.calibration_900, 10); % update these calibration paths as you get them...

lpc('1100') = LaserPowerControl(Output(DAQOutput(dq, 'port0/line4'), 'Shutter 1100'),...;
    ELL14(SerialInterface(sp), 2, 'Power 1100'),...
    power_calibration.calibration_1100, 10);

%calib_607 = 'C:\Users\holos\Documents\power-calibrations\240927_1030nm_100kHz_35AOM_uni_gate_calibration(607nm_first_order_suboptimal_dichroic).mat';
% calib_607 = 'C:\Users\holos\Documents\power-calibrations\241028_607nm_50kHz_35AOM_uni_gate_calibration.mat';
% laser_607 = LaserPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 607'),...
%     ELL14(SerialInterface(sp), 3, 'Power 607'), calib_607, 50); 

lpc('1030') = LaserPowerControl(Output(DAQOutput(dq, 'port0/line6'), 'Shutter 1030'),...
    LaserModulator(DAQOutput(dq, 'ao1'), 'Power 1030'),...
    power_calibration.calibration_1030A, 1250); % if we want to use this, make sure to do a more careful calibration on power...

%%
comm = HolochatInterface('daq'); % start here for preliminary info
calibration_wavelength = comm.read(10);


% create two more comms, one per power?
power_controllers = [];
for wv = calibration_wavelength
    comm(end+1) = HolochatInterface(wv); % 3 comms ... or more?
    power_controllers(end+1) = lpc(comm.id);
end

fprintf('Calibrating %dnm.\n', calibration_wavelength)
% initialize multiple of them, then dynamically switch between pwr as necessary (?)
% switch calibration_wavelength
%     case 900
%         pwr = laser_900;
%     case 607
%         pwr = laser_607;
%     case 1100
%         pwr = laser_1100;
%     case 1030
%         pwr = laser_1030;
% end

% pwr.initialize();

%% get data
disp('Waiting for Hologram info')
timeoutTime = 100000;

tstart = tic;
go =1;
while go;
    for c = comm
        % read the comms and determine which one to use
        invar = c.read(0.01);
        % we can get the correct one here by taking c's value
        if ~isempty(invar) && ~ischar(invar)
            pwr = power_controllers(c.id); % get the appropriate ID for the stuff
            fprintf('Update Power: ')
            PowerRequest = (invar(1)*invar(3))/invar(2);
            
            if PowerRequest == 0
                pwr.zero(); 
            else
                setpwr = min([(PowerRequest), pwr.max_pwr]);
                pwr.control.set(pwr.pwr_fun(setpwr));
                pwr.shutter.open()
            end
            
            fprintf('%0.8fmW\t', PowerRequest*1000);
            c.send('gotit', 'holo');
            
            fprintf(['Time since last run ' num2str(toc(tstart),2) 's\n']);
            tstart=tic;
        end
        
        if toc(tstart)>timeoutTime || strcmp(invar,'end');
            go=0;
            pwr.zero();
            c.send('kthx', 'holo');
        end
    end
end

for p = power_controllers
    p.zero()
end
disp('Ended')