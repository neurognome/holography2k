%% What is tthe rep rate?

addpath(genpath('C:\Users\holos\Documents\GitHub'))
addpath(genpath('C:\Users\holos\Documents\_code'))

clear
clc
close all;

% fill this in
wavelength = 1030; % 900, 1100, 1030
used_khz = 5000;
gate = 'uni';%'uni'; % or none or normal?

save_base = 'C:\Users\holos\Documents\power-calibrations\';

%% start visa thing (older matlabs)
% note: use instrhwinfo to find the correct dev
addpath(genpath('C:\Users\holos\Documents\GitHub\holodaq'));
addpath(genpath('C:\Users\holos\Documents\ThorlabsPowerMeter'));

meter_list = ThorlabsPowerMeter;
tpm = meter_list.connect(meter_list.listdevices, 1);


tpm.setWaveLength(wavelength);
tpm.setAverageTime(0.1);

%% Params
interStepPause = 1; % seconds (0.5 for the fast head)

%% initialize the thing

dq = daq('ni');

switch wavelength
    case 900
        dq.addoutput('Dev1', 'port0/line5', 'Digital'); % 0/5: 920, 0/6: 1100
        hwp = ELL14(SerialInterface(s), 1, 'hwp'); % 0: 900, 1:1100 % might need both? idk
    case 1100
        dq.addoutput('Dev1', 'port0/line4', 'Digital'); % 0/5: 920, 0/6: 1100
        hwp = ELL14(SerialInterface(s), 2, 'hwp'); % 0: 900, 1:1100 % might need both? idk
    case 1030
        dq.addoutput('Dev1', 'port0/line4', 'Digital');
        dq.addoutput('Dev1', 'ao1', 'Voltage');
end
disp('Devices connected.')

%% 
dq.write([0, 0]); % close shutter

disp('Set the laser rep rate to 25kHz')
disp('Running a full calibration across the whole range')
input('Press enter when ready: ')
%% initial search for low and high points
initial_search_queries = linspace(0, 5, 60); % 0 to 5v

initial_search_values = zeros(size(initial_search_queries));

% start
for ii = 1:numel(initial_search_queries)
    dq.write([1, initial_search_queries(ii)]);

    tpm.updateReading(interStepPause);
    initial_search_values(ii) = tpm.meterPowerReading * 1000;
    dq.write([0, 0]);
    disp(['Deg: ' num2str(initial_search_queries(ii)) ' Power (mW):  ' num2str(val)])
end
dq.write([0, 0]);
% first, we need to determine what is "safe", the limit of our laser head
% is 500mW, so let's keep it below that..
max_power = 0.3;
end_idx = find(initial_search_values * (used_khz/25) < max_power*1000, 1, 'last');

%%
disp(['Set the laser to ', num2str(used_khz), ' kHz'])
disp('Now we will collect the full power curve');
input('Press enter to continue: ')

% start
for ii = 1:numel(initial_search_queries(1:end_idx))
    dq.write([1, initial_search_queries(ii)]);
    pause(interStepPause);

    tpm.updateReading(interStepPause);
    initial_search_values(ii) = tpm.meterPowerReading * 1000;

    if val < 0
        val = 0;
    end

    full_power_values(ii) = val;
    dq.write([0, 0]);
    disp(['Deg: ' num2str(initial_search_queries(ii)) ' Power (mW):  ' num2str(val)])
end
dq.write([0, 0]);

%%
disp('Now we will collect the beginning of the curve in full power mode.')
input('Press enter to continue: ')

fine_power_queries = linspace(0, initial_search_queries(end_idx), 60);
% start
for ii = 1:numel(fine_power_queries)
    dq.write([1, fine_power_queries(ii)]);
    pause(interStepPause);

    fprintf(v, ['sense:average:count ', num2str(nsamplesPM)]);
    set(v, 'timeout', 3+1.1*nsamplesPM*3/1000)
    ret = query(v, 'read?');
    val = str2double(ret)*1000;
    if val < 0
        val = 0;
    end
    fine_power_values(ii) = val;
    dq.write([0, 0]);
    disp(['Deg: ' num2str(fine_power_queries(ii)) ' Power (mW):  ' num2str(val)])
end
dq.write([0, 0]);


%% now estimate
poly = polyfit(initial_search_values(1:end_idx), full_power_values(1:end_idx), 1);
estimated_power = polyval(poly, initial_search_values); % the whole thing now

figure; 
plot(initial_search_queries, estimated_power);
hold on
plot(fine_power_queries, fine_power_values);

%
combined_v = [fine_power_queries, initial_search_queries(end_idx:end)];
combined_power = [fine_power_values, estimated_power(end_idx:end)];

%% find hi-low
combined_power = combined_power/1000; % in W
%% information
calib = struct();
calib.max_deg = min(combined_v); % because we know the values
calib.min_deg = max(combined_v);
calib.max_power = max(combined_power);
calib.min_power = min(combined_power);

calib.degrees = combined_v;
calib.powers = combined_power;
calib.khz = used_khz;

%% save
fn = fullfile(save_base, sprintf('%s_%dnm_%dkHz_%s_gate_calibration.mat', datetime('now', 'Format', 'yyMMdd'), wavelength, used_khz, gate));
save(fn, 'calib')

fprintf('Saved in %s\n', fn)