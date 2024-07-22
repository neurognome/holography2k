%% What is tthe rep rate?

addpath(genpath('C:\Users\holos\Documents\GitHub'))
addpath(genpath('C:\Users\holos\Documents\_code'))

clear
clc
close all;

% fill this in
wavelength = 1030; % 870, 900, 1030
used_khz = 250;
% used_khz = 50;
% aom = 0.30;
% gate = 'uni';%'uni'; % or none or normal?

save_base = 'C:\Users\holos\Documents\power-calibrations\';

%% start visa thing (older matlabs)
% note: use instrhwinfo to find the correct dev
addpath(genpath('C:\Users\holos\Documents\GitHub\holodaq'));

instrreset()
vinfo = instrhwinfo('visa','ni');
v = eval(vinfo.ObjectConstructorName{1});
fopen(v);

%% Params
nsamplesPM = 1000; % counts (at 1000 Hz I think), produces an average
interStepPause = 5; % seconds (0.5 for the fast head, 5 for the thermal)

%% initialize the thing

dq = daq('ni');

switch wavelength
    case 900
        dq.addoutput('Dev1', 'port0/line5', 'Digital'); % 0/5: 920, 0/6: 1100
        dq.addoutput('Dev1', 'ao0', 'Voltage');
        max_v = 1.5;
    case 870
        dq.addoutput('Dev1', 'port0/line5', 'Digital'); % 0/5: 920, 0/6: 1100
        dq.addoutput('Dev1', 'ao0', 'Voltage');
        max_v = 1.5;
    case 1030
        dq.addoutput('Dev1', 'port0/line6', 'Digital');
        dq.addoutput('Dev1', 'ao1', 'Voltage');
        max_v = 5;
end
disp('Devices connected.')

%% 
dq.write([0, 0]); % start at 0

%%
% disp('Do this calibration with the z-block in place')
% disp('On the Holography Computer put a hologram near the zero order')
% disp('Put Laser Gate on Bypass') %added 3/18/21

%% initial search for low and high points
initial_search_queries = linspace(0, max_v, 60); % 0 to 70

initial_search_values = zeros(size(initial_search_queries));

% start
dq.write([1, 0]);
for ii = 1:numel(initial_search_queries)
    dq.write([1, initial_search_queries(ii)])
    pause(interStepPause);

    fprintf(v, ['sense:average:count ', num2str(nsamplesPM)]);
    set(v, 'timeout', 3+1.1*nsamplesPM*3/1000)
    ret = query(v, 'read?');
    val = str2double(ret)*1000;
    if val < 0
        val = 0;
    end
    initial_search_values(ii) = val;
    disp(['Deg: ' num2str(initial_search_queries(ii)) ' Power (mW):  ' num2str(val)])
end
dq.write([0, 0]);

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
fn = fullfile(save_base, sprintf('%s_%dCnm_%dkHz_calibration_MOD.mat', datetime('now', 'Format', 'yyMMdd'), wavelength, used_khz));
save(fn, 'calib')
fn2 = fullfile(save_base, sprintf('active_%daom.mat', wavelength));
fprintf('Saved in %s\n', fn)