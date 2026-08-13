%% Power calibration for a MODULATOR power path (EOM / Pockels / AOM)
%
% Produces the power->volts LUT that LaserPowerControl reads for a channel whose
% rig file declares `kind = 'eom'`. Without this file such a channel has no
% pwr_fun, and gated_waveform holds it at its resting voltage for the whole
% trial: the laser delivers NO light and you get a clean recording of nothing.
%
% This is NOT the same file as an old LaserPower.mat -- that feeds the older
% power-scaling code and has a different shape. See holodaq/rigs/README.md.
%
% Sibling scripts: AutoLaserPowerCalib_HWP.m (rotator/half-wave-plate path),
% AutoLaserPowerCalib_MOD.m.

% Put this checkout and its holodaq on the path. Replaces three hardcoded
% addpaths naming one machine's username; see holo_paths.
addpath(fileparts(mfilename('fullpath')), '-end');
holo_paths();

clear
clc
close all;

% fill this in
wavelength = 1030; % 900, 1100, 1030
used_khz = 1250;
gate = 'uni';%'uni'; % or none or normal?

% Which fpc module in the rig file describes this modulator. The EOM path is
% conventionally declared as `fpc_eom` (that is what holodaq's ExampleRig shows
% and what the Millennium Phoenix rig file uses), but a rig is free to name it
% per-wavelength like the rotator channels -- set that name here if so.
fpc_module = 'fpc_eom';

% Load the rig BEFORE anything reads it: a rig_get above this line silently
% returns its fallback.
try
    rig = load_rig();
catch
    rig = struct();
end

% Where the power->volts LUTs go. From the rig so another scope needs no edit
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
interStepPause = 5; % seconds (0.5 for the fast head)

% The command range to sweep, in VOLTS. A modulator biased below 0 V may need
% v_min set to its resting voltage rather than 0 -- check what your driver
% actually accepts before widening this.
v_min = 0;
v_max = 5;

% Laser-head safety limit, and the two rep rates this protocol steps between:
% the search runs at search_khz and the real curve at used_khz, so the search
% readings are scaled by their ratio to decide what will be safe at full rate.
max_power = 0.5;    % W at the objective -- the head's limit
search_khz = 100;

%% start visa thing (older matlabs)
% note: use instrhwinfo to find the correct dev
% (holodaq is already on the path via holo_paths above)
%
% NOTE this is still the raw-VISA meter path, deliberately left alone.
% AutoLaserPowerCalib_HWP.m was migrated to get_power_meter(), and that
% migration turned out to be subtle -- updateReading() hands back whatever the
% meter has already accumulated rather than triggering a fresh averaged read,
% so it needs an explicit wait the raw 'read?' query did not. Doing the same
% here is worth doing, but not blind: it wants the meter in front of you.
instrreset()
vinfo = instrhwinfo('visa','ni');
v = eval(vinfo.ObjectConstructorName{1});
fopen(v);

%% initialize the thing

dq = daq('ni');

% Device, shutter line and modulator output all come from the rig. This used to
% be a `switch wavelength` block that hardcoded 'Dev1' and drove **ao1** for the
% 1030 branch -- on Millennium Phoenix the stim modulator is ao3, and ao1 is the
% imaging-path Pockels cell, so the old literal calibrated the wrong device. The
% 900/1100 branches built an ELL14 from an undefined `s`, which cannot be right
% in a modulator script at all; they are gone.
dev = rig_get('daq.device', 'Dev1');
ao_line = rig_get(sprintf('modules.%s.output', fpc_module), 'ao1');

% Volts that mean DARK on this modulator. NOT 0: on a modulator biased to, say,
% -0.375 V, 0 V is a step away from dark and on some modulators it is fully
% OPEN. Everything below parks the line here between readings.
rest_v = rig_get(sprintf('modules.%s.rest', fpc_module), 0);

% A modulator rig may have no shutter at all -- the modulator itself is what
% gates the light, which is the whole point of `kind = 'eom'`. Only bind and
% drive a shutter line when the rig declares one.
has_shutter = rig_has(rig, sprintf('modules.%s.shutter', fpc_module));
if has_shutter
    dq.addoutput(dev, rig_get(sprintf('modules.%s.shutter', fpc_module)), 'Digital');
end
dq.addoutput(dev, ao_line, 'Voltage');

% Shutter column first, then the modulator, matching the addoutput order above.
if has_shutter
    write_dark  = @()  dq.write([0, rest_v]);
    write_probe = @(x) dq.write([1, x]);
else
    write_dark  = @()  dq.write(rest_v);
    write_probe = @(x) dq.write(x);
end

fprintf('Devices connected: %s, modulator on %s (rest %g V)', dev, ao_line, rest_v);
if has_shutter
    fprintf(', with a shutter.\n');
else
    fprintf(', no shutter -- the modulator gates.\n');
end

%%
write_dark();

fprintf('Set the laser rep rate to %dkHz\n', search_khz);
disp('Running a full calibration across the whole range')
input('Press enter when ready: ')
%% initial search for low and high points
initial_search_queries = linspace(v_min, v_max, 60);

initial_search_values = zeros(size(initial_search_queries));

% start
for ii = 1:numel(initial_search_queries)
    initial_search_values(ii) = read_power(v, write_probe, write_dark, ...
        initial_search_queries(ii), nsamplesPM, interStepPause);
end
write_dark();

% first, we need to determine what is "safe", the limit of our laser head
% is 500mW, so let's keep it below that..
end_idx = find(initial_search_values * (used_khz/search_khz) < max_power*1000, 1, 'last');
assert(~isempty(end_idx), 'AutoLaserPowerCalib_EOM:noSafePoint', ...
    ['Every point in the search already exceeds the %g W head limit once scaled ' ...
     'to %dkHz,\nso there is no safe range to calibrate. Check the meter units ' ...
     'and the sweep range.'], max_power, used_khz);

%%
fprintf('Set the laser to %dkHz\n', used_khz);
disp('Now we will collect the full power curve');
input('Press enter to continue: ')

full_power_queries = initial_search_queries(1:end_idx);
full_power_values = zeros(size(full_power_queries));

% start
for ii = 1:numel(full_power_queries)
    full_power_values(ii) = read_power(v, write_probe, write_dark, ...
        full_power_queries(ii), nsamplesPM, interStepPause);
end
write_dark();

%%
disp('Now we will collect the beginning of the curve in full power mode.')
input('Press enter to continue: ')

fine_power_queries = linspace(v_min, initial_search_queries(end_idx), 60);
fine_power_values = zeros(size(fine_power_queries));
% start
for ii = 1:numel(fine_power_queries)
    fine_power_values(ii) = read_power(v, write_probe, write_dark, ...
        fine_power_queries(ii), nsamplesPM, interStepPause);
end
write_dark();


%% now estimate
poly = polyfit(initial_search_values(1:end_idx), full_power_values, 1);
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
% The command column is VOLTS for a modulator. LaserPowerControl.get_pwr_fun
% prefers a 'volts' field and only falls back to 'degrees', so naming it
% correctly here is what keeps an EOM LUT from having to lie about its units.
calib.volts = combined_v;
calib.powers = combined_power;
calib.max_volts = max(combined_v);
calib.min_volts = min(combined_v);
calib.max_power = max(combined_power);
calib.min_power = min(combined_power);
calib.rest = rest_v;
calib.khz = used_khz;

%% save
% 'eom' in the name so a modulator LUT is never mistaken for a rotator one:
% they have the same shape but their command column means completely different
% things, and both live in this folder.
fn = fullfile(save_base, sprintf('%s_%dnm_%dkHz_%s_eom_gate_calibration.mat', datetime('now', 'Format', 'yyMMdd'), wavelength, used_khz, gate));
save(fn, 'calib')

fprintf('Saved in %s\n', fn)


function val = read_power(v, write_probe, write_dark, query_v, nsamplesPM, settle)
%READ_POWER Drive the modulator to query_v, read the meter, return to dark.
%   Factored out because the three sweeps below ran three copies of this, and
%   two of them printed a variable that did not exist (`full_power_queries` in
%   loops that had never assigned it), so the second sweep died on its first
%   iteration and this script could never have completed a calibration.
    write_probe(query_v);
    pause(settle);

    fprintf(v, ['sense:average:count ', num2str(nsamplesPM)]);
    set(v, 'timeout', 3+1.1*nsamplesPM*3/1000)
    ret = query(v, 'read?');
    val = str2double(ret)*1000;     % meter reports W; work in mW
    if val < 0
        val = 0;
    end
    write_dark();
    fprintf('Volts: %s  Power (mW): %s\n', num2str(query_v), num2str(val));
end
