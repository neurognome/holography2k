%clear; clc

% This repo, its holodaq, msocket and the SLM SDK. Replaces the two literals that
% used to sit further down (this checkout's path, and one machine's Desktop
% meadowlark folder) -- makePaths now takes the SDK folder from
% rig.paths.slm_sdk, so the SDK on the path is the same one
% start_holo_listener uses. It also cd's to the repo root, which is harmless
% here: nothing below opens a relative path.
%
% HOISTED to the top. It used to run after HolochatInterface and
% function_loadparameters2 had already been called, both of which need holodaq
% and this repo on the path -- so the script only worked in a session something
% else had set up, and pathing itself did nothing for its first two lines of real
% work. rig_remote_get below needs holodaq too.
makePaths();

wavelength = [1100, 900];%[1030];%[1030, 607]; %[900];% 607;%[1100, 900]; %[1100,  900];%[1100, 900];%[1100, 900]; % combinations: 900, 1030, 1100, 900/1100, 900/1030

test_z = false;
comm = HolochatInterface('holo');

timeout = 1700;

Setup = function_loadparameters2();
%Setup = function_loadparameters3(); %FOR YASAP IMAGING
% Hologram settings from the rig config the DAQ published, with the same fields
% and the same fallbacks start_holo_listener uses (2 = GSS, GPU on, 1700 ms), so
% the listener and this hand-run path cannot disagree about how holograms are
% compiled. The literals are the fallback, not the default.
try
    Setup.CGHMethod = rig_remote_get('holo.cgh_method', 2);
    Setup.useGPU    = double(logical(rig_remote_get('holo.use_gpu', true)));
    timeout         = rig_remote_get('holo.slm_timeout_ms', timeout);
catch
    Setup.CGHMethod = 2;      % GSS
    Setup.useGPU    = 1;
end
Setup.verbose = 0;
if wavelength == 1030
    %Setup.focal_SLM = .25;
    %Setup.focal_SLM = .2;
end

cycleiterations = 1; % Change this number to repeat the sequence N times instead of just once

%Overwrite delay duration
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
Setup.SLM.timeout_ms = timeout;     %No more than 2000 ms until time out
% Setup.calib = 'C:\Users\holos\Documents\calibs\ActiveCalib.mat'; % here we need to somehow feed multiple calibrations?

%%
% Calibrations now come from find_latest_calib, which is what this switch was
% replaced by: it takes the folder from rig.paths.calib_dir (via the config the
% DAQ published, since this machine has no rig file) and picks the NEWEST
% *_Calib_<nm>*.mat in it. Three things that used to be hand-maintained here are
% now automatic:
%   * the folder is no longer one machine's username;
%   * a fresh calibration is picked up without editing this file, which is the
%     whole point -- the dated filenames above had to be re-pasted by hand every
%     time anyone recalibrated;
%   * the 1030 rmfield of FitX/FitY/FitZ that was commented out here is done by
%     find_latest_calib's local_clean, so 1030 no longer needs a special case.
%     For 1030 this resolves to whichever *_Calib_1030*.mat is newest, which
%     includes the '_DE_calib' variant this switch used to name explicitly.
%
% 589 is the one case the wavelength->file rule cannot express: it has no
% calibration of its own and deliberately borrows 900's.
calib = [];
for w = wavelength
    lookup = w;
    if w == 589
        lookup = 900;   % no 589 calibration exists; 900's is used on purpose
    end
    calib = [calib, find_latest_calib(lookup)];
end

sequences = {};
for w = 1:numel(wavelength)
    fprintf('Waiting for holorequest for %dnm...\n', wavelength(w))
    % hololist = generate_holograms2D(comm, Setup, calib(w));

    %%%%%%%%%% test z-offset: start %%%%%%%%%%%%

    if test_z
        clear myloc;
        myloc = [  387   152     0
            210   484     0
            226   409     0
            256   274     0
            433   330     0
            313   394     0
            294   148     0];  %213   244     0
        if exist("myloc",'var')
            myOptionsIn = cell(1,2);
            myOptionsIn{1} = 'loc';
            myOptionsIn{2} = myloc;
        else
            myOptionsIn = {};
        end
        [hololist,locs] = generate_holograms_tuneZoffset(comm, Setup, calib(w),...
            'z offset',-11+12,'xy offset',[8-20+11,0-3,0],myOptionsIn{:}); % 'loc',[265,270,0],   ,'loc',myloc  test on Feb 28, 2025
        %                 Direction on Cam: [down, left]
        % -165
    else
        hololist = generate_holograms(comm, Setup, calib(w));

    end
    %%%%%%%%%% test z-offset: end %%%%%%%%%%%%


    sequences{end+1} = uint8(hololist);
end
if test_z
    xyc = [253,232,0];
    r2use = 24;
    flagBlock = vecnorm(locs - xyc,2,2) < r2use;
    if any(flagBlock)
        warning(['Possible blocking of (',num2str(sum(flagBlock)),') focus spot(s).']);
    end
end


% fprintf('Waiting for holorequest for 900nm...\n')
% hololist_900 = generate_holograms(control, Setup, CoC_900);
% fprintf('Waiting for holorequest for 1100nm...\n')
% hololist_1100 = generate_holograms(control, Setup, CoC_1100);
fprintf('All holorequests received.\n')

%totally remove sequences as a thing that exists basically
% sequences = {uint8(hololist_900), uint8(hololist_1100)}; %shouldn't change anything added 9/14/21
% slm = [get_slm(900), get_slm(1100)];
slm = [];
for w = wavelength
    slm = [slm, get_slm(w)];
end

% flushMSocket(masterSocket)
comm.flush();
%%
% slm = HoloeyePLUTO();
for s = slm
    s.stop();
    s.wait_for_trigger = 0; % set settintgs
    s.timeout_ms = timeout;

    s.start();
end

for ii = 1:numel(slm)
    slm(ii).feed(sequences{ii});
end
%
%
% orderBackup=[]; %Sequence list is archived in case the daq errors. normally disposed of after exp. 1/19/21
% c=1;
% while true
%     orderBackup{c} = ShootSequencesMsocket2K(slm, sequences, comm);
%     c=c+1;
% end
%

