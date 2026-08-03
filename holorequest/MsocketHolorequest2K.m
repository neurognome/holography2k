function MsocketHolorequest2K()
% choose wavelengths
%clear; clc
wavelength = [1100 900];%[1100, 900]; % combinations: 900, 1030, 1100, 900/1100, 900/1030

% This repo, its holodaq, msocket and the SLM SDK -- the SDK from
% rig.paths.slm_sdk rather than one machine's Desktop. Hoisted above
% HolochatInterface and function_loadparameters2 below, which both need holodaq
% and this repo already on the path; the two addpath literals this replaces sat
% AFTER them, so the function only ever worked in a session something else had set
% up. rig_remote_get below needs holodaq too.
makePaths();
%%

comm = HolochatInterface('holo');
    

timeout = 1700;

% Setup = function_loadparameters3(); % previously, we ran this one, but
% now let's try 2...

Setup = function_loadparameters2();
% Hologram settings from the rig config the DAQ published, same fields and same
% fallbacks start_holo_listener uses (2 = GSS, GPU on, 1700 ms), so the listener
% and this entry point cannot disagree about how holograms are compiled. The
% literals are the fallback, not the default.
try
    Setup.CGHMethod = rig_remote_get('holo.cgh_method', 2);
    Setup.useGPU    = double(logical(rig_remote_get('holo.use_gpu', true)));
    timeout         = rig_remote_get('holo.slm_timeout_ms', timeout);
catch
    Setup.CGHMethod = 2;      % GSS
    Setup.useGPU    = 1;
end

Setup.verbose = 0;
 
cycleiterations = 1; % Change this number to repeat the sequence N times instead of just once

%Overwrite delay duration
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
Setup.SLM.timeout_ms = timeout;     %No more than 2000 ms until time out
% Setup.calib = 'C:\Users\holos\Documents\calibs\ActiveCalib.mat'; % here we need to somehow feed multiple calibrations?

%% 
% Calibrations from find_latest_calib, which was written to replace exactly this
% switch (see its help) but was never actually wired in here. It takes the folder
% from rig.paths.calib_dir via the config the DAQ published, and picks the NEWEST
% *_Calib_<nm>*.mat in it, so a fresh calibration no longer means re-pasting a
% dated filename into this function. The 1030 rmfield of FitX/FitY/FitZ below is
% done by its local_clean, so 1030 needs no special case here either.
%
% 589 is the one case a wavelength->file rule cannot express: it has no
% calibration of its own and deliberately borrows 900's.
calib = [];
for w = wavelength
    lookup = w;
    if w == 589
        lookup = 900;   % no 589 calibration exists; 900's is used on purpose
    end
    calib = [calib, find_latest_calib(lookup)];
end

holograms = cell(1, numel(wavelength));
for w = 1:numel(wavelength)
    fprintf('Waiting for holorequest for %dnm...\n', wavelength(w))
    hololist = generate_holograms_new(comm, Setup, calib(w));
    % hololist = generate_holograms2D(comm, Setup, calib(w));
    holograms{w} = uint8(hololist);
end
% holograms contains the phase masks, 1 x n_slms, each is an Nx x Ny x
% patterns, each pattern is a specific hologram

fprintf('All holorequests received.\n')

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
    s.wait_for_trigger = 1; % set settintgs
    s.timeout_ms = timeout;

    s.start();
end

orderBackup=[]; %Sequence list is archived in case the daq errors. normally disposed of after exp. 1/19/21
c=1;
while true
    orderBackup{c} = ShootSequencesMsocket2K(slm, holograms, comm);
    c=c+1;
end

