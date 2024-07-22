% function MsocketHolorequest2K()
% choose wavelengths
clear; clc
wavelength = 1030;%[1100, 900]; % combinations: 900, 1030, 1100, 900/1100, 900/1030

comm = HolochatInterface('holo');

timeout = 1700;

addpath(genpath('C:\Users\holos\Documents\GitHub\holography2k'))
addpath(genpath('C:\Users\holos\Desktop\meadowlark'))

Setup = function_loadparameters3();
Setup.CGHMethod = 2; % now defaults to GSS
Setup.verbose = 0;
% Setup.useGPU = 1; % now defaults to GPU

cycleiterations = 1; % Change this number to repeat the sequence N times instead of just once

%Overwrite delay duration
Setup.TimeToPickSequence = 0.05;    %second window to select sequence ID
Setup.SLM.timeout_ms = timeout;     %No more than 2000 ms until time out
% Setup.calib = 'C:\Users\holos\Documents\calibs\ActiveCalib.mat'; % here we need to somehow feed multiple calibrations?

%% 
calib = [];
for w = wavelength
    switch w
        case 900
            c = importdata('C:\Users\holos\Documents\calibs\07-Nov-2023_Calib_900.mat');
        case 1100
            c = importdata('C:\Users\holos\Documents\calibs\08-Nov-2023_Calib_1100.mat');
        case 1030
            c = importdata ('C:\Users\holos\Documents\calibs\23-Apr-2024_Calib_1030.mat');
            % c = importdata('C:\Users\holos\Documents\calibs\06-Nov-2023_Calib_1030.mat');
    end
   calib = [calib, c]; 
end

sequences = {};
for w = 1:numel(wavelength)
    fprintf('Waiting for holorequest for %dnm...\n', wavelength(w))
    hololist = generate_holograms_new(comm, Setup, calib(w));
    % hololist = generate_holograms2D(comm, Setup, calib(w));
    sequences{end+1} = uint8(hololist);
end

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
    orderBackup{c} = ShootSequencesMsocket2K(slm, sequences, comm);
    c=c+1;
end

