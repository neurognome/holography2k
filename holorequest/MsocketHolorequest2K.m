% function MsocketHolorequest2K()
HoloPrepCode;

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
% calibID = 1;                        % Select the calibration ID (z1=1 but does not exist, Z1.5=2, Z1 sutter =3);
Setup.calib = 'C:\Users\holos\Documents\calibs\ActiveCalib.mat'; % here we need to somehow feed multiple calibrations?

disp('Loading current calibration...')
load(Setup.calib,'CoC');
disp(['Successfully loaded CoC from: ', Setup.calib])
%% now start msocket communication
CoC_900 = importdata(Setup.calib);
CoC_1100 = importdata(Setup.calib); 

hololist_900 = generate_holograms(control, Setup, CoC_900); 
hololist_1100 = generate_holograms(control, Setup, CoC_1100);

%totally remove sequences as a thing that exists basically
sequences = uint8(hololist); %shouldn't change anything added 9/14/21

% flushMSocket(masterSocket)
control.flush();

slm = MeadowlarkOneK();
% slm = HoloeyePLUTO();
slm.stop();
slm.wait_for_trigger = 1; % set settintgs
slm.timeout_ms = timeout;

slm.start();
% 
% [Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
% Setup.SLM.wait_For_Trigger = 1;
% Setup.SLM.timeout_ms = timeout;
% [ Setup.SLM ] = Function_Start_SLM( Setup.SLM );

%sendVar='B';
%mssend(masterSocket,sendVar)
orderBackup=[]; %Sequence list is archived in case the daq errors. normally disposed of after exp. 1/19/21
c=1;
while true
    orderBackup{c} = ShootSequencesMsocket2K(slm, sequences, masterSocket);
    c=c+1;
end

    
