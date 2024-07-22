function dirtyAlignment2D
%% Pathing
in = input('Run Calibration? (y/n)','s');
if ~strcmp(in,'y')
    disp('Did not detect ''y'' so not executing')
    return
end

clear;close all;clc
%%
tBegin = tic;
%
makePaths()
disp('done pathing')

castImg = @uint8;
castAs = 'uint8';
camMax = 255;

%% Setup Stuff
disp('Setting up stuff...');

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU = 0;

Setup.useThorCam =0;
Setup.maxFramesPerAcquire = 3; %set to 0 for unlimited (frames will return will be
Setup.camExposureTime = 10000;

calibration_wavelength = 1030;

if Setup.useGPU
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
end

delete(gcp('nocreate'));
parpool('IdleTimeout', 360);

% setup slms
slm = get_slm(calibration_wavelength);

slm.stop();
slm.wait_for_trigger = 0;
slm.start();

% setup sutter
sutter = sutterController();

disp('Ready')

%% Make mSocketConnections with DAQ and SI Computers
% RUN THIS FIRST, THEN DAQ CODE (alignCodeDAQ2K)

disp('Waiting for msocket communication From DAQ')
%then wait for a handshake
srvsock = mslisten(42130);
masterSocket = msaccept(srvsock,15);
msclose(srvsock);
sendVar = 'A';
mssend(masterSocket, sendVar);

invar = [];

while ~strcmp(invar,'B')
    invar = msrecv(masterSocket,.5);
end
disp('communication from Master To Holo Established');

% send the calibration parameters....
mssend(masterSocket, calibration_wavelength);

%% Connect with SI computer
% RUN THIS FIRST, THEN SI CODE (autoCalibSI)
disp('Waiting for msocket communication to ScanImage Computer')
%then wait for a handshake
srvsock2 = mslisten(42044);
SISocket = msaccept(srvsock2,15);
msclose(srvsock2);
sendVar = 'A';
mssend(SISocket, sendVar);

invar = [];

while ~strcmp(invar,'B');
    invar = msrecv(SISocket,.1);
end

disp('communication from Master To SI Established');

%% Set Power Levels
pwr = 40;
disp(['individual hologram power set to ' num2str(pwr) 'mW']);

%% ensure enough  power to burn slide a lil

disp('Find the spot and check if this is the right amount of power')
% this needs to change depending on which SLM... probably
slmCoords = [.4 .6 0 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords);

blankHolo = zeros([1024 1024]);

slm.feed(Holo);

mssend(masterSocket,[pwr/1000 1 1]); % again, check if this is mW or W

pause()
mssend(masterSocket,[0 1 1]);

%%  generate holograms
%% Create a random set of holograms or use flag to reload
disp('First step Acquire Holograms')
%ranges set by exploration moving holograms looking at z1 fov.
slmXrange = [0.35 0.85];%7/23/21 [.2 .9]; %[0.125 0.8]; %[0.5-RX 0.4+RX]; %you want to match these to the size of your imaging area
slmYrange = [0.2 0.85];%7/23/21 [.05 0.9];%9/19/19 [.01 .7];% [0.075 0.85];%[0.5-RY 0.5+RY];

% set Z range
slmZ = [0]; % no range = 2D

% make a grid...


x = linspace(slmXrange(1)+0.1, slmXrange(2)-0.1, 5);
y = linspace(slmYrange(1)+0.1, slmYrange(2)-0.1, 5);

[xx, yy] = meshgrid(x, y);

coords = [xx(:), yy(:)];
n_holos = size(coords, 1);

hololist = zeros(slm.Nx, slm.Ny, n_holos, 'uint8');
for n = 1:n_holos
    disp(n)
    hololist(:, :, n) = function_Make_3D_SHOT_Holos(Setup, [coords(n, :), slmZ, 1]);
end

%% hol burning

%confirm that SI computer in eval Mode
mssend(SISocket,'1+2');
invar=[];
while ~strcmp(num2str(invar),'3') %stupid cludge so that [] read as false
    invar = msrecv(SISocket,0.01);
end
disp('linked')

baseName = '''calib''';

mssend(SISocket,'hSI.hBeams.powers = 6;'); %power on SI laser. important no to use too much don't want to bleach

mssend(SISocket,['hSI.hStackManager.arbitraryZs = [0; 5; 10];']);
mssend(SISocket,['hSI.hStackManager.numVolumes = [5];']);
mssend(SISocket,'hSI.hStackManager.enable = 1 ;');

mssend(SISocket,'hSI.extTrigEnable = 0;'); %saassvign
mssend(SISocket,'hSI.hChannels.loggingEnable = 1;'); %savign
mssend(SISocket,'hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';');
mssend(SISocket,['hSI.hScan2D.logFileStem = ' baseName ';']);
mssend(SISocket,'hSI.hScan2D.logFileCounter = 1;');

mssend(SISocket, ['hSICtl.updateView;']);

%clear invar
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end

mssend(SISocket,'30+7');
invar=[];
while ~strcmp(num2str(invar),'37')
    invar = msrecv(SISocket,0.01);
end
disp('completed parameter set')

%AcquireBaseline
disp('Acquire Baseline')

mssend(SISocket,'hSI.startGrab()');
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end
wait = 1;
while wait
    mssend(SISocket,'hSI.acqState;');
    invar = msrecv(SISocket,0.01);
    while isempty(invar)
        invar = msrecv(SISocket,0.01);
    end

    if strcmp(invar,'idle')
        wait=0;
        disp(['Ready for Next'])
    else
        %             disp(invar)
    end
end

for i=1:size(hololist,3)%1:size(XYuse,2)
    t=tic;
    fprintf(['Blasting Hole ' num2str(i) '\n']);
    slm.feed(hololist(:, :, i));

  
    blastPower = pwr /1000;

    if blastPower>2 %cap for errors, now using a high divided mode so might be high
        blastPower =2;
    end

    stimT=tic;
    mssend(masterSocket,[blastPower 1 1]);
    while toc(stimT)<0.5
    end
    mssend(masterSocket,[0 1 1]);

    %flush masterSocket %flush and handshake added 9/20/19 by Ian
    invar='flush';
    while ~isempty(invar)
        invar = msrecv(masterSocket,0.01);
    end
    %re send 0
    mssend(masterSocket,[0 1 1]);
    %check for handshake
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = msrecv(masterSocket,0.01);
    end

    mssend(SISocket,'hSI.startGrab()');
    invar = msrecv(SISocket,0.01);
    while ~isempty(invar)
        invar = msrecv(SISocket,0.01);
    end

    wait = 1;
    while wait
        mssend(SISocket,'hSI.acqState');
        invar = msrecv(SISocket,0.01);
        while isempty(invar)
            invar = msrecv(SISocket,0.01);
        end

        if strcmp(invar,'idle')
            wait=0;
            %             disp(['Ready for Next'])
        else
            %             disp(invar)
        end
    end
    disp([' Took ' num2str(toc(t)) 's'])
end


disp('Done with Lasers and ScanImage now, you can turn it off')

%%
disp('Moving files')
tMov = tic;


%on ScanImage Computer
% destination = '''F:\frankenshare\FrankenscopeCalib''' ;
destination = '''K:\Calib\Temp''';
source = '''D:\Calib\Temp\calib*''';

%clear invar
invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = msrecv(SISocket,0.01);
end


mssend(SISocket,['movefile(' source ',' destination ')']);
invar = msrecv(SISocket,0.01);
while isempty(invar)
    invar = msrecv(SISocket,0.01);
end
disp(['Moved. Took ' num2str(toc(tMov)) 's']);
MovT= toc(tMov);


%% read/compute frame

mssend(SISocket,'end');

tLoad = tic;
pth = 'K:\Calib\Temp';
files = dir(sprintf('%s\\*.tif', pth));

baseN = eval(baseName);

[dummy fr] = bigread3(fullfile(pth,files(1).name) );


baseFr = mean(fr,3);

k=1;c=0; SIXYZ =[];
for i=2:numel(files)
    t = tic;
    fprintf(['Loading/Processing Frame ' num2str(i)]);
        [dummy fr] = bigread3(fullfile(pth,files(i).name) );
        % if c>=nBurnHoles
        %     k=k+1;
        %     c=0;
        %     nBurnHoles = size(XYtarg{k},2);
        % end
        c=c+1;

        Frame = mean(fr(:,:,:),3);%mean(fr(:,:,k:nOpto:end),3); %Probably more accurate to just do correct zoom, but sometimes having difficulty
        Frames{k}(:,:,c) = Frame;

        if c>1
            baseFrame = Frames{k}(:,:,c-1);

            %try to exclude those very bright spots
            maskFR = imgaussfilt(Frame,3) - imgaussfilt(Frame,16);
            mask = maskFR > mean(maskFR(:))+6*std(maskFR(:));

            %remove the low frequency slide illumination differences
            filtNum = 3;
            frameFilt = imgaussfilt(Frame,filtNum);
            baseFilt = imgaussfilt(baseFrame,filtNum);


            toCalc = (baseFrame-baseFilt) - (Frame-frameFilt);
            toCalc(mask)=0;

            %             testFr = Frames{k}(:,:,c-1) - Frame;
            [ x,y ] =function_findcenter(toCalc);

            figure(333)
            clf
            subplot(1,3,1)
            imagesc(Frame)

            subplot(1,3,2)
            imagesc(frameFilt)

            subplot(1,3,3)
            imagesc(toCalc)
            hold on
            scatter(y,x,[],'r')
            %             pause
        else
            x = 0;
            y=0;
        end


    SIXYZ(:,end+1) = [x,y];
    disp([' Took ' num2str(toc(t)) ' s']);
end


disp(['Done Loading/Processing SI files. Took ' num2str(toc(tLoad)) 's'])
loadT = toc(tLoad);
%% newfits
SIXYZbackup=SIXYZ;
coordsbackup = coords;

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;
coords = coords(~excl, :);
SIXYZ = SIXYZ(:, ~excl);

% fit tform
CoC = fitgeotform2d(coords, SIXYZ', 'polynomial', 4);
% transformPointsInverse(f, SIXYZ')


% 
% %%
% modelterms = [0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%     1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];
% modelterms = unique(modelterms(:, 1:2), 'rows');
% 
% 
% SItoSLM = function_2DCoCIterative(coords', SIXYZ', modelterms, 2.5, 0);
% 
% CoC.SItoSLM = SItoSLM;
pth = 'C:\Users\holos\Documents\calibs';
save([pth '\dirtycalib_240327.mat'],  'CoC', '-v7.3')
