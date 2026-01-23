function align_slm_to_camera_scope2k
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

calibration_wavelength = 1100; 

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

% setup basler camera
bas = bascam();
bas.start()

% setup sutter
sutter = sutterController();

disp('Ready')
%% look for objective in 1p

bas.preview()
% bas.preview_set_cmax(80)

%% Make mSocketConnections with DAQ and SI Computers
comm = HolochatInterface('holo');
comm.send(calibration_wavelength, 'daq');
%% Put all Manual Steps First so that it can be automated
%% Set Power Levels
pwr = 50; %in mW %0.25 normally
disp(['idividual hologra' ...
    'm; power set to ' num2str(pwr) 'mW']);


%%
disp('Find the spot and check if this is the right amount of power')
% this needs to change depending on which SLM... probably
slmCoords = [0.7 0.7 0.0 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );

blankHolo = zeros([1024 1024]);

slm.feed(Holo);

comm.send([pwr/1000, 1, 1], 'daq')

%bas.preview()
% bas.preview_set_cmax(80)
bas.preview();
comm.send([0, 1, 1], 'daq')
%% check power levels
comm.send([pwr/1000, 1, 1], 'daq')
data = bas.grab(5);
comm.send([0, 1, 1], 'daq')
frame = mean(data,3); % used to have a castImg, not sure why?
figure(668)
imagesc(frame)
colorbar
max(frame,[],'all')
%% Make Sure you're centered
disp('Find Focal Plane, Center and Zero the Sutter')
disp('Leave the focus at Zoom 1. at a power that is less likely to bleach (14% 25mW)') %25% 8/16/19
disp('Don''t forget to use Ultrasound Gel on the objective so it doesn''t evaporate')
comm.send([0, 1, 1], 'daq')

% function_Basler_Preview(Setup, 5);
bas.preview()

input('Turn off Focus and press any key to continue');
sutter.setRef()

comm.send([0, 0], 'si');

disp('Make Sure the DAQ computer is running testMultiTargetsDAQ. and the SI computer running autoCalibSI');
disp('also make those names better someday')
disp('Make sure both lasers are on and the shutters open')
disp('Scanimage should be idle, nearly in plane with focus. and with the gain set high enough to see most of the FOV without saturating')
disp('THE MOUSE MONITOR SHOULD BE TURNED OFF')

sutter.moveZ(100);
disp('testing the sutter double check that it moved to reference +100');
input('Ready to go (Press any key to continue)');

sutter.moveToRef();

%% Create a random set of holograms or use flag to reload
disp('First step Acquire Holograms')
reloadHolos = 0; % CHANGE THIS IF "RECALIFBRATION"
tSingleCompile = tic;
%ranges set by exploration moving holograms looking at z1 fov.
slmXrange = [0.18 0.9];%1030 or 1100
slmYrange = [0.0 0.9];%1030
% slmXrange = [0.15 0.85];%607 or 9001
% slmYrange = [0.15 0.85];%607
% slmXrange = [0.1 0.9];%607
% slmYrange = [0.1 0.9];%607


% set Z range
slmZrange = [-0.05 0.08];%1030
%slmZrange = [-0.10 0.03];%607

if ~reloadHolos
    create_holograms
else
    disp('Reloading old Holograms...')
    try
        load('tempHololist4_will.mat', 'out')
    catch
        [f, p] =uigetfile;
        load(fullfile(p,f),'out');
    end
    hololist = out.hololist;
    slmCoords = out.slmCoords;
    npts = size(slmCoords,2);
    figure(273);scatter3(slmCoords(1,:),slmCoords(2,:),slmCoords(3,:),'o')
    title('Pre Loaded Holograms in SLM space')
end

disp(['Done compiling holograms. Took ' num2str(toc(tSingleCompile)) 's']);

out.hololist=[];

%% Collect background frames for signal to noise testing
disp('Collecting Background Frames');
tSI=tic;

nBackgroundFrames = 10;

Bgdframe = bas.grab(nBackgroundFrames);
Bgd = castImg(mean(Bgdframe,3));
BGD=Bgd;
% BGD = mean(Bgdframe,3);
meanBgd = mean(single(Bgdframe(:)));
stdBgd =  std(single(Bgdframe(:)));

threshHold = meanBgd+3*stdBgd;

fprintf(['3\x03c3 above mean threshold ' num2str(threshHold,4) '\n'])

figure(328)
imagesc(Bgd)

%% Scan Image Planes Calibration
disp('Begining SI Depth calibration, we do this first incase spots burn holes with holograms')
clear SIdepthData

zsToUse = linspace(0, 90, 15);% %70 was about 125um on 3/11/21 %Newer optotune has more normal ranges 9/28/29; New Optotune has different range 9/19/19; [0:10:89]; %Scan Image Maxes out at 89

SIUZ = -15:5:130;% linspace(-120,200,SIpts);
SIpts = numel(SIUZ);

%generate xy grid
%this is used bc SI depths are not necessarily parallel to camera depths
%but SI just generates a sheet of illumination
%therefore, we want to check a grid of spots that we will get depth info
%for
sz = size(Bgd);
gridpts = 25;
xs = round(linspace(1,sz(1),gridpts+2));
ys = round(linspace(1,sz(2),gridpts+2));

xs([1 end])=[];
ys([1 end])=[];
range =15;

%frames to average for image (orig 6) %added 7/15/2020 -Ian
framesToAcquire = 6;

clear dimx dimy XYSI
c=0;
for i=1:gridpts
    for k=1:gridpts
        c=c+1;
        dimx(:,c) = xs(i)-range:xs(i)+range;
        dimy(:,c) = ys(k)-range:ys(k)+range;
        
        XYSI(:,c) = [xs(i) ys(k)];
    end
end

disp(['We will collect ' num2str(numel(zsToUse)) ' planes.'])

SIVals = zeros([SIpts c numel(zsToUse)]);

for k =1:numel(zsToUse)
    t=tic;
    z = zsToUse(k);
    fprintf(['Testing plane ' num2str(z) ': ']);
    
    comm.send([z, 1], 'si');
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = comm.read(0.01);
        % invar = comm.read(0.01);
    end
    
    
    dataUZ = zeros([sz SIpts]);
    for i = 1:numel(SIUZ)
        fprintf([num2str(round(SIUZ(i))) ' ']);
        
        sutter.moveZ(SIUZ(i));

        if i==1
            pause(3)
        else
            pause(0.3);
        end
        
        %change this part to change number of frames acquired
        frame = bas.grab(framesToAcquire);
        frame = castImg(mean(frame,3));
        frame =  max(frame-BGD,0);
        % frame = imgaussfilt(frame,2);
        dataUZ(:,:,i) =  frame;
    end
    sutter.moveToRef()
    pause(0.1)
    
    % mssend(SISocket,[z 0]);
    disp('Moving on...')
    comm.send([z, 0], 'si');
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = comm.read(0.01);
        disp(invar)
        % invar = comm.read(0.01);
    end
    
    figure(1212);
    for hbtmpi=1:numel(SIUZ)
        subplot(6,6,hbtmpi); colorbar;
        imagesc(dataUZ(:,:,hbtmpi));
        title(SIUZ(hbtmpi))
    end
    drawnow()
    
    for i =1:c
        SIVals(:,i,k) = squeeze(mean(mean(dataUZ(dimx(:,i),dimy(:,i),:))));
    end
    disp([' Took ' num2str(toc(t)) 's']);
end

% mssend(SISocket,'end');
comm.send('end', 'si');

disp(['Scanimage calibration done whole thing took ' num2str(toc(tSI)) 's']);
siT=toc(tSI);


%%
tFits=tic;
disp('No Longer Putting off the actual analysis until later, Just saving for now')
out.SIVals =SIVals;
out.XYSI =XYSI;
out.zsToUse =zsToUse;
out.SIUZ = SIUZ;

save('TempSIAlign.mat','out')

%% First Fits
%Extract data to fit OptotuneZ as a function of camera XYZ

disp('Fitting optotune to Camera... extracting optotune depths')

out.SIVals =SIVals;
out.XYSI = XYSI;
out.zsToUse = zsToUse;
out.SIUZ = SIUZ;
nGrids = size(SIVals,2);
nOpt = size(zsToUse,2);
fastWay = 1;

clear SIpeakVal SIpeakDepth
fprintf('Extracting point: ')
for i=1:nGrids
    for k=1:nOpt
        if fastWay
            [a, b] = max(SIVals(:,i,k));
            SIpeakVal(i,k)=a;
            SIpeakDepth(i,k) = SIUZ(b);
        else
            try
                ff = fit(SIUZ', SIVals(:,i,k), 'gauss1');
                SIpeakVal(i,k) = ff.a1;
                SIpeakDepth(i,k) = ff.b1;
            catch
                SIpeakVal(i,k) = nan;
                SIpeakDepth(i,k) = nan;
            end
        end
    end
    fprintf([num2str(i) ' '])
    if mod(i,25)==0
        disp(' ')
    end
end

fprintf('\ndone\n')
% b1 = SIpeakVal;
% b2 = SIpeakDepth;

% %%
% SIpeakVal = b1;
% SIpeakDepth = b2;
process_optotune; % script containing all the plotting code...

%% Coarse Data
tstart=tic;%Coarse Data Timer

disp('Begining Coarse Holo spot finding')
coarsePts = 9; %odd number please
coarseUZ = linspace(-25,150, coarsePts);
comm.send([0, 1, 1], 'daq');

invar='flush';
while ~isempty(invar)
    invar = comm.read(0.01);
    % invar = comm.read(0.01);
end

vals = nan(coarsePts,npts);
xyLoc = nan(2,npts);

sz = size(Bgd);
sizeFactor = 1;% changed bc camera size change 7/16/2020 by Ian %will have to manually test that this is scalable by 4
newSize = sz / sizeFactor;

% dataUZ2 = castImg(zeros([newSize  numel(coarseUZ) npts]));
dataUZ2 = zeros([newSize  numel(coarseUZ) npts], castAs);
maxProjections=castImg(zeros([newSize  npts]));


numFramesCoarseHolo = 10; %number of frames to collect here. added 7/16/2020 -Ian

figure(1234);
for i = 1:numel(coarseUZ)
    fprintf(['First Pass Holo, Depth: ' num2str(coarseUZ(i)) '. Holo : '])
    t = tic;
    sutter.moveZ(coarseUZ(i))
    
    if i==1
        pause(2)
    else
        pause(0.2);
    end
    
    for k=1:npts
        fprintf([num2str(k) ' ']);
        
        if mod(k,25)==0
            fprintf('\n')
        end
        
        slm.feed(hololist(:, :, k));
        comm.send([pwr/1000, 1, 1], 'daq')
        % mssend(masterSocket,[pwr/1000 1 1]);
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = comm.read(0.01);
            % invar = comm.read(0.01);
        end
        frame = bas.grab(numFramesCoarseHolo);
        frame = castImg(mean(frame, 3));

        comm.send([0, 1, 1], 'daq');

        invar=[];
        while ~strcmp(invar,'gotit')
            invar = comm.read(0.01);
        end
        frame =  max(frame-Bgd,0);
        frame = imgaussfilt(frame,2);
        frame = imresize(frame,newSize);
        dataUZ2(:,:,i,k) =  frame;
        imagesc(max(dataUZ2(:, :, i, :), [], 4));
        drawnow();
    end
    fprintf(['\nPlane Took ' num2str(toc(t)) ' seconds\n'])
    
end

sutter.moveToRef()
pause(0.1)

%%
disp('Calculating Depths and Vals')
range=ceil(5/sizeFactor);%Should be scaled by sizeFactor, also shrunk -Ian 7/16/2020 %Range for Hologram analysis window Changed to 5 9/16/19 by Ian
for k=1:npts
    dataUZ = dataUZ2(:,:,:,k);
    mxProj = max(dataUZ,[],3);
    [ x,y ] =function_findcenter(mxProj );
    xyLoc(:,k) = [x,y]*sizeFactor;
    
    maxProjections(:,:,k)=mxProj;
    
    dimx = max((x-range),1):min((x+range),size(mxProj,1));
    dimy =  max((y-range),1):min((y+range),size(mxProj,2));
    
    thisStack = squeeze(mean(mean(dataUZ(dimx,dimy,:))));
    vals(:,k) = thisStack;
    depthIndex = find(thisStack == max(thisStack),1);
    
    fprintf(['Spot ' num2str(k) ' centered at depth ' num2str(round(coarseUZ(depthIndex)))...
        'um. Value: ' num2str(round(vals(depthIndex,k))) '\n']);
end

fprintf(['All Done. Total Took ' num2str(toc(tstart)) 's\n']);

%% Second pass, multi-target version

% assign search params
finePts = 13; % odd number please
fineRange = 40;

disp('Begin multi-target z search...')
multi_time = tic;

% flush the socket
comm.flush();

% Generate multi-target holos based off coarse search data
create_multiholograms
% outputs planes, slmMultiCoordsIndiv, multiHolos, targListIndiv

%%
clear peakValue peakDepth peakFWHM
%%
box_range = 20; % 7/15/20 changed from 50 to 20 distance threshold is set to 50, this must be less to avoid trying to fit 2 holos
disp('shootin!')
% for every  plane

for i = 1:planes
    plane_time = tic;
    holos_this_plane = numel(slmMultiCoordsIndiv{i});
    
    disp(['Plane ' num2str(i) ' of ' num2str(planes)])
    
    % for every holo on that plane
    for j = 1:holos_this_plane
        holo_time = tic;
        disp(['Multi-target holo ' num2str(j) ' of ' num2str(holos_this_plane)])
        
        if size(slmMultiCoordsIndiv{i}{j},2) == 0 || size(slmMultiCoordsIndiv{i}{j},2) == 2 % or <3 ??
            continue
        end
        
        multi_pwr = size(slmMultiCoordsIndiv{i}{j},2) * pwr * 2; % remove *2 later
        slm.feed(multiHolos{i}(j, :, :));
        
        target_ref = targListIndiv{i}{j}(1);
        expected_z = xyzLoc(3,target_ref);

        fineUZ = linspace(expected_z-fineRange,expected_z+fineRange,finePts);
        dataUZ = castImg(nan([size(Bgdframe(:,:,1))  finePts]));
        
        fprintf('Depth: ')
        
        % for every sutter z plane
        for k = 1:finePts
            fprintf([num2str(round(fineUZ(k))) ' ']);
            sutter.moveZ(fineUZ(k));
            
            if i==1
                pause(1)
            else
                pause(0.1);
            end
            

            comm.send([multi_pwr/1000, 1, 1], 'daq');
            % requestPower(multi_pwr,masterSocket) % potentially check this but later
            
            % grab a frame, convert to uint8
            frame = bas.grab(10);
            frame = castImg(mean(frame, 3));
            
            % turn off the laser
            % requestPower(0,masterSocket)
            comm.send([0, 1, 1], 'daq');
            
            % subtract the background and filter
            frame =  max(frame-Bgd,0);
            % frame = imgaussfilt(frame,2);
            dataUZ(:,:,k) =  frame;
            figure(5)
            imagesc(frame)
            drawnow
        end
        
        % move sutter back to reference
        sutter.moveToRef()
        pause(0.1)
        
        % OK, now parse the basler data in expected holo spots
        for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)
            
            target_ref = targListIndiv{i}{j}(targ);
            expected_xyz = xyzLoc(:,target_ref);
            [x, y] = size(Bgd);
            
            targX = expected_xyz(1)-box_range:expected_xyz(1)+box_range;
            targY = expected_xyz(2)-box_range:expected_xyz(2)+box_range;
            
            if max(targX)>x
                targX = expected_xyz(1)-box_range:x;
            end
            if max(targY)>y
                targY = expected_xyz(2)-box_range:y;
            end
            if min(targX)<1
                targX = 1:expected_xyz(1)+box_range;
            end
            if min(targY)<1
                targY = 1:expected_xyz(2)+box_range;
            end
            
            
            % catch for small boxes!!!!
            try
                % method 1 - rely on XY from first step
                targ_stack = double(squeeze(max(max(dataUZ(targX,targY,:)))));
                mxProj = max(dataUZ(targX,targY,:),[],3);
                [ holo_x,holo_y ] = function_findcenter(mxProj );
                xyFine{i}{j}(:,targ) = [holo_x,holo_y];
            catch
                targ_stack = nan(finePts,1);
                xyFine{i}{j}(targ) = nan;
            end
            
            %
            
            try
                ff = fit(fineUZ', targ_stack, 'gauss1');
                peakValue{i}{j}(targ) = ff.a1;
                peakDepth{i}{j}(targ) = ff.b1;
                peakFWHM{i}{j}(targ) = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
            catch
                disp(['Error on fit! Holo: ', num2str(j), ' Target: ', num2str(targ)])
                peakValue{i}{j}(targ) = NaN;
                peakDepth{i}{j}(targ) = NaN;
                peakFWHM{i}{j}(targ) = NaN;
            end
            
            
        end
        
        fprintf('\n')
        disp(['Holo ' num2str(j) ' took ' num2str(toc(holo_time)) 's'])
    end
    
    disp(['Plane ' num2str(i) ' took ' num2str(toc(plane_time)) 's'])
end

fprintf(['All Done. Total Took ' num2str(toc(multi_time)) 's\n']);

process_holograms

%% 
nBurnGrid = 8; %number of points in the burn grid, deault 8
xpts = linspace(SImatchRangeX(1),SImatchRangeX(2),nBurnGrid);
ypts = linspace(SImatchRangeY(1),SImatchRangeY(2),nBurnGrid);


XYpts =[];
for i=1:nBurnGrid
    for k=1:nBurnGrid
        XYpts(:,end+1) = [xpts(i) ypts(k)];
    end
end

XYptsInBounds = inpolygon(XYpts(1,:),XYpts(2,:),SIxboundary,SIyboundary);
XYpts = XYpts(:,XYptsInBounds);

%figure out spacing of pts across zs
zsToBlast = zsToUse;% match to SI Calib %linspace(0,90,11);% Changed to account for newer optotune 9/28/20; Changed to account for new optotune range 9/19/19 by Ian 0:10:80; %OptoPlanes to Blast
interXdist = xpts(2)-xpts(1);
interYdist = ypts(2)-ypts(1);

gridSide = ceil(sqrt(numel(zsToBlast)));
xOff = round(interXdist/gridSide);
yOff = round(interYdist/gridSide);

%Turn into a more unique looking pattern
numPts = size(XYpts,2);
FractionOmit = 0.1; %changed down to 10% from 25% bc not really needed. 9/19/19 by Ian
XYpts(:,randperm(numPts,round(numPts*FractionOmit)))=[];
XYpts = reshape(XYpts,[2 numel(XYpts)/2]);

disp([num2str(size(XYpts,2)) ' points per plane selected. ' num2str(size(XYpts,2)*numel(zsToBlast)) ' total'])

intermediateFitsT = toc(tIntermediateFine);

%% Plot Hole Burn Stuff
tCompileBurn = tic;

%figure(6);
%scatter(XYpts(1,:),XYpts(2,:),'o');

figure(4); clf


clear XYtarg SLMtarg meanCamZ
for i = 1:numel(zsToBlast)
    a = mod(i-1,gridSide);
    b = floor((i-1)/gridSide);
    
    XYuse = bsxfun(@plus,XYpts,([xOff*a yOff*b])');
    optoZ = zsToBlast(i);
    
    zOptoPlane = ones([1 size(XYuse,2)])*optoZ;
    
    Ask = [XYuse; zOptoPlane];
    estCamZ = polyvaln(OptZToCam,Ask');
    meanCamZ(i) = nanmean(estCamZ); %for use by sutter
    Ask = [XYuse; estCamZ'];
    estSLM = function_Eval3DCoC(camToSLM,Ask');
    estPower = polyvaln(SLMtoPower,estSLM);
    
    % negative DE restrictions
    ExcludeBurns = ((estPower<=0) | (estSLM(:,1)<slmXrange(1)) | (estSLM(:,1)>slmXrange(2)) | (estSLM(:,2)<slmYrange(1)) | (estSLM(:,2)>slmYrange(2))); %don't shoot if you don't have the power
    estSLM(ExcludeBurns,:)=[];
    estPower(ExcludeBurns)=[];
    XYuse(:,ExcludeBurns)=[];
    zOptoPlane(ExcludeBurns)=[];
    estCamZ(ExcludeBurns)=[];
    
    XYtarg{i} = [XYuse; zOptoPlane];
    SLMtarg{i} = [estSLM estPower];
    
    subplot(1,2,1)
    scatter3(XYuse(1,:),XYuse(2,:),estCamZ,[],estPower,'filled')
    
    hold on
    subplot (1,2,2)
    scatter3(estSLM(:,1),estSLM(:,2),estSLM(:,3),[],estPower,'filled')
    
    disp([num2str(min(estSLM(:,3))) ' ' num2str(max(estSLM(:,3)))])
    hold on
end

subplot(1,2,1)
title('Targets in Camera Space')
zlabel('Depth \mum')
xlabel('X pixels')
ylabel('Y pixels')

subplot(1,2,2)
title('Targets in SLM space')
xlabel('X SLM')
ylabel('Y SLM')
zlabel('Z SLM')
c = colorbar;
c.Label.String = 'Estimated Power';

disp('Compiling Holos To Burn')

blankHolo = zeros(1024, 1024);
clear holos Diffraction
for k = 1:numel(zsToBlast)
    clear tempHololist
    for i=1:size(XYtarg{k},2)
        t=tic;
        fprintf(['Compiling Holo ' num2str(i) ' for depth ' num2str(k)]);
        subcoordinates =  [SLMtarg{k}(i,1:3) 1];
        %check to avoid out of range holos; added 11/1/19
        %12/5/19 now it allows negative zs
        if ~any(subcoordinates(1:2)>1 | subcoordinates(1:2) <0)
            DE(i) = SLMtarg{k}(i,4);
            [ Hologram,~,~ ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates );
        else
            DE(i)= 0;
            Hologram = blankHolo;
        end
        tempHololist(:,:,i)=Hologram;
        fprintf([' took ' num2str(toc(t)) 's\n']);
    end
    holos{k}=tempHololist;
    Diffraction{k}=DE;
end
disp(['Compiling Done took ' num2str(toc(tCompileBurn)) 's']);

%%
disp('Blasting Holes for SI to SLM alignment, this will take about an hour and take 25Gb of space')
tBurn = tic;

%confirm that SI computer in eval Mode
% mssend(SISocket,'1+2');
comm.send('1+2', 'si');
invar=[];
while ~strcmp(num2str(invar),'3') %stupid cludge so that [] read as false
    invar = comm.read(0.01);
end
disp('linked')

%setup acquisition

numVol = 5; %number of SI volumes to average
baseName = '''calib''';

% needed to turn it into a column for new scanimage?...
str = [];
for ii = 1:length(zsToBlast)
    str = strcat(str, num2str(zsToBlast(ii)), ';');
end

comm.send(['hSI.hStackManager.arbitraryZs = [' str '];'], 'si');
comm.send(['hSI.hStackManager.numVolumes = [' num2str(numVol) '];'], 'si');

comm.send('hSI.hStackManager.enable = 1 ;', 'si');

%comm.send('hSI.hBeams.pzAdjust = 0;', 'si');
comm.send('hSI.hBeams.powers = 10;', 'si'); %power on SI laser. important no to use too much don't want to bleach

comm.send('hSI.extTrigEnable = 0;', 'si'); %saassvign
comm.send('hSI.hChannels.loggingEnable = 1;', 'si'); %savign
comm.send('hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';', 'si');
comm.send(['hSI.hScan2D.logFileStem = ' baseName ';'], 'si');
comm.send('hSI.hScan2D.logFileCounter = 1;', 'si');

comm.send(['hSICtl.updateView;'], 'si');

%clear invar
invar = comm.read(0.01);
while ~isempty(invar)
    invar = comm.read(0.01);
end

comm.send('30+7', 'si');
invar=[];
while ~strcmp(num2str(invar),'37')
    invar = comm.read(0.01);
end
disp('completed parameter set')

%%Burn

%AcquireBaseline
disp('Acquire Baseline')

comm.send('hSI.startGrab()', 'si');
invar = comm.read(0.01);
while ~isempty(invar)
    invar = comm.read(0.01);
end
wait = 1;
while wait
    comm.send('hSI.acqState;', 'si');
    invar = comm.read(0.01);
    while isempty(invar)
        invar = comm.read(0.01);
    end
    
    if strcmp(invar,'idle')
        wait=0;
        disp(['Ready for Next'])
    else
        %             disp(invar)
    end
end

burnPowerMultiplier = 10;% 20 for 1030?%10;%10; % back to 10 bc better DE 12/29/22, WH %5; 10;%change to 10 3/11/21 %previously 5; added by Ian 9/20/19
burnTime = 10; %in seconds, very rough and not precise

disp('Now Burning')
fprintf('Expected number of images: %d\n', sum(cellfun(@(x) size(x, 3), holos)));

ct = 1;
for k=1:numel(zsToBlast)%1:numel(zsToBlast)
    
    offset = round(meanCamZ(k));
    sutter.moveZ(offset)
    if k==1
        pause(5)
    else
        pause(2);
    end
    
    tempHololist=holos{k};
    
    for i=1:size(XYtarg{k},2)%1:size(XYuse,2)
        t=tic;
        fprintf(['Blasting Hole ' num2str(i) '. Depth ' num2str(zsToBlast(k))]);
        slm.feed(tempHololist(:, :, i));
        
        DE = Diffraction{k}(i);
        if DE<.05 %if Diffraction efficiency too low just don't even burn %Ian 9/20/19
            DE=inf;
        end
        blastPower = pwr*burnPowerMultiplier /1000 /DE;
        
        if blastPower>2 %cap for errors, now using a high divided mode so might be high
            blastPower =2;
        end

        % blastPower = min(blastPower, 0.024); %  limit 0.024 for the 607
        
        stimT=tic;
        comm.send([blastPower 1 1], 'daq');
        while toc(stimT)<burnTime
        end
        comm.send([0, 1, 1], 'daq');
        
        %flush masterSocket %flush and handshake added 9/20/19 by Ian
        invar='flush';
        while ~isempty(invar)
            invar = comm.read(0.01);
        end
        %re send 0
        comm.send([0, 1, 1], 'daq');
        %check for handshake
        invar=[];
        while ~strcmp(invar,'gotit')
            invar = comm.read(0.01);
        end
        pause(10)
        comm.send('hSI.startGrab()', 'si');
        % invar = comm.read(0.01);
        % while ~isempty(invar)
        %     invar = comm.read(0.01);
        % end

        pause(5); % force a wait here...
        wait = 1;
        while wait
            comm.send('hSI.acqState', 'si');
            invar = comm.read(0.01);
            while isempty(invar)
                invar = comm.read(0.01);
            end
            
            if strcmp(invar,'idle')
                wait=0;
                %             disp(['Ready for Next'])
            else
                %             disp(invar)
            end
        end

        disp(ct)
        ct =ct + 1;
        disp([' Took ' num2str(toc(t)) 's'])
    end

end
fprintf('Expected number of images: %d\n', sum(cellfun(@(x) size(x, 3), holos)));

sutter.moveToRef();

burnT = toc(tBurn);
disp(['Done Burning. Took ' num2str(burnT) 's']);

disp('Done with Lasers and ScanImage now, you can turn it off')

%% Move file to modulation
process_holeburns

%% Save Output Function
pathToUse = 'C:\Users\holos\Documents\SLM_Management\Calib_Data';
disp('Saving...')
tSave = tic;

save(fullfile(pathToUse,[date '_Calib_' num2str(calibration_wavelength) '.mat']),'CoC')
save(fullfile(pathToUse,'ActiveCalib.mat'),'CoC')

pth = 'C:\Users\holos\Documents\calibs';
save(fullfile(pth,[date '_Calib_' num2str(calibration_wavelength) '.mat']),'CoC')
save(fullfile(pth,'ActiveCalib.mat'),'CoC')
save(fullfile(pth,['CalibWorkspace_' num2str(calibration_wavelength) '.mat']), '-v7.3');

save(fullfile(pathToUse,'CalibWorkspace_v73.mat'), '-v7.3');
% save(fullfile(pathToUse,['CalibWorkspace_will_' date '.mat']))
disp(['Saving took ' num2str(toc(tSave)) 's']);

disp(['All Done, total time from begining was ' num2str(toc(tBegin)) 's. Bye!']);

%% use this to run
%[SLMXYZP] = function_SItoSLM(SIXYZ,CoC);

slm.stop()
% [Setup.SLM ] = Function_Stop_SLM( Setup.SLM );
try; function_close_sutter( Sutter ); end
% try function_stopBasCam(Setup); end




