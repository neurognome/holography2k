function CoC = align_optotune_to_camera(CoC, zsToUse, SIUZ)
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

clear SIdepthData

% zsToUse = linspace(0, 90, 15);% %70 was about 125um on 3/11/21 %Newer optotune has more normal ranges 9/28/29; New Optotune has different range 9/19/19; [0:10:89]; %Scan Image Maxes out at 89

% SIUZ = -15:5:130;% linspace(-120,200,SIpts);
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
    
    devices.comm.send([z, 1], 'si');
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = devices.comm.read(0.01);
        % invar = comm.read(0.01);
    end
    
    
    dataUZ = zeros([sz SIpts]);
    for i = 1:numel(SIUZ)
        fprintf([num2str(round(SIUZ(i))) ' ']);
        
        devices.sutter.moveZ(SIUZ(i));

        if i==1
            pause(3)
        else
            pause(0.3);
        end
        
        %change this part to change number of frames acquired
        frame = devices.bas.grab(framesToAcquire);
        frame = castImg(mean(frame,3));
        frame =  max(frame-BGD,0);
        % frame = imgaussfilt(frame,2);
        dataUZ(:,:,i) =  frame;
    end
    devices.sutter.moveToRef()
    pause(0.1)
    
    % mssend(SISocket,[z 0]);
    disp('Moving on...')
    devices.comm.send([z, 0], 'si');
    invar=[];
    while ~strcmp(invar,'gotit')
        invar = devices.comm.read(0.01);
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
devices.comm.send('end', 'si');

disp(['Scanimage calibration done whole thing took ' num2str(toc(tSI)) 's']);
siT=toc(tSI);

%% First Fits
%Extract data to fit OptotuneZ as a function of camera XYZ

disp('Fitting optotune to Camera... extracting optotune depths')
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

SIThreshHold = 1.5 * stdBgd/sqrt(nBackgroundFrames + framesToAcquire);
excl = SIpeakVal< SIThreshHold;%changed to reflect difference better.\

disp([num2str(numel(SIpeakDepth)) ' points total before exclusions'])
disp([num2str(sum(excl(:))) ' points excluded b/c below threshold'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;

excl = SIpeakVal>255;
excl = SIpeakDepth<-15 | SIpeakDepth>110; %upper bound added 7/15/2020 -Ian

disp([num2str(sum(excl(:))) ' points excluded b/c too deep'])
SIpeakVal(excl)=nan;
SIpeakDepth(excl)=nan;
disp([num2str(sum(~isnan(SIpeakDepth), 'all')) ' points remaining'])


%%CamToOpt
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
    2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

camXYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
camXYZ(3,:) =  SIpeakDepth(:);

camPower = SIpeakVal(:);

optZ = repmat(zsToUse,[nGrids 1]);
optZ = optZ(:);

testSet = randperm(numel(optZ),50);

otherSet = ones([numel(optZ) 1]);
otherSet(testSet)=0;
otherSet = logical(otherSet);

refAsk = (camXYZ(1:3,otherSet))';
refGet = optZ(otherSet);

camToOpto =  polyfitn(refAsk,refGet,modelterms);


Ask = camXYZ(1:3,testSet)';
True = optZ(testSet);

Get = polyvaln(camToOpto,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);


%%fig

f201=figure(201);clf
f201.Units = 'normalized';
f201.Position = [0 0 0.25 0.45];
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],camPower,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Fluorescence intensity by space')
c = colorbar;
c.Label.String = 'Fluorescent Intensity';
axis square

f2 = figure(2);clf
f2.Units = 'normalized';
f2.Position = [0 0.4 0.4 0.4];
subplot(2,2,1)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],optZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,2)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)'),'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,3)
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],polyvaln(camToOpto,camXYZ(1:3,:)')-optZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

subplot(2,2,4)
c = sqrt((polyvaln(camToOpto,camXYZ(1:3,:)')-optZ).^2);
scatter3(camXYZ(1,:),camXYZ(2,:),camXYZ(3,:),[],c,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis \mum')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Optotune Depth';
axis square

%%optZtoCam
cam2XYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
cam2XYZ(3,:) =  optZ(:);
obsZ =  SIpeakDepth(:);

testSet = randperm(numel(obsZ),50);
otherSet = ones([numel(obsZ) 1]);
otherSet(testSet)=0;
otherSet = logical(otherSet);

refAsk = (cam2XYZ(1:3,otherSet))';
refGet = obsZ(otherSet);

OptZToCam =  polyfitn(refAsk,refGet,modelterms);


Ask = cam2XYZ(1:3,testSet)';
True = obsZ(testSet);

Get = polyvaln(OptZToCam,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp(['The mean error in Optotune depth prediction is : ' num2str(meanRMS) 'um']);
disp(['The Max error is: ' num2str(max(RMS)) 'um'])

CoC.OptZToCam= OptZToCam;

out.CoC=CoC;
out.SIfitModelTerms = modelterms;
%%fig
f3 = figure(3);clf
f3.Units = 'normalized';
f3.Position = [0.4 0.4 0.4 0.4];
subplot(2,2,1)
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],obsZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Measured Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,2)
% scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:), [], polyvaln(OptZToCam,cam2XYZ(1:3,:)'),'filled');
vals = polyvaln(OptZToCam,cam2XYZ(1:3,:)');
is_badval = vals < -130;
scatter3(cam2XYZ(1,~is_badval),cam2XYZ(2,~is_badval),cam2XYZ(3,~is_badval), [], vals(~is_badval), 'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Estimated Optotune Level (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,3)
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],polyvaln(OptZToCam,cam2XYZ(1:3,:)')-obsZ,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square

subplot(2,2,4)
c = sqrt((polyvaln(OptZToCam,cam2XYZ(1:3,:)')-obsZ).^2);
scatter3(cam2XYZ(1,:),cam2XYZ(2,:),cam2XYZ(3,:),[],c,'filled');
ylabel('Y Axis pixels')
xlabel('X axis pixels')
zlabel('Z axis Optotune Units')
title('Error RMS (A.U.)')
c = colorbar;
c.Label.String = 'Depth \mum';
axis square


CoC.camToOpto= camToOpto; % this is from the script, and is the important part...
