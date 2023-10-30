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

CoC.camToOpto= camToOpto;

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