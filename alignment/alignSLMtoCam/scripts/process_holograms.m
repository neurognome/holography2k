
%% reshape matrix of values into vectors
tIntermediateFine = tic;

peakValueList = peakValue(:);
peakDepthList = peakDepth(:);
peakFWHMList = peakFWHM(:);

[mx, mxi] = max(vals);
% current implementation does not threshold, should be included but more
% complicated here
c=0;
clear slmXYZ basXYZ1 basVal1
% for a = 1:npts %prob works but not always 750 pts bc of exlcusings
%     c = c+1;
for i=1:planes
    %     holos_this_plane = numel(slmMultiCoordsIndiv{i});
    
    %changed to account for when it skipped asked holos bc they were too small -Ian 7/16/2020
    if numel(peakDepth)<i
        holos_this_plane = 0;
    else
        holos_this_plane = numel(peakDepth{i});
    end
    
    for j=1:holos_this_plane
        for targ = 1:size(slmMultiCoordsIndiv{i}{j},2)
            c = c+1;
            slmXYZ(:,c) = slmMultiCoordsIndiv{i}{j}(:,targ);
            target_ref = targListIndiv{i}{j}(targ);
            
            % approach 1
            %             basXYZ1_fine(1:2,c) = xyFine{i}{j}(targ);
            basXYZ1(1:2,c) = xyzLoc(1:2,target_ref);
            basXYZ1(3,c) = peakDepth{i}{j}(targ);
            basVal1(c) = peakValue{i}{j}(targ);
            
            
            FWHMval(c) = peakFWHM{i}{j}(targ); %added 7/10/20 =Ian
        end
    end
end


%% Choose your favorite method of getting XYZ coords
approach = 1;
disp(['you chose approach ' num2str(approach)])
switch approach
    case 1
        basXYZ = basXYZ1;
        basVal = basVal1;
    case 2
        basXYZ = basXYZ2;
        basVal = basVal2;
end


%%

slmXYZBackup2 = slmXYZ;
basXYZBackup2 = basXYZ;
basValBackup2 = basVal;
FWHMBackup2   = FWHMval; %added 7/20/2020 -Ian

%% --- BEGIN INTERMEDIATE FITS ---- %%
%% exclude trials

slmXYZ = slmXYZBackup2;
basXYZ = basXYZBackup2;
basVal = basValBackup2;
FWHMVal = FWHMBackup2;%added 7/20/2020 -Ian

excludeTrials = all(basXYZ(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

excludeTrials = excludeTrials | basVal>camMax; %max of this camera is 255

basDimensions = size(Bgdframe);
excludeTrials = excludeTrials | basXYZ(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ(3,:)<-20; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ(3,:)>200;


excludeTrials = excludeTrials | any(isnan(basXYZ(:,:)));
excludeTrials = excludeTrials | basVal<2; %  02Nov2023 KS changedi from 10 --> 2
% excludeTrials = excludeTrials | basVal>(mean(basVal)+3*std(basVal)); %9/13/19 Ian Add
excludeTrials = excludeTrials | basVal>250;

slmXYZBackup = slmXYZ(:,~excludeTrials);
basXYZBackup = basXYZ(:,~excludeTrials);
basValBackup = basVal(:,~excludeTrials);
FWHMValBackup = FWHMVal(~excludeTrials); % added 7/20/2020 -Ian
disp(['num total holos: ' num2str(size(excludeTrials,2))])
disp(['num exlc holos: ' num2str(sum(excludeTrials))])
%%
f41=figure(41);
clf(41)
f41.Units = 'Normalized';
f41.Position = [0.05 0.4 0.5 0.5];
subplot(1,2,1)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 50, 'k', 'Filled', 'MarkerFaceAlpha',0.5)
hold on
scatter3(xyzLoc(1,:), xyzLoc(2,:), xyzLoc(3,:), 70,  'r', 'Filled', 'MarkerFaceAlpha', 0.7)
legend('Fine','Coarse')
title('Detected basXYZs')
subplot(1,2,2)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title('basXYZ and basVals (fine)')

f42=figure(42);
clf(42)
f42.Units = 'Normalized';
f42.Position = [0.05 0 0.40 0.5];
hold on
basx = 1:size(basVal,2);
plot(basx, basVal,'ro')
plot(basx(~excludeTrials), basVal(~excludeTrials), 'o')
legend('Excluded','Included')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')

%% fit SLM to Camera
%use model terms

basXYZ = basXYZBackup;
slmXYZ = slmXYZBackup;
basVal = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];
reOrder = randperm(size(slmXYZ,2));
slmXYZ = slmXYZ(:,reOrder);
basXYZ = basXYZ(:,reOrder);

holdback = 10;%50;

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = (basXYZ(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

errScalar = 2; %2.8;%2.5;
figure(1286);
clf
ax = gca();
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('SLM to Cam v1')

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
% legend('Measured targets', 'Estimated Targets');
title({'Reference Data'; 'SLM to Camera'})

refRMS = sqrt(sum((Get-True).^2,2));
subplot(1,2,2)
scatter3(True(:,1),True(:,2),True(:,3),[],refRMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title({'Reference Data'; 'RMS Error in position'})
caxis([0 30])


Ask = (slmXYZ(1:3,end-holdback:end))';
True = (basXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(101);clf
subplot(1,3,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')


ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
legend('Measured targets', 'Estimated Targets');
title('SLM to Camera')

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp('Error based on Holdback Data...')
disp(['The RMS error: ' num2str(meanRMS) ' pixels for SLM to Camera']);

% pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size
pxPerMu = 0.57723;

disp(['Thats approx ' num2str(meanRMS/pxPerMu) ' um']);

xErr = sqrt(sum((Get(:,1)-True(:,1)).^2,2));
yErr = sqrt(sum((Get(:,2)-True(:,2)).^2,2));
zErr = sqrt(sum((Get(:,3)-True(:,3)).^2,2));

disp('Mean:')
disp(['X: ' num2str(mean(xErr)/pxPerMu) 'um. Y: ' num2str(mean(yErr)/pxPerMu) 'um. Z: ' num2str(mean(zErr)) 'um.']);
disp('Max:')
disp(['X: ' num2str(max(xErr)/pxPerMu) 'um. Y: ' num2str(max(yErr)/pxPerMu) 'um. Z: ' num2str(max(zErr)) 'um.']);

subplot(1,3,2)
scatter3(True(:,1),True(:,2),True(:,3),[],RMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title('RMS Error in position')


refAsk = (basXYZ(1:3,1:end-holdback))';
refGet = (slmXYZ(1:3,1:end-holdback))';

camToSLM = function_3DCoC(refAsk,refGet,modelterms);

figure(1401);
clf
ax = gca();
[camToSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('cam to slm v1')

Ask = (basXYZ(1:3,end-holdback:end))';
True = (slmXYZ(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);


figure(101);
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

% RMS = sqrt(sum((Get-True).^2,2));
% meanRMS = nanmean(RMS);
%
% disp(['The RMS error: ' num2str(meanRMS) ' SLM units for Camera to SLM']);




CoC.camToSLM=camToSLM;
CoC.SLMtoCam = SLMtoCam;

out.CoC=CoC;
out.CoCmodelterms = modelterms;

rtXYZ = function_Eval3DCoC(SLMtoCam,function_Eval3DCoC(camToSLM,basXYZ(1:3,end-holdback:end)'));

err = sqrt(sum((rtXYZ - basXYZ(1:3,end-holdback:end)').^2,2));
meanRTerr = nanmean(err);
disp(['The Mean Round Trip RMS error: ' num2str(meanRTerr) ' pixels (' num2str(meanRTerr/pxPerMu) ' um) camera to SLM to camera']);

%% fit power as a function of SLM
disp('Fitting Power as a function of SLM')
%  modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%      1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%      2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;   ];  %XY spatial calibration model for Power interpolations
slmXYZ = slmXYZBackup;
basVal = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

intVal = basVal;
intVal = sqrt(intVal); %convert fluorescence intensity (2P) to 1P illumination intensity
intVal=intVal./max(intVal(:));

refAsk = (slmXYZ(1:3,1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ(1:3,end-holdback:end))';
True = intVal(end-holdback:end)';

Get = polyvaln(SLMtoPower,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

figure(1);clf
subplot(2,3,1)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
axis square

subplot(2,3,2)
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],basVal,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
axis square

subplot(2,3,6)
normVal = basVal./max(basVal(:));

c = (polyvaln(SLMtoPower,slmXYZ(1:3,:)').^2)-normVal';
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error 2P Power')
colorbar
axis square


disp(['The RMS error: ' num2str(meanRMS) ' A.U. Power Estimate']);
disp(['The Max power error: ' num2str(max(RMS)*100) '% of request']);

CoC.SLMtoPower = SLMtoPower;
out.CoC = CoC;
out.powerFitmodelTerms = modelterms;

%% THIS BEGINS THE NEW SECTION 3/15/21

%% Now using these CoC lets create holograms that shoot a pattern into a field of view
% disp('Picking Holes to Burn')


%module to restrict burn grid to most likely SI FOV.
%added 7/20/19 -Ian

binarySI = ~isnan(SIpeakVal);
SImatchProb = mean(binarySI');%probability that point was detected aka that was in SI range

SImatchXY = camXYZ(1:2,1:625); %location of points in XYZ

figure(8);clf;
s=scatter(SImatchXY(1,:),SImatchXY(2,:),[],SImatchProb,'filled');

SImatchThreshold = 0; % threshold for being in SI FOV (set to 0 to take whole range)

SIx = SImatchXY(1,SImatchProb>SImatchThreshold);
SIy = SImatchXY(2,SImatchProb>SImatchThreshold);

SIboundary = boundary(SIx',SIy');
SIxboundary = SIx(SIboundary);
SIyboundary = SIy(SIboundary);
hold on
p=plot(SIxboundary,SIyboundary);

sz = size(Bgd);
%for fine targets
bufferMargin = 0.1; %fraction of total area as buffer
SImatchRangeXforFine = [max(min(SIx(SIboundary))-sz(1)*bufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*bufferMargin,sz(1))];

SImatchRangeYforFine = [max(min(SIy(SIboundary))-sz(2)*bufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*bufferMargin,sz(2))];
r = rectangle('position',...
    [SImatchRangeXforFine(1) SImatchRangeYforFine(1) SImatchRangeXforFine(2)-SImatchRangeXforFine(1) SImatchRangeYforFine(2)-SImatchRangeYforFine(1)]);
r.EdgeColor='g';

%for hole burning
bufferMargin = 0.05; %fraction of total area as buffer
SImatchRangeX = [max(min(SIx(SIboundary))-sz(1)*bufferMargin,1) ...
    min(max(SIx(SIboundary))+sz(1)*bufferMargin,sz(1))];

SImatchRangeY = [max(min(SIy(SIboundary))-sz(2)*bufferMargin,1) ...
    min(max(SIy(SIboundary))+sz(2)*bufferMargin,sz(2))];

r = rectangle('position',...
    [SImatchRangeX(1) SImatchRangeY(1) SImatchRangeX(2)-SImatchRangeX(1) SImatchRangeY(2)-SImatchRangeY(1)]);
r.EdgeColor='r';

rline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','r');
gline = line(NaN,NaN,'LineWidth',1','LineStyle', '-','color','g');

colorbar
legend('Prob of SI FOV','Detected SI FOV','Burn Boundary Box','Calib Boundary Box')
xlabel('Camera X pixels')
ylabel('Camera Y Pixels')
% axis equal
axis image

nBurnGrid = 8; %number of points in the burn grid
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
%{
%% Simulate and create new Fine POints
% do a CoC to get more points to shoot

denseFineTimer = tic;

nSimulatedTargs = 10000;

% get basler targets to shoot in from basler range
% for X
a = min(basXYZBackup(1,:));
b = max(basXYZBackup(1,:));
% a = min(basXYZ(1,:));
% b = max(basXYZ(1,:));
% a = 100;
% b = 1000;
% a =SImatchRangeXforFine(1);
% b = SImatchRangeXforFine(2);
r = (b-a).*rand(nSimulatedTargs,1) + a;
rX = round(r);

% for Y
a = min(basXYZBackup(2,:));
b = max(basXYZBackup(2,:));
% a = min(basXYZ(2,:));
% b = max(basXYZ(2,:));
% a = 100;
% b = 1000;
r = (b-a).*rand(nSimulatedTargs,1) + a;
rY = round(r);

% for Z
a = min(basXYZBackup(3,:));
b = max(basXYZBackup(3,:));
% a = min(basXYZ(3,:));
% b = max(basXYZ(3,:));
% a = 5;
% b = 80;
r = (b-a).*rand(nSimulatedTargs,1) + a;
rZ = round(r);

bas2shoot = [rX rY rZ];

testSLM = function_Eval3DCoC(camToSLM, bas2shoot);
expectBas = function_Eval3DCoC(SLMtoCam, testSLM);

testSLM(:,4) = ones(size(testSLM,1),1);

% make sure the SLM vals are within range
excludeMe = testSLM(:,1) < 0 | testSLM(:,1) > 1;
excludeMe = excludeMe | testSLM(:,2) < 0 | testSLM(:,2) > 1;

testSLM = testSLM(~excludeMe,:);
expectBas = expectBas(~excludeMe,:);

% generate multi-target holos that are spread apart
multiholosize=20;
planes = 7;
holosperplane = 10;
ntotalPoints = multiholosize * planes * holosperplane;
disp(['Using ' num2str(ntotalPoints) ' points in round 2.'])

slm_coords = {};
bas_coords = {};
%c = 0;

[idx, cent] = kmeans(expectBas(:,3), planes);
% figure
% scatter3(expBas(:,1), expBas(:,2), expBas(:,3), [], categorical(idx))

for i=1:planes
    iter=0;
    h=0;
    while 1
        while h < holosperplane
            iter = iter+1;
            if iter > 10000
                disp(['****BAD WARNING! Exited hologram determination loop early. Could not find a suitable hologram for plane ' num2str(i) '.****'])
                break
            end
            
            targs_this_plane = find(idx==i);
            % choose rand holos
            holo_idxs = randperm(length(targs_this_plane), min(multiholosize, length(targs_this_plane)));
            dist = pdist2(expectBas(holo_idxs,:),expectBas(holo_idxs,:));
            temp = rand(size(dist,1));
            dist(find(diag(diag(temp))))=nan;
            if any(dist<100)
                continue
            end
            h = h+1;
            bas_coords{i}{h} = expectBas(targs_this_plane(holo_idxs),:);
            slm_coords{i}{h} = testSLM(targs_this_plane(holo_idxs),:);
            
            idx(targs_this_plane(holo_idxs)) = -i; %prevent shooting the same target twice. if there are too few in the simulation will error
        end
        break
    end
end



figure(1579)
clf
cmap = colormap(viridis(numel(slm_coords)*holosperplane));
c = 0;
for i = 1:numel(slm_coords)
    hold on
    for j = 1:numel(slm_coords{i})
        c = c + 1;
        hold on
        subplot(1,2,1)
        scatter3(bas_coords{i}{j}(:,1),bas_coords{i}{j}(:,2),bas_coords{i}{j}(:,3), [], cmap(c,:), 'filled')%, 'MarkerFaceAlpha',0.7)
        hold on
        title('Bas Coords')
        subplot(1,2,2)
        scatter3(slm_coords{i}{j}(:,1),slm_coords{i}{j}(:,2),slm_coords{i}{j}(:,3), [], cmap(c,:), 'filled')%, 'MarkerFaceAlpha',0.7)
        title('SLM Coords')
    end
end
%% compute holos
Setup.useGPU = 1;
c = 0;
for i=1:length(slm_coords)
    for j=1:length(slm_coords{i})
        c = c + 1;
        ht = tic;
        disp(['Compiling multi-target hologram ' num2str(c)])
        
        thisCoord = slm_coords{i}{j};
        thisCoord(:,3) = round(thisCoord(:,3),4); %Added 3/15/21 by Ian for faster compute times
        
        [ mtholo, Reconstruction, Masksg ] = function_Make_3D_SHOT_Holos(Setup,thisCoord);
        holos2shoot{i}{j} = mtholo;
        disp(['done in ' num2str(toc(ht)) 's.'])
    end
end

disp('now to shooting...')

%% now repeat multi-target search with new holograms
clear peakValue4 peakDepth4 peakFWHM4 dataUZ4

%%

background = bas.grab(20);
% background = function_BasGetFrame(Setup,20);
range = 6;
box_range = 100; % distance threshold is set to 100
disp('shooting!')

planes = numel(slm_coords);
for i = 1:planes
    holos_this_plane = numel(slm_coords{i});
    disp(['Plane ' num2str(i) ' of ' num2str(planes)])
    
    % find the mean z of the holo targets and set a range around it
    % rearanging so it doesn't move unnescessarily
    meanz = mean(cellfun(@(x) mean(x(:,3)),bas_coords{i}));
    %         minz = min(cellfun(@(x) min(x(:,3)),bas_coords{i}));
    %         maxz = min(cellfun(@(x) max(x(:,3)),bas_coords{i}));
    
    %      meanz = mean(bas_coords{i}{j}(:,3));
    fineUZ = linspace(meanz-fineRange, meanz+fineRange, finePts);
    
    % for every holo on that plane
    %     for j = 1:holos_this_plane
    dataUZPlane = castImg(nan([size(Bgdframe(:,:,1)) finePts holos_this_plane]));
    
    % set power
    %         multi_pwr = size(slm_coords{i}{j},1) * pwr;
    %         Function_Feed_SLM(Setup.SLM, holos2shoot{i}{j});
    %
    
    
    % for every sutter z plane
    fprintf('Depth: ')
    figure(4);clf;
    for k = 1:finePts
        fprintf([num2str(round(fineUZ(k))) ' ']);
        
        sutter.moveZ(-fineUZ(k));
        
        if i==1
            pause(1)
        else
            pause(0.1);
        end
        
        a = floor(sqrt(holos_this_plane));
        b = ceil(holos_this_plane/a);
        ro = min([a b]);
        co = max([a b]);
        
        
        for j= 1:holos_this_plane
            multi_pwr = size(slm_coords{i}{j},1) * pwr*2;
            % Function_Feed_SLM(Setup.SLM, holos2shoot{i}{j});
            slm.feed(holos2shoot{i}{j});
            
            requestPower(multi_pwr,masterSocket)
            
            % grab a frame, convert to uint8
            % frame = function_BasGetFrame(Setup,3);
            frame = bas.grab(3);
            frame = castImg(mean(frame,3));
            
            % turn off the laser
            requestPower(0,masterSocket)
            
            % subtract the background and filter
            frame =  max(frame-Bgd,0);
            frame = imgaussfilt(frame,2);
            % store into dataUZ(x,y,z-plane)
            dataUZPlane(:,:,k,j) =  frame;
            %             dataUZ4{i}{j} = dataUZ;
            
            figure(4);
            subplot(ro,co,j)
            imagesc(frame)
            colorbar
            caxis([0 10]);
            title({['Live Data. Depth ' num2str(round(fineUZ(k)))] ; ['Plane: ' num2str(i) '. Set ' num2str(j)]})
            drawnow
            
        end
    end
    
    for j= 1:holos_this_plane
        fineUZ4{i}{j} = fineUZ;
        dataUZ4{i}{j} = dataUZPlane(:,:,:,j); %brought out of for loop
        
        dataUZ = dataUZPlane(:,:,:,j);
        
        % move sutter back to reference
        % position = Sutter.Reference;
        % moveTime=moveTo(Sutter.obj,position);
        sutter.moveToRef();
        pause(0.1)
        
        % target parsing, might do later instead
        for targ = 1:size(slm_coords{i}{j},1)
            
            expected_xyz = bas_coords{i}{j}(targ,:);
            [x, y] = size(Bgd);
            
            targX = round(expected_xyz(1)-box_range:expected_xyz(1)+box_range);
            targY = round(expected_xyz(2)-box_range:expected_xyz(2)+box_range);
            
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
            
            try
                % method 1 - rely on XY from first step
                targ_stack = double(squeeze(max(max(dataUZ(targX,targY,:)))));
                mxProj = max(dataUZ(targX,targY,:),[],3);
     
                [ holo_x,holo_y ] =function_findcenter(mxProj );
                xyFine4{i}{j}(:,targ) = [holo_x+(min(targX)), holo_y+(min(targY))];
            catch
                targ_stack = nan(finePts,1);
                xyFine4{i}{j}(:,targ) =[nan, nan];
            end
            
            
            try
                ff = fit(fineUZ', targ_stack, 'gauss1');
                peakValue4{i}{j}(targ) = ff.a1;
                peakDepth4{i}{j}(targ) = ff.b1;
                peakFWHM4{i}{j}(targ) = 2*sqrt(2*log(2))*ff.c1/sqrt(2);
            catch
                disp(['Error on fit! Holo: ', num2str(j), ' Target: ', num2str(targ)])
                peakValue4{i}{j}(targ) = NaN;
                peakDepth4{i}{j}(targ) = NaN;
                peakFWHM4{i}{j}(targ) = NaN;
            end
        end
    end
end

fineT = toc(denseFineTimer);
disp(['Dense Fine Fits took ' num2str(fineT) 's']);

%% New Fits with denser Fine
denseFitsTimer = tic;

c=0;
slmXYZextra = [];
baxXYZextra =[];
basValextra=[];
FWHMValExtra = [];
peakDepthValExtra =[];

for i=1:planes
    for j=1:numel(slm_coords{i})
        for targ = 1:size(slm_coords{i}{j},1)
            c=c+1;
            slmXYZextra(c,:) = slm_coords{i}{j}(targ,:);
            baxXYZextra(c,:) = xyFine4{i}{j}(:,targ);
            basValextra(c) = peakValue4{i}{j}(targ);
            FWHMValExtra(c) = peakFWHM4{i}{j}(targ);
            peakDepthValExtra(c) = peakDepth4{i}{j}(targ);
        end
    end
end

%% exclude trials

slmXYZ4 = slmXYZextra';
basXYZ4 = [baxXYZextra peakDepthValExtra']';
basVal4 = basValextra;
FWHMVal4 = FWHMValExtra;%added 7/20/2020 -Ian

excludeTrials = all(basXYZ4(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

% excludeTrials = excludeTrials | basVal4>260; %max of this camera is 255

basDimensions = size(Bgdframe);
excludeTrials = excludeTrials | basXYZ4(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ4(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ4(3,:)<-25; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ4(3,:)>150;


excludeTrials = excludeTrials | any(isnan(basXYZ4(:,:)));
excludeTrials = excludeTrials | basVal4<5; %8/3 hayley add 5; Ian ammend to 1 9/13
% excludeTrials = excludeTrials | basVal4>(mean(basVal4)+2*std(basVal4)); %9/13/19 Ian Add
excludeTrials = excludeTrials | basVal4>245;

slmXYZBackup = slmXYZ4(:,~excludeTrials);
basXYZBackup = basXYZ4(:,~excludeTrials);
basValBackup = basVal4(:,~excludeTrials);
FWHMValBackup = FWHMVal4(~excludeTrials); % added 7/20/2020 -Ian
%basValBackup = basValBackup(:,1:386); % WH add to exlude trials with water loss
%slmXYZBackup = slmXYZBackup(:,1:386); % did this on 1/29/20
%basXYZBackup = basXYZBackup(:,1:386);

sum(excludeTrials)
%%
figure(1922)
clf
subplot(1,2,1)
% scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 50, 'k', 'Filled', 'MarkerFaceAlpha',0.5)
% hold on
% scatter3(xyzLoc(1,:), xyzLoc(2,:), xyzLoc(3,:), 70,  'r', 'Filled', 'MarkerFaceAlpha', 0.7)
% legend('Fine','Coarse')
% title('Detected basXYZs')
% subplot(1,2,2)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title({'Second Denser Fine';'basXYZ and basVals (fine)'})

% f42=figure(42);
% clf(42)
subplot(1,2,2)
hold on
basx = 1:size(basVal4,2);
plot(basx, basVal4,'ro')
plot(basx(~excludeTrials), basVal4(~excludeTrials), 'o')
legend('Excluded','Included')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')





%% fit SLM to Camera
%use model terms

basXYZ4 = basXYZBackup;
slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;

disp('Fitting SLM to Camera')
modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1 ; 1 1 1 ;...
    2 0 0; 0 2 0; 0 0 2;  ...
    2 0 1; 2 1 0; 0 2 1; 1 2 0; 0 1 2;  1 0 2; ... ];  %XY spatial calibration model for Power interpolations
    2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2;];
reOrder = randperm(size(slmXYZ4,2));
slmXYZ4 = slmXYZ4(:,reOrder);
basXYZ4 = basXYZ4(:,reOrder);

holdback = 300;%50

refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = (basXYZ4(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

errScalar = 3; %2.8;%2.5;
figure(1977);clf;subplot(1,2,1)
[SLMtoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('SLM to Cam v2')

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
% legend('Measured targets', 'Estimated Targets');
title({'Reference Data'; 'SLM to Camera'})

refRMS = sqrt(sum((Get-True).^2,2));
subplot(1,2,2)
scatter3(True(:,1),True(:,2),True(:,3),[],refRMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title({'Reference Data'; 'RMS Error in position'})
caxis([0 30])


Ask = (slmXYZ4(1:3,end-holdback:end))';
True = (basXYZ4(1:3,end-holdback:end))';
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(101);clf
subplot(1,2,1)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')


ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
legend('Measured targets', 'Estimated Targets');
title('SLM to Camera')

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);
disp('Error based on Holdback Data...')
disp(['The RMS error: ' num2str(meanRMS) ' pixels for SLM to Camera']);

% pxPerMu = size(frame,1) / 1000; %really rough approximate of imaging size
% pxPerMu = 0.57723;
pxPerMu = 1.3;

disp(['Thats approx ' num2str(meanRMS/pxPerMu) ' um']);

xErr = sqrt(sum((Get(:,1)-True(:,1)).^2,2));
yErr = sqrt(sum((Get(:,2)-True(:,2)).^2,2));
zErr = sqrt(sum((Get(:,3)-True(:,3)).^2,2));

disp('Mean:')
disp(['X: ' num2str(mean(xErr)/pxPerMu) 'um. Y: ' num2str(mean(yErr)/pxPerMu) 'um. Z: ' num2str(mean(zErr)) 'um.']);
disp('Max:')
disp(['X: ' num2str(max(xErr)/pxPerMu) 'um. Y: ' num2str(max(yErr)/pxPerMu) 'um. Z: ' num2str(max(zErr)) 'um.']);

subplot(1,3,2)
scatter3(True(:,1),True(:,2),True(:,3),[],RMS,'filled');
colorbar
ylabel('Y Axis Pixels')
xlabel('X axis Pixels')
zlabel('Depth \mum')
title('RMS Error in position')


refAsk = (basXYZ4(1:3,1:end-holdback))';
refGet = (slmXYZ4(1:3,1:end-holdback))';

%  camToSLM = function_3DCoC(refAsk,refGet,modelterms);

subplot(1,2,2)
[camToSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0);
title('Cam to SLM v2')

Ask = (basXYZ4(1:3,end-holdback:end))';
True = (slmXYZ4(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);


figure(101);
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

% RMS = sqrt(sum((Get-True).^2,2));
% meanRMS = nanmean(RMS);
%
% disp(['The RMS error: ' num2str(meanRMS) ' SLM units for Camera to SLM']);




CoC.camToSLM=camToSLM;
CoC.SLMtoCam = SLMtoCam;

out.CoC=CoC;
out.CoCmodelterms = modelterms;

rtXYZ = function_Eval3DCoC(SLMtoCam,function_Eval3DCoC(camToSLM,basXYZ4(1:3,end-holdback:end)'));

err = sqrt(sum((rtXYZ - basXYZ4(1:3,end-holdback:end)').^2,2));
meanRTerr = nanmean(err);
disp(['The Mean Round Trip RMS error: ' num2str(meanRTerr) ' pixels (' num2str(meanRTerr/pxPerMu) ' um) camera to SLM to camera']);

%% fit power as a function of SLM
disp('Fitting Power as a function of SLM')
%  modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%      1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%      2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;   ];  %XY spatial calibration model for Power interpolations
slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

intVal = basVal4;
intVal = sqrt(intVal); %convert fluorescence intensity (2P) to 1P illumination intensity
intVal=intVal./max(intVal(:));

refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = intVal(1:end-holdback);

SLMtoPower =  polyfitn(refAsk,refGet,modelterms);

Ask = (slmXYZ4(1:3,end-holdback:end))';
True = intVal(end-holdback:end)';

Get = polyvaln(SLMtoPower,Ask);

RMS = sqrt(sum((Get-True).^2,2));
meanRMS = nanmean(RMS);

figure(1);clf
subplot(2,3,1)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],intVal,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Measured Power (converted to 1p)')
colorbar
% axis square

subplot(2,3,2)
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],polyvaln(SLMtoPower,slmXYZ4(1:3,:)'),'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated Power Norm.')
colorbar
% axis square

subplot(2,3,4)
% scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],polyvaln(SLMtoPower,slmXYZ(1:3,:)')-intVal','filled');
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],basVal4,'filled');

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Raw Fluorescence')
colorbar
% axis square

subplot(2,3,5)
c = sqrt((polyvaln(SLMtoPower,slmXYZ4(1:3,:)')-intVal').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error RMS (A.U.)')
colorbar
% axis square

subplot(2,3,3)
c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2);
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Estimated 2P Power')
colorbar
% axis square

subplot(2,3,6)
normVal = basVal4./max(basVal4(:));

c = (polyvaln(SLMtoPower,slmXYZ4(1:3,:)').^2)-normVal';
scatter3(slmXYZ4(1,:),slmXYZ4(2,:),slmXYZ4(3,:),[],c,'filled');
ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Z axis SLM units')
title('Error 2P Power')
colorbar
% axis square


disp(['The RMS error: ' num2str(meanRMS) ' A.U. Power Estimate']);
disp(['The Max power error: ' num2str(max(RMS)*100) '% of request']);

CoC.SLMtoPower = SLMtoPower;
out.CoC = CoC;
out.powerFitmodelTerms = modelterms;

%% Plot FWHM
FWHM = FWHMValBackup;
depth = basXYZBackup(3,:);
slmXYZ = slmXYZBackup;

figure(1001); clf
subplot(1,2,1);
plot(FWHM,depth,'o')
% plot(FWHM,slmCoords(3,:),'o')

ylabel('Axial Depth \mum')
xlabel('FWHM \mum')
ylim([-25 125])
xlim([7.5 50])

refline(0,0)
refline(0,60)


subplot(1,2,2);
scatter3(slmXYZ(1,:),slmXYZ(2,:),slmXYZ(3,:),[],FWHM,'filled')
% scatter3(slmXYZ(1,:),slmXYZ(2,:),depth,[],FWHM,'filled')
caxis([10 50])
h= colorbar;
xlabel('SLM X')
ylabel('SLM Y')
zlabel('SLM Z')
set(get(h,'label'),'string','FWHM \mum')

fprintf(['FWHM in the typical useable volume (0 to 100um) is: ' num2str(mean(FWHM(depth>0 & depth<60))) 'um\n'])


finalFitsT = toc(denseFitsTimer);
%}