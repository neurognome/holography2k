
%% THIS BEGINS THE NEW SECTION 3/15/21

% we will use the model to generate a lot of points that we have a pretty
% good guess about where they are located, so we can optimize and multiplex
% the hologram search, using these CoC lets create holograms that shoot a
% pattern into a field of view

%% Simulate and create new Fine Points
% do a CoC to get more points to shoot

denseFineTimer = tic;

nSimulatedTargs = 10000;
multiholosize = 20;
planes = 7;
holosperplane = 10;

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
ntotalPoints = multiholosize * planes * holosperplane;
disp(['Using ' num2str(ntotalPoints) ' points in round 2.'])

slm_coords = {};
bas_coords = {};

% group
[idx, ~] = kmeans(expectBas(:,3), planes);

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
            holo_idxs = randperm(length(targs_this_plane),multiholosize);
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
        thisCoord(:,3) = round(thisCoord(:,3),4); % Added 3/15/21 by Ian for faster compute times

        [ mtholo, ~, ~ ] = function_Make_3D_SHOT_Holos(Setup,thisCoord);
        holos2shoot{i}{j} = mtholo;
        disp(['done in ' num2str(toc(ht)) 's.'])
    end
end

disp('now to shooting...')

%% now repeat multi-target search with new holograms
clear peakValue4 peakDepth4 peakFWHM4 dataUZ4

%% capture holograms

box_range = 20; % distance threshold is set to 100
nframesCapture = 5;

disp('shooting!')

planes = numel(slm_coords);
for i = 1:planes
    holos_this_plane = numel(slm_coords{i});
    disp(['Plane ' num2str(i) ' of ' num2str(planes)])

    % find the mean z of the holo targets and set a range around it
    % rearanging so it doesn't move unnescessarily
    meanz = mean(cellfun(@(x) mean(x(:,3)),bas_coords{i}));
    fineUZ = linspace(meanz-fineRange, meanz+fineRange, finePts);
    dataUZPlane = nan([size(Bgd(:,:,1)) finePts holos_this_plane]);

    % for every sutter z plane
    fprintf('Depth: ')
    figure(4);clf;

    for k = 1:finePts
        fprintf([num2str(round(fineUZ(k))) ' ']);

        % move the sutter
        sutter.moveZ(fineUZ(k))

        if i==1
            pause(mvLongWait)
        else
            pause(mvShortWait);
        end

        a = floor(sqrt(holos_this_plane));
        b = ceil(holos_this_plane/a);
        ro = min([a b]);
        co = max([a b]);

        for j= 1:holos_this_plane
            % 5X for 1100, 2X for 900
            multi_pwr = 5*(size(slm_coords{i}{j},1) * pwr)/1000;

            slm.feed(holos2shoot{i}{j})

            comm.send([multi_pwr, 1, 1], 'daq');
            frame = bas.grab(nframesCapture);
            comm.send([0, 1, 1], 'daq');


            % postprocess
            frame = castImg(mean(frame, 3));
            frame =  max(frame-Bgd,0);
            %frame = imgaussfilt(frame,2);

            % store into dataUZ(x,y,z-plane)
            dataUZPlane(:,:,k,j) =  frame;

            figure(4)
            subplot(ro,co,j)
            imagesc(frame)
            colorbar
            %caxis([0 15]);
            title({['Live Data. Depth ' num2str(round(fineUZ(k)))] ; ['Plane: ' num2str(i) '. Set ' num2str(j)]})
            drawnow

        end
    end
    fprintf('\n')

    for j= 1:holos_this_plane
        fineUZ4{i}{j} = fineUZ;
        dataUZ4{i}{j} = uint8(dataUZPlane(:,:,:,j)); %brought out of for loop

        dataUZ = dataUZPlane(:,:,:,j);

        % move sutter back to reference
        sutter.moveToRef()
        pause(mvShortWait)

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
                targ_stack = squeeze(max(max(dataUZ(targX,targY,:))));
                disp(targ_stack);
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
    for j=1:holos_this_plane
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

excludeTrials = all(basXYZ4(1:2,:)==[1 1]'); %if bas x and y are both one, exclude this trial

% excludeTrials = excludeTrials | basVal4>260; %max of this camera is 255

basDimensions = size(Bgd);
excludeTrials = excludeTrials | basXYZ4(1,:)>=basDimensions(1)-1;
excludeTrials = excludeTrials | basXYZ4(2,:)>=basDimensions(2)-1;
excludeTrials = excludeTrials | basXYZ4(3,:)<-25; %9/19/19 Ian Added to remove systematic low fits
excludeTrials = excludeTrials | basXYZ4(3,:)>150;


excludeTrials = excludeTrials | any(isnan(basXYZ4(:,:)));
excludeTrials = excludeTrials | basVal4<1; 
%excludeTrials = excludeTrials | basVal4>75;%hayley 2/5/24 switched 6 to 250 bc that doesnt make sense

slmXYZBackup = slmXYZ4(:,~excludeTrials);
basXYZBackup = basXYZ4(:,~excludeTrials);
basValBackup = basVal4(:,~excludeTrials);
FWHMValBackup = FWHMVal4(~excludeTrials);

disp(['Number of trials excluded: ' num2str(sum(excludeTrials))])

figure(1922)
clf
subplot(1,2,1)
scatter3(basXYZBackup(1,:),basXYZBackup(2,:),basXYZBackup(3,:), 75, basValBackup, 'Filled')
colorbar
colormap default
title({'Second Denser Fine';'basXYZ and basVals (fine)'})

subplot(1,2,2)
hold on
basx = 1:size(basVal4,2);
plot(basx, basVal4,'ko')
plot(basx(~excludeTrials), basVal4(~excludeTrials), 'o')
ylim([-1,260])
legend('Excluded','Included')
title('basVal by trial')
xlabel('time/holo/acq num')
ylabel('pixel intensity')

%% fit SLM to Camera
%use model terms

errScalar = 3;
holdback = 750;
pxPerMu = 1.0;


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


refAsk = (slmXYZ4(1:3,1:end-holdback))';
refGet = (basXYZ4(1:3,1:end-holdback))';

%  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

figure(1977)
clf
subplot(1,2,1)
ax = gca();
[SLMtoCam, ~] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,ax);
title('SLM to Cam v2')

Ask = refAsk;
True = refGet;
Get = function_Eval3DCoC(SLMtoCam,Ask);

figure(103);
clf
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

figure(101);
clf
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
figure(1977)
subplot(1,2,2)
ax = gca();
[camToSLM, ~] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0, ax);
title('Cam to SLM v2')

Ask = (basXYZ4(1:3,end-holdback:end))';
True = (slmXYZ4(1:3,end-holdback:end))';
Get = function_Eval3DCoC(camToSLM,Ask);

figure(101)
subplot(1,3,3)
scatter3(True(:,1),True(:,2),True(:,3),'*','k')
hold on
scatter3(Get(:,1), Get(:,2), Get(:,3),'o','r')

ylabel('Y Axis SLM units')
xlabel('X axis SLM units')
zlabel('Depth units')
legend('Measured targets', 'Estimated Targets');
title('Camera to SLM')

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

slmXYZ4 = slmXYZBackup;
basVal4 = basValBackup;


modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];

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

figure(1111);clf
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

fprintf(['FWHM in the typical useable volume (0 to 100um) is: ' num2str(median(FWHM(depth>0 & depth<60))) 'um\n'])


finalFitsT = toc(denseFitsTimer);