function [CoC, SImatchRangeX, SImatchRangeY] = align_slm_to_camera(CoC, devices, wavelength, hololist, slmCoords)
    % pass in the right devices
    comm = devices.comm;
    slm = devices.slm(wavelength);
    sutter = devices.sutter;

    tstart=tic;%Coarse Data Timer

    disp('Begining Coarse Holo spot finding')
    coarsePts = 9; %odd number please
    coarseUZ = linspace(-25,150, coarsePts);
    comm.send([0, 1, 1], wavelength);

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
            comm.send([pwr/1000, 1, 1], wv)
            % mssend(masterSocket,[pwr/1000 1 1]);
            invar=[];
            while ~strcmp(invar,'gotit')
                invar = comm.read(0.01);
                % invar = comm.read(0.01);
            end
            frame = bas.grab(numFramesCoarseHolo);
            frame = castImg(mean(frame, 3));

            comm.send([0, 1, 1], wv);

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
    targ_time = tic;

    [coarseVal, coarseZidx] =max(vals,[],1);
    zDepthVal = coarseUZ(coarseZidx);
    zdepths = unique(zDepthVal);
    n_planes = numel(zdepths);

    %{
    % changed for 2d...
    coarseInclusionThreshold = 2*stdBgd/sqrt(numFramesCoarseHolo + nBackgroundFrames); %inclusion threshold added based on frames acquired; more stringent then SI. Added 7/16/2020 -Ian
    zDepthVal(coarseVal<coarseInclusionThreshold)=NaN;
    %}

    xyzLoc = [xyLoc;zDepthVal]; %fix this later (?)

    clear slmMultiCoords basCoords targ_list targListIndiv slmMultiCoordsIndiv tempTargList

    for i=1:n_planes % this will be the number of holograms
        % index the z depth
        z = zdepths(i);
        targ_idx = find(xyzLoc(3,:)==z);
        slmMultiCoords{i} = slmCoords(:,targ_idx);
        basCoords{i} = xyzLoc(:,targ_idx);
        targ_list{i} = targ_idx;
    end

    for i=1:n_planes
        % get real and slm coords from coarse
        dist = pdist2(basCoords{i}',basCoords{i}');
        %             dist(find(diag(diag(dist))))=NaN;
        temp =rand(size(dist,1));
        dist(find(diag(diag(temp))))=nan;  %#ok<FNDSB>
        tempTargList = 1:numel(targ_list{i});
        iterCounter =0;
        multiHoloCounter = 0;
        keepGoing=1;
        iterationsBeforeStop =1000;
        distanceThreshold = 30; %changed from 50 on 7/15/20 bc new cam
        size_of_holo = 5;%hayley changed from 20 due to greater z spread of holograms %changed from 25 on 7/15/20 bc new cam
        doThisOnce =0;
        slmMultiCoordsIndiv{i} =[];
        targListIndiv{i}=[];
        
        while keepGoing
            iterCounter=iterCounter+1;
            if numel(tempTargList) <= size_of_holo
                testIdx = tempTargList;
                IdxofTempTargetList = 1:numel(tempTargList);
                keepGoing =0;
            else
                IdxofTempTargetList = randperm(numel(tempTargList),size_of_holo);
                testIdx = tempTargList(IdxofTempTargetList);
            end
            
            %test if good
            subDist = dist(testIdx,testIdx);
            if any(subDist(:)<distanceThreshold)
                good =0;
            else
                good =1;
            end
            
            if good
                multiHoloCounter=multiHoloCounter+1;
                slmMultiCoordsIndiv{i}{multiHoloCounter} = slmMultiCoords{i}(:,testIdx);
                targListIndiv{i}{multiHoloCounter} = targ_list{i}(testIdx) ;
                iterCounter=0;
                tempTargList(IdxofTempTargetList)=[];
            else
                if iterCounter>iterationsBeforeStop && doThisOnce
                    keepGoing=0;
                elseif iterCounter>iterationsBeforeStop
                    size_of_holo=max(round(size_of_holo/2),3);
                    iterCounter=0;
                    doThisOnce=1;
                end
            end
        end
    end

    disp('Setting up stuff for multi-targets...');
    Setup.CGHMethod=2;
    Setup.GSoffset=0;
    Setup.verbose =0;
    Setup.useGPU =1;

    cores=10; % check this

    if cores > 1
        p =gcp('nocreate');
        if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=cores
            delete(p);
            parpool(cores);
        end
    end

    % make the holos
    clear slmShootCoords
    holo_time = tic;
    disp('Compiling holograms...')
    planes = numel(slmMultiCoordsIndiv);

    for i=1:planes
        pt = tic;
        holos_this_plane = numel(slmMultiCoordsIndiv{i});
        
        mtholo_temp=[];
        for k=1:holos_this_plane
            ht = tic;
            [ mtholo, Reconstruction, Masksg ] = function_Make_3D_SHOT_Holos(Setup,slmMultiCoordsIndiv{i}{k}');
            
            mtholo_temp(k,:,:) = mtholo;
        end
        multiHolos{i} = mtholo_temp;
        disp(['Plane ' num2str(i) ' of ' num2str(planes) ' done!  Took ' num2str(toc(pt)) 's'])
    end
    disp(['Done. Took ' num2str(toc(holo_time)) 's'])
    disp(['took ' num2str(toc(targ_time)) 's to compile multi target holos'])
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
                

                comm.send([multi_pwr/1000, 1, 1], wv);
                % requestPower(multi_pwr,masterSocket) % potentially check this but later
                
                % grab a frame, convert to uint8
                frame = bas.grab(10);
                frame = castImg(mean(frame, 3));
                
                % turn off the laser
                % requestPower(0,masterSocket)
                comm.send([0, 1, 1], wv);
                
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

    % Below, we remove any of the "max" ones, because we have dead pixels that
    % always register as "max"
    excludeTrials = all(basXYZ(1:2,:)==[1 1]'); %hayley's understanding: if bas x and y are both one, exclude this trial

    % excludeTrials = excludeTrials | basVal>camMax; %max of this camera is 255

    basDimensions = size(Bgdframe);
    excludeTrials = excludeTrials | basXYZ(1,:)>=basDimensions(1)-1;
    excludeTrials = excludeTrials | basXYZ(2,:)>=basDimensions(2)-1;
    excludeTrials = excludeTrials | basXYZ(3,:)<-20; %9/19/19 Ian Added to remove systematic low fits
    excludeTrials = excludeTrials | basXYZ(3,:)>200;


    excludeTrials = excludeTrials | any(isnan(basXYZ(:,:)));
    excludeTrials = excludeTrials | basVal<10; %  02Nov2023 KS changedi from 10 --> 2
    % excludeTrials = excludeTrials | basVal>(mean(basVal)+3*std(basVal)); %9/13/19 Ian Add
    %excludeTrials = excludeTrials | basVal>250;

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

    holdback = 50;%50;

    refAsk = (slmXYZ(1:3,1:end-holdback))';
    refGet = (basXYZ(1:3,1:end-holdback))';

    %  SLMtoCam = function_3DCoC(refAsk,refGet,modelterms);

    errScalar = 2.5; %2.8;%2.5;
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

    % camToSLM = function_3DCoC(refAsk,refGet,modelterms);

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

    % out.CoC=CoC;
    % out.CoCmodelterms = modelterms;

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
    % out.CoC = CoC;
    % out.powerFitmodelTerms = modelterms;

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
    axis image
