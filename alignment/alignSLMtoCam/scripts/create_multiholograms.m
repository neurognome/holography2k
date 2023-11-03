targ_time = tic;

[coarseVal, coarseZidx] =max(vals,[],1);
zDepthVal = coarseUZ(coarseZidx);
zdepths = unique(zDepthVal);
n_planes = numel(zdepths);

coarseInclusionThreshold = 2*stdBgd/sqrt(numFramesCoarseHolo + nBackgroundFrames); %inclusion threshold added based on frames acquired; more stringent then SI. Added 7/16/2020 -Ian
zDepthVal(coarseVal<coarseInclusionThreshold)=NaN;

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