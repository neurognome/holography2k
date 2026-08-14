function hololist = generate_holograms_new(comm, Setup, CoC, preread)
% preread (optional): a holoRequest already read off the wire (e.g. by the
% config-primed listener's serve loop, which distinguishes a holoRequest struct
% from a firing-order cell). When given, skip the socket read and use it
% directly; the reply to the DAQ (below) is unchanged. Omit for the original
% behavior (block until a holoRequest arrives).

x = 1;
if nargin >= 4 && ~isempty(preread)
    HRin = preread;
else
    % msg/holo carries two different things: a holoRequest (STRUCT) and a
    % per-trial firing order (CELL). This read used to accept whatever turned up
    % and hand it straight to Pattern.from_struct, so a firing order left over
    % from an aborted run -- or posted early by a DAQ that got ahead -- was
    % consumed as if it were a holoRequest, died on holoRequest.patterns, and took
    % the real request down with it (msg is read-once: once swallowed, it is
    % gone). Discard anything that is not a struct and keep waiting.
    %
    % The wait is also bounded now. It used to be `while isempty(HRin)` with no
    % timeout, which meant a lost request parked the listener forever while the
    % DAQ sat out its own 900 s transfer_timeout and neither side said why.
    HR_WAIT = 900;   % s; matches the DAQ's HoloComputer.transfer_timeout
    HRin = [];
    t0 = tic;
    while isempty(HRin)
        HRin = comm.read(0.5);
        if ~isempty(HRin) && ~isstruct(HRin)
            warning('generate_holograms:notAHoloRequest', ...
                ['Discarded a %s on msg/holo while waiting for a holoRequest ' ...
                 '(a firing order\nfrom an aborted run looks like this). Still ' ...
                 'waiting for the real request.'], class(HRin));
            HRin = [];
        end
        if isempty(HRin) && toc(t0) > HR_WAIT
            error('generate_holograms:requestTimeout', ...
                ['No holoRequest arrived on msg/holo after %g s. The DAQ sends one ' ...
                 'per opto\nchannel from initialize_opto; if it errored or was ' ...
                 'stopped, nothing is coming.'], HR_WAIT);
        end
    end
end
disp('new File Detected - running HoloRequest')
holoRequest = HRin;


patterns = arrayfun(@Pattern.from_struct, holoRequest.patterns);
og_dims = size(patterns);
patterns = patterns(:)';
% KCZ modified to allow for scale; should do nothing if scale=1
scale = holoRequest.scale;
for p = patterns
    if scale ~= 1
        % center to 0,0, scale, then decenter
        scale_center = [.5, .5];

        p.targets(:,1)=scale*(p.targets(:,1)-scale_center(1))+holoRequest.xoffset+scale_center(:,1);
        p.targets(:,2)=scale*(p.targets(:,2)-scale_center(2))+holoRequest.yoffset+scale_center(:,2);
    else  % old behavior-
        p.targets(:,1)=p.targets(:,1)+holoRequest.xoffset;
        p.targets(:,2)=p.targets(:,2)+holoRequest.yoffset;
    end
end

%%quickly compute DEs and return them over msocket
% 
% if isfield(holoRequest,'roiWeights')
%     weightsToUse = holoRequest.roiWeights;
%     weightsToUse(isnan(weightsToUse))=1;
%     disp('Weighting Holograms based on roiWeights')
% else
%     weightsToUse = ones([1 LN]);
%     disp('NO weights detected using flat weight')
% end

%% new, 241029
% drop low DE's...
DE_floor = 0.05;
disp("DROPPING LOW DEs")
for p = patterns
    slm_coords = function_SItoSLM(p.targets, CoC);
    p.targets(slm_coords(:, end) < DE_floor, :) = [];
end


%% this is changed
% [AC, DE_list] = computeDEfromList(SICoordinates, holoRequest.rois, weightsToUse);
% ok how
% expose this function
% function [AttenuationCoeffs, DElist] = computeDEfromList(SICoordinates,ROIs,weights)

%load current CoC

% %catch targets that are below minimal diffraction efficiency 
% DEfloor = 0.05;
% 
% [SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';
% AttenuationCoeffs =SLMCoordinates(4,:);
% lowDE = AttenuationCoeffs<DEfloor;
% AttenuationCoeffs(lowDE)=DEfloor;
% disp([num2str(sum(lowDE)) ' Target(s) below Diffraction Efficiency floor (' num2str(DEfloor) ').']);


% if size(weights,1)~=1
%     weights=weights';
% end
% 
% for i=1:numel(ROIs)
%     ROIselection =ROIs{i};
%     myattenuation = AttenuationCoeffs(ROIselection);
%     energy = 1./myattenuation;
%     energy = energy.*weights(ROIselection);
%     energy = energy/sum(energy);
%     DElist(i) = sum(energy.*myattenuation);
% end


% send back the patterns
% control.io.send(DE_list)
% patterns

arrayfun(@(x) x.calculate_DE(CoC), patterns)
max_pattern_sz = max(arrayfun(@(x) size(x.targets, 1), patterns));
ct = 1;
for p = patterns
    p.id = ct;
    p.powerbias = p.powerbias * (length(p.powerbias)/max_pattern_sz);
    ct = ct + 1;
end

% here we should calculate any power biases before sending it back

% Deliberately BEFORE the compile loop below: the DAQ needs this channel's
% diffraction efficiencies to build its StimInfo, and letting it do that while we
% compile is what keeps a prime as fast as it was. What it does NOT mean is "this
% channel is finished" -- that used to be the only signal the DAQ had, which is
% precisely how it got ahead of us. The real per-channel barrier is
% HoloListener.wait_for_go, which runs after this function returns.
comm.send(arrayfun(@struct, reshape(patterns, og_dims)), 'daq'); % reshape patterns to the og format
disp('Sent patterns back to DAQ');

ct = 1;
hololist = zeros(Setup.Nx, Setup.Ny, numel(patterns));
for p = patterns
    % for each pattern, we generate a hologram
    fprintf('compiling hologram %d of %d\n', ct, numel(patterns));
    slm_coords = function_SItoSLM(p.targets, CoC); % this returns the 4th as an "attenuation", we want to inverse this, as in the original code

    slm_coords(:, 4) = 1./slm_coords(:, 4);

    % Apply the per-target power bias to the target amplitudes ALWAYS. This
    % multiply previously lived only inside the zero_order_dump branch below, so a
    % single-ensemble hologram (size == max_pattern_sz) silently ignored powerbias
    % -- making per-target / calibrated power control a no-op. After the 1/DE
    % inversion above, the delivered power at each target is proportional to
    % powerbias, which is exactly what the graded-power calibration relies on.
    % (The old struct path applied roiWeights unconditionally too; see
    % generate_holograms.m.) powerbias(:) forces a column to match slm_coords(:,4)
    % regardless of whether it was stored as a row or column vector.
    slm_coords(:, 4) = slm_coords(:, 4) .* p.powerbias(:);

    % Optionally dump the leftover (unrequested) power into the zero order so the
    % absolute laser power can stay fixed across differently-sized patterns. Clamp
    % the remainder at 0 so an over-unity powerbias sum can't create a negative
    % zero-order amplitude (which would become NaN in the GS phase retrieval).
    if p.zero_order_dump && (size(p.targets, 1) < max_pattern_sz)
        warning('Fixed laser power, dumping into 0 order...')
        slm_coords = cat(1, slm_coords, [0.5, 0.5, 0, max(0, 1 - sum(slm_coords(:, 4)))]);
    end
    disp(slm_coords)
    if holoRequest.spot_radius > 0
        hololist(:, :, ct) = function_Make_3D_SHOT_Holos_disks_KCZ(Setup, slm_coords, holoRequest.spot_radius);
    else  % old behavior
        hololist(:, :, ct) = function_Make_3D_SHOT_Holos(Setup, slm_coords);
    end

    ct = ct + 1;
end



%{
%%Compute SLM Coordinates
DEfloor = 0.05;

[SLMCoordinates] = function_SItoSLM(SICoordinates',CoC)';
% disp('BAD SUBTRACTION')
% SLMCoordinates([2], :) = 1 - SLMCoordinates([2], :);

AttenuationCoeffs =SLMCoordinates(4,:);
lowDE = AttenuationCoeffs<DEfloor;
AttenuationCoeffs(lowDE)=DEfloor;
SLMCoordinates(4,lowDE)=DEfloor;
disp([num2str(sum(lowDE)) ' Target(s) below Diffraction Efficiency floor (' num2str(DEfloor) ').']);

SLMCoordinates(4,:) = 1./SLMCoordinates(4,:);
SLMCoordinates(3,:) = round(SLMCoordinates(3,:),3); %Added 2/24/21 by Ian for faster compute times 

%%x
close
f = figure('units','normalized','outerposition',[0.125 0.5 0.75 0.5]);

subplot(1,3,1)
scatter3(SICoordinates(1,:),SICoordinates(2,:),SICoordinates(3,:),[],SLMCoordinates(4,:),'filled'); 
% colorbar;
xlabel('X, SI coordinates');ylabel('Y, SI coordinates'); zlabel('Z, SI coordinates'); title('Intensity Correction coefficients');
subplot(1,3,2)
scatter3(SLMCoordinates(1,:),SLMCoordinates(2,:),SLMCoordinates(3,:),[],SLMCoordinates(4,:),'filled'); 
% colorbar;
xlabel('X, SLM coordinates');ylabel('Y, SLM coordinates'); zlabel('Z, SLM coordinates'); title('Intensity Correction coefficients');
subplot(1,3,3)
hist(AttenuationCoeffs,20);
ylabel('Count')
xlabel('Single Target Diffraction Efficiency')
title('Histogram of Diffraction Efficiencies')
pause(1);



%%Sort Holograms

% holoRequest.rois;

numTargets = cellfun(@(x) numel(x),holoRequest.rois);

numSolo = sum(numTargets==1);
solos = find(numTargets==1);
numMid  = sum(numTargets>1 & numTargets < 25);
mids = find(numTargets>1 & numTargets < 25);
numLarge = sum(numTargets>=25 & numTargets < 50);
larges = find(numTargets>=25 & numTargets < 50);
numExtraLarge = sum( numTargets >= 50);
extraLarges =find( numTargets >= 50);


%%Compile Holograms
% optomized for our hologram computer that has a certain amount of gpu space
if size(weightsToUse,2)==1
    weightsToUse=weightsToUse';
end

SLMCoordinates(4,:)=SLMCoordinates(4,:).*weightsToUse; %I don't think this will do anything for single target holos but i'm not sure.

hololist=[];
%solos
if numSolo ==0
    disp('No Single Target Holos')
elseif numSolo<40
    disp('Less than 40 Single Holo Targets')
    for i=1:numSolo
        j=solos(i);
        disp(['Now compiling hologram ' int2str(i) ' of ' int2str(numel(solos))])
        ROIselection = holoRequest.rois{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; 
        energy = energy/sum(energy);
        DE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*DE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        if holoRequest.spot_radius > 0
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,subcoordinates', holoRequest.spot_radius );
            holoRequest.spot_radius
        else  % old behavior
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        end
        hololist(:,:,j) = Hologram;
    end
else 
    disp('More than 40 Single Holo Targets')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([solos]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=4
        delete(p);
        parpool(4);
    end
    %7/24 hayley switched to not parfor
    parfor j=1:numSolo
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(solos))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        if holoRequest.spot_radius > 0
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,subcoordinates', holoRequest.spot_radius );
        else  % old behavior
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        end
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numSolo
        j=solos(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end
end 


%Midsize holograms
if numMid ==0
    disp('No Midsized Holos')
else 
    disp('Computing Midsized Holos')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([mids]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=4
        delete(p);
        parpool(4);
    end
    parfor j=1:numMid
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(mids))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        if holoRequest.spot_radius > 0
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,subcoordinates', holoRequest.spot_radius );
            holoRequest.spot_radius
        else  % old behavior
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        end
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numMid
        j=mids(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end
end 


%Large holograms
if numLarge ==0
    disp('No Large Holos')
else 
    disp('Computing Large Holos')
    
    clear tempDE temphololist
    ROIs =  holoRequest.rois([larges]);
    p =gcp('nocreate');
    if isempty(p) || ~isprop(p,'NumWorkers') || p.NumWorkers ~=2;
        delete(p);
        parpool(2);
    end
    parfor j=1:numLarge
        disp(['Now compiling hologram ' int2str(j) ' of ' int2str(numel(larges))])
        ROIselection =ROIs{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        tempDE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*tempDE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        if holoRequest.spot_radius > 0
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,subcoordinates', holoRequest.spot_radius );
            holoRequest.spot_radius
        else  % old behavior
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        end
        temphololist(:,:,j) = Hologram;
    end
    for i = 1:numLarge
        j=larges(i);
        DE(j)=tempDE(i);
        hololist(:,:,j)=temphololist(:,:,i);
    end   
end

%Extra Large holograms
if numExtraLarge ==0
    disp('No Extra Large Holos')
else 
    disp('Computing Extra Large Holos')
        
    for i=1:numExtraLarge
        j=extraLarges(i);
        disp(['Now compiling hologram ' int2str(i) ' of ' int2str(numel(extraLarges))])
        ROIselection = holoRequest.rois{j};
        myattenuation = AttenuationCoeffs(ROIselection);
        energy = 1./myattenuation; energy = energy/sum(energy);
        DE(j) = sum(energy.*myattenuation);
        disp(['Diffraction efficiency of the hologram : ' int2str(100*DE(j)) '%']);
        subcoordinates = SLMCoordinates(:,ROIselection);
        if holoRequest.spot_radius > 0
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,subcoordinates', holoRequest.spot_radius );
            holoRequest.spot_radius
        else  % old behavior
            [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
        end
        hololist(:,:,j) = Hologram;
    end
end 
%}
disp('Done')