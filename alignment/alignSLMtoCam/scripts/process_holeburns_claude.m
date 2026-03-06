disp('Moving files')
tMov = tic;

%on ScanImage Computer
% destination = '''F:\frankenshare\FrankenscopeCalib''' ;
destination = '''K:\Calib\Temp''';
source = '''D:\Calib\Temp\calib*''';

%clear invar
invar = comm.read();
% invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = comm.read();
    % invar = msrecv(SISocket,0.01);
end


comm.send(['movefile(' source ',' destination ')'], 'si');
% mssend(SISocket,);
invar = comm.read();
% invar = msrecv(SISocket,0.01);
while ~isempty(invar)
    invar = comm.read();
    % invar = msrecv(SISocket,0.01);
end
disp(['Moved. Took ' num2str(toc(tMov)) 's']);
MovT= toc(tMov);


%% read/compute frame

% mssend(SISocket,'end');;
comm.send('end', 'si');

tLoad = tic;
pth = 'K:\Calib\Temp3';
files = dir(sprintf('%s\\*.tif', pth));

% can we guess the bad one?
for f = 2:numel(files)
    timediff(f-1) = files(f).datenum - files(f - 1).datenum;
end
%

baseN = eval(baseName);

[dummy fr] = bigread3(fullfile(pth,files(1).name) );

nOpto = numel(zsToBlast);
nBurnHoles = size(XYtarg{1}, 2);

baseFr = mean(fr(:,:,1:nOpto:end),3);

% -------------------------------------------------------------------------
% CONFIGURATION: hole-finding parameters
% -------------------------------------------------------------------------
SEARCH_RADIUS   = 40;   % px radius around expected XY to restrict search
FINE_GAUSS      = 1.5;  % fine-scale Gaussian for small feature preservation
COARSE_GAUSS    = 20;   % coarse-scale Gaussian for background estimation
DOG_SIGMA1      = 1.5;  % DoG small sigma (signal scale)
DOG_SIGMA2      = 5;    % DoG large sigma (surround suppression)
MIN_PEAK_ZSCORE = 2.5;  % minimum z-score of peak to accept detection
% -------------------------------------------------------------------------

k=1; c=0; SIXYZ=[]; Frames=[];

% Track which hole index within each group we are targeting so we never
% re-use the same expected position for consecutive frames.
expectedXY_history = [];   % [x, y] of all accepted detections so far

for i=2:numel(files) % start at 2 because the first frame is the "background"
    t = tic;
    fprintf(['Loading/Processing Frame ' files(i).name]);

    [dummy fr] = bigread3(fullfile(pth,files(i).name) );

    if c >= nBurnHoles
        k = k+1;
        c = 0;
        nBurnHoles = size(XYtarg{k}, 2);
    end
    c = c+1;

    Frame = mean(fr(:,:,:), 3);
    Frames{k}(:,:,c) = Frame;

    if c > 1
        baseFrame = Frames{k}(:,:,c-1);

        % ------------------------------------------------------------------
        % Step 1: Background-subtracted difference (low-freq illumination
        %         variation removed, same as before)
        % ------------------------------------------------------------------
        frameFilt  = imgaussfilt(Frame,     FINE_GAUSS);
        baseFilt   = imgaussfilt(baseFrame, FINE_GAUSS);
        toCalc_raw = baseFilt - frameFilt;   % bright where hole appeared

        % ------------------------------------------------------------------
        % Step 2: Difference-of-Gaussians (DoG) to enhance small holes.
        %         Suppresses broad background drift; amplifies spot-scale
        %         features. This is key for faint, small holes.
        % ------------------------------------------------------------------
        dog_frame = imgaussfilt(toCalc_raw, DOG_SIGMA1) - ...
                    imgaussfilt(toCalc_raw, DOG_SIGMA2);

        % ------------------------------------------------------------------
        % Step 3: Expected position from target list — constrain search so
        %         we never find the same spot twice.
        % ------------------------------------------------------------------
        % Expected XY for this hole (col=x, row=y in ScanImage convention)
        expX = XYtarg{k}(1, c);   % expected SI x
        expY = XYtarg{k}(2, c);   % expected SI y

        [nRows, nCols] = size(dog_frame);

        % Build a spatial mask centred on the expected position
        [CC, RR] = meshgrid(1:nCols, 1:nRows); % check this, is this correct or does this need to be swapped?
        searchMask = sqrt((CC - expY).^2 + (RR - expX).^2) <= SEARCH_RADIUS;

        % Also exclude a small exclusion zone around every PREVIOUSLY
        % detected hole to avoid re-detecting it.
        exclusionMask = false(nRows, nCols);
        EXCL_RADIUS = 12;  % px — tune to your hole spacing
        for prev = 1:size(expectedXY_history, 1)
            px = expectedXY_history(prev, 1);
            py = expectedXY_history(prev, 2);
            exclusionMask = exclusionMask | ...
                (sqrt((CC - py).^2 + (RR - px).^2) <= EXCL_RADIUS);
        end

        % Combined mask: inside search region AND not previously detected
        combinedMask = searchMask & ~exclusionMask;

        % ------------------------------------------------------------------
        % Step 4: Find peak within the masked DoG image.
        %         Use z-score gating: if the peak isn't strong enough,
        %         fall back to the expected position (safer than garbage).
        % ------------------------------------------------------------------
        masked_dog = dog_frame;
        masked_dog(~combinedMask) = NaN;

        validVals = masked_dog(combinedMask);
        peakVal   = max(validVals);
        mu        = mean(validVals, 'omitnan');
        sigma     = std(validVals,  0, 'omitnan');
        zScore    = (peakVal - mu) / (sigma + eps);

        if zScore >= MIN_PEAK_ZSCORE
            % Reliable detection — use centroid of top-N% region for
            % sub-pixel accuracy rather than pure max, which is noisy.
            topThresh = mu + 2.0 * sigma;
            hotMask   = combinedMask & (dog_frame >= topThresh);

            if sum(hotMask(:)) > 0
                weights = dog_frame .* hotMask;
                weights(weights < 0) = 0;
                wSum = sum(weights(:)) + eps;
                x = sum(RR(:) .* weights(:)) / wSum;   % row → SI-x
                y = sum(CC(:) .* weights(:)) / wSum;   % col → SI-y
            else
                % Fallback to simple max within mask
                [x, y] = function_findcenter(masked_dog);
            end
            detectionFlag = 'DETECTED';
        else
            % Weak / ambiguous — use expected position as best guess
            x = expX;
            y = expY;
            detectionFlag = 'FALLBACK(weak signal)';
            fprintf(' [%s, zScore=%.2f]', detectionFlag, zScore);
        end

        % Record this detection to exclude from future frames
        expectedXY_history(end+1, :) = [x, y]; %#ok<AGROW>

        % ------------------------------------------------------------------
        % Diagnostic figure
        % ------------------------------------------------------------------
        figure(333)
        clf
        subplot(2,2,1)
        imagesc(Frame);  title('Raw frame'); axis image off; colorbar

        subplot(2,2,2)
        imagesc(toCalc_raw); title('BG-sub diff'); axis image off; colorbar

        subplot(2,2,3)
        imagesc(dog_frame);  title('DoG enhanced'); axis image off; colorbar
        hold on
        % Show search region boundary
        theta = linspace(0, 2*pi, 200);
        plot(expY + SEARCH_RADIUS*cos(theta), expX + SEARCH_RADIUS*sin(theta), ...
             'y--', 'LineWidth', 1.2)
        scatter(y, x, 80, 'r', 'filled', 'MarkerEdgeColor', 'w')
        scatter(expY, expX, 80, 'y', '+', 'LineWidth', 2)

        subplot(2,2,4)
        % Show masked DoG so the operator can see what the algorithm used
        dispIm = masked_dog;
        dispIm(isnan(dispIm)) = min(dog_frame(:));
        imagesc(dispIm); title(sprintf('Masked (z=%.1f)', zScore)); 
        axis image off; colorbar
        hold on
        scatter(y, x, 80, 'r', 'filled', 'MarkerEdgeColor', 'w')
        drawnow

    else
        % First frame in a group: use expected position directly
        x = XYtarg{k}(1, c);
        y = XYtarg{k}(2, c);
        expectedXY_history(end+1, :) = [x, y]; %#ok<AGROW>
    end

    SIXYZ(:,end+1) = [x; y; zsToBlast(k)];
    disp([' Took ' num2str(toc(t)) ' s']);
end

SIXYZbackup = SIXYZ;
disp(['Done Loading/Processing SI files. Took ' num2str(toc(tLoad)) 's'])
loadT = toc(tLoad);


%% do non-cv SI to cam calculation
burnFitsTimer = tic;

modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];

cam3XYZ = [XYtarg{:};];
SIXYZ = SIXYZbackup;

cam3XYZ=cam3XYZ(:,1:size(SIXYZ,2));


figure(666)
clf
scatter3(SIXYZ(1,:), SIXYZ(2,:), SIXYZ(3,:), [], 'o')

excl = SIXYZ(1,:)<=9 | SIXYZ(1,:)>=503| SIXYZ(2,:)<=9| SIXYZ(2,:)>=503;
disp(['There were ' num2str(sum(excl)) ' of ' num2str(numel(excl)) ' points excluded.'])
cam3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];


refAsk = SIXYZ(1:3,:)';
refGet = (cam3XYZ(1:3,:))';
errScalar =2.2;

figure(2594)
clf

subplot(1,2,1)
aa = gca();
[SItoCam, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,1, aa);
title('SI to Cam')

subplot(1,2,2)
aa = gca();
[CamToSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,1, aa);
title('Cam to SI')

CoC.CamToSI = CamToSI;
CoC.SItoCam = SItoCam;
out.CoC=CoC;

%% alternate calculation
tempSLM = cellfun(@(x) x',SLMtarg,'UniformOutput',false);
slm3XYZ = [tempSLM{:}];
SIXYZ = SIXYZbackup;

slm3XYZ=slm3XYZ(1:3,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;

slm3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

refAsk = SIXYZ(1:3,:)';
refGet = (slm3XYZ(1:3,:))';
errScalar = 2.1;

figure(2616)
clf

subplot(1,2,1)
aaa = gca();

[SItoSLM, trialN] = function_3DCoCIterative(refAsk,refGet,modelterms,errScalar,0,aaa);
title('SI to SLM')
subplot(1,2,2)
aaa = gca();

[SLMtoSI, trialN] = function_3DCoCIterative(refGet,refAsk,modelterms,errScalar,0,aaa);
title('SLM to SI')

CoC.SItoSLM = SItoSLM;
CoC.SLMtoSI = SLMtoSI;


%% Calculate round trip errors
numTest = 10000;

rangeX = [0 511];
rangeY = [0 511];
rangeZ = [0 90];

clear test;
valX = round((rangeX(2)-rangeX(1)).*rand(numTest,1)+rangeX(1));
valY = round((rangeY(2)-rangeY(1)).*rand(numTest,1)+rangeY(1));
valZ = round((rangeZ(2)-rangeZ(1)).*rand(numTest,1)+rangeZ(1));

test = [valX valY valZ];

test2 = function_SLMtoSI(function_SItoSLM(test,CoC),CoC);
ER1xy = test2(:,1:2)-test(:,1:2);
RMSE1xy = sqrt(sum(ER1xy'.^2));

SIpxPerMu = 512/400;

ER1z = test2(:,3)-test(:,3);
RMSE1z = abs(ER1z);

meanE1rxy = mean(RMSE1xy);
meanE1rz = mean(RMSE1z);

figure(12);clf
subplot(4,2,1)
histogram(RMSE1xy/SIpxPerMu,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')
title({'4 Step CoC'; ['Mean RMS err: ' num2str(meanE1rxy) '\mum']})

subplot(4,2,2)
histogram(RMSE1z,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE1rz) ' optotune Units'])

estSLM = function_Eval3DCoC(CoC.SItoSLM,test);
test2 = function_Eval3DCoC(CoC.SLMtoSI,estSLM);
ER2xy = test2(:,1:2)-test(:,1:2);
RMSE2xy = sqrt(sum(ER2xy'.^2));

SIpxPerMu = 512/800;

ER2z = test2(:,3)-test(:,3);
RMSE2z = abs(ER2z);
meanE2rxy = mean(RMSE2xy);
meanE2rz = mean(RMSE2z);

subplot(4,2,3)
histogram(RMSE2xy/SIpxPerMu,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')
title({'1 Step CoC'; ['Mean RMS err: ' num2str(meanE2rxy) '\mum']})

subplot(4,2,4)
histogram(RMSE2z,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE2rz) ' optotune Units'])


estSLM = function_Eval3DCoC(CoC.SItoSLM,test);
estSIasym = function_SLMtoSI(estSLM,CoC);

ERA = estSIasym-test;
RMSErAxy = sqrt(sum(ERA(:,1:2)'.^2));
RMSErAz = abs(ERA(:,3));


subplot(4,2,5)
histogram(RMSErAxy,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')

meanE3rxy = mean(RMSErAxy);
meanE3rz = mean(RMSErAz);
title({'Asymetric CoC; 1S Forward, 4S Reverse'; ['Mean RMS err: ' num2str(meanE3rxy) '\mum']})

subplot(4,2,6)
histogram(RMSErAz,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE3rz) ' optotune Units'])



estSLM2 = function_SItoSLM(test,CoC);
estSLM2 = estSLM2(:,1:3);

estSIasym2 = function_Eval3DCoC(CoC.SLMtoSI,estSLM2);


ERA2 = estSIasym2-test;
RMSErAxy = sqrt(sum(ERA2(:,1:2)'.^2));
RMSErAz = abs(ERA2(:,3));


subplot(4,2,7)
histogram(RMSErAxy,0:0.1:12)
xlim([0 12])
xlabel('XY Error \mum')

meanE3rxy = mean(RMSErAxy);
meanE3rz = mean(RMSErAz);
title({'Asymetric CoC reverse. 4S Forward, 1S Reverse'; ['Mean RMS err: ' num2str(meanE3rxy) '\mum']})

subplot(4,2,8)
histogram(RMSErAz,0:0.1:12)
xlim([0 12])
xlabel('Z Error optoTuneUnits')
title(['Mean RMS err: ' num2str(meanE3rz) ' optotune Units'])

%%Plot scatter
N=10000;


figure(13);clf
subplot(1,2,1)
val=RMSErAxy;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated XY error, both methods')

subplot(1,2,2)
val=RMSErAz;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated Z error, both methods')

figure(600);clf
subplot(1,2,1)
val=RMSE1xy/SIpxPerMu;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated XY error, 1st methods')

subplot(1,2,2)
val=RMSE1z;
scatter3(test(1:N,1),test(1:N,2),test(1:N,3),[],val(1:N),'filled')
xlabel('SI X')
ylabel('SI Y')
zlabel('Opto Depth')
caxis([0 15])
colorbar
title('Simulated Z error, 1st methods')

burnFitsT = toc(burnFitsTimer);
