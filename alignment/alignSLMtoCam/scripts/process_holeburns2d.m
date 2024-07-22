
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

nOpto = numel(zsToBlast);
nBurnHoles = size(XYtarg{1},2);

baseFr = mean(fr(:,:,1:nOpto:end),3);%mean(fr(:,:,1:nOpto:end),3);%Probably more accurate to just do correct zoom, but sometimes having difficulty

k=1;c=0; SIXYZ =[];
for i=2:numel(files)
    t = tic;
    fprintf(['Loading/Processing Frame ' num2str(i)]);
    try
        [dummy fr] = bigread3(fullfile(pth,files(i).name) );
        if c>=nBurnHoles
            k=k+1;
            c=0;
            nBurnHoles = size(XYtarg{k},2);
        end
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
    catch
        fprintf('\nError in Hole analysis... probably loading.')
        x = 0;
        y=0;
    end
    
    
    SIXYZ(:,end+1) = [x,y,zsToBlast(k)];
    disp([' Took ' num2str(toc(t)) ' s']);
end

SIXYZbackup=SIXYZ;
disp(['Done Loading/Processing SI files. Took ' num2str(toc(tLoad)) 's'])
loadT = toc(tLoad);


%% do non-cv SI to cam calculation
burnFitsTimer = tic;

modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;];%...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

modelterms = unique(modelterms(:, 1:2), 'rows');
cam3XYZ = [XYtarg{:};];
SIXYZ = SIXYZbackup;

cam3XYZ=cam3XYZ(:,1:size(SIXYZ,2));

figure(666)
clf
% scatter3(cam3XYZ(1,:), cam3XYZ(2,:), cam3XYZ(3,:), [], '*')
scatter(SIXYZ(1,:), SIXYZ(2,:), [], 'o')

excl = SIXYZ(1,:)<=9 | SIXYZ(1,:)>=503| SIXYZ(2,:)<=9| SIXYZ(2,:)>=503;
disp(['There were ' num2str(sum(excl)) ' of ' num2str(numel(excl)) ' points excluded.'])
cam3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];


refAsk = SIXYZ(1:2,:)';
refGet = (cam3XYZ(1:2,:))';
errScalar = 2.5;

figure(2594)
clf

subplot(1,2,1)
aa = gca();
[SItoCam, trialN] = function_2DCoCIterative(refAsk,refGet,modelterms,errScalar,0, aa);
title('SI to Cam')

subplot(1,2,2)
aa = gca();
[CamToSI, trialN] = function_2DCoCIterative(refGet,refAsk,modelterms,errScalar,0, aa);
title('Cam to SI')

CoC.CamToSI = CamToSI;
CoC.SItoCam = SItoCam;
out.CoC=CoC;

%% alternate calculation
% modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
%     1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
%     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
%     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

tempSLM = cellfun(@(x) x',SLMtarg,'UniformOutput',false);
slm3XYZ = [tempSLM{:}];
SIXYZ = SIXYZbackup;

slm3XYZ=slm3XYZ(1:2,1:size(SIXYZ,2));

excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;

slm3XYZ(:,excl)=[];
SIXYZ(:,excl)=[];

refAsk = SIXYZ(1:2,:)';
refGet = (slm3XYZ(1:2,:))';
errScalar = 2;

figure(2616)
clf

subplot(1,2,1)
aaa = gca();

[SItoSLM, trialN] = function_2DCoCIterative(refAsk,refGet,modelterms,errScalar,0,aaa);
title('SI to SLM')
subplot(1,2,2)
aaa = gca();

[SLMtoSI, trialN] = function_2DCoCIterative(refGet,refAsk,modelterms,errScalar,0,aaa);
title('SLM to SI')

CoC.SItoSLM = SItoSLM;
CoC.SLMtoSI = SLMtoSI;


%% Calculate round trip errors
numTest = 10000;

rangeX = [0 511];%[0 511];
rangeY = [0 511];%[0 511];
rangeZ = [0 55];% Make Sure to match this to the correct range for this optotune;

clear test;
valX = round((rangeX(2)-rangeX(1)).*rand(numTest,1)+rangeX(1));
valY = round((rangeY(2)-rangeY(1)).*rand(numTest,1)+rangeY(1));
valZ = round((rangeZ(2)-rangeZ(1)).*rand(numTest,1)+rangeZ(1));

test = [valX valY];
%%display
test2 = function_SLMtoSI2D(function_SItoSLM2D(test,CoC),CoC);
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


subplot(4,2,7)%aysmetric reverse; foward with 4 chan
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