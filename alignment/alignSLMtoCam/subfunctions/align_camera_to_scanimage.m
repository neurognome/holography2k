function CoC = align_camera_to_scanimage(CoC, wavelength, devices, SImatchRangeX, SImatchRangeY, zsToUse, slmXrange, slmYrange, slmZrange, scanimage_power, burnPowerMultiplier, burnTime)
    sutter = devices.sutter;
    slm = devices.slm(wavelength);
    comm = devices.comm;

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
        estCamZ = polyvaln(CoC.OptZToCam,Ask');
        meanCamZ(i) = nanmean(estCamZ); %for use by sutter
        Ask = [XYuse; estCamZ'];
        estSLM = function_Eval3DCoC(CoC.camToSLM,Ask');
        estPower = polyvaln(CoC.SLMtoPower,estSLM);
        
        % negative DE restrictions
        %%%% CHECK CHECK CHECK CEHCKC -- how important is this????
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

    
    delay = 0.25;
    comm.send(['hSI.hStackManager.arbitraryZs = [' str '];'], 'si');
    pause(delay)
    comm.send(['hSI.hStackManager.numVolumes = [' num2str(numVol) '];'], 'si');
    pause(delay)

    comm.send('hSI.hStackManager.enable = 1 ;', 'si');
    pause(delay) 
    %comm.send('hSI.hBeams.pzAdjust = 0;', 'si');
    comm.send(sprintf('hSI.hBeams.powers = %d;', scanimage), 'si'); %power on SI laser. important no to use too much don't want to bleach
    pause(delay) 

    comm.send('hSI.extTrigEnable = 0;', 'si'); %saassvign
    pause(delay) 
    comm.send('hSI.hChannels.loggingEnable = 1;', 'si'); %savign
    pause(delay) 
    comm.send('hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';', 'si');
    pause(delay) 
    comm.send(['hSI.hScan2D.logFileStem = ' baseName ';'], 'si');
    pause(delay) 
    comm.send('hSI.hScan2D.logFileCounter = 1;', 'si');
    pause(delay) 

    comm.send(['hSICtl.updateView;'], 'si');
    pause(delay) 

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

    % burnPowerMultiplier = 10;% 20 for 1030?%10;%10; % back to 10 bc better DE 12/29/22, WH %5; 10;%change to 10 3/11/21 %previously 5; added by Ian 9/20/19
    % burnTime = 10; %in seconds, very rough and not precise

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
            comm.send([blastPower 1 1], wv);
            while toc(stimT)<burnTime
            end
            comm.send([0, 1, 1], wv);
            
            %flush masterSocket %flush and handshake added 9/20/19 by Ian
            invar='flush';
            while ~isempty(invar)
                invar = comm.read(0.01);
            end
            %re send 0
            comm.send([0, 1, 1], wv);
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

    %% now process the hole burns
    disp('Moving files')

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

    disp('Check if everything is moved, then press any key to continue...')
    pause


    %% read/compute frame
    % mssend(SISocket,'end');;
    comm.send('end', 'si');

    tLoad = tic;
    pth = 'K:\Calib\Temp';
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

    baseFr = mean(fr(:,:,1:nOpto:end),3);%mean(fr(:,:,1:nOpto:end),3);%Probably more accurate to just do correct zoom, but sometimes having difficulty

    k=1;c=0; SIXYZ =[]; Frames=[];
    for i=2:numel(files) % start at 2 because the first frame is the "background"
        t = tic;
        fprintf(['Loading/Processing Frame ' files(i).name]);
        % try
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
                % maskFR = imgaussfilt(Frame,3) - imgaussfilt(Frame,16);
                % mask = maskFR > mean(maskFR(:))+3*std(maskFR(:));
                % mask = maskFR < mean(maskFR(:))+6*std(maskFR(:));
                
                %remove the low frequency slide illumination differences
                filtNum = 2;
                frameFilt = imgaussfilt(Frame,filtNum);
                baseFilt = imgaussfilt(baseFrame,filtNum);
                
                
                % toCalc = (Frame-frameFilt) - (baseFrame-baseFilt);
            %    toCalc = (frameFilt-Frame) - (baseFilt-baseFrame);
                toCalc = baseFilt - frameFilt;
            % toCalc = -toCalc; % remove this later?
                % toCalc(mask)=0;
                
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
        % catch
        %     fprintf('\nError in Hole analysis... probably loading.')
        %     x = 0;
        %     y=0;
        % end
        
        
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

    cam3XYZ = [XYtarg{:};];
    SIXYZ = SIXYZbackup;

    cam3XYZ=cam3XYZ(:,1:size(SIXYZ,2));


    figure(666)
    clf
    % scatter3(cam3XYZ(1,:), cam3XYZ(2,:), cam3XYZ(3,:), [], '*')
    scatter3(SIXYZ(1,:), SIXYZ(2,:), SIXYZ(3,:), [], 'o')

    excl = SIXYZ(1,:)<=9 | SIXYZ(1,:)>=503| SIXYZ(2,:)<=9| SIXYZ(2,:)>=503;
    disp(['There were ' num2str(sum(excl)) ' of ' num2str(numel(excl)) ' points excluded.'])
    cam3XYZ(:,excl)=[];
    SIXYZ(:,excl)=[];


    refAsk = SIXYZ(1:3,:)'; % detected points from the hole burn
    refGet = (cam3XYZ(1:3,:))'; % expected points?
    errScalar = 3;%2.5

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

    %% alternate calculation
    % modelterms =[0 0 0; 1 0 0; 0 1 0; 0 0 1;...
    %     1 1 0; 1 0 1; 0 1 1; 1 1 1; 2 0 0; 0 2 0; 0 0 2;...
    %     2 0 1; 2 1 0; 0 2 1; 0 1 2; 1 2 0; 1 0 2;...
    %     2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2; ];  %XY spatial calibration model for Power interpolations

    tempSLM = cellfun(@(x) x',SLMtarg,'UniformOutput',false);
    slm3XYZ = [tempSLM{:}];
    SIXYZ = SIXYZbackup;

    slm3XYZ=slm3XYZ(1:3,1:size(SIXYZ,2));

    excl = SIXYZ(1,:)<=5 | SIXYZ(1,:)>=507| SIXYZ(2,:)<=5 | SIXYZ(2,:)>=507;
    %excl = SIXYZ(1,:)<=100 | SIXYZ(1,:)>=300| SIXYZ(2,:)<=100 | SIXYZ(2,:)>=300;

    slm3XYZ(:,excl)=[];
    SIXYZ(:,excl)=[];

    refAsk = SIXYZ(1:3,:)';
    refGet = (slm3XYZ(1:3,:))';
    errScalar = 2.5;

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

    rangeX = [0 511];%[0 511];
    rangeY = [0 511];%[0 511];
    rangeZ = [0 90];% Make Sure to match this to the correct range for this optotune;

    clear test;
    valX = round((rangeX(2)-rangeX(1)).*rand(numTest,1)+rangeX(1));
    valY = round((rangeY(2)-rangeY(1)).*rand(numTest,1)+rangeY(1));
    valZ = round((rangeZ(2)-rangeZ(1)).*rand(numTest,1)+rangeZ(1));

    test = [valX valY valZ];
    %%display
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