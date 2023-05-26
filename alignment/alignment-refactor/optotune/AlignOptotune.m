classdef AlignOptotune < handle
    properties
        optotunePlanes = 0:5:75
        sutterPlanes = -25:5:150
        nframesCapture = 5
        gridpts = 25
        
        SISocket
        masterSocket
        sutter
        bas
    end

    methods
        function obj = AlignOptotune(SISocket, masterSocket, sutter, bas)
            obj.SISocket = SISocket;
            obj.masterSocket = masterSocket;
            obj.sutter = sutter;
            obj.bas = bas;
        end

        function bgd = getBackground(obj, nFrames)
            bgd_frames = obj.bas.grab(nFrames);
            bgd = mean(bgd_frames, 3);
        end

        function SIfits = doFits(obj, SIdepthImages)
            range = 15;
            c=0;

            for i=1:obj.gridpts
                for k=1:obj.gridpts
                    c=c+1;
                    dimx(:,c) = xs(i)-range:xs(i)+range;
                    dimy(:,c) = ys(k)-range:ys(k)+range;
                end
            end

            SIVals = zeros([numel(obj.sutterPlanes), c, numel(obj.sutterPlanes)]);

            for k = 1:numel(SIdepthImages)
                for i =1:c
                    SIVals(:,i,k) = squeeze(mean(SIdepthImages{k}(dimx(:,i),dimy(:,i),:), [1, 2]));
                end
            end

            SIfits = extractScanImageFits(SIVals, obj.optotunePlanes, obj.sutterPlanes);
            SIfits = obj.excludeValues(SIfits);
        end

        function XYSI = getXYSI(obj, gridpts)
            if nargin < 2 || isempty(gridpts)
                gridpts = obj.gridpts;
            end
            xs = round(linspace(1,sz(1),gridpts+2));
            ys = round(linspace(1,sz(2),gridpts+2));
            xs([1 end])=[];
            ys([1 end])=[];

            [X,Y] = meshgrid(xs, ys);
            XYSI = [X(:) Y(:)]';
        end


        function getOpto2Cam(obj, SIfits)
            SIpeakVal = SIfits.SIpeakVal;
            SIpeakDepth = SIfits.SIpeakDepth;

            XYSI = obj.getXYSI();
            cam2XYZ(1:2,:) =  repmat(XYSI,[1 nOpt]);
            cam2XYZ(3,:) =  optZ(:);
            camPower = SIpeakVal(:);
            nGrids = size(SIVals, 2);

            optZ = repmat(app.optotunePlanes, [nGrids, 1]);
            optZ = optZ(:);

 
            obsZ =  SIpeakDepth(:);

            testSet = randperm(numel(obsZ),50);
            otherSet = ones([numel(obsZ) 1]);
            otherSet(testSet)=0;
            otherSet = logical(otherSet);

            refAsk = (cam2XYZ(1:3,otherSet))';
            refGet = obsZ(otherSet);

            OptZToCam =  polyfitn(refAsk,refGet,modelterms);
        end

        function getCam2Opto(obj, SIfits)
            SIpeakVal = SIfits.SIpeakVal;
            SIpeakDepth = SIfits.SIpeakDepth;

            XYSI = obj.getXYSI();
            camXYZ(1:2, :) = repmat(XYSI, [1 numel(app.optotunePlanes)]);
            camXYZ(3, :) = SIpeakDepth;

            camPower = SIpeakVal(:);
            nGrids = size(SIVals, 2);

            optZ = repmat(app.optotunePlanes, [nGrids, 1]);
            optZ = optZ(:);

            % split into test and train and fit
            testSet = randperm(numel(optZ),50);

            otherSet = ones([numel(optZ) 1]);
            otherSet(testSet)=0;
            otherSet = logical(otherSet);
            
            refAsk = (camXYZ(1:3,otherSet))';
            refGet = optZ(otherSet);
            
            camToOpto =  polyfitn(refAsk,refGet,modelterms);
            CoC.camToOpto = camToOpto;
        end

        function SIfits = excludeValues(SIfits)
            nBackgroundFrames = 10;
            bgd = obj.getBackground(nBackgroundFrames);

            SIThreshHold = 3*std(single(bgd(:)))/sqrt(nBackgroundFrames + obj.nframesCapture);
            isBelowThreshold = SIfits.SIpeakVal < SIThreshHold;
            isTooDeep = SIfits.SIpeakDepth < -25 | SIpeakDepth > 150;

            SIfits.SIpeakVal(isBelowThreshold | isTooDeep) = nan;
            SIfits.SIpeakDepth(isBelowThreshold | isTooDeep) = nan;

            fprintf('%d points total before exclusions\n', numel(SIpeakDepth));
            fprintf('%d points excluded b/c below threshold\n', sum(isBelowThreshold));
            fprintf('%d points excluded b/c too deep\n' sum(isTooDeep));
            fprintf('%d poinst remaining\n', sum(~isnan(SIpeakDepth)));
        end

        function optoSearch(obj)
            % main loop
            SIdepthArray = zeros([numel(obj.optotunePlanes), numel(obj.sutterPlanes), size(bas.grab(1))]);`
            SIdepthArrayLabels = {'opto plane', 'sutter plane', 'x', 'y'};
            SIdepthImages = cell(numel(obj.optotunePlanes), 1);
            for k = 1:numel(obj.optotunePlanes) % confirm proper size
                z = obj.optotunePlanes(k);
                fprintf('\n')
                fprintf(['Testing plane ' num2str(z) ': ']);
             
                mssend(obj.SISocket,[z 1]);
                invar=[];
                while ~strcmp(invar,'gotit')
                    invar = msrecv(obj.SISocket,0.01);
                end

                SIdepthImages{k} = obj.sutterSearch();

                mssend(obj.SISocket,[z 0]);
                invar=[];
                while ~strcmp(invar,'gotit')
                    invar = msrecv(obj.SISocket,0.01);
                end
            end
            mssend(obj.SISocket,'end');
        end

        function dataUZ = sutterSearch(obj, mvShortWait, mvLongWait)
            if nargin < 2 || isempty(mvShortWait)
                mvShortWait = 0.1;
            end

            if nargin < 3 || isempty(mvLongWait)
                mvLongWait = 1;
            end

            dataUZ = zeros([sz, numel(obj.sutterPlanes)]);
            for i = 1:numel(obj.sutterPlanes)
                fprintf([num2str(round(obj.sutterPlanes(i))) ' ']);
        
                obj.sutter.moveZ(obj.sutterPlanes(i))
        
                if i==1
                    pause(mvLongWait)
                else
                    pause(mvShortWait);
                end
                
                data = obj.bas.grab(obj.nframesCapture);
        
                data = mean(data, 3);
                frame = max(data-bgd, 0);
                frame = imgaussfilt(frame, 2);
        
                dataUZ(:,:,i) =  frame;
            end

            obj.sutter.moveToRef();
            pause(mvLongWait)
        end

    end

end
