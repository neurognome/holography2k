function [SLMXYZP] = function_SItoSLM(SIXYZ,CoC)
%Function takes coordinates in SI space XY and optotuneZ (n x 3 matrix)
%and CoC as created in alignSLMtoCam
%Outputs in SLM coordinates XYZ and power normalization (n x 4 Matrix)
approachMethod=1;

SItoCam = CoC.SItoCam;
OptZToCam = CoC.OptZToCam;
camToSLM = CoC.camToSLM;
SLMtoPower = CoC.SLMtoPower;

estCamXY = function_Eval3DCoC(SItoCam,SIXYZ);
estCamZ = polyvaln(OptZToCam,estCamXY);
estCamXYZ = [estCamXY(:,1:2), estCamZ];
estSLM = function_Eval3DCoC(camToSLM,estCamXYZ);

estPower = polyvaln(SLMtoPower,estSLM);


if approachMethod ==2
%shortcut approach
estSLM = function_Eval3DCoC(CoC.SItoSLM,SIXYZ);
% estSI = function_Eval3DCoC(CoC.SLMtoSI,SLMXYZ);
end

SLMXYZP = [estSLM estPower];


%% 12/19/24 KCZ edited to include power meter-based DE calibration:

if isfield(CoC, 'DE_calib')
    if length(CoC.DE_calib.slmZgrid) ~= 1
        error('DE compensation for now only implemented for 2D SLM coordinates at z=0.')
    else
        if CoC.DE_calib.slmZgrid ~= 0
            error('please redo DE calibration at SLM z=0 plane.')
        else
            x = CoC.DE_calib.slmXgrid;
            y = CoC.DE_calib.slmYgrid;
            P = CoC.DE_calib.powers;
            %P = P.^2;  % do we need to square cuz 2p?

            P = P / max(P(:));  % normalize
            
            center_row = ceil(size(P,1)/2);  % assumes odd
            center_col = ceil(size(P,2)/2);  % assumes odd
            P(center_row-1:center_row+1, center_col-1:center_col+1) = 1;  
            % ^set DE to 1 so that algo doesn't try to increase power near zero order!

            % get DE by interpolating the sampled coordinates:
            DE = interp2(x, y, P, SLMXYZP(:, 2), SLMXYZP(:, 1), 'cubic');  % for some reason, better results when swapping x and y ...
            
            % deal with low DE here; don't rely on generate_holograms:
            % set DE to 1 where low, so no power is redistributed to it
            DE_thresh = .1;
            DE(DE<DE_thresh) = 1;

            SLMXYZP(:, 4) = DE;

            disp(['Diffraction efficiencies: ', num2str(DE')])
        end
        disp('Using DE Calib')
    end
else
    % do nothing, accept the above-calculated SLMXYZP
end