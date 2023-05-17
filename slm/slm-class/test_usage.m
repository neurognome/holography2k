
makePaths()

clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

Setup = function_loadparameters2k();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU =0;

if Setup.useGPU
    parallel.gpu.enableCUDAForwardCompatibility(true) 
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
    
end

slm = MeadowlarkOneK();
slm.stop();
slm.start();
Setup = slm.add_slm(Setup);

disp('SLM Ready!')

%% Basler

bas = bascam();
bas.start()

disp('Basler Ready!')

%% Basler preview

% slmCoords = [0.45 0.45 0.052 1];
slmCoords = [0.65 0.65 0.052 1];

[ HoloTemp,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos(Setup, slmCoords);
% DEestimateTemp = DEfromSLMCoords(slmCoordsTemp); %
% disp(['Diffraction Estimate for this spot is: ' num2str(DEestimateTemp)])
slm.feed(HoloTemp);
disp('sent SLM')

figure(124)
clf
imagesc(HoloTemp)
title('Hologram sent to SLM')

bas.preview()
