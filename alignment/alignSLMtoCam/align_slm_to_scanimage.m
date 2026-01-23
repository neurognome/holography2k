function align_slm_to_camera_scope2k
%% Pathing
in = input('Run Calibration? (y/n)','s');
if ~strcmp(in,'y')
    disp('Did not detect ''y'' so not executing')
    return
end

clear;close all;clc
%%
tBegin = tic;
%
makePaths()
disp('done pathing')

castImg = @uint8;
castAs = 'uint8';
camMax = 255;

%% Setup Stuff
disp('Setting up stuff...');

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU = 0;

Setup.useThorCam =0;
Setup.maxFramesPerAcquire = 3; %set to 0 for unlimited (frames will return will be
Setup.camExposureTime = 10000;

calibration_wavelength = ['1100', '900']; % ENSURE 1100 is first or else.. bad..

if Setup.useGPU
    disp('Getting gpu...'); %this can sometimes take a while at initialization
    g= gpuDevice;
end

delete(gcp('nocreate'));
parpool('IdleTimeout', 360);

% setup slms
slm = [];
for w = calibration_wavelength
    slm(end+1) = get_slm(str2double(w)); % does this work? idk we'll see
end

for s = slm
    s.stop();
    s.wait_for_trigger = 0;
    s.start();
end

% setup basler camera
bas = bascam();
bas.start()

% setup sutter
sutter = sutterController();

% collect everything into devices
devices = struct();
devices.bas = bas;
devices.slm = slm;
devices.sutter = sutter;

% create an empty calibration structure ,which will store all the calibration information
% i think we can think about this a little, but let's not worry for now.
CoC = dictionary(calibration_wavelength, [[], []]); % each of them are stored independently?

disp('Ready')
%% look for objective in 1p

bas.preview()
% bas.preview_set_cmax(80)

%% Make mSocketConnections with DAQ and SI Computers
comm = HolochatInterface('holo');
comm.send(calibration_wavelength, 'daq');
devices.comm = comm;
%% Put all Manual Steps First so that it can be automated
%% Create a random set of holograms or use flag to reload
disp('First step Acquire Holograms')
tSingleCompile = tic;

% for 1100
x_1100 = [0.18 0.9];
y_1100 = [0.0 0.9];
z_1100 = [-0.05 0.08];

% for 900
x_900 = [0.18 0.9];
y_900 = [0.0 0.9];
z_900 = [-0.05 0.08];

slmXrange = dictionary(calibration_wavelength, [x_1100, x_900]);
slmYrange = dictionary(calibration_wavelength, [y_1100, y_900]);
slmZrange = dictionary(calibration_wavelength, [z_1100, z_900]);

hololist = dictionary();
slmCoords = dictionary();
for w = calibration_wavelength
    hololist(w), slmCoords(w) = create_holograms(slmXrange(w), slmYrange(w), slmZrange(w)); % loop through, i think this is the template..
end

disp(['Done compiling holograms. Took ' num2str(toc(tSingleCompile)) 's']);
%% Set Power Levels
p_1100 = 5; % in mW
p_900 = 5;
pwr = dictionary(calibration_wavelength, [p_1100, p_900]);
disp(['idividual hologram; power set to ' num2str(pwr) 'mW']);

%% test each
for wv = calibration_wavelength
    fprintf('Find the spot and check if this is the right amount of power for %s\n', w)
    % this needs to change depending on which SLM... probably
    slmCoords = [0.7 0.7 0.0 1];
    [Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup, slmCoords );

    devices.slm.feed(Holo);
    devices.comm.send([pwr(wv)/1000, 1, 1], wv) % be explicit

    devices.bas.preview();
    devices.comm.send([0, 1, 1], wv)
end

%% Make Sure you're centered
disp('Find Focal Plane, Center and Zero the Sutter')
disp('Leave the focus at Zoom 1. at a power that is less likely to bleach (14% 25mW)') %25% 8/16/19
disp('Don''t forget to use Ultrasound Gel on the objective so it doesn''t evaporate')
for wv = calibration_wavelength
    devices.comm.send([0, 1, 1], wv)
end

devices.bas.preview()

disp('Turn off Focus and press any key to continue');
pause
devices.sutter.setRef()

devices.comm.send([0, 0], 'si');

disp('Make Sure the DAQ computer is running testMultiTargetsDAQ. and the SI computer running autoCalibSI');
disp('also make those names better someday')
disp('Make sure both lasers are on and the shutters open')
disp('Scanimage should be idle, nearly in plane with focus. and with the gain set high enough to see most of the FOV without saturating')
disp('THE MOUSE MONITOR SHOULD BE TURNED OFF')

devices.sutter.moveZ(100);
disp('testing the sutter double check that it moved to reference +100');
disp('Ready to go (Press any key to continue)');
pause

devices.sutter.moveToRef();

%% Scan Image Planes Calibration
disp('Begining SI Depth calibration, we do this first incase spots burn holes with holograms')
zsToUse = linspace(0, 90, 15);% %70 was about 125um on 3/11/21 %Newer optotune has more normal ranges 9/28/29; New Optotune has different range 9/19/19; [0:10:89]; %Scan Image Maxes out at 89
SIUZ = -15:5:130;% linspace(-120,200,SIpts);
CoC_temp = align_optotune_to_camera(CoC, devices, zsToUse, SIUZ);
for wv = calibration_wavelength
    CoC(wv) = CoC_temp; % propagate to both
end

%% Coarse Data
% now the fun begins...
% let's just wrapt his all up for now..
% this will be run in a loop
for wv = calibration_wavelength % independent!!
    [CoC(wv), SImatchRangeX, SImatchRangeY] = align_slm_to_camera(CoC(wv), devices, wv, hololist(wv), slmCoords(wv));
end
%%%%%%%%%% hole burns start here
% we need some data from the slm to camera to pass into the next step, how can we pass this cleanly?
%%  align camera to scanimage, aka hole burns
scanimage_power = 10; % in percent (%)
burnPowerMultiplier = 10; %X
burnTime = 3; % in s
wv = '900'; % only perform hole burns on the 900 path, sad for 1100
CoC(wv) = align_slm_to_scanimage(CoC(wv), info, devices, wv, SImatchRangeX, SImatchRangeY, zsToUse, slmXrange(wv), slmYrange(wv), slmZrange(wv),...
                            scanimage_power, burnPowerMultiplier, burnTime);


% finally, propagate the CoC to both...
%% Save Output Function
pathToUse = 'C:\Users\holos\Documents\SLM_Management\Calib_Data';
disp('Saving...')
tSave = tic;

save(fullfile(pathToUse,[date '_Calib_' num2str(calibration_wavelength) '.mat']),'CoC')
save(fullfile(pathToUse,'ActiveCalib.mat'),'CoC')

pth = 'C:\Users\holos\Documents\calibs';
save(fullfile(pth,[date '_Calib_' num2str(calibration_wavelength) '.mat']),'CoC')
save(fullfile(pth,'ActiveCalib.mat'),'CoC')
save(fullfile(pth,['CalibWorkspace_' num2str(calibration_wavelength) '.mat']), '-v7.3');

save(fullfile(pathToUse,'CalibWorkspace_v73.mat'), '-v7.3');
% save(fullfile(pathToUse,['CalibWorkspace_will_' date '.mat']))
disp(['Saving took ' num2str(toc(tSave)) 's']);

disp(['All Done, total time from begining was ' num2str(toc(tBegin)) 's. Bye!']);

slm.stop()