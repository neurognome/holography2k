function tlCamera = thorSetupCam

% Adapted from Thorlabs
try
    thorStop;
    disp('stopped running cam')
catch
end
% clear
% close all

% Load TLCamera DotNet assembly. The assembly .dll is assumed to be in the
% same folder as the scripts.
oldLoc = cd;
% THIS folder -- the .dll's live beside this file (see the vendored
% Thorlabs.TSI.*.dll here, and 'Copy 64-bit DotNet managed libraries here.txt'),
% which is exactly what the comment above says. It used to be the literal
% 'C:\Users\Holography\Desktop\holography\ThorCam\': a different machine's
% username AND a copy outside any checkout, so on the current holography computer
% this cd already failed and the camera could not be opened at all.
%
% No rig field for this one on purpose: the assemblies are vendored in the repo,
% so the checkout's own location is the correct and fully portable answer.
cam_path = [fileparts(mfilename('fullpath')) filesep];
pause(0.1)
cd(cam_path);
NET.addAssembly([cam_path, 'Thorlabs.TSI.TLCamera.dll']);
% NET.addAssembly('C:\Users\ian\Dropbox\code\Scientific Camera Interfaces\Thorlabs.TSI.TLCamera.dll');

disp('Dot NET assembly loaded.');

tlCameraSDK = Thorlabs.TSI.TLCamera.TLCameraSDK.OpenTLCameraSDK;

% Get serial numbers of connected TLCameras.
serialNumbers = tlCameraSDK.DiscoverAvailableCameras;
disp([num2str(serialNumbers.Count), ' camera was discovered.']);

disp('Opening the first camera')
tlCamera = tlCameraSDK.OpenCamera(serialNumbers.Item(0), false);

% Set exposure time and gain of the camera.
tlCamera.ExposureTime_us = 1000000;

% Check if the camera supports setting "Gain"
gainRange = tlCamera.GainRange;
if (gainRange.Maximum > 0)
    tlCamera.Gain = 0;
end

% Set the FIFO frame buffer size. Default size is 1.
tlCamera.MaximumNumberOfFramesToQueue = 5;
% arm the camera
tlCamera.OperationMode = Thorlabs.TSI.TLCameraInterfaces.OperationMode.SoftwareTriggered;
tlCamera.FramesPerTrigger_zeroForUnlimited = 0;
tlCamera.Arm;
disp('Camera armed.')

other.serialNumbers = serialNumbers;
other.tlCameraSDK = tlCameraSDK;

global TLCAMSDK TLCAM SERIALNUMBERS
TLCAMSDK = tlCameraSDK;
TLCAM = tlCamera;
SERIALNUMBERS= serialNumbers;
cd(oldLoc)


