
clear
close all
clc

tBegin = tic;

disp('Setting up stuff...');

makePaths()
%% Setup sutter
sutter = sutterController();


%% connect to SI computer


%run this code first, then 'autoCalibSI'
disp('Waiting for msocket communication to ScanImage Computer')
%then wait for a handshake
srvsock2 = mslisten(42040);
SISocket = msaccept(srvsock2,15);
msclose(srvsock2);
sendVar = 'A';
mssend(SISocket, sendVar);
%MasterIP = '128.32.177.217';
%masterSocket = msconnect(MasterIP,3002);

invar = [];

while ~strcmp(invar,'B');
    invar = msrecv(SISocket,.1);
end;
disp('communication from Master To SI Established');

%% Setup
baseName = '''calib''';

% Where ScanImage stages these tiffs, ON THE SI COMPUTER. Derived from the rig's
% si_root (the SI tiff drive) instead of the literal 'D:\Calib\Temp' that was
% inlined twice below -- another scope changes one rig field rather than editing
% this script in two places. The literal is the fallback, so this rig is
% unchanged. Backslashes matter: this string is injected into a ScanImage command
% that runs on the Windows SI box.
try
    si_root = rig_remote_get('paths.si_root', 'D:');
catch
    si_root = 'D:';
end
si_calib_tmp = fullfile(si_root, 'Calib', 'Temp');

mssend(SISocket,'hSI.hStackManager.enable = 0 ;');

mssend(SISocket,'hSI.extTrigEnable = 0;'); %savign
mssend(SISocket,'hSI.hChannels.loggingEnable = 1;'); %savign
mssend(SISocket,['hSI.hScan2D.logFilePath = ''' si_calib_tmp ''';']);
mssend(SISocket,['hSI.hScan2D.logFileStem = ' baseName ';']);
mssend(SISocket,'hSI.hScan2D.logFileCounter = 1;');

mssend(SISocket,['hSICtl.updateView;']);

%% Collect PSF
sutter.setRef()

UZ= -50:1:50;%linspace(-100,100,21);

disp('Collecting PSF.')

for i=1:numel(UZ)
    fprintf('Plane %d/%d\n', i, numel(UZ))
    sutter.moveZ(UZ(i))

    if i==1
        pause(1)
    else
        pause(0.1);
    end

    mssend(SISocket,'hSI.startGrab()');
    invar = msrecv(SISocket,0.01);
    while ~isempty(invar)
        invar = msrecv(SISocket,0.01);
    end

    wait = 1;
    while wait
        mssend(SISocket,'hSI.acqState');
        invar = msrecv(SISocket,0.01);
        while isempty(invar)
            invar = msrecv(SISocket,0.01);
        end

        if strcmp(invar,'idle')
            wait=0;
        else
        end
    end

        
end

sutter.moveToRef()

%% Move onto server
disp('Moving files')
tMov = tic;

%on ScanImage Computer
% PER-SESSION value, deliberately still a literal: the folder name carries the
% date of the alignment it belongs to, so it is an operator choice per run, not
% rig configuration. Edit it when you run this. (A rig field here would suggest
% one fixed answer, which is what it must not be.)
destination = '''K:\objective_test\new_alignment_15Aug2023''';
% Same staging folder set above, from the rig -- not a second literal.
source = ['''' fullfile(si_calib_tmp, 'calib*') ''''];

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


%%
mssend(SISocket,'end');
disp('Done collecting PSF stack.')
