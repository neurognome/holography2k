
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
mssend(SISocket,'hSI.hStackManager.enable = 0 ;');

mssend(SISocket,'hSI.extTrigEnable = 0;'); %savign
mssend(SISocket,'hSI.hChannels.loggingEnable = 1;'); %savign
mssend(SISocket,'hSI.hScan2D.logFilePath = ''D:\Calib\Temp'';');
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
destination = '''K:\objective_test\new_alignment_15Aug2023''';
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


%%
mssend(SISocket,'end');
disp('Done collecting PSF stack.')
