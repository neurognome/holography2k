

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
% Setup.useGPU = 0;
%% choose one
wavelength = 589
slm = get_slm(wavelength);
blankHolo = zeros([1024 1024]);

%1030


%%
slm.stop();
slm.wait_for_trigger = 0;
slm.start();
%%
%slmCoords = [.475 .52 -.0 1]; % 0.
slmCoords = [513/1024 513/1024 0 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );


slm.feed(blankHolo);
% pwr= 3;
% mssend(masterSocket,[pwr/1000 1 1]); % again, check if this is mW or W
%  
% bas.preview()
% mssend(masterSocket,[0 1 1]);
% %1030/1100
% top right: 0.85, 0.08
% bottom right  0.11, 0.06
% bottom left 0.11, 0.93
% top left 0.85 0.93
% z range:  -0.01, 0.12, offset: 0


%900
% top right: 0.1, 0.1
% bottom right: 0.85, 0.1
% bottom left: 0.85. 0.85
% top left: 0.1 0.85 

%% make meshgrid (testing 589)
x = linspace(0.1, .9, 5);
[x,y]=meshgrid(x,x);
slmCoords = [x(:), y(:), x(:)*0, x(:)*0+1];
[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );


slm.feed(Holo);

%% test out making bigger spots

slmCoords = [.65 .69 0.0080 1]; % 0.

%[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,slmCoords,20 );

slm.feed(Holo);

%% make meshgrid (testing 589) of big spots
x = linspace(0.1, .9, 5);
[x,y]=meshgrid(x,x);
slmCoords = [x(:), y(:), x(:)*0, x(:)*0+1];
[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,slmCoords, 0);


slm.feed(Holo)

%% make big spots, but divert some power to zero order

slmCoords = [.65 .35 0 1]; % 0.
%slmCoords = [513/1024 513/1024 0 1]; 
%slmZero = [513/1024 513/1024 0 1]; 
%slmCoords = [slmCoords; slmZero];

%[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,slmCoords,0 );

slm.feed(Holo);