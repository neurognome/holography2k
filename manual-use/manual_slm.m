

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
% Setup.useGPU = 0;
%% choose one
wavelength = 1030;
slm = get_slm(wavelength);
blankHolo = zeros([1024 1024]);

% upper left: [.85 .9 0.00 1]
% lower left: [.15 .9 0.00 1]
% lower right: [.18 .07 0.00 1]
% upper right: [.88 .1 0.00 1]

% z lower: [-0.055 , 0.03]

% for 900:
% upper left: [.1 .88 0.00 1]
% lower left: [.77 .85 0.00 1]
% lower right: [.77 .05 0.00 1]
% upper right: [.1 .05 0.00 1]
% z lower: [-0.05 , 0.03]


%%
slm.stop();
slm.wait_for_trigger = 0;
slm.start();

slmCoords = [.4 .5 .00 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );


slm.feed(Holo);

bas.preview()