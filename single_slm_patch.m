[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
% Setup.useGPU = 0;
%% choose one
slm_1030 = get_slm(1030);
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


% 240326
% bottom left:  [1, 1]
% bottom right: [0.15, 1]
% top right: [0.15 ,0]
% top left: [1, 0.1]
%z 0.14

%%

slm_1030.stop();
slm_1030.wait_for_trigger = 0;
slm_1030.start();

%%
% slmCoords = [.44 .61 0 1]; % 0.
slmCoords = [0.4 0.4 0.02 1];
% slmCoords = [.64 .61 -.02 1]; % 0.
% slmCoords = [.64 .61 -.02 1;
%     .40 .64 -.02 1;];

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
% [Holo, ~, ~ ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,slmCoords, 5);
slm_1030.feed(Holo);
% slm_1030.feed(blankHolo);
bas.preview();

%%
imagesc(abs(fftshift(fft2(exp(1i*double(Holo)/255*2*pi)))))
shg