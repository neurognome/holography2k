%%
meadowlark2 = MeadowlarkOneK();
holoeye = HoloeyePLUTO();

[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
Setup.useGPU = 0;
%% choose one
slm2 = meadowlark2;
blankHolo = zeros([1024 1024]);

% upper left: [.85 .9 0.00 1]
% lower left: [.15 .9 0.00 1]
% lower right: [.18 .07 0.00 1]
% upper right: [.88 .1 0.00 1]

% z lower: [-0.055 , 0.03]

%%
slm2.stop();
slm2.wait_for_trigger = 0;
slm2.start();

slmCoords = [.4, .4 0.005 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );


slm.feed(Holo);
