%%
meadowlark = MeadowlarkOneK();
holoeye = HoloeyePLUTO();


%% choose one
slm = meadowlark;
blankHolo = zeros([1024 1024]);

% upper left: [.85 .9 0.00 1]
% lower left: [.15 .9 0.00 1]
% lower right: [.18 .07 0.00 1]
% upper right: [.88 .1 0.00 1]


%%
slm.stop();
slm.wait_for_trigger = 0;
slm.start();

slmCoords = [.4 .6 0.005 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );


slm.feed(Holo);
