[Setup ] = function_loadparameters2();
Setup.CGHMethod=2;
Setup.GSoffset=0;
Setup.verbose =0;
% Setup.useGPU = 0;
%% choose one
slm_1100 = get_slm(1030);
slm_900 = get_slm(900);
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
slm_1100.stop();
slm_1100.wait_for_trigger = 0;
slm_1100.start();

slm_900.stop();
slm_900.wait_for_trigger = 0;
slm_900.start();
%%
%slmCoords = [.54 .46 0 1];
slmCoords = [.6 .4 0 1];

% %
% 
% 0.
% dxy = .01;
% x = linspace(.485-dxy, .485+dxy, 2);
% y = linspace(.55-dxy, .55+dxy, 2);
% [x,y] = meshgrid(x, y);
% x = x(:);
% y = y(:);
% % x(5) = [];
% % y(5) = [];
% slmCoords = [x(:), y(:), x(:)*0-.00, x(:)*0+1];


% slmCoords = [.4 .6 -.0 1]; % 0.

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
%slm_1100.feed(zeros(1024, 1024));
slm_1100.feed(Holo);

%%
%slmCoords = [.475 .52 -.0 1]; % 0.
slmCoords = [.375 .37 0 1]; % 0.  
%slmCoords = [.4 .4 0 1]; % 0.  

% dxy = .02;
% x = linspace(.6650-dxy, .6650+dxy, 2);
% y = linspace(.34-dxy, .34+dxy, 2);
% [x,y] = meshgrid(x, y);>?
% slmCoords = [x(:), y(:), x(:)*0+.14, x(:)*0+1];

[Holo, ~, ~ ] = function_Make_3D_SHOT_Holos( Setup,slmCoords );
%slm_900.feed(zeros(1024, 1024));
slm_900.feed(Holo);%%

% bas.preview()