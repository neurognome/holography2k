function [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_tuneZ( Setup,Coordinates,zparams)
% Coordinates is a N points by 3 or by 4 map of x,y,z optional power
% coordinates.

[~,LP] = size(Coordinates);
if LP==4
    %Case where you provide adjusted intensity coefficients
    %disp('You have specified power in each target')
    power = Coordinates(:,4);
else
    power = Coordinates(:,3)-Coordinates(:,3);
    power = power +1;
    disp('Warning : You have not specified power distribution in targets, we ll use naive way')
end

xx = Coordinates(:,1);
yy = Coordinates(:,2);
zz = Coordinates(:,3);

badpoints =  double((double(xx<0) +double(yy<0)+ double(xx>1) +double(xx>1))>0);
if sum(badpoints)>0
    disp('Bad Warning : We are going to ignore some points that are out of SLM range') %to updat with an error idalog box
end
xx = xx(badpoints==0);
yy = yy(badpoints==0);
zz = zz(badpoints==0);
power = power(badpoints==0);

xx = floor(xx*(Setup.SLM.Nx-1)+1);
yy = floor(yy*(Setup.SLM.Ny-1)+1);

%Should you need to visualize your holograms uncomment the paragraph below
%f = figure(1);
%scatter3(xx,yy,zz,[],power,'filled'); xlabel('X, SLM coordinates');ylabel('X, SLM coordinates');zlabel('X, SLM coordinates');
%pause(2);close(f);

% Sorting points by depth
[a,b] = sort(zz);
xx = xx(b);
yy = yy(b);
zz = zz(b);
power = power(b);

Depths = unique(zz);
NZ = numel(Depths);
if Setup.verbose==1
disp(['Your hologram has ' int2str(NZ) ' distinct levels']);
end
Masks = zeros(Setup.SLM.Nx,Setup.SLM.Ny,NZ);

for i = 1:NZ
    select = (zz==Depths(i));
    sxx = xx(select);
    syy = yy(select);
    spower = power(select);
    for j = 1:numel(sxx)
        Masks(sxx(j),syy(j),i) = spower(j);
    end
end

Masksg = Masks;

if Setup.useGPU ==1
    Masks = gpuArray(Masks);
end

[ HStacks ] = function_Hstacks( Setup,Depths );

%%%%% COMPUTE HOLOGRAMS FROM MASKS HERE
if Setup.CGHMethod == 1 % Case Superposition
    [Holo] = function_Superposition( Setup, HStacks, Masks );
elseif  Setup.CGHMethod == 2 % Case Global GS
    [Holo] = function_globalGS(Setup, HStacks, Masks );
elseif  Setup.CGHMethod == 3 % NovoCGH
    [Holo] = function_NOVO_CGH_VarIEuclid( Setup, HStacks, Masks,Depths );
elseif  Setup.CGHMethod == 4 % Two photon NOVO CGH
    [Holo] = function_NOVO_CGH_TPEuclid( Setup, HStacks, sqrt(Masks),Depths );
else
    disp('There is no such CGH method')
end

%%%%%%%%% CLEAN UP


%%%%%%%%%%%%%%%% this manually add zOffset to the generated hologram
% system parameter
imsize      = Setup.Nx; % size of the entire FOV
wavelength  = Setup.lambda*1e9; % nm
NA          = 1;
mag         = 20;    
pixelSize   = Setup.psSLM*1e6; % micron
zDefocus    = zparams;   % micron

% calculate k-space parameter
k0 = 2*pi/wavelength;
kNA = 2*pi*NA/wavelength;


% convert to k-space grid
kMax = pi*mag/(pixelSize*1e3);
xc = floor(imsize/2+1);
yc = floor(imsize/2+1);
[Y,X] = meshgrid((1:imsize) - yc, (1:imsize) - xc);
kMap = abs(X + 1i*Y)*kMax/floor(imsize/2);

DefocusPhase = -(sqrt(k0.^2 - (kMap.*(kMap<kNA)).^2) - k0)*zDefocus*1e3;
Holo.phase = wrapToPi(Holo.phase + DefocusPhase);

%%%%%%%%%%%%%%% end: this manually add zOffset to the generated hologram


Hologram = uint8(floor((Setup.SLM.pixelmax*(Holo.phase+pi)/(2*pi))));
%f = figure(1); imagesc(Hologram); axis image; pause(2); close(f)

[ Reconstruction ] = function_VolumeIntensity( Setup,Holo.phase,HStacks);



end
