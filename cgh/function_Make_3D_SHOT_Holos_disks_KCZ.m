function [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos_disks_KCZ( Setup,Coordinates,disk_radius )
% Coordinates is a N points by 3 or by 4 map of x,y,z optional power
% coordinates.
% disk_diameter: radius in pixels; set to 0 if you just want dots

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

if disk_radius ~= 0
    n_sigma = 4;
    disk_window_size = disk_radius*2*n_sigma+1;
    x_disk = linspace(-n_sigma, n_sigma, disk_window_size);
    [x_disk, y_disk] = meshgrid(x_disk, x_disk);
    disk = exp(-(x_disk.^2 + y_disk.^2)/2) > 0.05;
    
    for i = 1:NZ
        select = (zz==Depths(i));
        sxx = xx(select);
        syy = yy(select);
        spower = power(select);
        for j = 1:numel(sxx)
            % r_start = sxx(j)-disk_radius*n_sigma;
            % r_end = sxx(j)+disk_radius*n_sigma;
            % c_start = syy(j)-disk_radius*n_sigma;
            % c_end = syy(j)+disk_radius*n_sigma;
            % Masks(r_start:r_end, c_start:c_end,i) = spower(j)*disk;

            % chatgpt code to handle edge cases:

            % Define the start and end indices
            r_start = sxx(j) - disk_radius * n_sigma;
            r_end = sxx(j) + disk_radius * n_sigma;
            c_start = syy(j) - disk_radius * n_sigma;
            c_end = syy(j) + disk_radius * n_sigma;
            
            % Get the size of the Masks array
            [mask_rows, mask_cols, ~] = size(Masks);
            
            % Check and adjust row indices
            if r_start < 1
                disk_r_start = 2 - r_start; % Adjust start index for disk
                r_start = 1;
            else
                disk_r_start = 1;
            end
            
            if r_end > mask_rows
                disk_r_end = size(disk, 1) - (r_end - mask_rows);
                r_end = mask_rows;
            else
                disk_r_end = size(disk, 1);
            end
            
            % Check and adjust column indices
            if c_start < 1
                disk_c_start = 2 - c_start; % Adjust start index for disk
                c_start = 1;
            else
                disk_c_start = 1;
            end
            
            if c_end > mask_cols
                disk_c_end = size(disk, 2) - (c_end - mask_cols);
                c_end = mask_cols;
            else
                disk_c_end = size(disk, 2);
            end
            
            % Place the cropped disk into Masks
            Masks(r_start:r_end, c_start:c_end, i) = spower(j) * disk(disk_r_start:disk_r_end, disk_c_start:disk_c_end);
        end
    end
    
else
    
    for i = 1:NZ
        select = (zz==Depths(i));
        sxx = xx(select);
        syy = yy(select);
        spower = power(select);
        for j = 1:numel(sxx)
            Masks(sxx(j),syy(j),i) = spower(j);
        end
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
Hologram = uint8(floor((Setup.SLM.pixelmax*(Holo.phase+pi)/(2*pi))));
%f = figure(1); imagesc(Hologram); axis image; pause(2); close(f)

[ Reconstruction ] = [];%function_VolumeIntensity( Setup,Holo.phase,HStacks);



end
