disp('Generating New Holograms...')
disp('Everything after this should be automated so sitback and enjoy')

npts = 250; %You can almost get through 750 with water before it evaporates.


% 12/29/22 WH - should be roughly +145 um (-0.04 SLM) to -30 um (0.025 SLM)

dummy = rand;

slmCoords=zeros(4,npts);
for i =1:npts
    slmCoords(:,i) = [...
        rand*(slmXrange(2)-slmXrange(1))+slmXrange(1),...
        rand*(slmYrange(2)-slmYrange(1))+slmYrange(1),...
        rand*(slmZrange(2)-slmZrange(1))+slmZrange(1),...
        1];
end

figure(205);scatter3(slmCoords(1,:),slmCoords(2,:),slmCoords(3,:),'o')
drawnow;
title('Holograms in SLM space')
%%compile random holograms

slmCoords(3,:) = round(slmCoords(3,:),3); %Added 3/15/21 by Ian for faster compute times
disp('Compiling Single Holograms')
parfor i =1:npts
    t=tic;
    fprintf(['Holo ' num2str(i)]);
    subcoordinates = slmCoords(:,i);
    
    [ Hologram,Reconstruction,Masksg ] = function_Make_3D_SHOT_Holos( Setup,subcoordinates' );
    hololist(:,:,i)=Hologram;
    fprintf([' took ' num2str(toc(t)) 's\n']);
end

out.hololist = hololist;
out.slmCoords = slmCoords;
save('tempHololist4_will.mat','out');