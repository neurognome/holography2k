function DE = DEfromSLMCoords(slmCoord)
% ActiveCalib.mat from the rig's calib folder -- the same folder
% align_slm_to_camera_scope2k saves ActiveCalib.mat into, and that
% find_latest_calib reads, so the read and write sides cannot drift apart.
%
% This was a literal naming user 'Holography' and the retired SLM_Management
% tree. That path does not exist on the current holography computer, and `load`
% of a missing file ERRORS -- so this function, which eleven other files here
% call, was already broken rather than merely unportable.
try
    calib_dir = rig_remote_get('paths.calib_dir', 'C:\Users\holos\Documents\calibs');
catch
    calib_dir = 'C:\Users\holos\Documents\calibs';
end
load(fullfile(calib_dir, 'ActiveCalib.mat'), 'CoC')

SIcoord = function_SLMtoSI(slmCoord,CoC);
DE = computeDEfromList(SIcoord',{1},1);
