function makePaths()
%MAKEPATHS Put this checkout, its holodaq, msocket and the SLM SDK on the path.
%   Called at the top of ~13 hand-run scripts here (alignment, PSF, tests). Kept
%   as the entry point they already call rather than editing all of them; new
%   code should call holo_paths() directly and add only what it needs.
%
%   WHAT THIS REPLACES. The previous version was four literals naming one
%   machine's username and folder layout:
%       rt = 'C:\Users\holos\Documents\GitHub\holography2k\'; cd(rt)
%       addpath(genpath('C:\Users\holos\Documents\GitHub\holography2k'));
%       addpath(genpath('C:\Users\holos\Documents\GitHub\msocket'));
%       addpath(genpath('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus'));
%       ... \SDK, and a version-pinned HOLOEYE folder
%   On any other holography computer every one of those was wrong, and the
%   failure was quiet: genpath of a missing folder returns '' and addpath('')
%   is a no-op, so the paths simply were not there and the first SLM call died
%   somewhere unrelated.
%
%   HOW EACH IS RESOLVED NOW:
%     this repo + holodaq : holo_paths() -- derived from this file's location,
%                           holodaq via HOLODAQ_HOME / sibling / error
%     msocket             : a sibling 'msocket' checkout, else the old literal
%     SLM SDK             : rig.paths.slm_sdk, via the config the DAQ published
%                           (this machine has no rig file of its own), else the
%                           old literal -- so this rig behaves identically
%     vendor runtimes     : still literals, but existence-checked. These are
%                           DEFAULT INSTALL locations under Program Files, not
%                           per-user paths, so they port far better than the
%                           above. A scope whose SLM lives elsewhere should set
%                           rig.paths.slm_sdk, which is the rig-driven hook.
%
%   Still cd's to the repo root, because callers rely on it for relative paths --
%   but to the DERIVED root, not a hardcoded one.
%
%   See also: holo_paths, start_holo_listener, find_latest_calib, rig_remote_get

    % This repo and holodaq. Errors with instructions if holodaq is not findable,
    % deliberately: without it there is no rig_remote_get, so nothing below could
    % be rig-driven and we would be back to guessing at literals.
    holo_paths();

    root = fileparts(mfilename('fullpath'));   % this file sits at the repo root
    cd(root)

    % ---- msocket -----------------------------------------------------------
    % Legacy transport (holochat replaces it), but the alignment and test scripts
    % here still use it. Resolved as a sibling checkout the way holo_paths finds
    % holodaq, so a machine that clones the two side by side needs no edit.
    ms = fullfile(fileparts(root), 'msocket');
    if ~isfolder(ms)
        ms = 'C:\Users\holos\Documents\GitHub\msocket';   % the old literal
    end
    local_add(ms, 'msocket');

    % ---- SLM SDK -----------------------------------------------------------
    % From the rig config the DAQ published; the literal is the fallback, not the
    % default. Same field and same resolution order start_holo_listener uses, so
    % the listener and these hand-run scripts cannot disagree about which SDK is
    % on the path.
    try
        slm_sdk = rig_remote_get('paths.slm_sdk', 'C:\Users\holos\Desktop\meadowlark');
    catch
        slm_sdk = 'C:\Users\holos\Desktop\meadowlark';
    end
    local_add(slm_sdk, 'slm_sdk (rig)');

    % ---- vendor runtimes ---------------------------------------------------
    % Default Program Files install locations. Existence-checked so a machine
    % without one does not silently accumulate a dead path entry. The HOLOEYE
    % entry pins a version in its folder name, which is why it is globbed.
    local_add('C:\Program Files\Meadowlark Optics\Blink OverDrive Plus', 'Meadowlark Blink');
    for h = local_glob('C:\Program Files\HOLOEYE Photonics', 'SLM Display SDK*')
        local_add(fullfile(h{1}, 'win64'), 'HOLOEYE SDK');
    end
end


% -------------------------------------------------------------------------
function local_add(p, what)
%LOCAL_ADD genpath+addpath one folder, saying so, and saying when it is missing.
%   The old code added these blind. genpath('') is '' and addpath('') is a no-op,
%   so a wrong path produced no error here and a confusing one much later.
    if isempty(p)
        return
    end
    if isfolder(p)
        addpath(genpath(p), '-end');   % append: prepending is how a stale tree shadows this one
        fprintf('makePaths: + %-18s %s\n', what, p);
    else
        fprintf('makePaths:   %-18s NOT FOUND, skipped: %s\n', what, p);
    end
end


function out = local_glob(parent, pattern)
%LOCAL_GLOB Subfolders of `parent` matching `pattern`, as a row cell of full paths.
%   Only the LAST component is wildcarded -- MATLAB's dir does not reliably glob a
%   middle path component, which is why the caller passes the parent separately
%   and appends any leaf itself.
    out = {};
    if ~isfolder(parent)
        return
    end
    d = dir(fullfile(parent, pattern));
    for i = 1:numel(d)
        if d(i).isdir && ~any(strcmp(d(i).name, {'.', '..'}))
            out{end+1} = fullfile(parent, d(i).name); %#ok<AGROW>
        end
    end
end
