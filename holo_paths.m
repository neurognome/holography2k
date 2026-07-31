function holo_paths()
%HOLO_PATHS Put this holography2k checkout and its holodaq on the MATLAB path.
%   Call once at the top of any hand-run script here that needs holodaq's modules
%   (FiberPowerControl, ELL14, Output/DAQOutput) or its rigs/ layer (load_rig,
%   rig_get, open_serial, rig_remote_get).
%
%   Replaces the hardcoded trio these scripts used to open with:
%       addpath(genpath('C:\Users\holos\Documents\GitHub'))
%       addpath(genpath('C:\Users\holos\Documents\_code'))
%       addpath(genpath('C:\Users\holos\Documents\GitHub\holodaq'))
%   which named one machine's username and folder layout, and pulled in a whole
%   code tree (_code) that none of these scripts actually needs -- every function
%   they call lives in holodaq.
%
%   Resolution, mirroring holoexpt/holodaq_root.m so there is ONE convention:
%     this repo : derived from this file's own location
%     holodaq   : 1) the HOLODAQ_HOME environment variable
%                 2) a sibling directory named 'holodaq'
%                 3) the only sibling holding default_setup.m and rigs/load_rig.m
%                 4) otherwise ERROR with instructions
%
%   It errors rather than guessing, deliberately. A wrong holodaq means a wrong
%   rigs/ layer, i.e. the wrong physical channel map -- silently continuing on a
%   stale checkout is how you drive the wrong hardware.
%
%   APPENDS (addpath ... '-end') rather than prepending. addpath prepends by
%   default, which is how a genpath'd tree can shadow the copy you are actually
%   editing -- a mistake that already cost a debugging session here, where
%   load_rig resolved to a different checkout than the one being changed.
%
%   Runs once per MATLAB session; later calls are no-ops.
%
%   See also: holodaq_root (holoexpt), load_rig, start_holo_listener

    persistent done
    if ~isempty(done) && done
        return
    end

    root = fileparts(mfilename('fullpath'));   % this file sits at the repo root
    addpath(genpath(root), '-end');

    [dq_root, source] = local_holodaq_root(root);
    addpath(genpath(dq_root), '-end');

    fprintf('holo_paths: holography2k %s\n            holodaq      %s (%s)\n', ...
        root, dq_root, source);
    done = true;
end


function [root, source] = local_holodaq_root(here)
    env = strtrim(getenv('HOLODAQ_HOME'));
    if ~isempty(env)
        root = env;
        source = 'HOLODAQ_HOME';
        local_check(root, source);
        return
    end

    parent = fileparts(here);

    named = fullfile(parent, 'holodaq');
    if local_is_holodaq(named)
        root = named;
        source = 'sibling ''holodaq''';
        return
    end

    % Any sibling that looks like a holodaq checkout, accepted only if unique.
    hits = {};
    d = dir(parent);
    for i = 1:numel(d)
        if ~d(i).isdir || any(strcmp(d(i).name, {'.', '..'}))
            continue
        end
        cand = fullfile(parent, d(i).name);
        if ~strcmp(cand, here) && local_is_holodaq(cand)
            hits{end+1} = cand; %#ok<AGROW>
        end
    end
    if numel(hits) == 1
        root = hits{1};
        source = 'the only holodaq-shaped sibling';
        return
    end

    error('holo_paths:noHolodaq', ...
        ['Could not find the holodaq checkout these scripts need.\n' ...
         'Looked for: $HOLODAQ_HOME, then a sibling of\n  %s\n' ...
         'holding default_setup.m and rigs/load_rig.m (found %d).\n\n' ...
         'Fix either way:\n' ...
         '  setenv(''HOLODAQ_HOME'', ''<path to holodaq>'')   %% this session\n' ...
         '  or clone holodaq beside this repo\n' ...
         'Set it permanently in this machine''s startup.m if the layout is ' ...
         'unusual.'], here, numel(hits));
end


function tf = local_is_holodaq(p)
    tf = isfolder(p) ...
        && isfile(fullfile(p, 'default_setup.m')) ...
        && isfile(fullfile(p, 'rigs', 'load_rig.m'));
end


function local_check(root, source)
    if ~local_is_holodaq(root)
        error('holo_paths:badHolodaq', ...
            ['%s points at\n  %s\nwhich does not look like a holodaq checkout ' ...
             '(no default_setup.m and rigs/load_rig.m).'], source, root);
    end
end
