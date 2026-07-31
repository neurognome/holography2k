function start_holo_listener(varargin)
%START_HOLO_LISTENER Config-primed persistent holography listener (remote).
%   Run once on the holo computer INSTEAD of MsocketHolorequest2K when driving
%   experiments remotely. It waits for a prime from the DAQ master
%   (config/holo, see prime_info.m), auto-loads the newest dated calibration
%   per wavelength (find_latest_calib), compiles holograms (via the DAQ's
%   transferHR), arms the SLM(s), acks to config/holo_status, then serves
%   firing-order sequences. It RE-PRIMES automatically for the next experiment
%   -- no restart, no hand-edited calibration paths.
%
%   The original MsocketHolorequest2K.m is unchanged and remains the manual
%   fallback.
%
%   Configuration comes from the rig. This machine has no rig file of its own, so
%   values are read from what the DAQ published to holochat via
%   publish_rig_config (rig_remote_get: config/rig -> a local rig -> a literal).
%   The listener prints each value and which of those three supplied it. If the
%   broker is cold it starts anyway on the old literals and says so, because a
%   listener that refuses to boot is worse than one running on last week's paths.
%
%   Options (an explicit value always outranks the rig):
%     'Wavelengths' (default: the rig's opto table) -- THIS MACHINE'S SLM
%         inventory. Pass a list only to open a deliberate subset. Changing the
%         set requires restarting this listener, because the boards are opened
%         once at startup.
%     'CalibDir'    (default rig.paths.calib_dir)
%     'Timeout'     (default rig.holo.slm_timeout_ms) -- SLM trigger timeout.
%
%   Each opto channel may pin its own slm_board and slm_lut in the rig file, which
%   is what a rig with two arms on one wavelength needs (the board cannot be
%   derived from the wavelength then). A channel that pins neither falls back to
%   get_slm's wavelength mapping, so an unpinned rig behaves exactly as before.
%
%   Paths: this checkout is added by deriving it from this file's own location,
%   so it no longer names one machine's username. The Meadowlark SDK is a
%   per-machine install rather than part of a checkout, so it comes from
%   rig.paths.slm_sdk. holodaq must already be on this machine's path -- it has to
%   be, since HolochatInterface is what reads the rig config in the first place.
%
%   How re-prime stays race-free: the serve loop reads msg/holo and switches on
%   type -- a firing-order CELL is played; a holoRequest STRUCT means the next
%   experiment has begun (its config prime is already posted), so it is handed
%   back to be compiled. Sequences and holoRequests therefore never collide.
%
%   See also: MsocketHolorequest2K (fallback), PlaySequence2K, find_latest_calib,
%             generate_holograms_new, prime_info.

    % '' means "ask the rig" -- resolved below, once holodaq is reachable. An
    % explicit value always wins, so every documented invocation still works.
    p = inputParser;
    % [] means "take the inventory from the rig's opto table" (local_slm_inventory).
    % An explicit list still wins, for a machine deliberately opening a subset.
    p.addParameter('Wavelengths', []);
    p.addParameter('CalibDir', '');
    p.addParameter('Timeout', []);
    p.parse(varargin{:});
    opt = p.Results;

    % This checkout, wherever it happens to live. Was
    % genpath('C:\Users\holos\Documents\GitHub\holography2k'), which named one
    % machine's username and folder layout; deriving it from this file's own
    % location does the same job on any machine. holodaq is NOT added here -- it is
    % already on this machine's path (HolochatInterface below has always been
    % resolved before these addpaths ran).
    addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

    comm = HolochatInterface('holo');
    comm.flush();   % drop any stale msg/holo from a previous session

    % ---- config -------------------------------------------------------------
    % This machine has no rig file, so everything comes from what the DAQ
    % published to holochat (rig_remote_get: config/rig -> local rig -> literal).
    % Reported explicitly, because the previous behaviour was a silent fallback to
    % whatever was hardcoded here.
    [calib_dir, calib_src] = local_cfg(opt.CalibDir, 'paths.calib_dir', ...
        'C:\Users\holos\Documents\calibs');
    % The Meadowlark SDK is genuinely per-machine (an install location, not part of
    % any checkout), so it stays a rig value rather than being derived.
    [slm_sdk, sdk_src] = local_cfg([], 'paths.slm_sdk', ...
        'C:\Users\holos\Desktop\meadowlark');
    if ~isempty(slm_sdk) && isfolder(slm_sdk)
        addpath(genpath(slm_sdk));
    elseif ~isempty(slm_sdk)
        warning('holo_listener:noSlmSdk', ...
            ['The Meadowlark SDK folder does not exist:\n  %s\nget_slm will fail ' ...
             'to open a board. Set paths.slm_sdk in the rig file and re-run\n' ...
             'publish_rig_config() on the DAQ.'], slm_sdk);
    end
    [cgh_method, cgh_src]  = local_cfg([], 'holo.cgh_method', 2);
    [use_gpu, gpu_src]     = local_cfg([], 'holo.use_gpu', true);
    [timeout_ms, to_src]   = local_cfg(opt.Timeout, 'holo.slm_timeout_ms', 1700);

    fprintf('--- holo listener config ---\n');
    fprintf('  calib dir   : %s   (%s)\n', calib_dir, calib_src);
    fprintf('  slm sdk     : %s   (%s)\n', slm_sdk, sdk_src);
    fprintf('  cgh method  : %d   (%s)\n', cgh_method, cgh_src);
    fprintf('  use gpu     : %d   (%s)\n', use_gpu, gpu_src);
    fprintf('  slm timeout : %d ms   (%s)\n', timeout_ms, to_src);
    if ~isfolder(calib_dir)
        warning('holo_listener:noCalibDir', ...
            ['The calibration folder does not exist:\n  %s\nfind_latest_calib ' ...
             'will fail when a prime arrives. Fix paths.calib_dir in the rig file ' ...
             'and\nre-run publish_rig_config() on the DAQ, then restart me.'], calib_dir);
    end

    Setup = function_loadparameters2();
    Setup.CGHMethod          = cgh_method;
    Setup.verbose            = 0;
    Setup.useGPU             = double(logical(use_gpu));
    Setup.TimeToPickSequence = 0.05;
    Setup.SLM.timeout_ms     = timeout_ms;

    % SLM objects are created ONCE, because a Meadowlark board cannot safely be
    % reopened per experiment. So this is THIS MACHINE'S SLM INVENTORY, not the
    % experiment's channel set: a prime's channels are resolved against it by
    % wavelength (local_resolve_slots), never by position. It is also exactly why
    % the inventory comes from PERSISTENT rig config rather than from the prime --
    % it has to exist before the first prime ever arrives.
    [slm_all, wl, inv_src] = local_slm_inventory(opt.Wavelengths);
    fprintf('  slm source  : %s\n', inv_src);

    % Two wavelengths can map to ONE physical board -- get_slm sends both 1100 and
    % 1030 to board 2. That would drive one SLM with two hologram stacks, each
    % overwriting the other, discovered only as garbage on the sample. Refuse at
    % startup instead. This is the check that makes pinning slm_board in the rig
    % file safe, so it applies to a rig-supplied inventory too.
    boards = [slm_all.board_id];
    assert(numel(unique(boards)) == numel(boards), 'holo_listener:sharedBoard', ...
        ['Wavelengths %s resolve to boards %s -- at least two share one board.\n' ...
         'That would drive one physical SLM with two hologram stacks. Give each ' ...
         'arm its own\nslm_board in the rig file (see ExampleRig''s split-arm ' ...
         'example), or restart with a\nwavelength set that maps 1:1.'], ...
        mat2str(wl), mat2str(boards));
    fprintf('SLM inventory: %s nm on boards %s.\n', mat2str(wl), mat2str(boards));

    last_seq = -inf;
    preread  = [];
    last_abort_seq = local_current_abort(comm);   % ignore any stale abort at startup
    % ignore a stale prime left on config/holo from a previous session
    try
        c0 = comm.scan_config('holo');
        if isstruct(c0) && isfield(c0, 'prime_seq') && ~isempty(c0.prime_seq)
            last_seq = c0.prime_seq;
        end
    catch
    end
    fprintf('Holo listener up; waiting for a prime from the DAQ...\n');
    while true
        prime = comm.get_config();
        if ~local_is_new(prime, last_seq)
            pause(0.25);
            continue
        end
        last_seq = prime.prime_seq;
        stem = local_stem(prime);

        % Which channels is the DAQ actually running, in which order? Taken from
        % the prime, NOT from this listener's launch parameter -- the old code
        % compiled hologram i against calib(i) of its OWN wavelength list, so a
        % listener started with the order reversed silently built every hologram
        % from the wrong wavelength's SLM calibration.
        [chans, source] = local_prime_channels(prime, wl);
        fprintf('Priming holo for %s (%s from %s)...\n', ...
            stem, local_describe(chans), source);

        try
            % 1) map each channel to an SLM slot BY WAVELENGTH, and load that
            %    channel's newest calibration. Position is never assumed.
            slot  = local_resolve_slots(chans, wl);
            slm   = slm_all(slot);
            % The boards were opened at startup; confirm they are still the ones
            % the DAQ declares, so a rig-file edit mid-session cannot quietly
            % point a channel at the wrong physical SLM.
            local_assert_boards(chans, slm);
            calib = [];
            cname = cell(1, numel(chans));
            for k = 1:numel(chans)
                [c, f] = local_calib_for(chans(k).wavelength, calib_dir);
                calib  = [calib, c]; %#ok<AGROW>
                cname{k} = f;
            end

            % Tell the DAQ what we resolved, BEFORE compiling. Its own topic, so
            % holo_status keeps meaning "primed and armed" and the launcher's
            % readiness lamp cannot go green early.
            local_declare(comm, prime, true, 'channels resolved', chans, source, cname);

            % 2) compile holograms. Each arriving holoRequest is matched to its
            %    channel by the tag transferHR stamps on it; an untagged request
            %    (a pre-opto DAQ) falls back to the next unfilled slot, which
            %    reproduces the old positional behaviour exactly.
            holograms = cell(1, numel(chans));
            got = false(1, numel(chans));
            for n = 1:numel(chans)
                HR = preread;
                preread = [];
                if isempty(HR)
                    HR = [];   % generate_holograms_new reads it itself
                end
                k = local_match_channel(HR, chans, got);
                if isempty(HR)
                    hololist = generate_holograms_new(comm, Setup, calib(k));
                else
                    hololist = generate_holograms_new(comm, Setup, calib(k), HR);
                end
                holograms{k} = uint8(hololist);
                got(k) = true;
            end
            assert(all(got), 'holo_listener:missingChannel', ...
                'No holoRequest arrived for channel(s): %s.', ...
                strjoin({chans(~got).name}, ', '));
            comm.flush();

            % 3) (re)arm the SLM(s) for triggered playback.
            for s = slm
                s.stop();
                s.wait_for_trigger = 1;
                s.timeout_ms = timeout_ms;
                s.start();
            end

            fprintf('Holo primed for %s; serving sequences.\n', stem);
            local_ack(comm, prime, true, 'holograms compiled');
        catch err
            fprintf('Holo prime FAILED: %s\n', err.message);
            % STOP THE SLMs FIRST. The arm loop is the last step above, so a throw
            % anywhere before it used to leave the boards armed with the PREVIOUS
            % experiment's holograms -- while the DAQ, having been told nothing,
            % went on to gate the laser against them.
            for s = slm_all
                try, s.stop(); catch, end
            end
            local_declare(comm, prime, false, err.message, chans, source, {});
            local_ack(comm, prime, false, err.message);
            preread = [];
            continue
        end

        % 4) serve until the next prime appears, the next holoRequest arrives, or
        %    an abort. Passing last_seq lets the serve loop notice a new prime on
        %    config/holo instead of only reacting to traffic on msg/holo.
        [preread, aborted, reprime] = local_serve(comm, slm, holograms, ...
            last_abort_seq, last_seq);
        if aborted
            for s = slm, try, s.stop(); catch, end, end   %#ok<AGROW>
            last_abort_seq = local_current_abort(comm);
            preread = [];
            fprintf('Holo priming aborted; waiting for next prime.\n');
        elseif reprime
            % Stop the SLMs before re-priming: they are still armed with the
            % finished experiment's holograms, and the next prime will load new
            % ones. Nothing is in flight, so this is the clean handover point.
            for s = slm, try, s.stop(); catch, end, end   %#ok<AGROW>
            preread = [];
            fprintf('New prime detected; re-priming.\n');
        end
    end
end

function [chans, source] = local_prime_channels(prime, inventory_wl)
%LOCAL_PRIME_CHANNELS The ordered channel table for this prime, and where it came from.
%   Three sources, best first, and the caller PRINTS which one was used so an
%   operator can see whether the DAQ is speaking the channel protocol or the
%   listener is guessing:
%     1) prime.opto        -- the rig's opto table (name + wavelength per channel)
%     2) prime.wavelengths -- ordered wavelengths only; names synthesised. The
%                             current DAQ already sends this (prime_info derives
%                             it), so order/count disagreement is caught TODAY,
%                             before any DAQ-side change.
%     3) this listener's own inventory -- the pre-existing behaviour, used only
%                             when the prime carries neither. Least safe, so it
%                             is labelled as a guess.
    % slm_board is in the field set so every branch produces the same shape --
    % local_signature reads it unconditionally.
    chans = struct('name', {}, 'wavelength', {}, 'slm_board', {});
    source = '';

    if isstruct(prime) && isfield(prime, 'opto') && ~isempty(prime.opto)
        t = prime.opto;
        for i = 1:numel(t)
            chans(i).name = char(local_field(t(i), 'name', sprintf('ch%d', i)));
            chans(i).wavelength = double(local_field(t(i), 'wavelength', NaN));
            % Carried through so local_signature can report the board the DAQ
            % DECLARED, which is what the two sides must agree on. Without this the
            % board never reached the signature and it could only ever say 'auto'.
            chans(i).slm_board = local_field(t(i), 'slm_board', []);
        end
        source = 'prime.opto';
    elseif isstruct(prime) && isfield(prime, 'wavelengths') && ~isempty(prime.wavelengths)
        w = reshape(double(prime.wavelengths), 1, []);
        for i = 1:numel(w)
            chans(i).name = sprintf('%dnm', w(i));
            chans(i).wavelength = w(i);
            chans(i).slm_board = [];   % not declared: 'auto', same as the DAQ sends
        end
        source = 'prime.wavelengths';
    else
        w = reshape(double(inventory_wl), 1, []);
        for i = 1:numel(w)
            chans(i).name = sprintf('%dnm', w(i));
            chans(i).wavelength = w(i);
            chans(i).slm_board = [];
        end
        source = 'this listener''s inventory (prime carried no channel info)';
    end

    assert(~isempty(chans), 'holo_listener:noChannels', ...
        'Prime carried no opto channels and this listener has no wavelengths.');
    bad = ~isfinite([chans.wavelength]);
    assert(~any(bad), 'holo_listener:badWavelength', ...
        'Primed channel(s) %s have no usable wavelength.', ...
        strjoin({chans(bad).name}, ', '));
end

function slot = local_resolve_slots(chans, inventory_wl)
%LOCAL_RESOLVE_SLOTS Channel -> SLM index, BY WAVELENGTH.
%   The old code paired the i-th holoRequest with calib(i) and slm(i) of this
%   listener's own list, so a listener launched as [900 1100] against a DAQ
%   sending 1100 first compiled each hologram against the other wavelength's
%   affine and returned the other's diffraction efficiency -- which then set the
%   laser power. Resolving by wavelength makes position irrelevant.
    inventory_wl = reshape(double(inventory_wl), 1, []);
    slot = zeros(1, numel(chans));
    for k = 1:numel(chans)
        hit = find(inventory_wl == chans(k).wavelength, 1);
        assert(~isempty(hit), 'holo_listener:unopenedWavelength', ...
            ['The DAQ primed channel ''%s'' at %d nm, but this listener opened ' ...
             'SLMs for %s.\nRestart it including that wavelength:\n    ' ...
             'start_holo_listener(''Wavelengths'', %s)'], ...
            chans(k).name, chans(k).wavelength, mat2str(inventory_wl), ...
            mat2str(unique([inventory_wl, chans(k).wavelength], 'stable')));
        slot(k) = hit;
    end
    assert(numel(unique(slot)) == numel(slot), 'holo_listener:duplicateWavelength', ...
        ['Two primed channels resolve to the same SLM: %s. Channels sharing a ' ...
         'wavelength need\nseparate boards, which this listener cannot infer from ' ...
         'the wavelength alone.'], local_describe(chans));
end

function k = local_match_channel(HR, chans, got)
%LOCAL_MATCH_CHANNEL Which channel does this holoRequest belong to?
%   Prefers the tag the DAQ stamps on the request (name, else wavelength). An
%   UNTAGGED request falls back to the next unfilled slot, which is exactly the
%   old positional behaviour -- so a pre-opto DAQ still works unchanged. A tag
%   naming an unknown channel, or one already filled, is refused rather than
%   silently attributed to a neighbour.
    k = [];
    if isstruct(HR) && ~isempty(HR)
        if isfield(HR, 'channel') && ~isempty(HR.channel)
            k = find(strcmp({chans.name}, char(HR.channel)), 1);
            assert(~isempty(k), 'holo_listener:unknownChannel', ...
                'holoRequest is tagged channel ''%s'', which was not primed (%s).', ...
                char(HR.channel), local_describe(chans));
        elseif isfield(HR, 'wavelength') && ~isempty(HR.wavelength)
            k = find([chans.wavelength] == double(HR.wavelength), 1);
            assert(~isempty(k), 'holo_listener:unknownChannel', ...
                'holoRequest is tagged %d nm, which was not primed (%s).', ...
                double(HR.wavelength), local_describe(chans));
        end
        if ~isempty(k)
            assert(~got(k), 'holo_listener:duplicateChannel', ...
                ['Two holoRequests both claim channel ''%s''. The extra one would ' ...
                 'otherwise become\nthe next experiment''s first request -- a ' ...
                 'persistent cross-experiment off-by-one.'], chans(k).name);
            return
        end
    end
    k = find(~got, 1);   % untagged: next unfilled slot (legacy DAQ)
    assert(~isempty(k), 'holo_listener:tooManyRequests', ...
        'Received more holoRequests than the %d channel(s) primed.', numel(chans));
end

function [c, fname] = local_calib_for(wavelength, calib_dir)
%LOCAL_CALIB_FOR Newest calibration for a wavelength, plus its filename.
%   The filename is reported to the DAQ so an operator can see WHICH calibration
%   each channel got, rather than trusting that the newest file is the right one.
    c = find_latest_calib(wavelength, calib_dir);
    fname = '';
    try
        if isstruct(c) && isfield(c, 'file') && ~isempty(c.file)
            [~, n, e] = fileparts(char(c.file));
            fname = [n e];
        end
    catch
    end
end

function local_declare(comm, prime, ok, message, chans, source, cnames)
%LOCAL_DECLARE Publish the resolved channel table on its OWN topic.
%   Deliberately not holo_status: that means "primed and armed", and the
%   launcher's readiness lamp reads it. Declaring resolution there would turn the
%   lamp green before a single hologram had been compiled.
    d = struct('who', 'holo', 'ok', logical(ok), 'message', message, ...
               'protocol', 'opto/1', 'source', source);
    if isstruct(prime) && isfield(prime, 'prime_seq')
        d.prime_seq = prime.prime_seq;
    end
    d.stem = local_stem(prime);
    d.n_channels = numel(chans);
    names = {}; wls = [];
    for k = 1:numel(chans)
        names{end+1} = chans(k).name;      %#ok<AGROW>
        wls(end+1)   = chans(k).wavelength; %#ok<AGROW>
    end
    d.names = names;
    d.wavelengths = wls;
    d.signature = local_signature(chans);
    d.calibrations = cnames;
    try
        comm.set_config(d, 'holo_channels');
    catch
    end
end

function [slm_all, inv_wl, src] = local_slm_inventory(explicit_wl)
%LOCAL_SLM_INVENTORY Open this machine's SLMs, once, from the rig's opto table.
%   Three sources, best first:
%     1) an explicit 'Wavelengths' argument -- legacy get_slm mapping, for a
%        machine deliberately opening a subset;
%     2) the rig's opto table from config/rig -- each channel's slm_board and
%        slm_lut used verbatim when it declares them, else get_slm for that
%        wavelength;
%     3) the historical [1100 900], when the broker is cold and no rig is around.
%
%   Board/LUT must come from PERSISTENT config, not from a prime: the boards are
%   opened here, once, before any prime exists, because a Meadowlark board cannot
%   safely be reopened per experiment.
%
%   get_slm is retained (not renamed) as the per-wavelength fallback -- it is
%   still the right answer for a channel that declares no board, and it has four
%   other callers (single_slm_patch, get_psf_v2, calibrate_DE_powermeter,
%   MsocketHolorequest2K) that have no reason to change.
    if ~isempty(explicit_wl)
        inv_wl = reshape(double(explicit_wl), 1, []);
        slm_all = local_open_by_wavelength(inv_wl);
        src = 'explicit Wavelengths argument (legacy get_slm mapping)';
        return
    end

    t = [];
    try, t = rig_remote_get('opto'); catch, end

    if isempty(t)
        inv_wl = [1100 900];
        slm_all = local_open_by_wavelength(inv_wl);
        src = 'historical default [1100 900] -- nothing published, run publish_rig_config on the DAQ';
        return
    end

    lut_dir = local_cfg([], 'paths.slm_lut_dir', ...
        'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files');

    inv_wl  = zeros(1, numel(t));
    slm_all = [];
    pinned  = 0;
    for i = 1:numel(t)
        wl_i  = double(t(i).wavelength);
        inv_wl(i) = wl_i;
        board = t(i).slm_board;
        lut   = char(t(i).slm_lut);

        if isempty(board) || isempty(lut)
            % Channel declares no board/LUT: derive both from the wavelength, the
            % behaviour every rig had before this was configurable.
            slm_all = [slm_all, local_open_by_wavelength(wl_i)]; %#ok<AGROW>
            continue
        end
        if ~local_is_abs(lut)
            lut = fullfile(lut_dir, lut);
        end
        if exist(lut, 'file') ~= 2
            error('holo_listener:noLut', ...
                ['Channel ''%s'' (%d nm) declares LUT\n  %s\nwhich does not ' ...
                 'exist. Fix slm_lut for that opto channel in the rig file and\n' ...
                 're-run publish_rig_config() on the DAQ.'], ...
                char(t(i).name), wl_i, lut);
        end
        slm_all = [slm_all, MeadowlarkOneK(double(board), lut)]; %#ok<AGROW>
        pinned  = pinned + 1;
        fprintf('Loaded %dnm SLM (board %d, pinned by the rig).\n', wl_i, board);
    end

    src = sprintf('config/rig opto table (%d of %d channels pin a board)', ...
        pinned, numel(t));
end


function slm = local_open_by_wavelength(wls)
%LOCAL_OPEN_BY_WAVELENGTH The legacy get_slm switch, with a usable error.
    slm = [];
    for w = reshape(double(wls), 1, [])
        try
            s = get_slm(w);
        catch err
            error('holo_listener:slmOpenFailed', ...
                ['get_slm(%d) failed: %s\nIt only knows 900, 1100, 1030, 589 and ' ...
                 '607. For any other wavelength, declare slm_board and slm_lut on ' ...
                 'that\nopto channel in the rig file instead of relying on this ' ...
                 'mapping.'], w, err.message);
        end
        assert(~isempty(s) && isscalar(s), 'holo_listener:slmOpenFailed', ...
            ['get_slm(%d) returned nothing. It only knows 900, 1100, 1030, 589 ' ...
             'and 607.\nDeclare slm_board and slm_lut on that opto channel in the ' ...
             'rig file instead.'], w);
        slm = [slm, s]; %#ok<AGROW>
    end
end


function tf = local_is_abs(p)
%LOCAL_IS_ABS Windows-style absolute path test (drive letter or UNC).
    p = char(p);
    tf = ~isempty(regexp(p, '^([A-Za-z]:[\\/]|\\\\)', 'once'));
end


function local_assert_boards(chans, slm)
%LOCAL_ASSERT_BOARDS The opened board must be the one the DAQ declared.
%   Only checks channels that actually declare a board. A mismatch means the rig
%   file changed after this listener opened its boards, so the fix is a restart --
%   the alternative is driving the wrong physical SLM.
    for k = 1:numel(chans)
        want = chans(k).slm_board;
        if isempty(want)
            continue
        end
        got = slm(k).board_id;
        assert(double(want) == double(got), 'holo_listener:boardMismatch', ...
            ['Channel ''%s'' (%d nm) is declared on SLM board %d, but this ' ...
             'listener opened board %d\nfor that wavelength. The rig file changed ' ...
             'since startup; restart this listener.\nRefusing to run rather than ' ...
             'drive the wrong SLM.'], ...
            chans(k).name, chans(k).wavelength, double(want), double(got));
    end
end


function [v, src] = local_cfg(explicit, rig_path, fallback)
%LOCAL_CFG One config value: an explicit argument, else the rig, else a literal.
%   Returns the value and a short string naming where it came from, so the
%   listener's startup banner can say which config it is actually running on.
%   rig_remote_get does the config/rig -> local rig -> fallback resolution; this
%   only adds the "caller passed it explicitly" tier on top.
    if ~isempty(explicit)
        v = explicit;
        src = 'explicit argument';
        return
    end
    try
        [v, src] = rig_remote_get(rig_path, fallback);
    catch err
        % holodaq not on this machine's path, or the broker unreachable in a way
        % rig_remote_get did not absorb. Still start, on the old literal.
        v = fallback;
        src = sprintf('literal (rig lookup failed: %s)', err.message);
    end
end


function sig = local_signature(chans)
%LOCAL_SIGNATURE Must match holodaq's opto_signature for the same table.
%   Format: name@wavelength#board, joined by '|'.
%
%   The board reported is the one the DAQ DECLARED (the prime's slm_board), not
%   the one this machine happened to open. That distinction is the whole point:
%   the signature is an agreement about the DECLARATION, so both sides compute it
%   from the same input and get the same string. Whether the opened hardware
%   matches the declaration is a separate check (local_assert_boards).
%
%   This used to hardcode '#auto'. That was correct only while no rig pinned a
%   board -- the moment one did, opto_signature emitted '#<board>' here while this
%   side still said '#auto', the exact-match test in
%   Experiment.confirm_opto_agreement failed, and every run fell through to the
%   weak wavelength-only branch printing "the signatures differ by name only".
%   That silently downgraded the check whose entire job is to stop a beam being
%   steered by the wrong phase mask.
%
%   Kept duplicated rather than calling opto_signature because this file must run
%   on a machine that has no rig file; the two implementations have to stay in
%   step, so the 'auto' convention is spelled out identically in both.
    if isempty(chans)
        sig = 'none';
        return
    end
    parts = cell(1, numel(chans));
    for i = 1:numel(chans)
        board = 'auto';
        if isfield(chans, 'slm_board') && ~isempty(chans(i).slm_board)
            board = sprintf('%d', double(chans(i).slm_board));
        end
        parts{i} = sprintf('%s@%d#%s', chans(i).name, chans(i).wavelength, board);
    end
    sig = strjoin(parts, '|');
end

function s = local_describe(chans)
    if isempty(chans)
        s = '<no channels>';
        return
    end
    p = cell(1, numel(chans));
    for i = 1:numel(chans)
        p{i} = sprintf('%s@%dnm', chans(i).name, chans(i).wavelength);
    end
    s = strjoin(p, ', ');
end

function v = local_field(s, name, default)
    v = default;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    end
end

function seq = local_current_abort(comm)
    seq = -inf;
    try
        a = comm.scan_config('abort');
        if isstruct(a) && isfield(a, 'abort_seq') && ~isempty(a.abort_seq)
            seq = a.abort_seq;
        end
    catch
    end
end

function tf = local_is_new(prime, last_seq)
    tf = isstruct(prime) && isfield(prime, 'prime_seq') && ~isempty(prime.prime_seq) ...
        && prime.prime_seq > last_seq;
end

function s = local_stem(prime)
    if isstruct(prime) && isfield(prime, 'stem') && ~isempty(prime.stem)
        s = prime.stem;
    else
        s = '?';
    end
end

function local_ack(comm, prime, ok, message)
    a = struct('who', 'holo', 'ok', logical(ok), 'message', message);
    if isstruct(prime) && isfield(prime, 'prime_seq')
        a.prime_seq = prime.prime_seq;
    end
    a.stem = local_stem(prime);
    try
        comm.set_config(a, 'holo_status');
    catch
    end
end

function [preread, aborted, reprime] = local_serve(comm, slm, holograms, abort_baseline, serving_seq)
    % Poll msg/holo with a short timeout so the shared config/abort signal is
    % caught promptly. A firing-order CELL -> play it. A STRUCT is the next
    % experiment's holoRequest -> hand it back so the outer loop re-primes.
    %
    % ALSO watches config/holo for a NEWER prime_seq than the one being served.
    % Without that, the only ways out of here were an abort or a holoRequest
    % arriving on msg/holo -- so between experiments the listener was not
    % "waiting for the next prime" at all, it was waiting for a sequence. The DAQ
    % writes the next prime to config/holo and then, since this listener had not
    % declared its channels for it, sat in confirm_opto_agreement until
    % HoloAckTimeout (120 s) before warning and pressing on. A circular wait: the
    % DAQ would not send a holoRequest until the listener declared, and the
    % listener would not look at the prime until a holoRequest arrived. Every
    % experiment after the first paid two minutes and a spurious warning.
    preread = []; aborted = false; reprime = false;
    while true
        a = [];
        try, a = comm.scan_config('abort'); catch, end
        if isstruct(a) && isfield(a, 'abort_seq') && ~isempty(a.abort_seq) ...
                && a.abort_seq > abort_baseline
            aborted = true;
            return
        end

        % A new prime means the DAQ has moved on to the next experiment.
        c = [];
        try, c = comm.scan_config('holo'); catch, end
        if isstruct(c) && isfield(c, 'prime_seq') && ~isempty(c.prime_seq) ...
                && c.prime_seq > serving_seq
            reprime = true;
            return
        end

        msg = comm.read(2);          % short timeout: re-check abort between waits
        if isempty(msg)
            continue                 % idle; keep waiting (no crash, unlike ShootSequences)
        end
        if isstruct(msg)
            preread = msg;
            return
        end
        PlaySequence2K(slm, holograms, msg);
    end
end
