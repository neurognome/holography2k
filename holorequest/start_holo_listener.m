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
%   Options:
%     'Wavelengths' (default [1100 900]) -- must match the order the DAQ
%         transfers holoRequests (FullExperiment.initialize: hr1100 then hr900).
%         Changing the wavelength SET requires restarting this listener.
%     'CalibDir'    (default C:\Users\holos\Documents\calibs)
%     'Timeout'     (default 1700 ms) -- SLM trigger timeout.
%
%   How re-prime stays race-free: the serve loop reads msg/holo and switches on
%   type -- a firing-order CELL is played; a holoRequest STRUCT means the next
%   experiment has begun (its config prime is already posted), so it is handed
%   back to be compiled. Sequences and holoRequests therefore never collide.
%
%   See also: MsocketHolorequest2K (fallback), PlaySequence2K, find_latest_calib,
%             generate_holograms_new, prime_info.

    p = inputParser;
    p.addParameter('Wavelengths', [1100 900]);
    p.addParameter('CalibDir', 'C:\Users\holos\Documents\calibs');
    p.addParameter('Timeout', 1700);
    p.parse(varargin{:});
    opt = p.Results;
    wl  = opt.Wavelengths;

    addpath(genpath('C:\Users\holos\Documents\GitHub\holography2k'))
    addpath(genpath('C:\Users\holos\Desktop\meadowlark'))

    comm = HolochatInterface('holo');
    comm.flush();   % drop any stale msg/holo from a previous session

    Setup = function_loadparameters2();
    Setup.CGHMethod          = 2;      % GSS
    Setup.verbose            = 0;
    Setup.useGPU             = 1;
    Setup.TimeToPickSequence = 0.05;
    Setup.SLM.timeout_ms     = opt.Timeout;

    % SLM objects are created ONCE, because a Meadowlark board cannot safely be
    % reopened per experiment. So 'Wavelengths' is THIS MACHINE'S SLM INVENTORY,
    % not the experiment's channel set: a prime's channels are resolved against
    % it by wavelength (local_resolve_slots), never by position.
    slm_all = [];
    for w = wl
        slm_all = [slm_all, get_slm(w)]; %#ok<AGROW>
    end
    assert(numel(slm_all) == numel(wl), 'holo_listener:slmOpenFailed', ...
        ['Opened %d SLM(s) for %d wavelength(s). get_slm returns nothing for a ' ...
         'wavelength it\ndoes not know (it handles 900, 1100, 1030, 589, 607), so ' ...
         'check the set: %s.'], numel(slm_all), numel(wl), mat2str(wl));

    % Two wavelengths can map to ONE physical board -- get_slm sends both 1100 and
    % 1030 to board 2. That would drive one SLM with two hologram stacks, each
    % overwriting the other, discovered only as garbage on the sample. Refuse at
    % startup instead.
    boards = [slm_all.board_id];
    assert(numel(unique(boards)) == numel(boards), 'holo_listener:sharedBoard', ...
        ['Wavelengths %s resolve to boards %s -- at least two share one board.\n' ...
         'get_slm maps both 1100 and 1030 to board 2, so that set would drive one ' ...
         'physical\nSLM with two hologram stacks. Restart with a set that maps 1:1.'], ...
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
            calib = [];
            cname = cell(1, numel(chans));
            for k = 1:numel(chans)
                [c, f] = local_calib_for(chans(k).wavelength, opt.CalibDir);
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
                s.timeout_ms = opt.Timeout;
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

        % 4) serve until the next experiment's holoRequest arrives (or abort).
        [preread, aborted] = local_serve(comm, slm, holograms, last_abort_seq);
        if aborted
            for s = slm, try, s.stop(); catch, end, end   %#ok<AGROW>
            last_abort_seq = local_current_abort(comm);
            preread = [];
            fprintf('Holo priming aborted; waiting for next prime.\n');
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
    chans = struct('name', {}, 'wavelength', {});
    source = '';

    if isstruct(prime) && isfield(prime, 'opto') && ~isempty(prime.opto)
        t = prime.opto;
        for i = 1:numel(t)
            chans(i).name = char(local_field(t(i), 'name', sprintf('ch%d', i)));
            chans(i).wavelength = double(local_field(t(i), 'wavelength', NaN));
        end
        source = 'prime.opto';
    elseif isstruct(prime) && isfield(prime, 'wavelengths') && ~isempty(prime.wavelengths)
        w = reshape(double(prime.wavelengths), 1, []);
        for i = 1:numel(w)
            chans(i).name = sprintf('%dnm', w(i));
            chans(i).wavelength = w(i);
        end
        source = 'prime.wavelengths';
    else
        w = reshape(double(inventory_wl), 1, []);
        for i = 1:numel(w)
            chans(i).name = sprintf('%dnm', w(i));
            chans(i).wavelength = w(i);
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

function sig = local_signature(chans)
%LOCAL_SIGNATURE Must match holodaq's opto_signature for the same table.
%   Duplicated rather than shared because this file runs on the holography
%   computer, which does not have holodaq's rigs/ on its path. Format:
%   name@wavelength#board, joined by '|'. Board is 'auto' here: this side
%   resolves the board from the wavelength, so it has nothing else to report.
    if isempty(chans)
        sig = 'none';
        return
    end
    parts = cell(1, numel(chans));
    for i = 1:numel(chans)
        parts{i} = sprintf('%s@%d#auto', chans(i).name, chans(i).wavelength);
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

function [preread, aborted] = local_serve(comm, slm, holograms, abort_baseline)
    % Poll msg/holo with a short timeout so the shared config/abort signal is
    % caught promptly. A firing-order CELL -> play it. A STRUCT is the next
    % experiment's holoRequest -> hand it back so the outer loop re-primes.
    preread = []; aborted = false;
    while true
        a = [];
        try, a = comm.scan_config('abort'); catch, end
        if isstruct(a) && isfield(a, 'abort_seq') && ~isempty(a.abort_seq) ...
                && a.abort_seq > abort_baseline
            aborted = true;
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
