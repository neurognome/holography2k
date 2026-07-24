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

    Setup = function_loadparameters2();
    Setup.CGHMethod          = 2;      % GSS
    Setup.verbose            = 0;
    Setup.useGPU             = 1;
    Setup.TimeToPickSequence = 0.05;
    Setup.SLM.timeout_ms     = opt.Timeout;

    % SLM objects are created ONCE for the fixed wavelength set (as in the
    % original script); only the compiled holograms change per experiment.
    slm = [];
    for w = wl
        slm = [slm, get_slm(w)]; %#ok<AGROW>
    end

    last_seq = -inf;
    preread  = [];
    last_abort_seq = local_current_abort(comm);   % ignore any stale abort at startup
    fprintf('Holo listener up; waiting for a prime from the DAQ...\n');
    while true
        prime = comm.get_config();
        if ~local_is_new(prime, last_seq)
            pause(0.25);
            continue
        end
        last_seq = prime.prime_seq;
        stem = local_stem(prime);
        fprintf('Priming holo for %s (wavelengths %s)...\n', stem, mat2str(wl));

        try
            % 1) newest calibration per wavelength.
            calib = [];
            for w = wl
                calib = [calib, find_latest_calib(w, opt.CalibDir)]; %#ok<AGROW>
            end

            % 2) compile holograms (blocks on the DAQ transferHR per wavelength,
            %    in the same order the DAQ sends them). wl(1)'s holoRequest may
            %    already be in hand from the serve loop (preread).
            holograms = cell(1, numel(wl));
            for wi = 1:numel(wl)
                if wi == 1 && ~isempty(preread)
                    hololist = generate_holograms_new(comm, Setup, calib(wi), preread);
                    preread  = [];
                else
                    hololist = generate_holograms_new(comm, Setup, calib(wi));
                end
                holograms{wi} = uint8(hololist);
            end
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
