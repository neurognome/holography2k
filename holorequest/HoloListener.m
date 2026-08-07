classdef HoloListener < handle
    %HOLOLISTENER The config-primed holography listener, as a tickable state machine.
    %   This is start_holo_listener's loop with the `while true` taken out. The
    %   behaviour is unchanged -- wait for a prime, resolve channels, compile
    %   holograms, arm the SLM(s), serve firing orders, re-prime -- but it is now
    %   driven one step at a time, which is what makes an async listener (and so a
    %   Stop button and a status lamp) possible at all.
    %
    %   THREE PHASES, one per tick:
    %       waiting     polling config/holo for a NEW prime_seq
    %       compiling   resolving channels, loading calibrations, compiling
    %                   holograms, arming the SLMs. ONE tick, and a long one --
    %                   generate_holograms_new blocks until the DAQ's holoRequests
    %                   arrive, with no timeout. The prompt and the status window
    %                   are frozen for its whole duration; this is inherent, not a
    %                   defect, and the phase is reported so it does not look hung.
    %       serving     playing firing orders until an abort, a new prime, or the
    %                   next experiment's holoRequest arrives
    %
    %   Usage (on the holo computer -- normally via start_holo_listener):
    %       l = HoloListener();          % opens the SLMs, resolves rig config
    %       l.listen_async();            % timer-driven; the prompt stays usable
    %       l.status()                   % am I actually listening?
    %       l.stop()                     % halt
    %       l.listen();                  % blocking, the pre-2026-08 behaviour
    %
    %   TIMING. listen_async polls every poll_period (0.05 s) and, while serving,
    %   spends serve_read (0.4 s) of each tick inside a single blocking read. So a
    %   firing order is picked up at most poll_period late rather than immediately,
    %   which is the one behavioural cost of going async. That 50 ms sits inside the
    %   pause(0.1) HolochatInterface.send already performs on the DAQ AFTER posting
    %   the sequence and BEFORE the trial can start, so it cannot make a sequence
    %   land after its trigger. listen() (blocking) does not pause while serving, so
    %   it keeps the old continuous read exactly.
    %
    %   See also start_holo_listener, HoloListenerMonitor, PlaySequence2K,
    %   generate_holograms_new, find_latest_calib, Receiver (holodaq, the same
    %   listen/listen_async/stop/status shape for the SI and PTB boxes).

    properties
        name = 'holo'
        comm                 % HolochatInterface('holo'); [] when Offline
        Setup                % function_loadparameters2 struct, CGH settings applied
        slm_all              % THIS MACHINE'S SLM inventory, opened once at startup
        wl                   % inventory wavelengths, aligned with slm_all
        calib_dir
        timeout_ms

        % Tick spacing for listen_async. Small, because it bounds how late a firing
        % order can be picked up -- see the TIMING note above.
        poll_period = 0.05
        % How long each serving tick spends in one blocking msg read. Sized so the
        % listener is reading ~90% of the time while serving (0.4 / (0.4 + 0.05)),
        % which keeps the broker load and the responsiveness of the old read(2).
        serve_read = 0.4
        % Config topics (prime / abort) are polled at most this often, whatever the
        % tick rate. The old loop's pause(0.25); one webread per topic per 0.25 s is
        % all a config channel needs, and hammering it is what scan_config exists to
        % avoid.
        config_period = 0.25

        poll = []            % async timer ([] unless listen_async is running)

        last_seq = -inf
        last_abort_seq = -inf
        preread = []
        prime = []
        chans = struct('name', {}, 'wavelength', {}, 'slm_board', {})
        source = ''
        slm = []             % the subset of slm_all this prime resolved to
        holograms = {}
        sequences = 0        % firing orders played since this prime

        last_ok = true
        last_message = 'not started'

        % Called with no arguments whenever the phase changes. HoloListenerMonitor
        % sets it so the lamp repaints BEFORE a long blocking phase begins rather
        % than after it ends -- a window that only says COMPILING once compiling is
        % over is worse than no window.
        on_phase = []
    end

    properties (SetAccess = private)
        % 'waiting' | 'compiling' | 'serving'. Read freely; write only through
        % set_phase, which is what fires on_phase -- a phase changed behind the
        % hook's back is a lamp showing the previous phase, which is the one thing
        % this window exists not to do.
        phase = 'waiting'
    end

    properties (Access = private)
        last_config_scan = []     % tic id; [] means "scan on the next tick"
        consecutive_errors = 0
        ERROR_REPORT_EVERY = 200  % ~10 s at the default poll_period
    end

    properties (Constant)
        % Named so an orphaned poll timer can be found and killed by name, the way
        % holodaq's Receiver timers can.
        TIMER_NAME = 'holoListenerPoll'
    end

    methods
        function obj = HoloListener(varargin)
            %HOLOLISTENER Resolve rig config and open this machine's SLMs, once.
            %   Options are start_holo_listener's, plus 'Offline'.
            p = inputParser;
            % [] means "take the inventory from the rig's opto table"
            % (local_slm_inventory). An explicit list still wins, for a machine
            % deliberately opening a subset.
            p.addParameter('Wavelengths', []);
            p.addParameter('CalibDir', '');
            p.addParameter('Timeout', []);
            % Offline: build the object without opening SLM hardware or a broker
            % connection. It never primes anything -- it exists so the status
            % window, the phase machine and the argument handling can be exercised
            % off the rig, where there is neither a Meadowlark board nor a broker.
            p.addParameter('Offline', false, @(v) islogical(v) || isnumeric(v));
            p.parse(varargin{:});
            opt = p.Results;

            % This checkout, wherever it happens to live. holodaq is NOT added here
            % -- it is already on this machine's path (HolochatInterface below has
            % always been resolved before these addpaths ran).
            addpath(genpath(fileparts(fileparts(mfilename('fullpath')))));

            if logical(opt.Offline)
                obj.comm = [];
                obj.slm_all = [];
                obj.wl = [];
                obj.calib_dir = '';
                obj.timeout_ms = NaN;
                obj.last_message = 'offline (no SLMs, no broker)';
                return
            end

            obj.comm = HolochatInterface('holo');
            obj.comm.flush();   % drop any stale msg/holo from a previous session

            % ---- config ------------------------------------------------------
            % This machine has no rig file, so everything comes from what the DAQ
            % published to holochat (rig_remote_get: config/rig -> local rig ->
            % literal). Reported explicitly, because the previous behaviour was a
            % silent fallback to whatever was hardcoded here.
            [obj.calib_dir, calib_src] = local_cfg(opt.CalibDir, 'paths.calib_dir', ...
                'C:\Users\holos\Documents\calibs');
            % The Meadowlark SDK is genuinely per-machine (an install location, not
            % part of any checkout), so it stays a rig value rather than being derived.
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
            [cgh_method, cgh_src]   = local_cfg([], 'holo.cgh_method', 2);
            [use_gpu, gpu_src]      = local_cfg([], 'holo.use_gpu', true);
            [obj.timeout_ms, to_src] = local_cfg(opt.Timeout, 'holo.slm_timeout_ms', 1700);

            fprintf('--- holo listener config ---\n');
            fprintf('  calib dir   : %s   (%s)\n', obj.calib_dir, calib_src);
            fprintf('  slm sdk     : %s   (%s)\n', slm_sdk, sdk_src);
            fprintf('  cgh method  : %d   (%s)\n', cgh_method, cgh_src);
            fprintf('  use gpu     : %d   (%s)\n', use_gpu, gpu_src);
            fprintf('  slm timeout : %d ms   (%s)\n', obj.timeout_ms, to_src);
            if ~isfolder(obj.calib_dir)
                warning('holo_listener:noCalibDir', ...
                    ['The calibration folder does not exist:\n  %s\nfind_latest_calib ' ...
                     'will fail when a prime arrives. Fix paths.calib_dir in the rig file ' ...
                     'and\nre-run publish_rig_config() on the DAQ, then restart me.'], ...
                    obj.calib_dir);
            end

            obj.Setup = function_loadparameters2();
            obj.Setup.CGHMethod          = cgh_method;
            obj.Setup.verbose            = 0;
            obj.Setup.useGPU             = double(logical(use_gpu));
            obj.Setup.TimeToPickSequence = 0.05;
            obj.Setup.SLM.timeout_ms     = obj.timeout_ms;

            % SLM objects are created ONCE, because a Meadowlark board cannot safely
            % be reopened per experiment. So this is THIS MACHINE'S SLM INVENTORY,
            % not the experiment's channel set: a prime's channels are resolved
            % against it by wavelength (local_resolve_slots), never by position. It
            % is also exactly why the inventory comes from PERSISTENT rig config
            % rather than from the prime -- it has to exist before the first prime
            % ever arrives.
            [obj.slm_all, obj.wl, inv_src] = local_slm_inventory(opt.Wavelengths);
            fprintf('  slm source  : %s\n', inv_src);

            % Two wavelengths can map to ONE physical board -- get_slm sends both
            % 1100 and 1030 to board 2. That would drive one SLM with two hologram
            % stacks, each overwriting the other, discovered only as garbage on the
            % sample. Refuse at startup instead. This is the check that makes pinning
            % slm_board in the rig file safe, so it applies to a rig-supplied
            % inventory too.
            boards = [obj.slm_all.board_id];
            assert(numel(unique(boards)) == numel(boards), 'holo_listener:sharedBoard', ...
                ['Wavelengths %s resolve to boards %s -- at least two share one board.\n' ...
                 'That would drive one physical SLM with two hologram stacks. Give each ' ...
                 'arm its own\nslm_board in the rig file (see ExampleRig''s split-arm ' ...
                 'example), or restart with a\nwavelength set that maps 1:1.'], ...
                mat2str(obj.wl), mat2str(boards));
            fprintf('SLM inventory: %s nm on boards %s.\n', mat2str(obj.wl), mat2str(boards));
            obj.last_message = 'waiting for a prime';
        end

        function listen(obj)
            %LISTEN Blocking listener. Owns this MATLAB session until Ctrl-C.
            %   The pre-2026-08 behaviour of start_holo_listener, and still the
            %   lowest-latency one: while serving it does NOT pause between reads,
            %   so firing-order delivery is continuous exactly as it always was.
            %   No status window is possible here -- nothing can click a button in a
            %   loop that never yields.
            obj.flush_stale();
            fprintf(['[holo] listener up (BLOCKING; Ctrl-C to stop); ' ...
                     'waiting for a prime from the DAQ...\n']);
            while true
                obj.tick();
                if ~strcmp(obj.phase, 'serving')
                    % Serving already spends serve_read inside the read; pausing on
                    % top of that would open a gap the old loop never had.
                    pause(obj.poll_period);
                end
            end
        end

        function t = listen_async(obj)
            %LISTEN_ASYNC The same listener on a timer, leaving the prompt free.
            %   t = obj.listen_async()   % also stored in obj.poll
            %
            %   TWO THINGS TO KNOW, both inherent to timers rather than to this code:
            %
            %   1. The event queue is serviced only when MATLAB is idle or hits
            %      pause/drawnow/waitfor, so a long blocking operation at the prompt
            %      DELAYS ticks. status() reports the achieved period so you can see
            %      when that is happening.
            %   2. The COMPILING phase is one long blocking tick (minutes are normal
            %      while the DAQ sends holoRequests). Nothing else in this session
            %      runs during it, including the Stop button. That is the one place
            %      the window genuinely freezes, and it says COMPILING while it does.
            obj.stop();          % never run two poll timers
            obj.flush_stale();
            obj.poll = timer('Name', HoloListener.TIMER_NAME, ...
                'ExecutionMode', 'fixedSpacing', ...  % space AFTER the callback: a slow
                'BusyMode', 'drop', ...               % tick must not build a backlog
                'Period', obj.poll_period, ...
                'TimerFcn', @(~, ~) obj.tick(), ...
                'ErrorFcn', @(src, evt) obj.on_timer_error(src, evt));
            start(obj.poll);
            fprintf(['[holo] listener up (ASYNC, every %.2f s; l.stop() to halt, ' ...
                     'l.status() to check); waiting for a prime...\n'], obj.poll_period);
            t = obj.poll;
        end

        function stop(obj)
            %STOP Halt the async listener. No-op if it is not running.
            %   Deliberately does NOT stop the SLMs: a listener stopped between
            %   experiments has nothing armed, and a listener stopped mid-experiment
            %   should not blank a pattern the DAQ is still gating a laser against.
            %   Use l.stop_slms() for that.
            if ~isempty(obj.poll) && isvalid(obj.poll)
                try, stop(obj.poll);   catch, end
                try, delete(obj.poll); catch, end
                fprintf('[holo] async listener stopped.\n');
            end
            obj.poll = [];
        end

        function tf = is_listening(obj)
            tf = ~isempty(obj.poll) && isvalid(obj.poll) && strcmp(obj.poll.Running, 'on');
        end

        function s = status(obj)
            %STATUS Am I actually listening? One call instead of a timerfindall hunt.
            s = struct('mode', 'not running', 'running', false, 'phase', obj.phase, ...
                       'period', obj.poll_period, 'expected_period', obj.expected_period(), ...
                       'achieved_period', NaN, 'last_prime_seq', obj.last_seq, ...
                       'stem', obj.stem(), 'sequences', obj.sequences, ...
                       'ok', obj.last_ok, 'message', obj.last_message);
            if ~isempty(obj.poll) && isvalid(obj.poll)
                s.mode = 'async';
                s.running = strcmp(obj.poll.Running, 'on');
                s.achieved_period = obj.poll.AveragePeriod;
            end
            fprintf('[holo] %s', s.mode);
            if s.running
                fprintf(', running, phase %s', s.phase);
            elseif strcmp(s.mode, 'async')
                fprintf(', STOPPED (timer exists but is not running -- see any error above)');
            end
            fprintf('; last prime seq %g, stem %s, %d sequence(s) played\n', ...
                s.last_prime_seq, s.stem, s.sequences);
            fprintf('       %s: %s\n', local_okword(s.ok), s.message);
        end

        function e = expected_period(obj)
            %EXPECTED_PERIOD What one tick SHOULD cost in the current phase.
            %   The starvation test compares against this, not against poll_period:
            %   a serving tick deliberately blocks for serve_read, so measuring it
            %   against the bare timer period would report STARVED on every healthy
            %   experiment. NaN while compiling -- a tick that blocks for minutes is
            %   the phase working as designed, not a starved one.
            switch obj.phase
                case 'serving'
                    e = obj.poll_period + obj.serve_read;
                case 'compiling'
                    e = NaN;
                otherwise
                    % The throttled config webread, amortised over the ticks between
                    % scans, plus room for ordinary timer jitter.
                    e = obj.poll_period + 0.02;
            end
        end

        function tick(obj)
            %TICK One step of the listener: the body both drivers call.
            %   Never throws. Under a timer an escaping error STOPS the timer, which
            %   would leave a listener that looks up but is dead; in the blocking
            %   loop it would kill the loop outright. Either way a broker hiccup, or
            %   a malformed prime, must not end the session.
            %
            %   This is a real change from the old loop, where local_prime_channels'
            %   asserts ran OUTSIDE the per-prime try/catch: a prime naming an
            %   unopened wavelength threw straight out of `while true` and killed the
            %   listener. Now it is reported, acked as a failure so the DAQ is not
            %   left waiting out its 300 s transfer timeout, and the listener stays up.
            try
                switch obj.phase
                    case 'waiting',   obj.step_waiting();
                    case 'compiling', obj.step_compiling();
                    case 'serving',   obj.step_serving();
                end
                obj.note_ok();
            catch err
                obj.note_error(err);
            end
        end

        function set_phase(obj, p)
            %SET_PHASE Move to a new phase and tell whoever is watching.
            %   The ONLY way phase changes, so on_phase can be relied on to fire.
            %   Public because a test has to be able to drive the phase machine
            %   without SLM hardware; on the rig only the step_* methods call it.
            assert(any(strcmp(p, {'waiting', 'compiling', 'serving'})), ...
                'HoloListener:badPhase', ...
                'Phase must be waiting / compiling / serving, got ''%s''.', char(p));
            if strcmp(obj.phase, p), return; end
            obj.phase = p;
            % The new phase polls immediately rather than waiting out the previous
            % phase's throttle window -- otherwise a re-prime detected while serving
            % sits idle for config_period before anyone looks at it.
            obj.last_config_scan = [];
            if ~isempty(obj.on_phase)
                % A status window must never be able to cost you the listener: this
                % runs inside tick(), where an escaping error would stop the timer.
                try, obj.on_phase(); catch, end
            end
        end

        function stop_slms(obj, which)
            %STOP_SLMS Stop the SLM(s), swallowing per-board failures.
            if nargin < 2, which = obj.slm_all; end
            for s = which
                try, s.stop(); catch, end
            end
        end

        function flush_stale(obj)
            %FLUSH_STALE Adopt the current prime/abort as baselines, so a message
            %   left over from a previous session does not fire on startup.
            %
            %   Called by BOTH listen() and listen_async(), so restarting from the
            %   Stop/Start button re-baselines: a prime posted while the listener was
            %   stopped is not picked up when it comes back. That is holodaq
            %   Receiver's rule too, and the same reason -- a prime the DAQ sent to
            %   a stopped box has already had its ack time out.
            if isempty(obj.comm), return; end
            try
                c0 = obj.comm.scan_config('holo');
                if isstruct(c0) && isfield(c0, 'prime_seq') && ~isempty(c0.prime_seq)
                    obj.last_seq = c0.prime_seq;
                end
            catch
            end
            obj.last_abort_seq = local_current_abort(obj.comm);
            obj.set_phase('waiting');
        end

        function s = stem(obj)
            s = local_stem(obj.prime);
        end

        function on_timer_error(obj, src, ~)
            % ErrorFcn: a timer whose callback throws STOPS. tick() already catches
            % everything, so reaching here means something escaped it -- and the
            % listener is now dead while still looking like it was started.
            %
            % DELETE the timer, do not just drop the handle: MATLAB stops an errored
            % timer but leaves the object alive, so clearing obj.poll on its own
            % orphans it -- still in timerfindall, invisible to status(), and a later
            % listen_async() cannot stop what it can no longer see.
            for h = {src, obj.poll}
                t = h{1};
                if ~isempty(t) && isa(t, 'timer') && isvalid(t)
                    try, stop(t);   catch, end
                    try, delete(t); catch, end
                end
            end
            obj.poll = [];
            obj.set_message(false, 'poll timer errored -- listener is dead');
            warning('HoloListener:asyncStopped', ...
                ['[holo] ASYNC LISTENER STOPPED -- the poll timer errored, so this box ' ...
                 'is no longer\nbeing primed. Restart it with l.listen_async(). ' ...
                 '(status() reports "not running" from here on.)']);
        end
    end

    methods (Static)
        function closeStaleTimers()
            %CLOSESTALETIMERS Delete orphaned poll timers left by a dead listener.
            ts = timerfindall('Name', HoloListener.TIMER_NAME);
            for i = 1:numel(ts)
                t = ts(i);
                if isvalid(t)
                    try, stop(t);   catch, end
                    try, delete(t); catch, end
                end
            end
        end
    end

    methods (Access = private)
        function set_message(obj, ok, msg)
            obj.last_ok = logical(ok);
            obj.last_message = msg;
        end

        function tf = due_for_config(obj)
            %DUE_FOR_CONFIG Rate-limit config polling to config_period.
            %   The tick rate is set by how promptly a firing order must be picked
            %   up; the config channels do not need it and should not pay for it.
            if isempty(obj.last_config_scan) || toc(obj.last_config_scan) >= obj.config_period
                obj.last_config_scan = tic;
                tf = true;
            else
                tf = false;
            end
        end

        function step_waiting(obj)
            if isempty(obj.comm) || ~obj.due_for_config(), return; end

            p = [];
            try
                p = obj.comm.scan_config('holo');
            catch err
                obj.note_error(err);
                return
            end
            if ~local_is_new(p, obj.last_seq), return; end
            obj.last_seq = p.prime_seq;
            obj.prime = p;
            obj.sequences = 0;

            try
                [c, src] = local_prime_channels(p, obj.wl);
            catch err
                % Was fatal: these asserts sat outside the old loop's try/catch.
                fprintf('Holo prime REFUSED: %s\n', err.message);
                obj.chans = struct('name', {}, 'wavelength', {}, 'slm_board', {});
                obj.source = 'refused';
                local_ack(obj.comm, p, false, err.message);
                obj.set_message(false, err.message);
                return
            end
            obj.chans = c;
            obj.source = src;

            % A prime that speaks the channel protocol and declares ZERO channels is
            % a vis-only experiment. There is nothing to compile, and it must NOT be
            % mistaken for "the prime carried no channel information".
            if isempty(c)
                fprintf('Prime %s declares no opto channels (%s) -- nothing to compile.\n', ...
                    local_stem(p), src);
                obj.set_message(true, sprintf('%s: no opto channels (%s)', local_stem(p), src));
                return
            end

            fprintf('Priming holo for %s (%s from %s)...\n', ...
                local_stem(p), local_describe(c), src);
            obj.set_message(true, sprintf('%s: compiling %s', local_stem(p), local_describe(c)));
            % Paint the lamp BEFORE handing the session to the compile, which does
            % not yield: the phase change is the only warning the operator gets.
            obj.set_phase('compiling');
        end

        function step_compiling(obj)
            if isempty(obj.comm)
                % Offline: there is nothing to compile against and no DAQ to tell.
                % Say so and fall back, rather than failing deep inside the compile
                % with a dot-index error on an empty broker handle.
                obj.set_message(false, 'offline: nothing to compile');
                obj.set_phase('waiting');
                return
            end
            p = obj.prime;
            c = obj.chans;
            try
                % 1) map each channel to an SLM slot BY WAVELENGTH, and load that
                %    channel's newest calibration. Position is never assumed.
                slot = local_resolve_slots(c, obj.wl);
                active = obj.slm_all(slot);
                % The boards were opened at startup; confirm they are still the ones
                % the DAQ declares, so a rig-file edit mid-session cannot quietly
                % point a channel at the wrong physical SLM.
                local_assert_boards(c, active);
                calib = [];
                cname = cell(1, numel(c));
                for k = 1:numel(c)
                    [cc, f] = local_calib_for(c(k).wavelength, obj.calib_dir);
                    calib = [calib, cc]; %#ok<AGROW>
                    cname{k} = f;
                end

                % Tell the DAQ what we resolved, BEFORE compiling. Its own topic, so
                % holo_status keeps meaning "primed and armed" and the launcher's
                % readiness lamp cannot go green early.
                local_declare(obj.comm, p, true, 'channels resolved', c, obj.source, cname);

                % 2) compile holograms. Each arriving holoRequest is matched to its
                %    channel by the tag transferHR stamps on it; an untagged request
                %    (a pre-opto DAQ) falls back to the next unfilled slot, which
                %    reproduces the old positional behaviour exactly.
                holos = cell(1, numel(c));
                got = false(1, numel(c));
                for n = 1:numel(c)
                    HR = obj.preread;
                    obj.preread = [];
                    if isempty(HR)
                        HR = [];   % generate_holograms_new reads it itself
                    end
                    k = local_match_channel(HR, c, got);
                    if isempty(HR)
                        hololist = generate_holograms_new(obj.comm, obj.Setup, calib(k));
                    else
                        hololist = generate_holograms_new(obj.comm, obj.Setup, calib(k), HR);
                    end
                    holos{k} = uint8(hololist);
                    got(k) = true;
                end
                assert(all(got), 'holo_listener:missingChannel', ...
                    'No holoRequest arrived for channel(s): %s.', ...
                    strjoin({c(~got).name}, ', '));
                obj.comm.flush();

                % 3) (re)arm the SLM(s) for triggered playback.
                for s = active
                    s.stop();
                    s.wait_for_trigger = 1;
                    s.timeout_ms = obj.timeout_ms;
                    s.start();
                end
                obj.slm = active;
                obj.holograms = holos;

                fprintf('Holo primed for %s; serving sequences.\n', local_stem(p));
                local_ack(obj.comm, p, true, 'holograms compiled');
                obj.set_message(true, sprintf('%s primed (%s); serving', ...
                    local_stem(p), local_describe(c)));
                obj.set_phase('serving');
            catch err
                fprintf('Holo prime FAILED: %s\n', err.message);
                % STOP THE SLMs FIRST. The arm loop is the last step above, so a
                % throw anywhere before it used to leave the boards armed with the
                % PREVIOUS experiment's holograms -- while the DAQ, having been told
                % nothing, went on to gate the laser against them.
                obj.stop_slms(obj.slm_all);
                local_declare(obj.comm, p, false, err.message, c, obj.source, {});
                local_ack(obj.comm, p, false, err.message);
                obj.preread = [];
                obj.slm = [];
                obj.holograms = {};
                obj.set_message(false, err.message);
                obj.set_phase('waiting');
            end
        end

        function step_serving(obj)
            if isempty(obj.comm), return; end   % Offline: nothing to serve from
            % Abort and re-prime are read from config, so they are rate-limited;
            % firing orders are read from msg on every tick.
            if obj.due_for_config()
                a = [];
                try, a = obj.comm.scan_config('abort'); catch, end
                if isstruct(a) && isfield(a, 'abort_seq') && ~isempty(a.abort_seq) ...
                        && a.abort_seq > obj.last_abort_seq
                    obj.stop_slms(obj.slm);
                    obj.last_abort_seq = local_current_abort(obj.comm);
                    obj.preread = [];
                    fprintf('Holo priming aborted; waiting for next prime.\n');
                    obj.set_message(true, 'aborted by the DAQ; waiting for the next prime');
                    obj.set_phase('waiting');
                    return
                end

                % A new prime means the DAQ has moved on to the next experiment.
                % Without this check the only ways out of serving were an abort or a
                % holoRequest, so between experiments the listener was not "waiting
                % for the next prime" at all -- it was waiting for a sequence, while
                % the DAQ sat in confirm_opto_agreement burning its HoloAckTimeout.
                c = [];
                try, c = obj.comm.scan_config('holo'); catch, end
                if isstruct(c) && isfield(c, 'prime_seq') && ~isempty(c.prime_seq) ...
                        && c.prime_seq > obj.last_seq
                    % Stop the SLMs before re-priming: they are still armed with the
                    % finished experiment's holograms, and the next prime will load
                    % new ones. Nothing is in flight, so this is the clean handover.
                    obj.stop_slms(obj.slm);
                    obj.preread = [];
                    fprintf('New prime detected; re-priming.\n');
                    obj.set_phase('waiting');   % step_waiting picks it up next tick
                    return
                end
            end

            msg = obj.comm.read(obj.serve_read);
            if isempty(msg)
                return                       % idle; keep serving
            end
            if isstruct(msg)
                % The next experiment's holoRequest, whose config prime is already
                % posted. Hand it to the next compile rather than playing it.
                obj.preread = msg;
                obj.set_phase('waiting');
                return
            end
            obj.sequences = obj.sequences + 1;
            PlaySequence2K(obj.slm, obj.holograms, msg);
        end

        function note_error(obj, err)
            % Report a polling failure without spamming: the first one, then every
            % ERROR_REPORT_EVERY.
            obj.consecutive_errors = obj.consecutive_errors + 1;
            n = obj.consecutive_errors;
            if n == 1 || mod(n, obj.ERROR_REPORT_EVERY) == 0
                fprintf('[holo] poll error (%d in a row, still trying): %s\n', n, err.message);
            end
            obj.set_message(false, err.message);
        end

        function note_ok(obj)
            if obj.consecutive_errors > 0
                fprintf('[holo] polling recovered after %d failed attempt(s).\n', ...
                    obj.consecutive_errors);
                obj.consecutive_errors = 0;
            end
        end
    end
end

% =========================================================================
% Everything below is start_holo_listener's helper set, moved here unchanged.
% They are file-local functions rather than methods because they are pure
% functions of a prime / a channel table, and keeping them identical is what
% makes this a restructure rather than a rewrite.
% =========================================================================

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
    elseif local_declares_zero(prime)
        % The prime SPEAKS the protocol and says there are no opto channels: a
        % vis-only experiment. That is information, not the absence of it. Return
        % an empty table; the caller skips the prime entirely.
        %
        % This branch must come BEFORE the inventory fallback. Without it, a
        % vis prime (prime_info sends opto empty, n_channels 0, wavelengths [])
        % failed both tests above and fell through to "no channel info", so the
        % listener synthesised its whole inventory, loaded every calibration,
        % declared itself ready for channels nobody asked about, and then blocked
        % forever in generate_holograms_new -- which has no timeout and no abort
        % check.
        source = sprintf('prime declares 0 channels, protocol %s', ...
            char(local_field(prime, 'opto_protocol', '<unversioned>')));
    else
        w = reshape(double(inventory_wl), 1, []);
        for i = 1:numel(w)
            chans(i).name = sprintf('%dnm', w(i));
            chans(i).wavelength = w(i);
            chans(i).slm_board = [];
        end
        source = 'this listener''s inventory (prime carried no channel info)';
    end

    % Empty is legitimate ONLY via the zero-declaration branch above; anywhere else
    % it means this listener has no wavelengths either.
    assert(~isempty(chans) || ~isempty(source), 'holo_listener:noChannels', ...
        'Prime carried no opto channels and this listener has no wavelengths.');
    if isempty(chans)
        return
    end
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
    if isempty(comm), return; end
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

function tf = local_declares_zero(prime)
%LOCAL_DECLARES_ZERO True when the prime explicitly says "no opto channels".
%   Distinguishes a vis-only prime from a legacy prime that carries no channel
%   information at all. Either marker is enough: n_channels present and 0, or the
%   protocol tag present with an empty opto table. prime_info sends both.
    tf = false;
    if ~isstruct(prime)
        return
    end
    if isfield(prime, 'n_channels') && ~isempty(prime.n_channels) ...
            && double(prime.n_channels) == 0
        tf = true;
        return
    end
    % Protocol tag present => this DAQ knows how to declare channels, so an EMPTY
    % opto table is a statement rather than an omission. The isempty test is
    % load-bearing: without it a full multi-channel prime satisfies this too, and
    % while the caller's branch order happens to catch prime.opto first today, a
    % predicate that answers "declares zero channels" for a 2-channel prime is a
    % trap for whoever reorders those branches next.
    tf = isfield(prime, 'opto_protocol') && ~isempty(prime.opto_protocol) ...
        && isfield(prime, 'opto') && isempty(prime.opto);
end

function [v, src] = local_cfg(explicit, rig_path, fallback)
%LOCAL_CFG One config value: an explicit argument, else the rig, else a literal.
%   Returns the value and a short string naming where it came from, so the
%   listener's startup banner can say which config it is actually running on.
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
%   from the same input and get the same string.
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
    if isempty(comm), return; end
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
    if isempty(comm), return; end
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

function w = local_okword(ok)
    if ok, w = 'ok'; else, w = 'FAILED'; end
end
