function l = start_holo_listener(varargin)
%START_HOLO_LISTENER Config-primed persistent holography listener (remote).
%   Run once on the holo computer INSTEAD of MsocketHolorequest2K when driving
%   experiments remotely. It waits for a prime from the DAQ master
%   (config/holo, see prime_info.m), auto-loads the newest dated calibration
%   per wavelength (find_latest_calib), compiles holograms (via the DAQ's
%   transferHR), arms the SLM(s), acks to config/holo_status, then serves
%   firing-order sequences. It RE-PRIMES automatically for the next experiment
%   -- no restart, no hand-edited calibration paths.
%
%   ASYNC BY DEFAULT (changed 2026-08-06). The listener runs on a timer and the
%   MATLAB prompt stays usable, matching start_si_listener on the ScanImage box.
%
%   A small status window opens with it (HoloListenerMonitor): a lamp saying
%   whether the listener is waiting / compiling / serving / starved / stopped, and
%   one Start/Stop button. CLOSING THAT WINDOW STOPS THE LISTENER (it does not stop
%   the SLMs); 'nogui' skips it entirely, which is what a headless invocation wants.
%
%   Usage (on the holo computer):
%       l = start_holo_listener                          % async, with the window
%       l = start_holo_listener('nogui')                 % async, no window
%       l = start_holo_listener('blocking')              % the pre-2026-08 loop, Ctrl-C
%       l = start_holo_listener('Wavelengths', [1100 900])
%       l = start_holo_listener('Timeout', 1700, 'nogui')
%
%       l.status()         % am I actually listening? also flags a starved timer
%       l.stop()           % halt the async listener
%       l.listen_async()   % start it again
%       HoloListenerMonitor(l)   % re-open the status window
%
%   WHAT ASYNC DOES NOT GIVE YOU: an uninterrupted session, only a usable one.
%   Timer callbacks run on MATLAB's event queue, serviced when MATLAB is idle or at
%   pause/drawnow/waitfor, so a long blocking operation at the prompt DELAYS the
%   listener -- l.status() reports the achieved period and says when it is starved.
%   And the COMPILING phase is one long blocking tick (generate_holograms_new waits
%   on the DAQ's holoRequests with no timeout), during which nothing else in this
%   session runs, including the Stop button. The lamp says COMPILING while that is
%   happening; it is not hung.
%
%   Firing-order latency: an async serving tick picks up a sequence at most
%   poll_period (0.05 s) after it lands, where the blocking loop read continuously.
%   That sits inside the pause(0.1) the DAQ's HolochatInterface.send already
%   performs after posting a sequence and before the trial can start, so it cannot
%   make a sequence arrive after its trigger. Use 'blocking' if you want the old
%   continuous read regardless; MsocketHolorequest2K remains the manual fallback.
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
%     'Offline'     (default false) -- build the listener without opening SLMs or
%         a broker connection. It never primes anything; it exists so the window
%         and the argument handling can be exercised off the rig.
%
%   Each opto channel may pin its own slm_board and slm_lut in the rig file, which
%   is what a rig with two arms on one wavelength needs (the board cannot be
%   derived from the wavelength then). A channel that pins neither falls back to
%   get_slm's wavelength mapping, so an unpinned rig behaves exactly as before.
%
%   Paths: this checkout is added by deriving it from HoloListener's own location,
%   so it does not name one machine's username. The Meadowlark SDK is a per-machine
%   install rather than part of a checkout, so it comes from rig.paths.slm_sdk.
%   holodaq must already be on this machine's path -- it has to be, since
%   HolochatInterface is what reads the rig config in the first place.
%
%   How re-prime stays race-free: the serving phase reads msg/holo and switches on
%   type -- a firing-order CELL is played; a holoRequest STRUCT means the next
%   experiment has begun (its config prime is already posted), so it is handed
%   back to be compiled. Sequences and holoRequests therefore never collide.
%
%   See also: HoloListener (the loop itself), HoloListenerMonitor (the window),
%             MsocketHolorequest2K (manual fallback), PlaySequence2K,
%             find_latest_calib, generate_holograms_new, prime_info,
%             start_si_listener (holodaq -- the same shape for the SI box).

    % Name/value options stay as they were; 'blocking'/'async'/'gui'/'nogui' are
    % keywords in any position. Split out into its own function so it is testable --
    % the constructor below addpaths the whole repo and opens SLM hardware.
    [opts, mode, gui] = parse_holo_listener_args(varargin{:});

    % Make sure this checkout is on the path before naming HoloListener, so the
    % function works when only holorequest/ happens to be reachable.
    here = fileparts(mfilename('fullpath'));   % holography2k/holorequest
    addpath(genpath(fileparts(here)));

    l = HoloListener(opts{:});

    if strcmp(mode, 'blocking')
        if gui
            % Blocking mode wins over the window: listen() is a while-true loop, so
            % a Stop button could never be clicked -- only Ctrl-C can interrupt it.
            fprintf('HoloListener: no status window in blocking mode (Ctrl-C stops it).\n');
        end
        l.listen();        % owns this session until Ctrl-C
    else
        % Listener FIRST, window second: a GUI that fails to build must never cost
        % you the listener, which is the thing that actually arms the SLMs.
        l.listen_async();  % returns immediately; l.stop() to halt
        if gui
            try
                HoloListenerMonitor(l);
            catch err
                fprintf(['[holo] status window failed to open (%s); the listener is ' ...
                         'running -- use l.status() / l.stop().\n'], err.message);
            end
        end
    end

    if nargout == 0
        % The old signature returned nothing, so `start_holo_listener` at the prompt
        % printed no ans. Keep that, but the handle is still reachable from the
        % window: findall(0, 'Name', 'Holo Listener').UserData.Listener
        clear l
    end
end
