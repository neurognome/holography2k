classdef HoloListenerMonitor < handle
    %HOLOLISTENERMONITOR Status window for the async holo listener: lamp + Start/Stop.
    %   start_holo_listener pops one of these up. It exists because the async
    %   listener is INVISIBLE at the prompt: several states look identical unless
    %   you remember to type l.status() (and still hold l) --
    %
    %       green  WAITING     polling for a prime; nothing armed
    %       green  SERVING     primed and armed, playing firing orders
    %       blue   COMPILING   loading calibrations and compiling holograms. This
    %                          phase BLOCKS the session, so the window is frozen
    %                          while it shows -- including the Stop button. Minutes
    %                          are normal; generate_holograms_new waits on the DAQ.
    %       amber  STARVED     running, but ticks are late, so something is blocking
    %                          the prompt and firing orders are being picked up late
    %       red    STOPPED     not polling -- either you stopped it, or the poll
    %                          timer errored (HoloListener.on_timer_error) and the
    %                          listener died while still looking started
    %       grey   NO LISTENER / DISPLAY STALLED -- the window itself is not live
    %
    %   CLOSING THE WINDOW STOPS THE LISTENER. This is the on/off control, not a
    %   passive monitor, so the X leaves this box un-primed for the next experiment;
    %   delete() prints a line saying so, because a silently dead holo listener is
    %   the expensive failure here (the DAQ waits out a 300 s transfer timeout and
    %   then runs anyway).
    %
    %   It does NOT stop the SLMs. A window closed mid-experiment must not blank a
    %   pattern the DAQ is still gating a laser against; the listener's own abort
    %   and re-prime paths are what stop the boards.
    %
    %   Usage:
    %       start_holo_listener                    % window opens with the listener
    %       start_holo_listener('nogui')           % listener only, no window
    %       HoloListenerMonitor(l)                 % re-open over a listener you hold
    %       HoloListenerMonitor(l, 'Visible', false)   % off-screen, for tests
    %
    %   See also start_holo_listener, HoloListener, HoloListener/listen_async,
    %   HoloListener/status, SIListenerMonitor (holodaq -- the same window for the
    %   SI box), test_holo_listener_monitor.

    properties (SetAccess = private)
        Listener          % the HoloListener being watched
        UIFigure
        Lamp
        StateLabel
        ToggleButton
        PollLabel
        PrimeLabel
        ChanLabel
        SlmLabel
        MsgLabel
        Refresh           % display timer ([] once stopped)
    end

    properties (Access = private)
        % Baseline for the achieved-period estimate: the poll timer's tick count
        % and a stopwatch id from the previous refresh. [] means "no baseline yet".
        LastTasks = []
        LastClock = []
    end

    properties (Constant)
        % Lamp colours: green safe / red live is the convention
        % PowerControllerCalibrated and SIListenerMonitor already use, so a rig
        % operator reads the same colours the same way on every window.
        COLOR_OK      = [0.30 0.65 0.40]
        COLOR_BUSY    = [0.20 0.55 0.95]
        COLOR_STARVED = [0.93 0.69 0.13]
        COLOR_STOPPED = [0.85 0.16 0.16]
        COLOR_NONE    = [0.60 0.60 0.60]
        % Button colours from the launcher's Start/Abort pair.
        BTN_STOP  = [0.85 0.33 0.31]
        BTN_START = [0.20 0.55 0.95]
        DIM_TEXT  = [0.45 0.45 0.45]
        % 1 s, deliberately NOT the listener's poll period: a faster lamp tells you
        % nothing, and this timer shares the event queue whose congestion is exactly
        % what STARVED reports.
        REFRESH_PERIOD = 1
        TIMER_NAME = 'holoListenerMonitor'
        FIG_NAME   = 'Holo Listener'
    end

    methods
        function app = HoloListenerMonitor(l, varargin)
            assert(nargin >= 1 && ~isempty(l) && isa(l, 'HoloListener') && isvalid(l), ...
                'HoloListenerMonitor:badListener', ...
                ['HoloListenerMonitor needs the HoloListener handle start_holo_listener ' ...
                 'returned: HoloListenerMonitor(l).']);

            p = inputParser();
            p.addParameter('Visible', true, @(v) islogical(v) || isnumeric(v));
            p.parse(varargin{:});

            app.Listener = l;
            % Take over from any previous window BEFORE building ours, so the
            % sweep cannot eat the window we are about to create.
            HoloListenerMonitor.closeStaleMonitors();
            app.buildUI(logical(p.Results.Visible));
            % Repaint the moment the listener changes phase, rather than up to a
            % second later. This is what puts COMPILING on the lamp BEFORE the
            % compile blocks the session, instead of after it finishes.
            l.on_phase = @() app.refresh(true);
            app.refresh();            % paint a real state before the first tick
            % Realizing a uifigure blocks the event queue for a second or two, which
            % genuinely starves the listener -- but blaming the listener for a stall
            % this window just caused is a false amber on every single launch. Let
            % the figure finish, then start timing from a quiet prompt.
            drawnow
            app.resetTickBaseline();
            app.startRefreshTimer();
        end

        function delete(app)
            % Closing the window stops the listener. Deliberate, and announced:
            % this leaves the holo box un-primed, and the whole point of the window
            % is that the listener's state is never a surprise.
            app.stopRefreshTimer();
            l = app.Listener;
            stopped = false;
            if ~isempty(l) && isvalid(l)
                if ~isempty(l.poll) && isvalid(l.poll)
                    try
                        l.stop();
                        stopped = true;
                    catch
                    end
                end
                % Drop the hook so a dead window cannot be called back into.
                try, l.on_phase = []; catch, end
            end
            if stopped
                fprintf(['[holo] status window closed -> async listener STOPPED. ' ...
                         'Restart with start_holo_listener.\n' ...
                         '       The SLMs were left as they were.\n']);
            end
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                app.UIFigure.CloseRequestFcn = '';   % already tearing down
                delete(app.UIFigure);
            end
        end

        function detach(app)
            %DETACH Tear this window down WITHOUT stopping the listener.
            %   Only for handing over to a replacement window; a normal close
            %   goes through delete(), which does stop it.
            app.stopRefreshTimer();
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                app.UIFigure.CloseRequestFcn = '';
                delete(app.UIFigure);
            end
            % Drop the listener so a later destructor call on this dead app cannot
            % reach in and stop the listener the new window owns.
            app.Listener = [];
        end

        function refresh(app, force)
            %REFRESH Re-read the poll timer and repaint. Never throws.
            %   Under a timer an escaping error stops the timer, which would freeze
            %   the lamp on a stale reading -- a window still saying SERVING after
            %   the listener died is worse than no window.
            %
            %   force=true is the phase-change hook: it uses a full drawnow so the
            %   new phase is actually on screen before a blocking phase starts,
            %   where drawnow limitrate would be free to skip the paint.
            if nargin < 2, force = false; end
            if isempty(app.UIFigure) || ~isvalid(app.UIFigure), return; end
            try
                [state, period, achieved] = app.readListener();
                app.applyState(state, period, achieved);
            catch err
                app.Lamp.Color = app.COLOR_NONE;
                app.StateLabel.Text = 'DISPLAY ERROR';
                app.PollLabel.Text = err.message;
            end
            if force
                drawnow
            else
                drawnow limitrate
            end
        end

        function tf = isListening(app)
            l = app.Listener;
            tf = ~isempty(l) && isvalid(l) && l.is_listening();
        end

        function onToggle(app)
            l = app.Listener;
            if isempty(l) || ~isvalid(l), return; end
            try
                if app.isListening()
                    l.stop();
                else
                    l.listen_async();
                end
            catch err
                % A callback must not throw into the event queue; show it instead.
                app.StateLabel.Text = 'ERROR';
                app.PollLabel.Text = err.message;
                return
            end
            app.resetTickBaseline();   % a restart must not be timed against the old run
            app.refresh();             % the lamp agrees with the click now, not in a second
        end

        function onRefreshError(app, src)
            %ONREFRESHERROR ErrorFcn for the display timer. Public like
            %   HoloListener.on_timer_error, and for the same reason: it is the only
            %   way to exercise the dead-display path from a test.
            %
            %   refresh() catches everything, so reaching here means the display is
            %   dead. DELETE both handles rather than dropping them -- an orphaned
            %   timer is the zombie on_timer_error exists to prevent, one level down.
            %   Then say the window is stale: a frozen lamp still reading SERVING is
            %   the one outcome worse than no lamp at all.
            for h = {src, app.Refresh}
                t = h{1};
                if ~isempty(t) && isa(t, 'timer') && isvalid(t)
                    try, stop(t);   catch, end
                    try, delete(t); catch, end
                end
            end
            app.Refresh = [];
            if ~isempty(app.Lamp) && isvalid(app.Lamp)
                app.Lamp.Color = app.COLOR_NONE;
                app.StateLabel.Text = 'DISPLAY STALLED';
                app.PollLabel.Text = 'refresh timer died; this window is no longer live';
            end
        end
    end

    methods (Static)
        function state = stateFrom(alive, running, phase, expected, achieved)
            %STATEFROM Lamp state from the numbers the listener exposes.
            %   Static and pure so the truth table is testable off the rig: a
            %   starved listener cannot be staged by driving a real timer, and the
            %   compiling phase cannot be staged without SLM hardware.
            %
            %   EXPECTED, not the bare timer period, is what starvation is measured
            %   against: a serving tick deliberately blocks for serve_read, so
            %   comparing it to the timer period would paint amber on every healthy
            %   experiment. And COMPILING is never starved -- a tick that blocks for
            %   minutes is that phase working as designed.
            if ~alive
                state = 'none';
            elseif ~running
                state = 'stopped';
            elseif strcmp(phase, 'compiling')
                state = 'compiling';
            elseif ~isnan(expected) && ~isnan(achieved) && achieved > 3 * expected
                % The same 3x rule holodaq's Receiver.status prints STARVED for --
                % one threshold, not two that can disagree with each other.
                state = 'starved';
            else
                state = phase;   % 'waiting' | 'serving'
            end
        end

        function closeStaleMonitors()
            %CLOSESTALEMONITORS One window, one display timer.
            %   The rule listen_async enforces for poll timers, applied to the
            %   display: running start_holo_listener twice must not leave a dead
            %   window ticking beside the live one.
            figs = findall(groot, 'Type', 'figure', 'Name', HoloListenerMonitor.FIG_NAME);
            for i = 1:numel(figs)
                f = figs(i);
                old = f.UserData;
                if isa(old, 'HoloListenerMonitor') && isvalid(old)
                    try, old.detach(); catch, end   % detach, NOT delete: keep the listener up
                end
                if isvalid(f)
                    f.CloseRequestFcn = '';
                    delete(f);
                end
            end
            % Orphans whose app object is already gone: the timer is all that is
            % left of them, and it is invisible except through timerfindall.
            ts = timerfindall('Name', HoloListenerMonitor.TIMER_NAME);
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
        function buildUI(app, visible)
            app.UIFigure = uifigure('Name', app.FIG_NAME, ...
                'Position', [100 100 460 200], 'Resize', 'off', ...
                'Visible', local_onoff(visible));
            app.UIFigure.CloseRequestFcn = @(s, e) delete(app);
            % Makes the window (and so the app) findable again from the prompt:
            % findall(0, 'Name', 'Holo Listener').UserData
            app.UIFigure.UserData = app;

            main = uigridlayout(app.UIFigure, [2 1]);
            main.RowHeight   = {'fit', 'fit'};
            main.ColumnWidth = {'1x'};
            main.RowSpacing  = 12;
            main.Padding     = [12 12 12 12];

            % Row 1: lamp, state word, one toggle button
            top = uigridlayout(main, [1 3]);
            top.ColumnWidth   = {24, '1x', 90};
            top.ColumnSpacing = 8;
            top.Padding       = [0 0 0 0];
            app.Lamp = uilamp(top, 'Color', app.COLOR_NONE);
            app.StateLabel = uilabel(top, 'Text', '', 'FontWeight', 'bold');
            app.ToggleButton = uibutton(top, 'Text', 'Stop', ...
                'FontColor', [1 1 1], 'FontWeight', 'bold', ...
                'BackgroundColor', app.BTN_STOP, ...
                'ButtonPushedFcn', @(s, e) app.onToggle());

            % Row 2: the detail lines
            det = uigridlayout(main, [5 2]);
            det.RowHeight     = repmat({'fit'}, 1, 5);
            det.ColumnWidth   = {52, '1x'};
            det.RowSpacing    = 2;
            det.ColumnSpacing = 6;
            det.Padding       = [0 0 0 0];
            app.PollLabel  = app.detailRow(det, 'poll');
            app.PrimeLabel = app.detailRow(det, 'prime');
            app.ChanLabel  = app.detailRow(det, 'channels');
            app.SlmLabel   = app.detailRow(det, 'slm');
            app.MsgLabel   = app.detailRow(det, 'last');
        end

        function lbl = detailRow(app, parent, name)
            uilabel(parent, 'Text', name, 'FontColor', app.DIM_TEXT, 'FontSize', 11);
            lbl = uilabel(parent, 'Text', '', 'FontColor', app.DIM_TEXT, 'FontSize', 11);
        end

        function [state, expected, achieved] = readListener(app)
            l = app.Listener;
            alive = false; running = false; ph = 'waiting';
            expected = NaN; achieved = NaN;
            if ~isempty(l) && isvalid(l)
                alive = true;
                ph = l.phase;
                expected = l.expected_period();
                % Read the timer directly, NOT l.status(): status() prints two lines
                % on every call, and at 1 Hz that buries the command window the
                % listener also reports primes and errors on.
                if ~isempty(l.poll) && isvalid(l.poll)
                    running  = strcmp(l.poll.Running, 'on');
                    achieved = app.tickRate(l);
                end
            end
            if ~running
                app.resetTickBaseline();   % a restart must not be timed against this
            end
            state = HoloListenerMonitor.stateFrom(alive, running, ph, expected, achieved);
        end

        function resetTickBaseline(app)
            app.LastTasks = [];
            app.LastClock = [];
        end

        function achieved = tickRate(app, l)
            %TICKRATE Achieved tick period over the ticks SINCE THE LAST REFRESH.
            %   Deliberately not timer.AveragePeriod: that average is cumulative
            %   over the timer's whole life, so one long COMPILING tick biases it for
            %   good and leaves an amber lamp for the rest of the session -- while a
            %   long healthy run washes out a stall happening right now. A lamp has
            %   to mean now, so measure the recent interval instead.
            achieved = NaN;
            n = l.poll.TasksExecuted;
            if ~isempty(app.LastClock)
                elapsed = toc(app.LastClock);
                delta   = n - app.LastTasks;
                % delta < 0 means the timer was restarted under us: re-baseline
                % rather than report nonsense. And give it at least a couple of
                % expected periods' worth of elapsed time before passing a verdict.
                floor_s = 2 * max(l.poll_period, 0.01);
                if ~isnan(l.expected_period())
                    floor_s = 2 * l.expected_period();
                end
                if delta >= 0 && elapsed >= floor_s
                    % delta == 0 -> Inf: not one tick in that whole window, which is
                    % starvation in its most severe form, not a missing reading.
                    achieved = elapsed / delta;
                end
            end
            app.LastTasks = n;
            app.LastClock = tic;
        end

        function applyState(app, state, expected, achieved)
            switch state
                case 'waiting'
                    app.Lamp.Color = app.COLOR_OK;
                    app.StateLabel.Text = 'WAITING FOR A PRIME';
                case 'serving'
                    app.Lamp.Color = app.COLOR_OK;
                    app.StateLabel.Text = 'SERVING';
                case 'compiling'
                    app.Lamp.Color = app.COLOR_BUSY;
                    app.StateLabel.Text = 'COMPILING (window frozen)';
                case 'starved'
                    app.Lamp.Color = app.COLOR_STARVED;
                    app.StateLabel.Text = 'STARVED';
                case 'stopped'
                    app.Lamp.Color = app.COLOR_STOPPED;
                    app.StateLabel.Text = 'STOPPED';
                case 'none'
                    app.Lamp.Color = app.COLOR_NONE;
                    app.StateLabel.Text = 'NO LISTENER';
            end

            if any(strcmp(state, {'waiting', 'serving', 'compiling', 'starved'}))
                app.ToggleButton.Text = 'Stop';
                app.ToggleButton.BackgroundColor = app.BTN_STOP;
            else
                app.ToggleButton.Text = 'Start';
                app.ToggleButton.BackgroundColor = app.BTN_START;
            end
            app.ToggleButton.Enable = local_onoff(~strcmp(state, 'none'));

            app.PollLabel.Text  = app.pollText(state, expected, achieved);
            app.PrimeLabel.Text = app.primeText();
            app.ChanLabel.Text  = app.chanText();
            app.SlmLabel.Text   = app.slmText();
            app.MsgLabel.Text   = app.msgText();
        end

        function s = pollText(app, state, expected, achieved)
            if strcmp(state, 'compiling')
                s = 'blocked in the compile -- no ticks until it finishes';
                return
            end
            if isnan(expected)
                s = '—';
                return
            end
            if isnan(achieved)
                s = sprintf('%.2f s expected / — achieved', expected);
            elseif isinf(achieved)
                s = sprintf('%.2f s expected / NOT TICKING', expected);
            else
                s = sprintf('%.2f s expected / %.2f s achieved', expected, achieved);
            end
            if strcmp(state, 'starved')
                s = [s '   something is blocking the prompt'];
            end
        end

        function s = primeText(app)
            l = app.Listener;
            s = '—';
            if isempty(l) || ~isvalid(l), return; end
            st = l.stem();
            if strcmp(st, '?')
                % No stem but a finite seq is the normal state at startup:
                % flush_stale adopts whatever prime is already on the channel as its
                % baseline WITHOUT running it, so that seq is a prime this listener
                % has deliberately not served. Say that, rather than showing a bare
                % '?' next to a number and letting it read as a prime that happened.
                if isinf(l.last_seq)
                    s = 'nothing yet';
                else
                    s = sprintf('nothing yet (channel at seq %g)', l.last_seq);
                end
                return
            end
            s = sprintf('%s   (seq %g, %d sequence(s) played)', st, l.last_seq, l.sequences);
        end

        function s = chanText(app)
            l = app.Listener;
            s = '—';
            if isempty(l) || ~isvalid(l) || isempty(l.chans), return; end
            parts = cell(1, numel(l.chans));
            for i = 1:numel(l.chans)
                parts{i} = sprintf('%s@%dnm', l.chans(i).name, l.chans(i).wavelength);
            end
            s = sprintf('%s   (from %s)', strjoin(parts, ', '), l.source);
        end

        function s = slmText(app)
            l = app.Listener;
            s = '—';
            if isempty(l) || ~isvalid(l) || isempty(l.slm_all), return; end
            boards = [l.slm_all.board_id];
            % "armed" is this listener's own record of having started the boards it
            % resolved -- not a read back off the hardware, which the SLM objects do
            % not expose. It says what the listener did, not what the board is doing.
            if isempty(l.slm)
                armed = 'none armed';
            else
                armed = sprintf('%d armed', numel(l.slm));
            end
            s = sprintf('%s nm on boards %s   (%s)', mat2str(l.wl), mat2str(boards), armed);
        end

        function s = msgText(app)
            l = app.Listener;
            s = '—';
            if isempty(l) || ~isvalid(l), return; end
            if l.last_ok
                s = l.last_message;
            else
                s = ['FAILED: ' l.last_message];
            end
        end

        function startRefreshTimer(app)
            app.Refresh = timer('Name', app.TIMER_NAME, ...
                'ExecutionMode', 'fixedSpacing', ...  % space AFTER the callback, as the
                'BusyMode', 'drop', ...               % listener's own poll timer does
                'Period', app.REFRESH_PERIOD, ...
                'StartDelay', app.REFRESH_PERIOD, ...  % start() fires at once otherwise,
                ...                                    % sampling an interval of zero
                'TimerFcn', @(~, ~) app.refresh(), ...
                'ErrorFcn', @(src, ~) app.onRefreshError(src));
            start(app.Refresh);
        end

        function stopRefreshTimer(app)
            if ~isempty(app.Refresh) && isvalid(app.Refresh)
                try, stop(app.Refresh);   catch, end
                try, delete(app.Refresh); catch, end
            end
            app.Refresh = [];
        end
    end
end

function s = local_onoff(tf)
    % Char 'on'/'off' rather than a logical: accepted by every release that has
    % uifigure, unlike the logical shorthand.
    if tf, s = 'on'; else, s = 'off'; end
end
