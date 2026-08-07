function test_holo_listener_monitor()
%TEST_HOLO_LISTENER_MONITOR  The holo listener status window. Offline, no rig, no broker.
%
%   >> addpath(genpath(<holography2k>)); addpath(genpath(<holodaq>));
%   >> test_holo_listener_monitor
%
%   Windows are built with 'Visible', false so this runs unattended, and every
%   listener is built with 'Offline', true so nothing opens an SLM board or reaches
%   for the broker. The cases are the ways a status light can lie: showing SERVING
%   for a listener that stopped, showing amber for a listener that is only doing
%   what the COMPILING phase is supposed to do, showing green for one that really
%   is starved, a second window ticking beside the live one, and a torn-down window
%   leaving an orphaned timer behind.
%
%   Sections 4 onwards need a JVM (uifigure). They are SKIPPED, loudly, on a
%   headless MATLAB -- which is the only way the pure sections can be run on a
%   machine that is not the rig.
%
%   See also HoloListenerMonitor, HoloListener, parse_holo_listener_args,
%   test_si_listener_monitor (holodaq -- the same suite for the SI box).

    % --- 1. 'gui' / 'nogui' / 'blocking' argument forms ------------------------
    % Every previously valid invocation must keep meaning what it meant, so the
    % default is a window and the name/value options pass through untouched.
    cases = { ...
        {},                                {},                    'async',    true  ; ...
        {'nogui'},                         {},                    'async',    false ; ...
        {'gui'},                           {},                    'async',    true  ; ...
        {'blocking'},                      {},                    'blocking', true  ; ...
        {'NOGUI'},                         {},                    'async',    false ; ...  % case-insensitive
        {'  nogui  '},                     {},                    'async',    false ; ...  % whitespace-tolerant
        {"nogui"},                         {},                    'async',    false ; ...  % string, not char
        {'Timeout', 1700},                 {'Timeout', 1700},     'async',    true  ; ...
        {'nogui', 'Timeout', 1700},        {'Timeout', 1700},     'async',    false ; ...  % order-insensitive
        {'Timeout', 1700, 'nogui'},        {'Timeout', 1700},     'async',    false ; ...
        {'CalibDir', 'C:\c', 'blocking'},  {'CalibDir', 'C:\c'},  'blocking', true  };

    for i = 1:size(cases, 1)
        [opts, mode, gui] = parse_holo_listener_args(cases{i, 1}{:});
        want = cases{i, 2};
        assert(numel(opts) == numel(want), 'case %d: %d opts, want %d', ...
            i, numel(opts), numel(want));
        for k = 1:numel(want)
            assert(isequal(opts{k}, want{k}), 'case %d: opt %d differs', i, k);
        end
        assert(strcmp(mode, cases{i, 3}), 'case %d: mode "%s", want "%s"', ...
            i, mode, cases{i, 3});
        assert(islogical(gui) && gui == cases{i, 4}, ...
            'case %d: gui %d, want %d', i, gui, cases{i, 4});
    end

    % last keyword wins, so a scripted override can append one
    [~, ~, gui] = parse_holo_listener_args('nogui', 'gui');
    assert(gui, 'the last gui keyword should win');
    [~, mode] = parse_holo_listener_args('blocking', 'async');
    assert(strcmp(mode, 'async'), 'the last mode keyword should win');

    % Stripping a keyword must never leave a name without its value. This is the
    % failure that would otherwise reach inputParser as a generic complaint.
    threw = false;
    try
        parse_holo_listener_args('Timeout');
    catch err
        threw = strcmp(err.identifier, 'start_holo_listener:oddArgs');
    end
    assert(threw, 'a dangling option name should be refused here, with a usable message');

    % --- 2. stateFrom truth table ---------------------------------------------
    % A pure static function precisely because the interesting states cannot be
    % staged: timer.AveragePeriod is read-only, and COMPILING needs SLM hardware.
    st = @HoloListenerMonitor.stateFrom;
    assert(strcmp(st(false, false, 'waiting', NaN,  NaN),  'none'), ...
        'no listener -> none');
    assert(strcmp(st(true,  false, 'serving', 0.45, NaN),  'stopped'), ...
        'not running -> stopped, whatever phase it stopped in');
    assert(strcmp(st(true,  true,  'waiting', 0.07, 0.07), 'waiting'), ...
        'on schedule while idle -> waiting');
    assert(strcmp(st(true,  true,  'serving', 0.45, 0.46), 'serving'), ...
        'on schedule while primed -> serving');
    assert(strcmp(st(true,  true,  'waiting', 0.07, NaN),  'waiting'), ...
        'a timer that has not yet been sampled must not read as starved');

    % The whole reason expected_period exists: a serving tick blocks for
    % serve_read, so measuring it against the bare timer period would paint amber
    % on every healthy experiment.
    assert(strcmp(st(true, true, 'serving', 0.45, 0.46), 'serving'), ...
        'a serving tick at its expected cost is NOT starved');
    assert(strcmp(st(true, true, 'serving', 0.05, 0.46), 'starved'), ...
        'the same reading IS starved when measured against the bare poll period');

    assert(strcmp(st(true, true, 'waiting', 0.07, 0.30), 'starved'), ...
        '4x the expected period -> starved');
    assert(strcmp(st(true, true, 'waiting', 0.10, 0.29), 'waiting'), ...
        'just under 3x must stay green (same threshold as Receiver.status)');
    assert(strcmp(st(true, true, 'waiting', 0.10, 0.31), 'starved'), ...
        'just over 3x -> starved');
    assert(strcmp(st(true, true, 'waiting', 0.07, Inf), 'starved'), ...
        'Inf is "not one tick in the whole window" -- the worst starvation, not a gap');

    % COMPILING outranks starvation: it blocks the session by design, for minutes,
    % and calling that starved would make the one phase that always looks stalled
    % also always look broken.
    assert(strcmp(st(true, true, 'compiling', NaN, NaN), 'compiling'), ...
        'compiling -> compiling');
    assert(strcmp(st(true, true, 'compiling', 0.07, 600), 'compiling'), ...
        'a compiling tick is never starved, however long it took');
    % ...but a stopped listener still reads STOPPED even if it died compiling.
    assert(strcmp(st(true, false, 'compiling', NaN, NaN), 'stopped'), ...
        'stopped outranks compiling -- a dead listener is not busy');

    % --- 3. the offline listener's phase machine ------------------------------
    % No SLMs, no broker: ticking must be a safe no-op, and expected_period must
    % track the phase, since that is what the lamp measures starvation against.
    l = HoloListener('Offline', true);
    cleanL = onCleanup(@() l.stop());
    assert(isempty(l.comm) && isempty(l.slm_all), 'Offline must open nothing');
    assert(strcmp(l.phase, 'waiting'), 'a fresh listener waits');
    assert(~l.is_listening(), 'a fresh listener is not listening');

    l.tick();   % must not throw with no broker
    l.tick();
    assert(strcmp(l.phase, 'waiting'), 'ticking offline must not change phase');

    e_wait = l.expected_period();
    l.set_phase('serving');
    e_serve = l.expected_period();
    l.set_phase('compiling');
    e_comp = l.expected_period();
    assert(e_serve > e_wait, ...
        'a serving tick must be expected to cost more than an idle one (it blocks in the read)');
    assert(abs(e_serve - (l.poll_period + l.serve_read)) < 1e-9, ...
        'serving cost must be the poll period plus the blocking read');
    assert(isnan(e_comp), 'compiling has no expected period -- it is never starved');
    l.set_phase('waiting');

    % A serving tick with no broker must not throw -- the Offline guard, which is
    % also what stops a test from hammering a nonexistent holochat.
    l.set_phase('serving');
    l.tick();
    assert(strcmp(l.phase, 'serving'), 'an offline serving tick is a no-op, not a fault');
    l.set_phase('waiting');

    % phase is read-only from outside: set_phase is the only door in, because it is
    % the only door that fires on_phase.
    threw = false;
    try
        l.phase = 'serving'; %#ok<NASGU>
    catch
        threw = true;
    end
    assert(threw, 'phase must not be settable behind set_phase''s back');
    threw = false;
    try
        l.set_phase('nonsense');
    catch err
        threw = strcmp(err.identifier, 'HoloListener:badPhase');
    end
    assert(threw, 'an unknown phase must be refused');

    % --- 3b. the on_phase hook -------------------------------------------------
    % This is what paints COMPILING on the lamp BEFORE the compile blocks the
    % session for minutes. A Map because a plain counter would be captured by value.
    hits = containers.Map('KeyType', 'char', 'ValueType', 'double');
    hits('n') = 0;
    l.on_phase = @() local_bump(hits);
    l.set_phase('serving');
    assert(hits('n') == 1, 'a phase change must fire the hook');
    l.set_phase('serving');
    assert(hits('n') == 1, 'setting the SAME phase must not fire the hook');
    l.set_phase('waiting');
    assert(hits('n') == 2, 'the hook fires on every real change');

    % A hook that THROWS must not take the listener with it: it runs inside tick(),
    % where an escaping error stops the poll timer and kills the listener outright.
    l.on_phase = @() error('holo:testHook', 'boom');
    l.set_phase('serving');
    assert(strcmp(l.phase, 'serving'), 'a throwing hook must not block the phase change');
    l.tick();
    assert(strcmp(l.phase, 'serving'), 'a throwing hook must not derail the listener');
    l.on_phase = [];
    l.set_phase('waiting');

    % --- 4. GUI sections ------------------------------------------------------
    if ~usejava('jvm') || ~usejava('awt')
        fprintf(['test_holo_listener_monitor: sections 1-3 PASSED; sections 4+ ' ...
                 'SKIPPED (no JVM,\n    so uifigure is unavailable). Run this on a ' ...
                 'MATLAB with a desktop to cover the window.\n']);
        return
    end

    % 4a. the lamp tracks a live listener
    l2 = HoloListener('Offline', true);
    l2.poll_period = 0.05;
    cleanL2 = onCleanup(@() l2.stop());
    l2.listen_async();
    app = HoloListenerMonitor(l2, 'Visible', false);
    cleanApp = onCleanup(@() delete(app));
    assert(app.isListening(), 'a running listener must read as listening');
    pause(0.5); app.refresh();
    assert(strcmp(app.StateLabel.Text, 'WAITING FOR A PRIME'), ...
        'an idle listener should read WAITING, got "%s"', app.StateLabel.Text);
    assert(strcmp(app.ToggleButton.Text, 'Stop'), 'a running listener offers Stop');

    % 4b. the monitor installs the phase hook, so the lamp repaints on the phase
    %     change itself rather than up to a second later on the display timer.
    %     The listener has to stay RUNNING for this (a stopped one reads STOPPED
    %     whatever phase it is in), but its ticks must not land between the phase
    %     change and the assert -- an offline compile bounces straight back to
    %     waiting -- so slow the timer right down for the duration.
    l2.stop();
    l2.poll_period = 30;
    l2.listen_async();
    l2.set_phase('compiling');
    assert(strcmp(app.StateLabel.Text, 'COMPILING (window frozen)'), ...
        'the phase hook must repaint the lamp immediately, got "%s"', app.StateLabel.Text);
    l2.set_phase('waiting');
    assert(strcmp(app.StateLabel.Text, 'WAITING FOR A PRIME'), ...
        'and again on the way back, got "%s"', app.StateLabel.Text);
    l2.stop();
    l2.poll_period = 0.05;
    l2.listen_async();

    % 4c. the Stop button really stops, and the lamp says so
    app.onToggle();
    assert(~l2.is_listening(), 'Stop must halt the poll timer');
    assert(strcmp(app.StateLabel.Text, 'STOPPED'), ...
        'a stopped listener must read STOPPED, got "%s"', app.StateLabel.Text);
    assert(strcmp(app.ToggleButton.Text, 'Start'), 'a stopped listener offers Start');

    % ...and Start restarts it
    app.onToggle();
    assert(l2.is_listening(), 'Start must restart the poll timer');
    app.refresh();
    assert(~strcmp(app.StateLabel.Text, 'STOPPED'), 'the lamp must follow the restart');

    % 4d. one window, one display timer: a second monitor retires the first
    app2 = HoloListenerMonitor(l2, 'Visible', false);
    cleanApp2 = onCleanup(@() delete(app2));
    assert(l2.is_listening(), ...
        'opening a replacement window must NOT stop the listener the old one watched');
    assert(numel(timerfindall('Name', HoloListenerMonitor.TIMER_NAME)) == 1, ...
        'a retired window must not leave its display timer ticking');
    figs = findall(groot, 'Type', 'figure', 'Name', HoloListenerMonitor.FIG_NAME);
    assert(numel(figs) == 1, 'only one status window may exist at a time');

    % 4e. closing the window stops the listener -- the documented, deliberate
    %     behaviour, and the one most worth pinning down
    delete(app2);
    assert(~l2.is_listening(), 'closing the status window must stop the listener');
    assert(isempty(timerfindall('Name', HoloListenerMonitor.TIMER_NAME)), ...
        'a closed window must not leave an orphaned display timer');
    assert(isempty(l2.on_phase), ...
        'a closed window must drop its hook, or a dead app gets called back into');

    % 4f. a dead display says so rather than freezing on a stale reading
    l3 = HoloListener('Offline', true);
    cleanL3 = onCleanup(@() l3.stop());
    l3.listen_async();
    app3 = HoloListenerMonitor(l3, 'Visible', false);
    cleanApp3 = onCleanup(@() delete(app3));
    app3.onRefreshError(app3.Refresh);
    assert(strcmp(app3.StateLabel.Text, 'DISPLAY STALLED'), ...
        'a dead refresh timer must be announced, not left showing the last good state');
    assert(isempty(timerfindall('Name', HoloListenerMonitor.TIMER_NAME)), ...
        'onRefreshError must delete the timer, not just drop the handle');
    assert(l3.is_listening(), 'a dead DISPLAY must not take the listener with it');

    delete(cleanApp3); delete(cleanApp2); delete(cleanApp);
    delete(cleanL3);   delete(cleanL2);   delete(cleanL);
    fprintf('test_holo_listener_monitor: PASSED\n');
end

function local_bump(m)
    m('n') = m('n') + 1;
end
