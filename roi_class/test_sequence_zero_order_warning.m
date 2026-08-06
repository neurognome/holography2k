function test_sequence_zero_order_warning()
%TEST_SEQUENCE_ZERO_ORDER_WARNING  The ragged/zero-order guard fires when it should.
%
%   >> addpath(pwd); test_sequence_zero_order_warning
%
%   ONE laser power is commanded per sequence (calculate_power returns a scalar and
%   LaserPowerControl.set applies it once per trial), so a pattern with fewer targets
%   than the largest in the same sequence concentrates that fixed power into fewer
%   spots. zero_order_dump corrects it and DEFAULTS TO FALSE, so the mistake is silent
%   -- which is what this warning exists to stop.
%
%   Each case runs in a fresh MATLAB-free way by clearing the persistent guard, since
%   the warning is deliberately once-per-session.

    a3 = [100 100 0; 120 120 0; 140 140 0];      % 3 targets
    b1 = [200 200 40];                            % 1 target
    c3 = [300 300 60; 320 320 60; 340 340 60];   % 3 targets

    %% 1. ragged + multiple + dump OFF -> warns, and names the overdrive ------------
    id = local_capture(@() Sequence([Pattern(a3, [], false), Pattern(b1, [], false)]));
    assert(strcmp(id, 'Sequence:raggedWithoutZeroOrderDump'), ...
        'ragged sequence without the dump must warn, got ''%s''', id);

    %% 2. ragged + multiple + dump ON  -> silent -----------------------------------
    id = local_capture(@() Sequence([Pattern(a3, [], true), Pattern(b1, [], true)]));
    assert(isempty(id), 'the dump is the fix; it must not warn (got ''%s'')', id);

    %% 3. UNIFORM sizes + dump off -> silent ---------------------------------------
    % Nothing to correct: every pattern gets the same share of the trial's power.
    id = local_capture(@() Sequence([Pattern(a3, [], false), Pattern(c3, [], false)]));
    assert(isempty(id), 'uniform sizes need no dump (got ''%s'')', id);

    %% 4. a SINGLE pattern, any size, dump off -> silent ---------------------------
    % One pattern per trial is the other way to be correct: power_per_cell is per pool
    % entry, so each trial can command its own power.
    id = local_capture(@() Sequence(Pattern(b1, [], false)));
    assert(isempty(id), 'a single-pattern sequence needs no dump (got ''%s'')', id);

    %% 5. only the SMALL pattern needs the dump ------------------------------------
    % The largest pattern has nothing to dump (size == max), so leaving it false is
    % harmless -- the guard must not fire on that alone.
    id = local_capture(@() Sequence([Pattern(a3, [], false), Pattern(b1, [], true)]));
    assert(isempty(id), ...
        'only sub-max patterns need the dump; the max one being false is fine (got ''%s'')', id);

    % ...but the reverse must still warn
    id = local_capture(@() Sequence([Pattern(a3, [], true), Pattern(b1, [], false)]));
    assert(strcmp(id, 'Sequence:raggedWithoutZeroOrderDump'), ...
        'the SMALL pattern is the one that matters, got ''%s''', id);

    disp('PASS test_sequence_zero_order_warning');
end

function id = local_capture(fn)
%LOCAL_CAPTURE  'Sequence:raggedWithoutZeroOrderDump' if fn raises it, '' otherwise.
%   Clears Sequence's persistent guard first, because the warning is once-per-session
%   by design and every case here would otherwise be swallowed after the first.
%
%   Pattern.to_struct calls struct() on a handle object, which raises
%   MATLAB:structOnObject on every Sequence construction. That is pre-existing noise
%   unrelated to this guard, so silence it -- otherwise lastwarn returns IT rather than
%   the warning under test and every case looks like a pass/fail at random.
    clear Sequence
    lastwarn('', '');
    state = warning('off', 'Sequence:raggedWithoutZeroOrderDump');
    state2 = warning('off', 'MATLAB:structOnObject');
    state3 = warning('off', 'Pattern:noZeroOrderDump');
    cleanup = onCleanup(@() [warning(state), warning(state2), warning(state3)]);
    fn();
    [~, last_id] = lastwarn();
    if strcmp(last_id, 'Sequence:raggedWithoutZeroOrderDump')
        id = last_id;
    else
        id = '';
    end
end
