function test_firing_order()
%TEST_FIRING_ORDER  The firing order survives, and refuses to go out wrong.
%
%   >> addpath(pwd); addpath(fullfile(pwd,'roi_class')); addpath(fullfile(pwd,'holorequest'))
%   >> test_firing_order
%
%   Two silent failures, both of which sent the SLM to the wrong holograms while the
%   laser gate fired on time -- so the light looked correct and nothing in the saved
%   stim record disagreed:
%
%     1. PlaySequence2K bounded the feed on `counter` (position in the order) while
%        indexing with order(counter). A firing order LONGER than the hologram stack
%        therefore stopped being fed partway, and the SLM held its last frame.
%     2. Sequence.ids used [patterns.id], which silently drops empty ids.

    warning('off', 'MATLAB:structOnObject');
    warning('off', 'Pattern:noZeroOrderDump');
    warning('off', 'Sequence:raggedWithoutZeroOrderDump');
    warning('off', 'StimInfo:slmTriggersMerge');

    %% 1. ids round-trip in the requested (permuted) order ------------------------
    pats = local_patterns(4);
    local_assign_ids(pats);
    order = [3 1 4 2];
    seq = Sequence(pats(order));
    assert(isequal(seq.ids, order), ...
        'firing order not preserved: wanted %s, got %s', ...
        mat2str(order), mat2str(seq.ids));

    %% 2. missing ids REFUSE rather than shortening the order ----------------------
    none = local_patterns(3);                       % never round-tripped
    local_must_error(@() local_ids(Sequence(none)), 'Sequence:missingPatternIds');

    partial = local_patterns(3);
    partial(1).id = 1; partial(3).id = 3;           % middle one empty
    local_must_error(@() local_ids(Sequence(partial)), 'Sequence:missingPatternIds');

    % a non-scalar id would flatten into the wrong length
    weird = local_patterns(2);
    weird(1).id = 1; weird(2).id = [2 3];
    local_must_error(@() local_ids(Sequence(weird)), 'Sequence:badPatternIds');

    %% 3. PlaySequence2K feeds EVERY entry, including orders longer than the stack --
    % This is the spiking_holography_clicked shape: one pattern id reused for a burst
    % of spikes, so 20 pulses over 4 holograms.
    cases = {
        'order == stack',                [3 1 4 2],             4
        'order LONGER than stack',       repelem([3 1 4 2], 5), 4
        'single pattern, many pulses',   ones(1, 8),            4
        'order shorter than stack',      [2 1],                 4
    };
    for k = 1:size(cases, 1)
        ord = cases{k, 2};
        fed = local_feed(ord, cases{k, 3});
        assert(isequal(fed, ord), ...
            '%s: fed %s, wanted %s', cases{k, 1}, mat2str(fed), mat2str(ord));
    end

    %% 4. an out-of-range id is still refused, not fed ----------------------------
    % PlaySequence2K's max(order) check blanks and returns before the feed loop; the
    % relaxed bound must not have opened a hole there.
    fed = local_feed([1 5 2], 4);       % 5 > 4 holograms
    assert(isempty(fed), 'an out-of-range id must not be fed, got %s', mat2str(fed));

    disp('PASS test_firing_order');
end

function pats = local_patterns(n)
    pats = Pattern.empty(1, 0);
    for k = 1:n
        pats(k) = Pattern([10*k 10*k 0], [], true);
        pats(k).diffraction_efficiency = 0.5;
    end
end

function local_assign_ids(pats)
    for k = 1:numel(pats)
        pats(k).id = k;
    end
end

function out = local_ids(seq)
    out = seq.ids;
end

function fed = local_feed(order, n_holo)
%LOCAL_FEED  PlaySequence2K's feed logic, with its pre-loop range check.
%   Mirrors the real thing rather than calling it, since that needs SLM handles.
    fed = [];
    if max([order(:); 0]) > n_holo
        return                                   % blanks and returns
    end
    counter = 1;
    while counter <= numel(order)
        if counter <= numel(order) && n_holo >= order(counter)
            fed(end+1) = order(counter); %#ok<AGROW>
            counter = counter + 1;
        else
            break
        end
    end
end

function local_must_error(fn, id)
    try
        fn();
    catch err
        assert(strcmp(err.identifier, id), ...
            'expected %s, got %s (%s)', id, err.identifier, err.message);
        return
    end
    error('expected %s, but nothing was thrown', id);
end
