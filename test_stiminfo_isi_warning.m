function test_stiminfo_isi_warning()
%TEST_STIMINFO_ISI_WARNING  The merged-SLM-trigger guard fires when it should.
%
%   >> addpath(pwd); addpath(fullfile(pwd,'roi_class')); test_stiminfo_isi_warning
%
%   SLMComm requests [pulse_start - 0.025, +0.025] per pulse, so pulses spaced <= 25 ms
%   produce ONE rising edge instead of one each, and PlaySequence2K -- which advances a
%   frame per edge -- slides out of step. Verified against the real PulseGenerator:
%   5 pulses at 0.025 s give 1 edge, at 0.026 s give 5.

    a = [10 10 0];
    b = [20 20 0];
    c = [30 30 40];

    %% 1. distinct patterns closer than 25 ms -> warns ------------------------------
    id = local_capture(@() local_stim({a, b, c}, [0.20 0.22 0.24]));   % 20 ms
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        'distinct patterns at 20 ms must warn, got ''%s''', id);

    % exactly at the boundary: the windows abut, the line never drops
    id = local_capture(@() local_stim({a, b}, [0.50 0.525]));          % 25 ms
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        '25 ms is the boundary and must warn, got ''%s''', id);

    %% 2. just over the boundary -> silent -----------------------------------------
    id = local_capture(@() local_stim({a, b}, [0.50 0.526]));          % 26 ms
    assert(isempty(id), '26 ms clears the boundary (got ''%s'')', id);

    %% 3. comfortably spaced -> silent ---------------------------------------------
    id = local_capture(@() local_stim({a, b, c}, [0.2 0.4 0.6]));
    assert(isempty(id), 'a 200 ms burst must not warn (got ''%s'')', id);

    %% 4. a single pulse -> silent -------------------------------------------------
    id = local_capture(@() local_stim({a}, 0.5));
    assert(isempty(id), 'one pulse has no gap to check (got ''%s'')', id);

    %% 5. the SAME pattern repeated fast -> silent ---------------------------------
    % Legitimate: a burst of spikes on one target. The frame left on screen is the
    % right one, so a merged trigger changes nothing.
    id = local_capture(@() local_stim({a, a, a, a, a}, 0.5 + (0:4)*0.01));
    assert(isempty(id), ...
        'repeating one pattern is the documented safe case (got ''%s'')', id);

    %% 6. same pattern fast, then a DIFFERENT one still close -> warns -------------
    id = local_capture(@() local_stim({a, a, b}, [0.50 0.51 0.52]));
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        'a different pattern inside the merged run must still warn, got ''%s''', id);

    disp('PASS test_stiminfo_isi_warning');
end

function local_stim(target_list, starts)
%LOCAL_STIM  A StimInfo whose Sequence holds one Pattern per entry, ids 1..n.
    pats = Pattern.empty(1, 0);
    for k = 1:numel(target_list)
        pats(k) = Pattern(target_list{k}, [], true);
    end
    % ids are what the guard compares to decide whether neighbours share a pattern;
    % generate_holograms assigns them on the rig, so assign them here the same way --
    % identical targets get the SAME id, which is the case the guard forgives.
    [~, ~, uid] = unique(cellfun(@(t) mat2str(t), target_list, 'uni', 0));
    for k = 1:numel(pats)
        pats(k).id = uid(k);
        pats(k).diffraction_efficiency = 1;
    end
    StimInfo(Sequence(pats), 0.001, starts, 0.01 * ones(1, numel(starts)));
end

function id = local_capture(fn)
%LOCAL_CAPTURE  'StimInfo:slmTriggersMerge' if fn raises it, '' otherwise.
%   Clears the persistent guard, which is once-per-session by design. Silences the
%   unrelated warnings a Sequence construction emits so lastwarn reports the one under
%   test rather than whichever fired last.
    clear StimInfo
    lastwarn('', '');
    s1 = warning('off', 'StimInfo:slmTriggersMerge');
    s2 = warning('off', 'MATLAB:structOnObject');
    s3 = warning('off', 'Pattern:noZeroOrderDump');
    s4 = warning('off', 'Sequence:raggedWithoutZeroOrderDump');
    cleanup = onCleanup(@() [warning(s1), warning(s2), warning(s3), warning(s4)]);
    fn();
    [~, last_id] = lastwarn();
    if strcmp(last_id, 'StimInfo:slmTriggersMerge')
        id = last_id;
    else
        id = '';
    end
end
