function test_stiminfo_isi_warning()
%TEST_STIMINFO_ISI_WARNING  The merged-SLM-trigger guard fires when it should.
%
%   >> addpath(pwd); addpath(fullfile(pwd,'roi_class')); test_stiminfo_isi_warning
%
%   SLMComm requests [pulse_start - pre_delay, + trigger_width] per pulse, so pulses
%   spaced <= trigger_width produce ONE rising edge instead of one each, and
%   PlaySequence2K -- which advances a frame per edge -- slides out of step.
%
%   The threshold is the WIDTH, not pre_delay: pre_delay shifts both windows equally
%   and cancels. They were both 0.025 until the width dropped to 0.005, so the cases
%   below are written against whatever SLMComm currently declares rather than a literal.

    w = local_width();
    fprintf('trigger width in force: %g s\n', w);

    a = [10 10 0];
    b = [20 20 0];
    c = [30 30 40];

    tight = w * 0.8;        % inside the window -> merges
    edge  = w;              % exactly abutting -> merges
    clear_gap = w * 2;      % comfortably apart -> distinct edges

    %% 1. distinct patterns closer than the width -> warns --------------------------
    id = local_capture(@() local_stim({a, b, c}, 0.2 + (0:2)*tight));
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        'distinct patterns at %g s must warn, got ''%s''', tight, id);

    % exactly at the boundary: the windows abut, the line never drops
    id = local_capture(@() local_stim({a, b}, [0.5, 0.5 + edge]));
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        '%g s is the boundary and must warn, got ''%s''', edge, id);

    %% 2. clear of the boundary -> silent ------------------------------------------
    id = local_capture(@() local_stim({a, b}, [0.5, 0.5 + clear_gap]));
    assert(isempty(id), '%g s clears the boundary (got ''%s'')', clear_gap, id);

    id = local_capture(@() local_stim({a, b, c}, [0.2 0.4 0.6]));
    assert(isempty(id), 'a 200 ms burst must not warn (got ''%s'')', id);

    %% 3. a single pulse -> silent -------------------------------------------------
    id = local_capture(@() local_stim({a}, 0.5));
    assert(isempty(id), 'one pulse has no gap to check (got ''%s'')', id);

    %% 4. the SAME pattern repeated fast -> silent ---------------------------------
    % Legitimate: a burst of spikes on one target. The frame left on screen is the
    % right one, so a merged trigger changes nothing.
    id = local_capture(@() local_stim({a, a, a, a, a}, 0.5 + (0:4)*tight));
    assert(isempty(id), ...
        'repeating one pattern is the documented safe case (got ''%s'')', id);

    %% 5. same pattern fast, then a DIFFERENT one still close -> warns -------------
    id = local_capture(@() local_stim({a, a, b}, 0.5 + (0:2)*tight));
    assert(strcmp(id, 'StimInfo:slmTriggersMerge'), ...
        'a different pattern inside the merged run must still warn, got ''%s''', id);

    %% 6. the guard tracks SLMComm rather than a literal ---------------------------
    % If SLMComm is reachable its constant must be what the guard used; a stale literal
    % here is exactly the drift the read-through exists to prevent.
    if exist('SLMComm', 'class') == 8
        assert(abs(w - SLMComm.trigger_width) < eps, ...
            'guard threshold %g does not match SLMComm.trigger_width %g', ...
            w, SLMComm.trigger_width);
        fprintf('guard threshold matches SLMComm.trigger_width\n');
    else
        fprintf('SLMComm not on the path; guard used its fallback\n');
    end

    disp('PASS test_stiminfo_isi_warning');
end

function w = local_width()
    w = 0.005;
    if exist('SLMComm', 'class') == 8
        try
            w = SLMComm.trigger_width;
        catch
            w = 0.025;
        end
    end
end

function local_stim(target_list, starts)
%LOCAL_STIM  A StimInfo whose Sequence holds one Pattern per entry.
    pats = Pattern.empty(1, 0);
    for k = 1:numel(target_list)
        pats(k) = Pattern(target_list{k}, [], true);
    end
    % ids decide whether the guard treats neighbours as the same pattern; on the rig
    % generate_holograms assigns them, so mirror that -- identical targets share an id,
    % which is the case the guard forgives.
    [~, ~, uid] = unique(cellfun(@(t) mat2str(t), target_list, 'uni', 0));
    for k = 1:numel(pats)
        pats(k).id = uid(k);
        pats(k).diffraction_efficiency = 1;
    end
    StimInfo(Sequence(pats), 0.001, starts, 0.001 * ones(1, numel(starts)));
end

function id = local_capture(fn)
%LOCAL_CAPTURE  'StimInfo:slmTriggersMerge' if fn raises it, '' otherwise.
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
