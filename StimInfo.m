classdef StimInfo < matlab.mixin.Copyable
    properties
        sequence
        power
        pulse_start
        pulse_duration
       
    end


    methods
        function obj = StimInfo(sequence, power, pulse_start, pulse_duration)
            if numel(power) > 2
                disp('Can''t do multiple powers rn')
                return
            end
            obj.sequence = sequence; % expects a sequence object
            % checks
            obj.pulse_start = obj.check_and_equalize(pulse_start);
            obj.pulse_duration = obj.check_and_equalize(pulse_duration);
            obj.power = power;
            obj.warn_if_slm_triggers_merge();
        end

        function warn_if_slm_triggers_merge(obj)
            %WARN_IF_SLM_TRIGGERS_MERGE  Pulses closer than the SLM trigger width.
            %
            %   SLMComm asks for a trigger of
            %   [pulse_start - pre_delay, + SLMComm.trigger_width] per pulse. So two
            %   pulses spaced dt apart give windows that ABUT at dt == trigger_width and
            %   OVERLAP below it -- pre_delay shifts both equally and cancels. Either
            %   way the digital line never returns low, so the SLM sees ONE rising edge
            %   instead of two.
            %
            %   That matters because PlaySequence2K advances one frame per trigger: with
            %   n pulses merged into one edge it feeds ONE hologram and the rest of the
            %   sequence slides out of step, so later pulses land on whatever pattern is
            %   still displayed. The laser gate is unaffected -- LaserGate writes every
            %   pulse window independently -- so the light still fires on time, at the
            %   WRONG targets. Nothing about that is visible in the saved stim record.
            %
            %   Measured with the real PulseGenerator back when trigger_width was
            %   0.025: 5 pulses at 0.025 s -> 1 edge; 10 at 0.020 s -> 1 edge; 5 at
            %   0.026 s -> 5 edges; 8 at 0.2 s -> 8. At the current 0.005 the same
            %   arithmetic applies one fifth of the way down.
            %
            %   Harmless in one case, which is why this warns rather than errors: if the
            %   merged pulses all point at the SAME pattern, the frame that stays on
            %   screen is the right one. Repeating one pattern_id for a burst of spikes
            %   on one target is legitimate. Distinct patterns closer than 25 ms are not.
            starts = sort(obj.pulse_start(:))';
            if numel(starts) < 2
                return
            end

            % The floor is the trigger WIDTH, not pre_delay. The trigger for a pulse at
            % t spans [t - pre_delay, t - pre_delay + width], so shifting both by
            % pre_delay cancels it: two pulses dt apart abut at dt == width. The two
            % constants were both 0.025 until the width dropped to 0.005, which is
            % exactly the sort of coincidence that hides a wrong threshold -- so read
            % the real one off SLMComm rather than restating it here.
            width = local_slm_trigger_width();

            % Compared with a tolerance, not exactly. The sweep is quantised --
            % Generator.to_samples is round(s * sample_rate) -- so a gap that exceeds
            % the width by less than a sample period still lands on contiguous samples
            % and merges. It also absorbs float noise: a nominal 25 ms gap written as
            % 0.525 - 0.5 is 0.025000000000000022 in double, which failed a bare <=
            % while the real PulseGenerator sweep at 20 kHz still produced ONE rising
            % edge (measured). One sample at 10 kHz is 100 us, so that is the bound.
            tol = 1e-4;
            gaps = diff(starts);
            merged = gaps <= width + tol;
            if ~any(merged)
                return
            end

            % Only a real problem when the merged neighbours are DIFFERENT patterns.
            ids = [];
            if ~isempty(obj.sequence) && ~isempty(obj.sequence.patterns)
                ids = reshape([obj.sequence.patterns.id], 1, []);
            end
            if numel(ids) == numel(starts)
                differs = ids(1:end-1) ~= ids(2:end);
                if ~any(merged & differs)
                    return                     % same pattern either side: harmless
                end
            end

            persistent warned
            if ~isempty(warned)
                return
            end
            warned = true;

            warning('StimInfo:slmTriggersMerge', ...
                ['%d of %d pulse gaps are <= %g s (min %g s), so their SLM triggers ' ...
                 'merge into a\nsingle rising edge -- SLMComm holds each trigger ' ...
                 '%g s. PlaySequence2K advances one\nframe per edge, so the sequence ' ...
                 'slides out of step and later pulses fire on the\nwrong pattern. The ' ...
                 'laser gate is NOT affected, so the light still looks correct.\n' ...
                 '  fix:  space distinct patterns more than %g s apart, or\n' ...
                 '  fix:  reuse ONE pattern_id for pulses meant to hit the same target.\n' ...
                 'Warned once per session; silence with ' ...
                 'warning(''off'', ''StimInfo:slmTriggersMerge'').'], ...
                sum(merged), numel(gaps), width, min(gaps), width, width);
        end

        function out = check_and_equalize(obj, to_check)
            % desired_length = obj.sequence.length;
            % if length(to_check) == 1
            %     out = repmat(to_check, 1, desired_length);
            %     return
            % end
            
            % if length(to_check) ~= obj.sequence.N
            %     error('Parameters much match sequence length (%d patterns)', obj.sequence.N)
            % end

            out = to_check;

        end

        function out = to_struct(obj)
            out = copy(obj);
            out.sequence.patterns = struct(out.sequence.patterns);
            out.sequence = struct(out.sequence);
            out = struct(out);
        end

        % function out = trial_length(obj)
        %     out = sum(obj.total_stimulation_time); % sum of all durations is the trial length
        % end
        % 
        % function out = N(obj)
        %     out = numel(obj.firing_order);
        % end
    end
end

function w = local_slm_trigger_width()
%LOCAL_SLM_TRIGGER_WIDTH  SLMComm.trigger_width, or the value it currently holds.
%   Read from the class when it is on the path -- which it is on any machine that can
%   actually run a stim -- so shortening the trigger there cannot leave this guard
%   checking a stale threshold. holography2k does not otherwise depend on holodaq, so
%   fall back rather than error when it is absent (offline analysis, tests).
%
%   Keep the fallback equal to SLMComm.trigger_width.
    persistent cached
    if ~isempty(cached)
        w = cached;
        return
    end
    w = 0.005;
    try
        if exist('SLMComm', 'class') == 8
            w = SLMComm.trigger_width;
        end
    catch
        % older SLMComm without the constant: its width was the hardcoded 0.025
        w = 0.025;
    end
    cached = w;
end