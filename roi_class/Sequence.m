% classdef Sequence < handle
%     properties
%         patterns
%     end
% 
%     methods
%         function obj = Sequence(patterns)
%             obj.patterns = copy(patterns);
%             % obj.equalize_patterns();
%         end
% 
%         function equalize_patterns(obj)
%             warning('Pattern equalization enabled, modulating power with SLM.')  
%             for p = 1:numel(obj.patterns)
%                 obj.patterns(p).powerbias = obj.patterns(p).powerbias * (length(obj.patterns(p).powerbias)/obj.max_pattern_N);
%             end
%         end
% 
%         function out = ids(obj)            
%             out = [obj.patterns.id];
%         end
% 
%         function out = average_de(obj)
%             % right now we can't really "flip" fast enough to change power
%             % that carefully... so let's just choose an average lol
%             out = sum([obj.patterns.diffraction_efficiency])/numel(obj.patterns);
%         end
% 
%         function out = N(obj)
%             out = numel(obj.patterns);
%         end
% 
%         function out = max_pattern_N(obj)
%                     out = max(arrayfun(@(x) size(x.targets, 1), obj.patterns));
%         end
%     end
% end

classdef Sequence < handle
    properties
        patterns
    end

    methods
        function obj = Sequence(patterns)
            obj.patterns = arrayfun(@to_struct, patterns);
            % obj.equalize_patterns();
        end

        function equalize_patterns(obj)
            warning('Pattern equalization enabled, modulating power with SLM.')
            for p = 1:numel(obj.patterns)
                obj.patterns(p).powerbias = obj.patterns(p).powerbias * ...
                    (length(obj.patterns(p).powerbias) / obj.max_pattern_N);
            end
        end

        function out = ids(obj)
            out = [obj.patterns.id];
        end

        function out = calculate_power(obj, power_per_cell)
            sz = zeros(numel(obj.patterns), 1);

            for k = 1:numel(obj.patterns)
                sz(k) = size(obj.patterns(k).targets, 1);
            end

            de = obj.average_de;
            % Safety guard: diffraction_efficiency is only populated after the
            % holoRequest is round-tripped through the SLM (transferHR/calculate_DE).
            % If it is missing, average_de is 0/NaN and out -> Inf, which
            % LaserPowerControl then CLAMPS to the laser's MAX power. Fail loudly
            % here instead of silently commanding max power.
            if isempty(de) || ~isfinite(de) || de <= 0
                error('Sequence:calculate_power:badDE', ...
                    ['average_de = %g -- patterns have no valid diffraction_efficiency. ', ...
                     'Run transferHR (or calculate_DE) before calculate_power.'], de);
            end

            out = (power_per_cell * max(sz)) / de;
        end

        function out = average_de(obj)
            out = sum([obj.patterns.diffraction_efficiency]) / numel(obj.patterns);
        end

        function out = N(obj)
            out = numel(obj.patterns);
        end

        function out = max_pattern_N(obj)
            out = max(arrayfun(@(x) size(x.targets, 1), obj.patterns));
        end
    end
end