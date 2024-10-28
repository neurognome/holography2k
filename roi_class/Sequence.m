classdef Sequence < handle
    properties
        patterns
    end

    methods
        function obj = Sequence(patterns)
            obj.patterns = copy(patterns);
            % obj.equalize_patterns();
        end

        function equalize_patterns(obj)
            warning('Pattern equalization enabled, modulating power with SLM.')  
            for p = 1:numel(obj.patterns)
                obj.patterns(p).powerbias = obj.patterns(p).powerbias * (length(obj.patterns(p).powerbias)/obj.max_pattern_N);
            end
        end
        
        function out = ids(obj)            
            out = [obj.patterns.id];
        end

        function out = average_de(obj)
            % right now we can't really "flip" fast enough to change power
            % that carefully... so let's just choose an average lol
            out = sum([obj.patterns.diffraction_efficiency])/numel(obj.patterns);
        end

        function out = N(obj)
            out = numel(obj.patterns);
        end

        function out = max_pattern_N(obj)
                    out = max(arrayfun(@(x) size(x.targets, 1), obj.patterns));
        end
    end
end