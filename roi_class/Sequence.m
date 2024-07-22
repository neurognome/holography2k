classdef Sequence < handle
    properties
        patterns
    end

    methods
        function obj = Sequence(patterns)
            obj.patterns = patterns;
        end

        function equalize_patterns(obj)
        % this is a bandaid fix because of lack of power modulation
            warning('Pattern equalization enabled, modulating power with SLM.')
            
        end
        
        function out = ids(obj)            
            out = [obj.patterns.id];
        end

        function out = average_de(obj)
            % right now we can't really "flip" fast enough to change power
            % that carefully... so let's just choose an average lol
            out = sum([obj.patterns.diffraction_efficiency])/numel(obj.patterns);
        end
    end
end