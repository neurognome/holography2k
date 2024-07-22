classdef Sequence < handle
    properties
        patterns
    end

    methods
        function obj = Sequence(patterns)
            obj.patterns = patterns;
        end
        
        function out = generate_sequence(obj)            
            out = zeros(1, length(obj.patterns));
            ct = 1;
            for p = obj.patterns
                out(ct) = p.id;
                ct = ct + 1;
            end
        end
    end
end