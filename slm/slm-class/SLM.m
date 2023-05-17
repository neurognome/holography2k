classdef SLM < handle
properties
    % lut_file
    % reg_lut
    % true_frames
    % pixelmax
end

methods
    function obj = SLM()
    end

    function start(obj)
    end

    function stop(obj)
    end

    function feed(obj)
    end

    function out =  get_slm(obj)
        props = properties(obj);
        out = struct();
        for p = props
            out.(p{:}) = obj.(p{:});
        end
    end
end
end