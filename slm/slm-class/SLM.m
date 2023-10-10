classdef SLM < handle
properties
    Nx
    Ny
    psSLM
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

    function blank(obj)
        obj.feed(zeros(obj.Nx, obj.Ny));
    end

    function out =  get_slm_parameters(obj)
        props = properties(obj);
        out = struct();
        for p = props'
            out.(p{:}) = obj.(p{:});
        end
    end

    function Setup = add_slm(obj, Setup)
        Setup.SLM = obj.get_slm_parameters();
        Setup.psSLM = obj.psSLM;
        Setup.Nx = obj.Nx;
        Setup.Ny = obj.Ny;
        Setup.intensity = 1;
        Setup.source = sqrt(Setup.intensity)*(1/(Setup.Nx* Setup.Ny))*ones(Setup.Nx, Setup.Ny);
    end
end
end