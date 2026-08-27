classdef SLM < handle
properties
    Nx
    Ny
    psSLM
    % Static wavefront correction for THIS panel, in radians, Nx-by-Ny, laid out
    % (x, y) like a hologram. Empty means none. add_slm copies every property into
    % Setup.SLM, so setting this on the object is what puts Setup.SLM.correction in
    % front of phase_to_frame. It lives per-SLM rather than in Setup because two
    % boards on one machine (900 and 1100 nm) have different aberrations, and a
    % single shared Setup would silently give both the same correction.
    correction = []
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