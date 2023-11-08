classdef Holomaker2 < handle
    properties
        holo_request
    end

    methods
        function obj = Holomaker2(holo_request, slm_trigger, shutter_trigger, power_)
            obj.holo_request = holo_request;
        end

        function set_rois(obj, rois)
            if ~iscell(rois)
                disp('rois need to be a cell array')
                return
            end
            obj.holo_request.rois = rois;
        end

    end
end