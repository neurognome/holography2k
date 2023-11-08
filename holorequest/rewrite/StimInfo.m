classdef StimInfo < handle
    properties
        firing_order

        hz
        power
        on_time
        duration
    end


    methods
        function obj = StimInfo(firing_order, hz, power, on_time, duration)
        obj.firing_order = firing_order;
        obj.hz = hz;
        obj.power = power;
        obj.on_time = on_time;
        obj.duration = duration;     
        end

        function out = trial_length(obj)
            out = sum(obj.duration); % sum of all durations is the trial length
        end

        function out = N(obj)
            out = numel(obj.firing_order);
        end
    end
end