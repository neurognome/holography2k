classdef StimInfo < handle
    properties
        firing_order

        hz
        power
        total_stimulation_time
        pulse_duration
    end


    methods
        function obj = StimInfo(firing_order, hz, power, total_stimulation_time, pulse_duration)
            if numel(power) > 2
                disp('Can''t do multiple powers rn')
                return
            end
            obj.firing_order = firing_order;
            % checks
            obj.hz = obj.check_and_equalize(hz);
            obj.power = power;
            obj.total_stimulation_time = obj.check_and_equalize(total_stimulation_time);
            obj.pulse_duration= obj.check_and_equalize(pulse_duration);  
        end

        function out = check_and_equalize(obj, to_check)
            desired_length = numel(obj.firing_order);
            if length(to_check) == 1
                out = repmat(to_check, 1, desired_length);
                return
            end
            
            if length(to_check) ~= desired_length
                error('Either give me something with 1 length or matching')
            end

            out = to_check;

        end

        function out = trial_length(obj)
            out = sum(obj.total_stimulation_time); % sum of all durations is the trial length
        end

        function out = N(obj)
            out = numel(obj.firing_order);
        end
    end
end