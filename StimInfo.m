classdef StimInfo < handle
    properties
        sequence
        power
        pulse_start
        pulse_duration
       
    end


    methods
        function obj = StimInfo(sequence, power, pulse_start, pulse_duration)
            if numel(power) > 2
                disp('Can''t do multiple powers rn')
                return
            end
            obj.sequence = sequence; % expects a sequence object
            % checks
            obj.pulse_start = obj.check_and_equalize(pulse_start);
            obj.pulse_duration = obj.check_and_equalize(pulse_duration);
            obj.power = power;
        end

        function out = check_and_equalize(obj, to_check)
            % desired_length = obj.sequence.length;
            % if length(to_check) == 1
            %     out = repmat(to_check, 1, desired_length);
            %     return
            % end
            
            % if length(to_check) ~= obj.sequence.N
            %     error('Parameters much match sequence length (%d patterns)', obj.sequence.N)
            % end

            out = to_check;

        end

        % function out = trial_length(obj)
        %     out = sum(obj.total_stimulation_time); % sum of all durations is the trial length
        % end
        % 
        % function out = N(obj)
        %     out = numel(obj.firing_order);
        % end
    end
end