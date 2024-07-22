classdef Pattern < handle
    properties
        diffraction_efficiency % de of this pattern
        targets % targets present in this pattern
        powerbias % bias of power to each target
        zero_order_dump
        id
    end

    methods
        function obj = Pattern(targets, powerbias, zero_order_dump)      
            obj.validate(targets);

            if nargin < 2 || isempty(powerbias)
                powerbias = ones(1, size(targets, 1));
            end

            if nargin < 3 || isempty(zero_order_dump)
                warning('Dumping power into zero order!!')
                zero_order_dump = true;
            end
            obj.zero_order_dump = zero_order_dump;
            obj.targets = targets;
            obj.powerbias = powerbias;
        end

        function validate(obj, targets)
            if size(targets, 2) ~= 3 % x y z
                error('Targets are not correct')
            end
            fprintf('Input %d targets\n', size(targets, 1))
        end

        function out = to_struct(obj)
            out = struct(obj);
        end

        function out = get.id(obj)
            out = obj.id;
        end

        function set.id(obj, id)
            obj.id = id;
        end

        function calculate_DE(obj, CoC, de_floor)
            if nargin < 3 || isempty(de_floor)
                de_floor = 0.05;
            end
            % make sure everything is pathed, or else you won't find this..
            slm_coords = function_SItoSLM(obj.targets, CoC);
            attenuation_coeffs = max(slm_coords(:, 4), de_floor);

            energy = 1./attenuation_coeffs;
            energy = energy.*obj.powerbias;
            energy = energy/sum(energy);
            obj.diffraction_efficiency = sum(energy.*attenuation_coeffs);
        end
    end

    methods (Static = true)
        function out = from_struct(input_struct)
            obj = Pattern(input_struct.targets, input_struct.powerbias);
            obj.diffraction_efficiency = input_struct.diffraction_efficiency;
            if isfield(input_struct, 'id')
                obj.id = input_struct.id;
            end
            out = obj;
        end
    end
end