classdef Pattern < matlab.mixin.Copyable
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
                % DEFAULT IS FALSE, and that is only safe when every pattern in a
                % sequence is the same size: one laser power is commanded per sequence,
                % so a smaller pattern would otherwise concentrate it into fewer targets.
                % Sequence checks the real condition (ragged AND more than one pattern)
                % and warns with the actual overdrive factor -- this one only records
                % that the default was taken.
                %
                % Warned once per session and given an identifier. It used to be a bare
                % warning() on every construction, which printed hundreds of times per
                % experiment build, could not be suppressed selectively, and said
                % nothing about when it actually matters -- so it read as noise.
                persistent told
                if isempty(told)
                    told = true;
                    warning('Pattern:noZeroOrderDump', ...
                        ['zero_order_dump defaults to FALSE. Fine when every pattern in ' ...
                         'a sequence is the\nsame size; pass Pattern(targets, powerbias, ' ...
                         'true) for a RAGGED sequence, or the\nsmaller patterns are ' ...
                         'overdriven. Warned once per session.']);
                end
                zero_order_dump = false;
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

            % powerbias(:) forces a COLUMN, matching slm_coords(:,4) -- the same guard
            % generate_holograms_new.m already applies when it builds the hologram.
            % Without it, a row powerbias (which is what calibrated_powerbias,
            % channel_powerbias and the Pattern constructor's default all return) met a
            % column `energy` and implicitly expanded into an Nsomething x N OUTER
            % PRODUCT; the next line's `/sum(...)` then became mrdivide (a least-squares
            % solve) rather than elementwise. The result was numerically the DE for a
            % FLAT bias -- i.e. powerbias was silently ignored here, while the hologram
            % itself applied it correctly. Sequence.calculate_power divides the laser
            % command by this DE, so calibrated patterns were under-driven: measured
            % 1.96% at X=1.5 and 3.36% at X=2 on a 5-cell ensemble with attenuation
            % spread 0.9..0.2. Flat-bias patterns were unaffected, which is why every
            % uncalibrated experiment looked correct.
            assert(numel(obj.powerbias) == size(obj.targets, 1), ...
                'Pattern:powerbiasSize', ...
                ['powerbias has %d entries for %d targets. They must match -- a ' ...
                 'mismatch used to expand into an outer product and silently return a ' ...
                 'meaningless diffraction efficiency.'], ...
                numel(obj.powerbias), size(obj.targets, 1));

            energy = 1./attenuation_coeffs;
            energy = energy.*obj.powerbias(:);
            energy = energy./sum(energy);
            obj.diffraction_efficiency = min(sum(energy.*attenuation_coeffs), 1);
        end
    end

    methods (Static = true)
        function out = from_struct(input_struct)
            obj = Pattern(input_struct.targets, input_struct.powerbias);
            obj.diffraction_efficiency = input_struct.diffraction_efficiency;
            if isfield(input_struct, 'id')
                obj.id = input_struct.id;
            end
            obj.zero_order_dump = input_struct.zero_order_dump;
            out = obj;
        end
    end
end