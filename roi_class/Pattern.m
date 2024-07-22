classdef Pattern < handle
    properties
        diffraction_efficiency % de of this pattern
        targets % targets present in this pattern
        powerbias % bias of power to each target
    end

    properties (Access = protected)
        id % unique id for proper indexing of patterns across computers
    end

    methods
        function obj = Pattern(targets, powerbias)
            obj.validate(targets);

            if nargin < 2 || isempty(powerbias)
                powerbias = ones(1, size(targets, 1));
            end
            
            obj.targets = targets;
            obj.powerbias = powerbias;
        end

        function validate(obj, targets)
            if size(targets, 2) ~= 3 % x y z
                error('Targets are not correct')
            end
            fprintf('Input %d targets\n', size(targets, 1))
        end

        function out = get.id(obj)
            out = obj.id;
        end

        function set.id(obj, id)
            obj.id = id;
        end
    end
end