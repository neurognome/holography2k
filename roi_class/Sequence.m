% classdef Sequence < handle
%     properties
%         patterns
%     end
% 
%     methods
%         function obj = Sequence(patterns)
%             obj.patterns = copy(patterns);
%             % obj.equalize_patterns();
%         end
% 
%         function equalize_patterns(obj)
%             warning('Pattern equalization enabled, modulating power with SLM.')  
%             for p = 1:numel(obj.patterns)
%                 obj.patterns(p).powerbias = obj.patterns(p).powerbias * (length(obj.patterns(p).powerbias)/obj.max_pattern_N);
%             end
%         end
% 
%         function out = ids(obj)            
%             out = [obj.patterns.id];
%         end
% 
%         function out = average_de(obj)
%             % right now we can't really "flip" fast enough to change power
%             % that carefully... so let's just choose an average lol
%             out = sum([obj.patterns.diffraction_efficiency])/numel(obj.patterns);
%         end
% 
%         function out = N(obj)
%             out = numel(obj.patterns);
%         end
% 
%         function out = max_pattern_N(obj)
%                     out = max(arrayfun(@(x) size(x.targets, 1), obj.patterns));
%         end
%     end
% end

classdef Sequence < handle
    properties
        patterns
    end

    methods
        function obj = Sequence(patterns)
            obj.patterns = arrayfun(@to_struct, patterns);
            obj.warn_if_ragged_without_dump();
            % obj.equalize_patterns();
        end

        function warn_if_ragged_without_dump(obj)
            %WARN_IF_RAGGED_WITHOUT_DUMP  Ragged multi-pattern sequence, dump off = overdrive.
            %
            %   ONE laser power is commanded for a whole sequence: calculate_power below
            %   returns a single scalar and LaserPowerControl.set applies it once per
            %   trial. So every pattern in the sequence shares it, and a pattern with
            %   FEWER targets than the largest concentrates that fixed power into fewer
            %   spots -- overdriving each of them by max_size/this_size.
            %
            %   Pattern.zero_order_dump is what corrects it: generate_holograms_new.m
            %   parks the unrequested remainder in the zero order "so the absolute laser
            %   power can stay fixed across differently-sized patterns". It DEFAULTS TO
            %   FALSE (see the Pattern constructor), so a ragged sequence built without
            %   asking for it is silently mis-dosed.
            %
            %   Only a warning: it is legitimate to run this way if the zero order is not
            %   blocked and you would rather not park power there, or if the sizes differ
            %   so little that you do not care. But it should never be a surprise.
            %
            %   Not an issue BETWEEN sequences -- power_per_cell is per pool entry, so
            %   each trial can command its own power. Uniform sizes are also fine.
            if numel(obj.patterns) < 2
                return
            end

            sz = arrayfun(@(x) size(x.targets, 1), obj.patterns);
            if all(sz == sz(1))
                return                      % uniform: nothing to correct
            end

            dumped = arrayfun(@(x) isscalar(x.zero_order_dump) && x.zero_order_dump, ...
                obj.patterns);
            small = sz < max(sz);           % only patterns below the max are affected
            if all(dumped(small))
                return
            end

            % Once per session: build_opto_stims constructs a Sequence per trial per
            % channel, so an unguarded warning would print thousands of times.
            persistent warned
            if ~isempty(warned)
                return
            end
            warned = true;

            bad = find(small & ~dumped);
            warning('Sequence:raggedWithoutZeroOrderDump', ...
                ['This sequence mixes ensemble sizes %s and shares ONE laser power, but ' ...
                 '%d of its\npatterns have zero_order_dump = FALSE (the Pattern default). ' ...
                 'Those patterns will be\nOVERDRIVEN by up to %.1fx per cell, because the ' ...
                 'fixed trial power lands on fewer\ntargets.\n' ...
                 '  fix:  Pattern(targets, powerbias, true)   %% dump the remainder to zero order\n' ...
                 '  or:   keep every pattern in a sequence the same size,\n' ...
                 '  or:   one pattern per trial and vary power_per_cell between trials.\n' ...
                 'Warned once per session; silence with ' ...
                 'warning(''off'', ''Sequence:raggedWithoutZeroOrderDump'').'], ...
                mat2str(sz), numel(bad), max(sz) / min(sz(small)));
        end

        function equalize_patterns(obj)
            warning('Pattern equalization enabled, modulating power with SLM.')
            for p = 1:numel(obj.patterns)
                obj.patterns(p).powerbias = obj.patterns(p).powerbias * ...
                    (length(obj.patterns(p).powerbias) / obj.max_pattern_N);
            end
        end

        function out = ids(obj)
            %IDS  The firing order for this sequence: one pattern id per pulse.
            %
            %   Checked rather than trusted, because [obj.patterns.id] SILENTLY DROPS
            %   empty ids -- [1 [] 3] concatenates to [1 3]. Ids are assigned by
            %   generate_holograms_new on the holo computer and come back through
            %   transferHR, so an empty one means the holoRequest was never
            %   round-tripped. The consequences were both invisible:
            %
            %     all empty  -> ids is [], PlaySequence2K's loop never runs, nothing is
            %                   fed, and the SLM holds the PREVIOUS trial's hologram
            %                   while the laser gate fires this trial's pulses.
            %     some empty -> a short order that no longer lines up with
            %                   delay/duration, so the sequence plays offset.
            %
            %   Either way the light looks right and the targets are wrong, with
            %   nothing in the saved stim record to show it. Fail here instead.
            raw = {obj.patterns.id};
            missing = cellfun(@isempty, raw);
            assert(~any(missing), 'Sequence:missingPatternIds', ...
                ['%d of %d patterns have no id, so the firing order would be sent ' ...
                 'SHORT (or\nempty) and the SLM would play the wrong holograms. Ids ' ...
                 'are assigned by\ngenerate_holograms_new and returned by transferHR ' ...
                 '-- an empty one means this\nholoRequest was never round-tripped ' ...
                 'through the holography computer.'], ...
                sum(missing), numel(raw));

            bad = ~cellfun(@(x) isscalar(x) && isnumeric(x), raw);
            assert(~any(bad), 'Sequence:badPatternIds', ...
                ['%d of %d pattern ids are not numeric scalars. [patterns.id] would ' ...
                 'flatten\nthem into a firing order of the wrong length.'], ...
                sum(bad), numel(raw));

            out = [obj.patterns.id];
        end

        function out = calculate_power(obj, power_per_cell)
            sz = zeros(numel(obj.patterns), 1);

            for k = 1:numel(obj.patterns)
                sz(k) = size(obj.patterns(k).targets, 1);
            end

            de = obj.average_de;
            % Safety guard: diffraction_efficiency is only populated after the
            % holoRequest is round-tripped through the SLM (transferHR/calculate_DE).
            % If it is missing, average_de is 0/NaN and out -> Inf, which
            % LaserPowerControl then CLAMPS to the laser's MAX power. Fail loudly
            % here instead of silently commanding max power.
            if isempty(de) || ~isfinite(de) || de <= 0
                error('Sequence:calculate_power:badDE', ...
                    ['average_de = %g -- patterns have no valid diffraction_efficiency. ', ...
                     'Run transferHR (or calculate_DE) before calculate_power.'], de);
            end

            out = (power_per_cell * max(sz)) / de;
        end

        function out = average_de(obj)
            out = sum([obj.patterns.diffraction_efficiency]) / numel(obj.patterns);
        end

        function out = N(obj)
            out = numel(obj.patterns);
        end

        function out = max_pattern_N(obj)
            out = max(arrayfun(@(x) size(x.targets, 1), obj.patterns));
        end
    end
end