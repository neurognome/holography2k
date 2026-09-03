function [hr, scale_center] = apply_holo_adjustments(hr)
%APPLY_HOLO_ADJUSTMENTS Prefer the measured SI->SLM offsets over the stamped ones.
%   [hr, scale_center] = APPLY_HOLO_ADJUSTMENTS(hr) looks up the shared offsets
%   store for this holoRequest's wavelength and, if there is a usable record,
%   REPLACES hr.xoffset / hr.yoffset / hr.scale with the measured values. It also
%   returns the pivot that scale should be taken about.
%
%   REPLACE, NEVER ADD. The ScanImage computer stamps the same request from the
%   same store when it builds it (get_holo_adjustments -> saveAllToHoloRequest),
%   so adding here would apply the correction twice. Only this function's values
%   ever reach a hologram, which is what keeps the two sides from compounding.
%
%   Falls back on every failure -- no store, no wavelength tag, a zoom mismatch,
%   holodaq missing from the path -- to whatever the request already carried,
%   which is the behaviour that predates this function. It never throws: a prime
%   must not fail because a drive is unmounted.
%
%   Every path prints what it decided, including the values it did NOT use.
%   "Which offsets did that session actually run with" is the question the whole
%   store exists to answer, and a value that arrives silently is the same problem
%   in a new place.
%
%   It lives in its own file rather than as a local function of
%   generate_holograms_new because that file's main function is not terminated
%   with `end`, so it cannot host one -- and because this is worth testing on its
%   own.
%
%   See also: generate_holograms_new, holo_adjustments, get_holo_adjustments

    LEGACY_CENTER = [.5, .5];   % what the code did before; see below

    % An explicit opt-out, for a request that IS the calibration. holeburn_grid
    % sets it: a measure pass must command exactly 0/0/1 so the burn reports the
    % raw mapping, and a verify pass must command the CANDIDATE values -- the
    % ones about to be written -- not whatever the store currently holds. Either
    % way, looking the store up here would measure the wrong thing.
    if isfield(hr, 'ignore_adjustments') && ~isempty(hr.ignore_adjustments) ...
            && logical(hr.ignore_adjustments)
        scale_center = LEGACY_CENTER;
        if isfield(hr, 'scale_center') && ~isempty(hr.scale_center)
            scale_center = reshape(double(hr.scale_center), 1, []);
        end
        fprintf(['Adjustments: this holoRequest asked for its own values ' ...
                 '(ignore_adjustments).\n']);
        fprintf('   applying             dx %+.2f  dy %+.2f  scale %.4f, pivot [%g %g]\n', ...
            hr.xoffset, hr.yoffset, hr.scale, scale_center(1), scale_center(2));
        return
    end

    wl = [];
    if isfield(hr, 'wavelength') && ~isempty(hr.wavelength)
        wl = double(hr.wavelength);
    end
    zoom = [];
    if isfield(hr, 'zoom') && ~isempty(hr.zoom)
        zoom = double(hr.zoom);
    end

    adj  = [];
    note = 'holoRequest carries no wavelength tag, so no record can be looked up';
    if ~isempty(wl)
        try
            [adj, note] = holo_adjustments(wl, zoom);
        catch err
            % holodaq missing from the path, K: unmounted, whatever it is: this
            % must never be the reason a prime fails. Fall back and say so.
            adj  = [];
            note = sprintf('adjustments lookup failed (%s)', err.message);
        end
    end

    fprintf('Adjustments: %s\n', note);
    fprintf('   holoRequest carries  dx %+.2f  dy %+.2f  scale %.4f\n', ...
        hr.xoffset, hr.yoffset, hr.scale);

    if isempty(adj)
        scale_center = LEGACY_CENTER;
        if isfield(hr, 'scale_center') && ~isempty(hr.scale_center)
            scale_center = reshape(double(hr.scale_center), 1, []);
        elseif hr.scale ~= 1
            % Only warn when it actually bites. Saying this on every unit-scale
            % request would train everyone to ignore it.
            warning('generate_holograms:legacyScaleCenter', ...
                ['scale = %.4f with no measured record, so the LEGACY pivot ' ...
                 '[0.5 0.5] is being\nused -- that is the frame CORNER in pixel ' ...
                 'coordinates, so this adds roughly\n%+.1f px of translation on ' ...
                 'top of the offsets. It is what this rig has always\ndone, and ' ...
                 'the offsets in use were tuned through it, so nothing is being ' ...
                 'changed\nhere. Run the hole burn calibration for this ' ...
                 'wavelength to replace both with\nmeasured values about the ' ...
                 'frame centre.'], hr.scale, (hr.scale - 1) * 256);
        end
        fprintf('   applying the holoRequest''s own values, pivot [%g %g]\n', ...
            scale_center(1), scale_center(2));
        return
    end

    hr.xoffset = adj.xoffset;
    hr.yoffset = adj.yoffset;
    hr.scale   = adj.scale;
    if isfield(adj, 'zMap') && ~isempty(adj.zMap)
        hr.zMap = adj.zMap;
    end
    scale_center = reshape(double(adj.scale_center), 1, []);
    hr.scale_center = scale_center;

    fprintf('   applying             dx %+.2f  dy %+.2f  scale %.4f, pivot [%g %g]\n', ...
        hr.xoffset, hr.yoffset, hr.scale, scale_center(1), scale_center(2));
end
