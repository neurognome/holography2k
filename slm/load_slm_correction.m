function [base, info] = load_slm_correction(src, varargin)
%LOAD_SLM_CORRECTION Load a static SLM wavefront correction as a phase map in radians.
%   base = LOAD_SLM_CORRECTION(src, 'Size', [Nx Ny]) returns an Nx-by-Ny array of
%   RADIANS ready to assign to Setup.SLM.correction, which phase_to_frame then adds
%   to every hologram before wrapping (AO bead protocol, Phase 4 item 18).
%
%   SRC may be:
%     * a path to an image  (.bmp .png .tif .tiff) -- the Meadowlark
%       wavefront-correction form, one 8-bit count per pixel;
%     * a path to a .mat    -- either a single numeric variable holding a phase map,
%       or a struct with a 'coefficients' field (see below);
%     * a numeric matrix    -- used as-is, after units/orientation handling;
%     * a struct with a 'coefficients' field (+ optional 'modes', 'convention',
%       'normalized', 'pupil'), rendered by zernike_phase.
%
%   Name/value options:
%     'Size'       [Nx Ny] of the panel. Required for a coefficient source, and
%                  checked against the data for every other source.
%     'Units'      'auto' (default) | 'dn' | 'radians' | 'waves'.
%     'PixelMax'   255 (default). The DN that spans a full 2*pi -- Setup.SLM.pixelmax.
%     'Sign'       +1 (default) | -1. Phase 4 item 19 says to test both signs and
%                  keep whichever improves the spot, because the convention between
%                  pupil-phase estimation and SLM display is genuinely ambiguous.
%     'Transpose'  For an image source the DEFAULT IS TRUE; see ORIENTATION below.
%     'FlipX'      false (default). Flip along SLM x (dimension 1) after transpose.
%     'FlipY'      false (default). Flip along SLM y (dimension 2) after transpose.
%     'Verbose'    true (default). Print exactly what was applied.
%
%   ORIENTATION -- read this before trusting a result.
%   imread returns an image as (row, column) = (y, x). A hologram in this repo is
%   indexed (x, y): function_Make_3D_SHOT_Holos builds Masks(sxx, syy) with x first,
%   and that array is handed to Write_image unchanged. An image file therefore needs
%   a TRANSPOSE, which is the default for image sources here.
%
%   On a square panel -- and the scope2k Meadowlark is 1024x1024 -- NO size check
%   can catch a wrong transpose or flip. There is no way to settle it from the file,
%   so this function prints what it did and you settle the rest on the rig: feed a
%   known tilt as the correction and confirm the spot moves the predicted direction,
%   which is the same measurement Phase 4 item 19 already requires for the sign.
%
%   UNITS. 'auto' resolves as:
%     * integer-valued data inside [0, PixelMax]        -> 'dn'
%     * anything reaching beyond +-3.2 rad (~half turn) -> 'radians'
%     * otherwise                                       -> ERROR, because a small
%       correction is numerically indistinguishable in radians and in waves (a
%       0.19-wave map peaks near 0.7 in waves and 4.4 in radians, but a 0.1-wave one
%       is ambiguous), and guessing wrong scales the whole correction by 2*pi.
%
%   DN converts as phase = DN / PixelMax * 2*pi. That is the inverse of the forward
%   map in phase_to_frame up to a constant piston, which is dropped: a uniform phase
%   offset cannot change a far-field intensity.
%
%   Example -- a Meadowlark-style .bmp on the 900 nm board:
%       slm  = get_slm(900);
%       Setup.SLM.correction = load_slm_correction( ...
%           'C:\slm\corrections\ao_900nm.bmp', 'Size', [slm.Nx slm.Ny]);
%
%   See also: phase_to_frame, zernike_phase
%
%   Existing static correction on this rig is the Meadowlark hardware LUT, which is
%   a per-pixel voltage response curve and NOT a wavefront correction; the two are
%   independent and both apply.

    p = inputParser;
    p.FunctionName = 'load_slm_correction';
    p.addParameter('Size', []);
    p.addParameter('Units', 'auto');
    p.addParameter('PixelMax', 255);
    p.addParameter('Sign', 1);
    p.addParameter('Transpose', []);     % [] = per-source default
    p.addParameter('FlipX', false);
    p.addParameter('FlipY', false);
    p.addParameter('Verbose', true);
    p.parse(varargin{:});
    o = p.Results;

    assert(ismember(o.Sign, [1 -1]), 'load_slm_correction:badSign', ...
        'Sign must be +1 or -1, got %s.', mat2str(o.Sign));
    assert(isnumeric(o.PixelMax) && isscalar(o.PixelMax) && o.PixelMax > 0, ...
        'load_slm_correction:badPixelMax', 'PixelMax must be a positive scalar.');

    info = struct('source', '', 'kind', '', 'units', '', 'sign', o.Sign, ...
                  'transpose', false, 'flipx', logical(o.FlipX), ...
                  'flipy', logical(o.FlipY), 'zernike', []);

    % ---- 1. get raw numbers, and a default orientation for that source ----------
    [raw, kind, srcname, tr_default] = local_read(src, o);
    info.source = srcname;
    info.kind   = kind;

    if isempty(o.Transpose)
        do_transpose = tr_default;
    else
        do_transpose = logical(o.Transpose);
    end

    % ---- 2. units -> radians ----------------------------------------------------
    if strcmp(kind, 'zernike')
        % zernike_phase already returns radians on the requested grid, and it built
        % the grid in (x, y), so no transpose is owed.
        units = 'radians';
        do_transpose = false;
        phase = raw;
    else
        units = local_units(raw, o);
        switch units
            case 'dn',      phase = double(raw) / o.PixelMax * 2*pi;
            case 'radians', phase = double(raw);
            case 'waves',   phase = double(raw) * 2*pi;
        end
    end
    info.units = units;

    % ---- 3. orientation and sign ------------------------------------------------
    if do_transpose, phase = phase.'; end
    if o.FlipX,      phase = flip(phase, 1); end
    if o.FlipY,      phase = flip(phase, 2); end
    info.transpose = do_transpose;

    phase = o.Sign * phase;

    % A piston is unobservable in the far field, so remove it. This also makes the
    % reported RMS/PV comparable between a DN file (which sits around +pi) and a
    % zero-mean coefficient rendering.
    phase = phase - mean(phase(:));

    % ---- 4. validate against the panel -----------------------------------------
    assert(ismatrix(phase) && ~isempty(phase), 'load_slm_correction:notAMatrix', ...
        'The correction must be a non-empty 2-D array, got %s.', local_dims(phase));
    assert(all(isfinite(phase(:))), 'load_slm_correction:nonFinite', ...
        ['The correction contains %d non-finite value(s). phase_to_frame refuses ' ...
         'those\nrather than let uint8(NaN) == 0 become a silent flat spot.'], ...
        sum(~isfinite(phase(:))));

    if ~isempty(o.Size)
        want = double(o.Size(:))';
        assert(numel(want) == 2, 'load_slm_correction:badSize', ...
            'Size must be [Nx Ny], got %s.', mat2str(o.Size));
        assert(isequal(size(phase), want), 'load_slm_correction:sizeMismatch', ...
            ['The correction is %s but the panel is %dx%d.\nIf the two are ' ...
             'transposes of each other, set ''Transpose'' explicitly -- note that ' ...
             'an\nimage source already transposes by default, so passing an ' ...
             'already-(x,y) array\nfrom a file needs ''Transpose'', false.'], ...
            local_dims(phase), want(1), want(2));
    end

    base = phase;

    info.rms_waves = std(base(:)) / (2*pi);
    info.pv_waves  = (max(base(:)) - min(base(:))) / (2*pi);
    info.summary   = sprintf(['%s: %s, %s -> radians, sign %+d, transpose %d, ' ...
                              'flip [%d %d], RMS %.3f waves, PV %.3f waves'], ...
        srcname, local_dims(base), units, o.Sign, do_transpose, ...
        o.FlipX, o.FlipY, info.rms_waves, info.pv_waves);

    % A DN file is stored ALREADY WRAPPED, so the RMS and PV above describe the
    % wrapped map rather than the wavefront. That is harmless for the correction --
    % adding the wrapped phase is identical modulo 2*pi to adding the unwrapped one
    % -- but the numbers must not be read as a wavefront error. For a map wrapping
    % more than once they say almost nothing: PV tends to one wave and RMS to the
    % 0.289 waves of a uniform distribution over a turn, whatever the aberration is.
    info.wrapped = strcmp(units, 'dn');

    if o.Verbose
        fprintf('SLM correction loaded -- %s\n', info.summary);
        if info.wrapped
            fprintf(['  NOTE these are the WRAPPED map''s RMS/PV, not the ' ...
                     'wavefront''s: a DN file arrives\n  already wrapped. Unwrap ' ...
                     'it before quoting an aberration magnitude.\n']);
        end
        if ~strcmp(kind, 'zernike')
            fprintf(['  orientation is NOT verifiable from the file on a square ' ...
                     'panel. Confirm on the rig\n  (feed a known tilt; see Phase 4 ' ...
                     'item 19, which already asks you to test both signs).\n']);
        end
    end
end

% -------------------------------------------------------------------------
function [raw, kind, srcname, tr_default] = local_read(src, o)
%LOCAL_READ Resolve SRC to raw numbers plus the natural transpose for that source.

    if isstruct(src)
        raw = local_zernike(src, o);
        kind = 'zernike'; srcname = 'zernike coefficients'; tr_default = false;
        return
    end

    if isnumeric(src)
        raw = src;
        % A matrix handed in from the workspace is assumed to already be in the
        % hologram's (x, y) layout -- the caller had to choose one, and this is the
        % layout everything else in this file talks about.
        kind = 'array'; srcname = 'numeric array'; tr_default = false;
        return
    end

    assert(ischar(src) || isstring(src), 'load_slm_correction:badSource', ...
        ['SRC must be a file path, a numeric matrix, or a struct with a ' ...
         '''coefficients'' field;\ngot %s.'], class(src));
    f = char(src);
    assert(exist(f, 'file') == 2, 'load_slm_correction:noFile', ...
        'No such correction file:\n  %s', f);
    srcname = f;
    [~, ~, ext] = fileparts(f);

    switch lower(ext)
        case {'.bmp', '.png', '.tif', '.tiff', '.jpg', '.jpeg'}
            a = imread(f);
            if ndims(a) == 3
                % An RGB export of a grayscale map is common. Accept it only when
                % the channels agree, rather than quietly keeping the red one.
                same = all(a(:,:,1) == a(:,:,2), 'all') && all(a(:,:,1) == a(:,:,3), 'all');
                assert(same, 'load_slm_correction:colourImage', ...
                    ['%s has %d colour channels that do NOT agree, so it is not a ' ...
                     'grayscale\nphase map saved as RGB. Export it as single-channel ' ...
                     'grayscale.'], f, size(a, 3));
                a = a(:,:,1);
            end
            raw = a;
            kind = 'image';
            tr_default = true;   % imread gives (y, x); the hologram is (x, y)

        case '.mat'
            s = load(f);
            fn = fieldnames(s);
            if numel(fn) == 1 && isstruct(s.(fn{1})) && isfield(s.(fn{1}), 'coefficients')
                raw = local_zernike(s.(fn{1}), o);
                kind = 'zernike'; tr_default = false;
                return
            end
            num = fn(cellfun(@(k) isnumeric(s.(k)) && ~isvector(s.(k)), fn));
            assert(numel(num) == 1, 'load_slm_correction:ambiguousMat', ...
                ['%s must hold exactly one 2-D numeric variable (or one struct ' ...
                 'with a\n''coefficients'' field); it holds %d candidate(s): %s.'], ...
                f, numel(num), strjoin(fn', ', '));
            raw = s.(num{1});
            kind = 'mat';
            % A .mat written by MATLAB code in this repo is already (x, y).
            tr_default = false;

        otherwise
            error('load_slm_correction:badExt', ...
                ['Unsupported correction file type ''%s'' (%s).\nUse an image ' ...
                 '(.bmp .png .tif), a .mat, or pass Zernike coefficients.'], ext, f);
    end
end

% -------------------------------------------------------------------------
function phase = local_zernike(s, o)
%LOCAL_ZERNIKE Render a coefficient struct via zernike_phase.
    assert(isfield(s, 'coefficients'), 'load_slm_correction:noCoefficients', ...
        'A struct source needs a ''coefficients'' field.');
    assert(~isempty(o.Size) && numel(o.Size) == 2, 'load_slm_correction:sizeRequired', ...
        ['Rendering Zernike coefficients needs the panel size: pass ' ...
         '''Size'', [Nx Ny].']);

    args = {};
    if isfield(s, 'modes'),      args = [args, {'Modes',      s.modes}];      end
    if isfield(s, 'convention'), args = [args, {'Convention', s.convention}]; end
    if isfield(s, 'normalized'), args = [args, {'Normalized', s.normalized}]; end
    if isfield(s, 'pupil'),      args = [args, {'Pupil',      s.pupil}];      end

    [phase, zinfo] = zernike_phase(double(o.Size(1)), double(o.Size(2)), ...
                                   s.coefficients, args{:});
    if o.Verbose
        fprintf('  zernike: %s\n', zinfo.summary);
    end
end

% -------------------------------------------------------------------------
function units = local_units(raw, o)
%LOCAL_UNITS Resolve 'auto', or validate an explicit choice.
    want = lower(strtrim(char(o.Units)));
    if ~strcmp(want, 'auto')
        assert(ismember(want, {'dn', 'radians', 'waves'}), ...
            'load_slm_correction:badUnits', ...
            'Units must be auto / dn / radians / waves, got ''%s''.', want);
        units = want;
        return
    end

    v = double(raw(:));
    is_int  = all(v == round(v));
    in_dn   = min(v) >= 0 && max(v) <= o.PixelMax;

    if is_int && in_dn
        units = 'dn';
        return
    end
    if max(abs(v)) > 3.2
        units = 'radians';
        return
    end

    error('load_slm_correction:ambiguousUnits', ...
        ['Cannot tell radians from waves for this data (range %.4g to %.4g).\n' ...
         'A small correction looks the same in both, and guessing wrong scales it ' ...
         'by 2*pi,\nso say which: ''Units'', ''radians'' or ''Units'', ''waves''.'], ...
        min(v), max(v));
end

% -------------------------------------------------------------------------
function s = local_dims(a)
    s = strjoin(arrayfun(@(n) sprintf('%d', n), size(a), 'UniformOutput', false), 'x');
end
