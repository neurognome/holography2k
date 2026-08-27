function [phase, info] = zernike_phase(nx, ny, coeffs, varargin)
%ZERNIKE_PHASE Render Zernike coefficients into a phase map on the SLM grid.
%   phase = ZERNIKE_PHASE(nx, ny, coeffs) returns an nx-by-ny map of RADIANS,
%   laid out the way the CGH indexes a hologram: dimension 1 is SLM x,
%   dimension 2 is SLM y.
%
%   Phase 3 item 17 of the AO bead protocol asks for the correction to come back as
%   "Zernike coefficients with the convention stated", and warns that an unlabeled
%   phase map costs a day of reconciliation. This function is the other half of that
%   bargain: the convention is named in the call, in the rig file, and in the info
%   struct, so it stays auditable instead of being frozen into a rendered image.
%
%   COEFFS is always a VECTOR OF VALUES. What those values mean depends on whether
%   you pass 'Modes':
%     * without 'Modes', they are consecutive single-index modes in the chosen
%       convention (OSA/ANSI starts at j = 0 with piston, Noll starts at j = 1);
%     * with 'Modes', they pair up with an m-by-2 list of [n m] orders, which needs
%       no ordering convention at all and is the form to prefer when you can get it.
%
%   There is deliberately no "m-by-3 [n m value] matrix" form. It reads well but
%   collides with the vector form at exactly one input -- a single [n m value] row is
%   also a 3-element vector -- and the two readings differ silently: [2 0 1] is
%   either "1 radian of defocus" or "2 radians of piston plus 1 of x-tilt". Naming
%   the modes separately makes that impossible to get wrong.
%
%   Name/value options:
%     'Modes'       m-by-2 [n m] rows, same length as COEFFS. Omit to use the
%                   single-index convention below.
%     'Convention'  'osa' (default) | 'noll'   -- only used when 'Modes' is omitted.
%     'Normalized'  false (default) | true     -- true applies the orthonormal
%                   sqrt(2(n+1)/(1+(m==0))) factor. Getting this wrong rescales
%                   every term by a mode-dependent factor of order 2-3, so it is
%                   worth confirming rather than assuming.
%     'Pupil'       [cx cy r] in PIXELS, defaulting to the circle inscribed in the
%                   grid. Pass the values measured by ao_bead_calibration's
%                   slm_pupil_mapping (Phase 1) -- the inscribed default is a
%                   placeholder and every coefficient's meaning depends on it.
%     'Mask'        true (default) zeroes the phase outside the pupil.
%                   Zernikes are polynomials: they are defined outside r = 1 but
%                   not orthogonal there, and the r^4 terms reach 4x their
%                   pupil-edge value in the corners of a square grid (r = sqrt 2).
%                   Extrapolating would add waves of meaningless phase where the
%                   beam mostly is not, so the default zeroes it.
%     'Piston'      false (default) drops the n = 0 term. A uniform phase offset
%                   cannot change a far-field intensity, so it is dropped to keep
%                   the reported RMS honest.
%
%   INFO reports the modes used, the pupil, and the RMS / peak-to-valley in waves
%   over the pupil.
%
%   See also: load_slm_correction, phase_to_frame

    p = inputParser;
    p.FunctionName = 'zernike_phase';
    p.addParameter('Modes', []);
    p.addParameter('Convention', 'osa');
    p.addParameter('Normalized', false);
    p.addParameter('Pupil', []);
    p.addParameter('Mask', true);
    p.addParameter('Piston', false);
    p.parse(varargin{:});
    o = p.Results;

    conv = lower(strtrim(char(o.Convention)));
    assert(ismember(conv, {'osa', 'ansi', 'noll'}), 'zernike_phase:badConvention', ...
        'Convention must be ''osa'' (= ''ansi'') or ''noll'', got ''%s''.', conv);
    if strcmp(conv, 'ansi'), conv = 'osa'; end

    modes = local_modes(coeffs, o.Modes, conv);

    if ~o.Piston
        modes = modes(modes(:,1) > 0, :);
    end
    assert(~isempty(modes), 'zernike_phase:noModes', ...
        ['No modes left to render. Either COEFFS was empty, or it held only a ' ...
         'piston term\nand ''Piston'' is false (a uniform offset cannot change ' ...
         'the far field).']);

    pupil = o.Pupil;
    if isempty(pupil)
        pupil = [(nx-1)/2, (ny-1)/2, (min(nx, ny) - 1)/2];
    end
    assert(numel(pupil) == 3 && all(isfinite(pupil)) && pupil(3) > 0, ...
        'zernike_phase:badPupil', ...
        'Pupil must be [cx cy r] in pixels with r > 0, got %s.', mat2str(pupil));

    % ndgrid, not meshgrid: dimension 1 must be x, to match the Masks(x,y) layout
    % that function_Make_3D_SHOT_Holos builds and phase_to_frame consumes.
    [XX, YY] = ndgrid(0:(nx-1), 0:(ny-1));
    X = (XX - pupil(1)) / pupil(3);
    Y = (YY - pupil(2)) / pupil(3);
    R = hypot(X, Y);
    TH = atan2(Y, X);

    phase = zeros(nx, ny);
    for i = 1:size(modes, 1)
        n = modes(i,1); m = modes(i,2); c = modes(i,3);
        if c == 0, continue; end
        Z = local_radial(n, abs(m), R);
        if m >= 0
            Z = Z .* cos(abs(m) * TH);
        else
            Z = Z .* sin(abs(m) * TH);
        end
        if o.Normalized
            Z = Z * sqrt(2*(n+1) / (1 + double(m == 0)));
        end
        phase = phase + c * Z;
    end

    inside = R <= 1;
    if o.Mask
        phase(~inside) = 0;
    end

    info = struct();
    info.modes       = modes;
    info.convention  = conv;
    info.normalized  = logical(o.Normalized);
    info.pupil       = pupil(:)';
    info.masked      = logical(o.Mask);
    info.rms_waves   = std(phase(inside)) / (2*pi);
    info.pv_waves    = (max(phase(inside)) - min(phase(inside))) / (2*pi);
    info.summary     = sprintf(['%d mode(s), %s%s, pupil [%.1f %.1f r=%.1f]px, ' ...
                                'RMS %.3f waves, PV %.3f waves'], ...
        size(modes, 1), upper(conv), local_iff(o.Normalized, ' normalized', ' unnormalized'), ...
        pupil(1), pupil(2), pupil(3), info.rms_waves, info.pv_waves);
end

% -------------------------------------------------------------------------
function modes = local_modes(coeffs, given, conv)
%LOCAL_MODES Pair COEFFS with mode orders, giving [n m value] rows.
    assert(isnumeric(coeffs) && ~isempty(coeffs) && isvector(coeffs), ...
        'zernike_phase:badCoeffs', ...
        ['COEFFS must be a non-empty numeric VECTOR of values. To name the modes ' ...
         'explicitly,\npass ''Modes'', [n m; ...] alongside it -- there is no ' ...
         '[n m value] matrix form (see help).']);
    c = double(coeffs(:));

    if ~isempty(given)
        assert(size(given, 2) == 2, 'zernike_phase:badModes', ...
            'Modes must be an m-by-2 list of [n m] rows, got %d column(s).', ...
            size(given, 2));
        assert(size(given, 1) == numel(c), 'zernike_phase:modeCountMismatch', ...
            ['Got %d mode(s) but %d coefficient(s). They pair up one-to-one, so a ' ...
             'mismatch\nwould silently drop or invent terms.'], ...
            size(given, 1), numel(c));
        n = double(given(:,1)); m = double(given(:,2));
        assert(all(n >= 0 & mod(n,1) == 0), 'zernike_phase:badN', ...
            'Zernike radial order n must be a non-negative integer.');
        bad = ~(abs(m) <= n & mod(n - abs(m), 2) == 0);
        assert(~any(bad), 'zernike_phase:badM', ...
            ['Zernike azimuthal order m must satisfy |m| <= n and n-|m| even. ' ...
             'Offending row(s): %s.'], mat2str(find(bad)'));
        modes = [n, m, c];
        return
    end

    modes = zeros(numel(c), 3);
    for j = 1:numel(c)
        if strcmp(conv, 'osa')
            [n, m] = local_osa(j - 1);   % OSA/ANSI single index starts at 0
        else
            [n, m] = local_noll(j);      % Noll single index starts at 1
        end
        modes(j,:) = [n, m, c(j)];
    end
end

function [n, m] = local_osa(j)
%LOCAL_OSA OSA/ANSI single index -> (n, m). j = 0 is piston.
    n = floor((-1 + sqrt(1 + 8*j)) / 2);
    m = 2*(j - n*(n+1)/2) - n;
end

function [n, m] = local_noll(j)
%LOCAL_NOLL Noll single index -> (n, m). j = 1 is piston.
%   Noll's rules, and the two that are easy to get backwards:
%     * within one radial order n, |m| ASCENDS;
%     * m > 0 (cosine) for even j, m < 0 (sine) for odd j.
%   The slot arithmetic is the subtle part: |m| = 0 exists only for even n and
%   occupies ONE slot, while every other |m| occupies TWO (the sine/cosine pair).
%   Treating them all as pairs shifts every term of even n by one index, which
%   silently relabels defocus as astigmatism.
%
%   Reference table this reproduces:
%     j : 1     2     3      4     5      6     7      8     9      10
%         (0,0) (1,1) (1,-1) (2,0) (2,-2) (2,2) (3,-1) (3,1) (3,-3) (3,3)
%     j : 11    12    13     14    15
%         (4,0) (4,2) (4,-2) (4,4) (4,-4)
    n = 0;
    while (n+1)*(n+2)/2 < j
        n = n + 1;
    end
    idx = j - n*(n+1)/2;          % 1-based position within this radial order
    if mod(n, 2) == 0
        % even n: slot 1 is |m| = 0, then pairs 2,2, 4,4, ...
        if idx == 1
            mm = 0;
        else
            mm = 2 * ceil((idx - 1) / 2);
        end
    else
        % odd n: no |m| = 0; pairs 1,1, 3,3, ...
        mm = 2 * ceil(idx / 2) - 1;
    end
    if mm == 0
        m = 0;
    elseif mod(j, 2) == 0
        m = mm;
    else
        m = -mm;
    end
end

function R = local_radial(n, m, rho)
%LOCAL_RADIAL Zernike radial polynomial R_n^m(rho).
    R = zeros(size(rho));
    for k = 0:((n - m)/2)
        c = (-1)^k * factorial(n - k) / ...
            (factorial(k) * factorial((n + m)/2 - k) * factorial((n - m)/2 - k));
        R = R + c * rho.^(n - 2*k);
    end
end

function s = local_iff(tf, a, b)
    if tf, s = a; else, s = b; end
end
