function fit = fit_holeburn_offsets(commanded, measured, found, varargin)
%FIT_HOLEBURN_OFFSETS Turn a burned grid into holoRequest xoffset/yoffset/scale.
%   fit = FIT_HOLEBURN_OFFSETS(commanded, measured, found)
%
%   commanded  Nx2 [row col] positions that were ASKED for
%   measured   Nx2 [row col] positions the holes actually landed at
%   found      Nx1 logical from detect_holeburns
%
%   Returns a struct carrying the three numbers a holoRequest takes -- xoffset,
%   yoffset, scale -- plus the pivot they are defined about and enough
%   provenance for the offsets store.
%
%   THE MODEL, and why the answer is not simply "minus the residual".
%   In measure mode the request carries 0, 0, 1, so what the burn measures is
%   the raw delivery mapping:
%
%       land(q) = s*(q - c) + c + t
%
%   Apply an adjustment and the request's target q is rewritten before delivery
%   (this is generate_holograms_new's arithmetic, exactly):
%
%       q' = s_a*(q - c) + c + off      then      land(q') = s*(q' - c) + c + t
%
%   Setting land(q') = q for all q gives
%
%       s_a = 1/s          off = -t/s
%
%   So the scale is INVERTED and the translation is divided by it. Writing -t
%   would be right only when s == 1, and it is the sort of wrong that looks
%   nearly right on a well-aligned rig and drifts as soon as the scale moves.
%   The pivot c cancels, so long as the fit and generate_holograms_new use the
%   same one -- which is why scale_center travels with the record instead of
%   being assumed at each end.
%
%   AXIS ORDER. Column 1 is the row (ScanImage Y) and becomes xoffset; column 2
%   is the column (ScanImage X) and becomes yoffset. That looks transposed and
%   is not: a holoRequest's targets are [y x z] because saveAllToHoloRequest
%   fliplr's scanfields.centerXY, and xoffset is added to targets(:,1). The
%   comments in the historical get_holo_adjustments -- "xoffset ... moves holes
%   down", "yoffset ... more left" -- describe exactly this and are correct.
%
%   TWO RESIDUALS, and they answer different questions:
%     rms_before_px  rms|measured - commanded|. In a VERIFY pass this is the
%                    number that matters: it says whether the offsets already in
%                    use put the light where it was asked for.
%     rms_after_px   rms of what is left once s and t are fitted out. This says
%                    whether a translation-plus-uniform-scale model can describe
%                    the error at all. If it stays high while rms_before is
%                    high, the problem is not an offset.
%
%   rotation_deg is a DIAGNOSTIC only, from a full affine fit that is otherwise
%   discarded: the holoRequest has no field for rotation, so this exists so that
%   leftover rotation shows up as a number instead of hiding inside a slightly
%   wrong scale.
%
%   Name-value options:
%     'ScaleCenter'  pivot, [row col]        (default [256 256])
%     'MinFraction'  detection floor         (default 0.7)
%     'Sigma'        outlier rejection       (default 2.5)
%     'Iterations'   robust refits           (default 2)
%
%   See also: detect_holeburns, holo_adjustments_write, apply_holo_adjustments

    p = inputParser;
    p.addParameter('ScaleCenter', [256 256]);
    p.addParameter('MinFraction', 0.7);
    p.addParameter('Sigma',       2.5);
    p.addParameter('Iterations',  2);
    p.parse(varargin{:});
    o = p.Results;

    assert(size(commanded, 2) == 2 && size(measured, 2) == 2, ...
        'fit_holeburn_offsets:shape', 'commanded and measured must be Nx2.');
    assert(size(commanded, 1) == size(measured, 1), ...
        'fit_holeburn_offsets:count', ...
        'commanded has %d rows but measured has %d.', ...
        size(commanded, 1), size(measured, 1));
    if nargin < 3 || isempty(found)
        found = all(isfinite(measured), 2);
    end
    found = logical(found(:)) & all(isfinite(measured), 2);

    c = reshape(double(o.ScaleCenter), 1, 2);
    n = size(commanded, 1);

    fit = struct();
    fit.ok            = false;
    fit.why           = '';
    fit.scale_center  = c;
    fit.n_commanded   = n;
    fit.n_detected    = sum(found);
    fit.n_used        = 0;
    fit.used          = false(n, 1);
    fit.residual      = nan(n, 2);
    fit.rms_before_px = NaN;
    fit.rms_after_px  = NaN;
    fit.rotation_deg  = NaN;
    fit.xoffset = NaN; fit.yoffset = NaN; fit.scale = NaN;
    fit.s = NaN; fit.t = [NaN NaN];

    % Refuse early. A fit to a handful of holes converges happily and produces a
    % confident, wrong offset -- which then steers the beam. Better to say the
    % burn did not work.
    if fit.n_detected < 3
        fit.why = sprintf('only %d hole(s) detected; need at least 3.', fit.n_detected);
        return
    end
    frac = fit.n_detected / n;
    if frac < o.MinFraction
        fit.why = sprintf( ...
            ['only %d of %d holes detected (%.0f%%, floor %.0f%%). Too few to ' ...
             'trust a fit -- raise the burn power or dwell, check focus, or ' ...
             'confirm the PMT was on for the post image.'], ...
            fit.n_detected, n, 100*frac, 100*o.MinFraction);
        return
    end

    fit.rms_before_px = local_rms(measured(found, :) - commanded(found, :));

    % --- robust fit of measured = s*(commanded - c) + c + t ------------------
    use = found;
    s = NaN; t = [NaN NaN];
    for it = 1:max(1, o.Iterations)
        [s, t] = local_solve(commanded(use, :), measured(use, :), c);
        pred = s .* (commanded - c) + c + t;
        res  = measured - pred;
        rn   = sqrt(sum(res.^2, 2));
        rn(~found) = NaN;

        med = median(rn(found));
        mad = 1.4826 * median(abs(rn(found) - med));
        if ~isfinite(mad) || mad <= 0
            break   % a perfect fit; nothing to reject
        end
        keep = found & (rn <= med + o.Sigma * mad);
        if sum(keep) < 3 || isequal(keep, use)
            break
        end
        use = keep;
    end

    pred = s .* (commanded - c) + c + t;
    fit.residual = measured - pred;
    fit.residual(~found, :) = NaN;
    fit.used   = use;
    fit.n_used = sum(use);
    fit.rms_after_px = local_rms(fit.residual(use, :));

    fit.s = s;
    fit.t = t;

    % --- invert into what a holoRequest takes --------------------------------
    if ~isfinite(s) || abs(s) < 1e-6
        fit.why = sprintf('degenerate scale (%g); the grid may have collapsed.', s);
        return
    end
    adj = -t / s;
    fit.scale   = 1 / s;
    fit.xoffset = adj(1);    % targets(:,1) -- ScanImage Y
    fit.yoffset = adj(2);    % targets(:,2) -- ScanImage X

    % --- diagnostic: is there rotation the 3-dof model cannot express? -------
    fit.rotation_deg = local_rotation(commanded(use, :), measured(use, :), c);

    fit.ok  = true;
    % n_detected and n_used are different numbers and the difference matters:
    % detected is how many holes the burn actually made and the detector found,
    % used is how many survived outlier rejection. A big gap between them means
    % some holes landed nowhere near a translation-plus-scale prediction, which
    % is a different problem from not burning at all.
    fit.why = sprintf('%d found, %d used of %d, rms %.2f -> %.2f px', ...
        fit.n_detected, fit.n_used, n, fit.rms_before_px, fit.rms_after_px);
end

% -------------------------------------------------------------------------
function [s, t] = local_solve(q, m, c)
%LOCAL_SOLVE Least squares for one shared scale and a per-axis translation.
%   Both axes share s (that is what "uniform" means) but get their own offset,
%   so the design matrix stacks the two axes into one system rather than
%   fitting them separately -- fitting separately would give two scales.
    u = q - c;
    v = m - c;
    N = size(q, 1);
    A = [u(:,1), ones(N,1), zeros(N,1)
         u(:,2), zeros(N,1), ones(N,1)];
    b = [v(:,1); v(:,2)];
    x = A \ b;
    s = x(1);
    t = [x(2), x(3)];
end

function ang = local_rotation(q, m, c)
%LOCAL_ROTATION Rotation implied by a full affine fit, in degrees.
%   Discarded except as a printed diagnostic: a holoRequest cannot carry a
%   rotation, so the only useful thing to do with it is show it.
    u = q - c;
    v = m - c;
    N = size(q, 1);
    if N < 3
        ang = NaN;
        return
    end
    X = [u, ones(N,1)];
    A1 = X \ v(:,1);      % [a11 a12 b1]
    A2 = X \ v(:,2);      % [a21 a22 b2]
    a11 = A1(1); a12 = A1(2);
    a21 = A2(1); a22 = A2(2);
    ang = atan2d(a21 - a12, a11 + a22);
end

function r = local_rms(v)
    v = v(all(isfinite(v), 2), :);
    if isempty(v)
        r = NaN;
    else
        r = sqrt(mean(sum(v.^2, 2)));
    end
end
