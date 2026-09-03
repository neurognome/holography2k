function [measured, found, det] = detect_holeburns(base, post, commanded, varargin)
%DETECT_HOLEBURNS Find each commanded hole in a before/after image pair.
%   [measured, found, det] = DETECT_HOLEBURNS(base, post, commanded)
%
%   base, post  the SAME field, imaged before and after the burn (2D doubles)
%   commanded   Nx2 target positions, [row col] -- see the note below
%   measured    Nx2 sub-pixel positions, NaN where nothing was found
%   found       Nx1 logical
%   det         per-target diagnostics: z-score, peak value, search radius used
%
%   COORDINATES. commanded is [row col], i.e. [ScanImage Y, ScanImage X], which
%   is the same order a holoRequest's targets are in. That is not a coincidence
%   and not a choice made here: saveAllToHoloRequest builds targets by
%   fliplr'ing scanfields.centerXY, so targets(:,1) is Y and targets(:,2) is X.
%   MATLAB indexes images (row, col) = (Y, X), so the target columns line up
%   with image indices with no swap anywhere. Keeping that alignment is why this
%   function takes [row col] rather than the more natural-looking [x y]: the one
%   place a swap could hide is the place it would be hardest to find.
%
%   METHOD, lifted from process_holeburns_claude and simplified for a known grid:
%     1. difference base - post, so a burned hole is POSITIVE (fluorescence lost)
%     2. difference of Gaussians, to suppress illumination drift and enhance
%        spot-scale features
%     3. per target, search a window around where the hole was ASKED for
%     4. accept the peak only if it stands out from that window robustly
%     5. sub-pixel centroid on the accepted peak
%
%   Searching per target, rather than finding all the blobs and matching them
%   afterwards, is what a known grid buys: two holes can never be confused for
%   each other as long as the search radius is under half the grid pitch, which
%   is the default.
%
%   Name-value options:
%     'Radius'      search radius in px (default: 0.4 * the grid pitch)
%     'FineSigma'   pre-smoothing sigma            (default 1)
%     'DoGSigma'    [small large] DoG sigmas       (default [1.5 5])
%     'MinZ'        minimum robust z-score to accept (default 2.5)
%     'CentroidR'   radius of the centroid patch   (default 4)
%
%   See also: read_holeburn_tiff, fit_holeburn_offsets

    p = inputParser;
    p.addParameter('Radius',    []);
    p.addParameter('FineSigma', 1);
    p.addParameter('DoGSigma',  [1.5 5]);
    p.addParameter('MinZ',      2.5);
    p.addParameter('CentroidR', 4);
    p.parse(varargin{:});
    o = p.Results;

    base = double(base);
    post = double(post);
    assert(isequal(size(base), size(post)), 'detect_holeburns:size', ...
        'base is %s but post is %s -- they must be the same field.', ...
        mat2str(size(base)), mat2str(size(post)));
    assert(size(commanded, 2) == 2, 'detect_holeburns:commanded', ...
        'commanded must be Nx2 ([row col]).');

    [H, W] = size(base);
    n = size(commanded, 1);

    if isempty(o.Radius)
        o.Radius = 0.4 * local_pitch(commanded);
    end
    o.Radius = max(3, o.Radius);

    % 1-2. bright where a hole appeared, then DoG.
    d   = local_gauss(base, o.FineSigma) - local_gauss(post, o.FineSigma);
    dog = local_gauss(d, o.DoGSigma(1)) - local_gauss(d, o.DoGSigma(2));

    measured = nan(n, 2);
    found    = false(n, 1);
    det = struct('z', num2cell(nan(n,1)), 'peak', num2cell(nan(n,1)), ...
                 'radius', num2cell(repmat(o.Radius, n, 1)));

    r = round(o.Radius);
    for i = 1:n
        r0 = round(commanded(i, 1));
        c0 = round(commanded(i, 2));

        rows = max(1, r0 - r) : min(H, r0 + r);
        cols = max(1, c0 - r) : min(W, c0 + r);
        if isempty(rows) || isempty(cols)
            continue    % commanded outside the frame; not a detection failure
        end

        win = dog(rows, cols);
        [pk, k] = max(win(:));
        [pr, pc] = ind2sub(size(win), k);

        % Robust z against the window itself. Median/MAD rather than mean/std
        % because the window CONTAINS the thing we are measuring -- a real hole
        % would inflate a standard deviation and hide itself.
        med = median(win(:));
        sd  = 1.4826 * median(abs(win(:) - med));
        if sd <= 0
            z = 0;
        else
            z = (pk - med) / sd;
        end

        det(i).z    = z;
        det(i).peak = pk;
        if z < o.MinZ
            continue
        end

        % 5. sub-pixel centroid on a small patch around the peak, over the part
        % of it above half the peak height. Half-max keeps the centroid on the
        % hole rather than dragging it toward whatever the local background does.
        gr = rows(pr);
        gc = cols(pc);
        pr_rows = max(1, gr - o.CentroidR) : min(H, gr + o.CentroidR);
        pr_cols = max(1, gc - o.CentroidR) : min(W, gc + o.CentroidR);
        patch = dog(pr_rows, pr_cols);
        wgt = patch - max(0, med);
        wgt(wgt < 0.5 * (pk - med)) = 0;
        if ~any(wgt(:) > 0)
            continue
        end
        [cc, rr] = meshgrid(pr_cols, pr_rows);
        tot = sum(wgt(:));
        measured(i, 1) = sum(rr(:) .* wgt(:)) / tot;
        measured(i, 2) = sum(cc(:) .* wgt(:)) / tot;
        found(i) = true;
    end
end

% -------------------------------------------------------------------------
function pitch = local_pitch(commanded)
%LOCAL_PITCH Nearest-neighbour spacing of the commanded grid.
%   Used only to size the default search radius, so the median is the right
%   summary: a couple of dropped grid points must not widen every window.
    n = size(commanded, 1);
    if n < 2
        pitch = 20;
        return
    end
    d = inf(n, 1);
    for i = 1:n
        v = commanded - commanded(i, :);
        v(i, :) = [];
        d(i) = min(sqrt(sum(v.^2, 2)));
    end
    pitch = median(d);
    if ~isfinite(pitch) || pitch <= 0
        pitch = 20;
    end
end

function out = local_gauss(img, sigma)
%LOCAL_GAUSS imgaussfilt when the Image Processing Toolbox is there, else conv2.
%   The fallback is not a nicety: it means the hole detector, and therefore the
%   tests for it, run on a machine without IPT. A separable Gaussian with
%   symmetric edge padding is what imgaussfilt does by default anyway.
    persistent has_ipt
    if isempty(has_ipt)
        has_ipt = exist('imgaussfilt', 'file') == 2 || exist('imgaussfilt', 'builtin') == 5;
    end
    if sigma <= 0
        out = img;
        return
    end
    if has_ipt
        out = imgaussfilt(img, sigma);
        return
    end
    r = max(1, ceil(3 * sigma));
    x = -r:r;
    k = exp(-(x.^2) / (2 * sigma^2));
    k = k / sum(k);
    pad = padarray_sym(img, r);
    out = conv2(conv2(pad, k, 'same'), k', 'same');
    out = out(r+1:end-r, r+1:end-r);
end

function out = padarray_sym(img, r)
%PADARRAY_SYM Symmetric padding, without needing the Image Processing Toolbox.
    ri = [r+1:-1:2, 1:size(img,1), size(img,1)-1:-1:size(img,1)-r];
    ci = [r+1:-1:2, 1:size(img,2), size(img,2)-1:-1:size(img,2)-r];
    ri = min(max(ri, 1), size(img,1));
    ci = min(max(ci, 1), size(img,2));
    out = img(ri, ci);
end
