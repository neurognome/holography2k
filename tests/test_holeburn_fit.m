function test_holeburn_fit()
%TEST_HOLEBURN_FIT Detect a synthetic burn, fit it, and prove the fix works.
%   Run from the repo root (no hardware, no images needed):
%       matlab -batch "addpath(genpath(pwd)); test_holeburn_fit"
%   Errors on the first failure; prints PASS on success.
%
%   THE CLAIM THIS TEST EXISTS TO CHECK. The adjustment to write is not "minus
%   the residual". Because the delivery mapping has a scale in it, the correct
%   inversion is scale = 1/s and offset = -t/s. On a rig with s = 1 the two are
%   indistinguishable, so a sign or inversion error here would pass every casual
%   check and then drift as soon as the scale moved off unity. So the test does
%   what the operator's verify pass does, in silico: fit the offsets, apply them
%   through the SAME arithmetic generate_holograms_new uses, and require the
%   residual to collapse.
%
%   It also pins the axis convention. commanded is [row col] = [ScanImage Y,
%   ScanImage X], xoffset goes to column 1 and yoffset to column 2, matching a
%   holoRequest's [y x z] targets. The simulated truth here is deliberately
%   ASYMMETRIC (different offset per axis) so a transposition cannot pass.
%
%   See also: detect_holeburns, fit_holeburn_offsets, apply_holo_adjustments

    root = fileparts(fileparts(mfilename('fullpath')));
    p = strsplit(genpath(root), pathsep);
    p = p(~cellfun(@isempty, p));
    % Skip dot-directories BELOW root (.git, nested .claude worktrees, which
    % hold stale copies that would shadow the real files) -- but judge them
    % RELATIVE to root, so a root that itself sits under a dot-directory is not
    % refused outright. That is the case whenever this runs from a worktree.
    rel = cellfun(@(x) x(min(numel(root)+1, numel(x)+1):end), p, 'UniformOutput', false);
    p = p(~contains(rel, [filesep '.']));
    addpath(strjoin(p, pathsep));

    rng(7);

    SZ = 512;
    C  = [256 256];          % the pivot everything is defined about
    S_TRUE = 1.018;          % delivery is 1.8% too big
    T_TRUE = [-13.2, 11.4];  % [row col] -- deliberately different per axis

    % ---- a 6x6 commanded grid, inset from the edge ------------------------
    lo = 70; hi = 442; n = 6;
    v = linspace(lo, hi, n);
    [gc, gr] = meshgrid(v, v);
    commanded = [gr(:), gc(:)];          % [row col]

    % ---- where the light actually lands -----------------------------------
    land = S_TRUE .* (commanded - C) + C + T_TRUE;

    [base, post] = local_scene(SZ, land);

    % ---- detect ------------------------------------------------------------
    [measured, found] = detect_holeburns(base, post, commanded);
    assert(sum(found) >= 34, 'only detected %d of 36 holes', sum(found));

    err = sqrt(sum((measured(found,:) - land(found,:)).^2, 2));
    assert(max(err) < 1.5, ...
        'detection is off by up to %.2f px from the simulated truth', max(err));

    % ---- fit ---------------------------------------------------------------
    fit = fit_holeburn_offsets(commanded, measured, found, 'ScaleCenter', C);
    assert(fit.ok, 'fit refused: %s', fit.why);

    assert(abs(fit.s - S_TRUE) < 5e-4, ...
        'recovered s %.5f, expected %.5f', fit.s, S_TRUE);
    assert(norm(fit.t - T_TRUE) < 0.5, ...
        'recovered t %s, expected %s', mat2str(fit.t, 4), mat2str(T_TRUE));

    % The inversion: scale = 1/s, offset = -t/s. Checked as VALUES, so a lazy
    % "offset = -t" cannot pass -- at s = 1.018 the two differ by ~0.2 px, which
    % is small but systematically wrong and grows with the scale error.
    assert(abs(fit.scale - 1/S_TRUE) < 1e-4, ...
        'scale should be 1/s = %.5f, got %.5f', 1/S_TRUE, fit.scale);
    expect_off = -T_TRUE / S_TRUE;
    assert(abs(fit.xoffset - expect_off(1)) < 0.5, ...
        'xoffset should be %.3f (row axis), got %.3f', expect_off(1), fit.xoffset);
    assert(abs(fit.yoffset - expect_off(2)) < 0.5, ...
        'yoffset should be %.3f (col axis), got %.3f', expect_off(2), fit.yoffset);

    % A transposed answer must NOT pass: the two offsets differ enough that
    % swapping them is a >20 px error.
    assert(abs(fit.xoffset - expect_off(2)) > 5, ...
        'the test grid is too symmetric to catch a row/col transposition');

    assert(fit.rms_before_px > 5, ...
        'the simulated misalignment should be obvious, rms_before was %.2f', ...
        fit.rms_before_px);
    assert(fit.rms_after_px < 1.0, ...
        'a translation+scale model should describe this exactly, rms_after %.2f', ...
        fit.rms_after_px);
    assert(abs(fit.rotation_deg) < 0.2, ...
        'no rotation was simulated but %.3f deg was reported', fit.rotation_deg);

    % ---- the verify pass, in silico ----------------------------------------
    % Rewrite the targets exactly as generate_holograms_new does, then deliver
    % them through the same physical mapping. This is the acceptance test the
    % operator runs on the rig, and it is the only thing that actually proves
    % the sign and the inversion.
    q2 = local_apply(commanded, fit.xoffset, fit.yoffset, fit.scale, C);
    land2 = S_TRUE .* (q2 - C) + C + T_TRUE;
    verify_rms = sqrt(mean(sum((land2 - commanded).^2, 2)));
    % 0.1 px, not 0: the adjustment is fitted from sub-pixel centroids of noisy
    % spots, so it inherits their scatter. Measured at ~0.02 px here, which is
    % two orders below the misalignment being corrected and far below the ~5 px
    % the naive inversion leaves behind (checked next) -- so this bound
    % discriminates between a right and a wrong answer without pretending the
    % fit is exact.
    assert(verify_rms < 0.1, ...
        ['applying the fitted adjustment should put the light where it was ' ...
         'asked for, but the residual is %.3f px'], verify_rms);

    % ---- and a wrong answer must fail that same check -----------------------
    % "offset = -t, scale = 1" is the plausible mistake. It must not pass.
    q3 = local_apply(commanded, -T_TRUE(1), -T_TRUE(2), 1, C);
    land3 = S_TRUE .* (q3 - C) + C + T_TRUE;
    naive_rms = sqrt(mean(sum((land3 - commanded).^2, 2)));
    assert(naive_rms > 10 * verify_rms && naive_rms > 1, ...
        ['the naive inversion happens to work on this test case, so the test ' ...
         'does not discriminate (rms %.3f)'], naive_rms);

    % ---- refusal when too little was detected -------------------------------
    few = found; few(4:end) = false;
    bad = fit_holeburn_offsets(commanded, measured, few, 'ScaleCenter', C);
    assert(~bad.ok, 'a 3-of-36 detection must be refused');
    assert(contains(bad.why, 'detected'), 'refusal should explain: %s', bad.why);

    fprintf('PASS test_holeburn_fit  (s %.4f, t %s, verify rms %.4f px, naive %.2f px)\n', ...
        fit.s, mat2str(fit.t, 4), verify_rms, naive_rms);
end

% -------------------------------------------------------------------------
function q = local_apply(q, xoffset, yoffset, scale, c)
%LOCAL_APPLY generate_holograms_new's target rewrite, verbatim in shape.
    q(:,1) = scale * (q(:,1) - c(1)) + c(1) + xoffset;
    q(:,2) = scale * (q(:,2) - c(2)) + c(2) + yoffset;
end

function [base, post] = local_scene(sz, land)
%LOCAL_SCENE A fluorescent slide before and after burning holes at LAND.
%   Uneven illumination and shot-ish noise, because a detector that only works
%   on a flat field is not one worth testing.
    [xx, yy] = meshgrid(1:sz, 1:sz);
    illum = 800 + 180 * exp(-((xx - 180).^2 + (yy - 300).^2) / (2 * 260^2));
    base = illum + 12 * randn(sz);

    post = illum;
    sigma = 2.2;
    depth = 260;            % how much fluorescence a hole removes
    for i = 1:size(land, 1)
        r0 = land(i,1); c0 = land(i,2);
        post = post - depth * exp(-((yy - r0).^2 + (xx - c0).^2) / (2 * sigma^2));
    end
    post = post + 12 * randn(sz);
end
