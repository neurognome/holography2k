function test_pattern_de()
%TEST_PATTERN_DE  Pattern.calculate_DE must honour powerbias, row or column. No rig.
%
%   >> test_pattern_de
%
%   Regression for the row/column broadcast bug: powerbias arrived as a ROW (which is
%   what calibrated_powerbias, channel_powerbias and Pattern's own default all produce)
%   while slm_coords(:,4) is a COLUMN, so `energy .* obj.powerbias` implicitly expanded
%   into an outer product and `energy/sum(energy)` became mrdivide. calculate_DE then
%   returned the FLAT-bias DE for every input -- powerbias was silently ignored.
%
%   Sequence.calculate_power divides the laser command by this DE, so the whole ensemble
%   was under-driven by a few percent on calibrated (non-flat) patterns. The invariant
%   below is what makes the delivered power come out equal to the requested power:
%
%       DE = sum(powerbias) / sum(powerbias ./ attenuation)
%
%   which for mean(powerbias) == 1 over M targets is M / sum(pb./att). Combined with
%   Sequence.calculate_power = ppc*M/DE and the hologram's own powerbias/att weighting,
%   the attenuation cancels and each target receives exactly ppc*powerbias(k).
%
%   Uses a local function_SItoSLM stub on the path (see stub_dir below) so no alignment
%   or SLM is needed. See also: test_powerbias_delivery.

    here = fileparts(mfilename('fullpath'));
    stub_dir = fullfile(here, 'stubs');
    addpath(fullfile(here, '..', 'roi_class'));
    addpath(stub_dir);
    assert(exist('function_SItoSLM', 'file') == 2, ...
        'need the function_SItoSLM stub in %s', stub_dir);

    att = [0.90 0.75 0.50 0.35 0.20]';
    closed_form = @(pb) min(sum(pb(:)) / sum(pb(:) ./ max(att, 0.05)), 1);

    % --- 1. a non-flat bias must give a DIFFERENT DE from a flat one ----------------
    pb_flat = ones(1, 5);
    pb_bias = [0.3991 1.1653 0.3991 2.4557 0.5807];   % calibrated_powerbias at X=2
    de_flat = de_of(pb_flat, att);
    de_bias = de_of(pb_bias, att);
    assert(abs(de_flat - de_bias) > 1e-6, ...
        ['DE is identical for a flat and a strongly biased pattern (%.10f) -- ', ...
         'powerbias is being ignored, which is exactly the old bug.'], de_flat);

    % --- 2. both must match the closed form ----------------------------------------
    for pb = {pb_flat, pb_bias}
        got = de_of(pb{1}, att); want = closed_form(pb{1});
        assert(abs(got - want) < 1e-12, 'DE %.12f != closed form %.12f', got, want);
    end

    % --- 3. ROW and COLUMN powerbias must agree ------------------------------------
    assert(abs(de_of(pb_bias, att) - de_of(pb_bias(:), att)) < 1e-12, ...
        'row and column powerbias give different DE -- the broadcast bug is back');

    % --- 4. scale-invariance: generate_holograms_new calls calculate_DE BEFORE it
    %        rescales powerbias by length/max_pattern_sz, so DE must not care ---------
    for c = [0.5 2 3.7]
        assert(abs(de_of(pb_bias * c, att) - de_of(pb_bias, att)) < 1e-12, ...
            'DE changed under a uniform powerbias scale of %g', c);
    end

    % --- 5. the delivered-power identity this all exists for ------------------------
    % laser = ppc*M/DE, hologram weights = pb./att normalised, delivered = laser*frac.*att
    %
    % The cancellation is exact only when mean(powerbias) == 1, i.e. sum(pb) == M --
    % Sequence.calculate_power multiplies by max(sz) = M, and DE's numerator is sum(pb),
    % so the two agree only under that normalisation. calibrated_powerbias guarantees it
    % by construction (powerbias = p/mean(p)); normalise here so the fixture does too,
    % rather than loosening the tolerance and hiding a real dependency.
    pb_norm = pb_bias / mean(pb_bias);
    assert(abs(mean(pb_norm) - 1) < eps(4), 'fixture must have mean(powerbias) == 1');
    ppc = 0.626356e-3;  M = numel(pb_norm);
    DE = de_of(pb_norm, att);
    laser = ppc * M / DE;
    w = pb_norm(:) ./ att;  frac = w / sum(w);
    delivered = laser * frac .* att;
    requested = ppc * pb_norm(:);
    err = max(abs(delivered - requested));
    assert(err < 1e-15, ...
        'delivered power differs from requested by %.3e W (should cancel exactly)', err);

    % --- 6. a mismatched powerbias length must ERROR, not outer-product silently ----
    p = Pattern([zeros(5,2) (1:5)'], ones(1,4), false);
    try
        p.calculate_DE(struct('att', att));
        error('test_pattern_de:noAssert', 'a 4-entry powerbias on 5 targets must error');
    catch e
        assert(strcmp(e.identifier, 'Pattern:powerbiasSize'), ...
            'wrong error for a size mismatch: %s', e.identifier);
    end

    fprintf(['PASS test_pattern_de (DE honours powerbias, row==col, scale-invariant, ', ...
             'delivery exact, size mismatch caught).\n']);
end

function de = de_of(pb, att)
    p = Pattern([zeros(numel(pb),2) (1:numel(pb))'], pb, false);
    p.calculate_DE(struct('att', att));
    de = p.diffraction_efficiency;
end
