function test_phase_to_frame()
%TEST_PHASE_TO_FRAME  The phase->frame wrap, with and without a static correction. No rig.
%
%   >> test_phase_to_frame
%
%   phase_to_frame replaced a line that was copy-pasted into four CGH functions, and
%   it is now the only place a static wavefront correction can enter a hologram. Two
%   properties have to hold or this change is not safe to put on the rig:
%
%   1. WITHOUT a correction it must be BYTE-IDENTICAL to the expression it replaced.
%      Every hologram this repo has ever compiled came out of that line, so any drift
%      would silently change alignment, PSF and diffraction-efficiency measurements
%      as well as stimulation -- and the DE numbers feed the power model, so the
%      change would show up as a power error nobody would connect to a refactor.
%
%   2. WITH a correction the addition must WRAP, not saturate. MATLAB uint8
%      arithmetic saturates (uint8(200)+uint8(100) == 255), so the obvious
%      client-side implementation -- add the correction to the compiled uint8 stack
%      -- silently clips every pixel whose sum passes 255 into a flat maximum-phase
%      region. That is a diffuse haze in the reconstruction, not an error. The test
%      below reproduces the saturating version to show the two genuinely differ.
%
%   See also: phase_to_frame, load_slm_correction, test_slm_correction

    Setup.SLM.pixelmax = 255;
    Setup.SLM.Nx = 16;
    Setup.SLM.Ny = 16;

    rng(7);   % deterministic; this test must not be flaky on the rig
    phase = (rand(16, 16) * 2 - 1) * pi;

    total = 0; failed = 0;

    % ---- 1. no correction == the legacy expression, exactly --------------------
    legacy = uint8(floor((Setup.SLM.pixelmax*(phase+pi)/(2*pi))));
    [total, failed] = check(total, failed, 'empty correction is byte-identical', ...
        isequal(phase_to_frame(Setup, phase), legacy));

    % Absent field, empty field, and empty double all mean "none".
    S2 = Setup; S2.SLM = rmfield(S2.SLM, 'Nx');       % unrelated field missing
    [total, failed] = check(total, failed, 'no correction field at all', ...
        isequal(phase_to_frame(S2, phase), legacy));
    S3 = Setup; S3.SLM.correction = [];
    [total, failed] = check(total, failed, 'empty correction field', ...
        isequal(phase_to_frame(S3, phase), legacy));

    % The endpoints are where a mod() rewrite would have differed, so pin them.
    edge = [-pi, 0, pi];
    [total, failed] = check(total, failed, 'endpoints match the legacy expression', ...
        isequal(phase_to_frame(Setup, edge), ...
                uint8(floor((Setup.SLM.pixelmax*(edge+pi)/(2*pi))))));

    % ---- 2. a full turn of correction is a no-op ------------------------------
    Sc = Setup; Sc.SLM.correction = 2*pi*ones(16, 16);
    [total, failed] = check(total, failed, 'a 2*pi correction changes nothing', ...
        isequal(phase_to_frame(Sc, phase), phase_to_frame(Setup, phase)));

    Sc.SLM.correction = -4*pi*ones(16, 16);
    [total, failed] = check(total, failed, 'a -4*pi correction changes nothing', ...
        isequal(phase_to_frame(Sc, phase), phase_to_frame(Setup, phase)));

    % ---- 3. correction is equivalent to pre-adding it to the phase -------------
    corr = (rand(16, 16) * 2 - 1) * 3*pi;      % deliberately beyond one turn
    Sc.SLM.correction = corr;
    wrapped_phase = mod(phase + corr + pi, 2*pi) - pi;
    [total, failed] = check(total, failed, 'correction == pre-added, pre-wrapped phase', ...
        isequal(phase_to_frame(Sc, phase), phase_to_frame(Setup, wrapped_phase)));

    % ---- 4. wrapping, not saturating ------------------------------------------
    % Push most pixels past the top of the range and confirm the result is spread
    % across the range rather than piled at 255.
    Sc.SLM.correction = 0.9 * 2*pi * ones(16, 16);
    got = phase_to_frame(Sc, phase);
    naive = uint8(double(phase_to_frame(Setup, phase)) + 0.9*255);   % saturating add
    [total, failed] = check(total, failed, 'wrapped result is not pinned at pixelmax', ...
        mean(got(:) == 255) < 0.02);
    [total, failed] = check(total, failed, 'the saturating add really does clip', ...
        mean(naive(:) == 255) > 0.5);
    [total, failed] = check(total, failed, 'wrap and saturating add differ', ...
        ~isequal(got, naive));

    % ---- 5. a linear tilt produces a linear ramp of the right period -----------
    % One full turn of tilt across the panel must give exactly one DN ramp, which is
    % what makes the on-rig tilt test (does the spot move by the predicted amount?)
    % interpretable.
    nx = Setup.SLM.Nx;
    tilt = 2*pi * ((0:(nx-1))' / nx) * ones(1, Setup.SLM.Ny);
    St = Setup; St.SLM.correction = tilt;
    ramp = phase_to_frame(St, zeros(nx, Setup.SLM.Ny));
    step = diff(double(ramp(:,1)));
    % The ramp wraps EXACTLY ONCE, and not at the end: the forward map is
    % floor(pixelmax*mod(phase+pi,2*pi)/(2*pi)), so the +pi offset puts the wrap in
    % the middle of a tilt that starts at zero. Every step is therefore the same
    % size modulo pixelmax, with one negative raw step where it turns over.
    [total, failed] = check(total, failed, 'one turn of tilt steps uniformly mod pixelmax', ...
        all(abs(mod(step, 255) - 255/nx) < 1.5));
    [total, failed] = check(total, failed, 'one turn of tilt wraps exactly once', ...
        sum(step < 0) == 1);
    [total, failed] = check(total, failed, 'tilt varies along x only', ...
        all(all(diff(double(ramp), 1, 2) == 0)));

    % ---- 6. the guards fire ---------------------------------------------------
    Sb = Setup; Sb.SLM.correction = zeros(8, 8);
    [total, failed] = check(total, failed, 'size mismatch errors', ...
        threw(@() phase_to_frame(Sb, phase), 'phase_to_frame:sizeMismatch'));

    Sn = Setup; Sn.SLM.correction = zeros(16, 16); Sn.SLM.correction(3,4) = NaN;
    [total, failed] = check(total, failed, 'NaN in the correction errors', ...
        threw(@() phase_to_frame(Sn, phase), 'phase_to_frame:nonFiniteCorrection'));

    fprintf('\ntest_phase_to_frame: %d/%d passed\n', total - failed, total);
    if failed > 0
        error('test_phase_to_frame:failed', '%d check(s) failed.', failed);
    end
end

% -------------------------------------------------------------------------
function [total, failed] = check(total, failed, name, tf)
    total = total + 1;
    if tf
        fprintf('  ok    %s\n', name);
    else
        failed = failed + 1;
        fprintf('  FAIL  %s\n', name);
    end
end

function tf = threw(fn, id)
    tf = false;
    try
        fn();
    catch err
        tf = strcmp(err.identifier, id);
        if ~tf
            fprintf('        (threw %s, wanted %s)\n', err.identifier, id);
        end
    end
end
