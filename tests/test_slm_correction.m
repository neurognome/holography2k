function test_slm_correction()
%TEST_SLM_CORRECTION  load_slm_correction and zernike_phase, offline. No rig, no SLM.
%
%   >> test_slm_correction
%
%   Covers the three things about a correction file that are wrong-by-default rather
%   than wrong-by-accident:
%
%   ORIENTATION. imread returns (row, column) = (y, x); a hologram in this repo is
%   indexed (x, y), because function_Make_3D_SHOT_Holos builds Masks(sxx, syy) and
%   hands that array straight to Write_image. An image source therefore has to be
%   transposed. On the scope2k panel (1024x1024) no size check can catch it, so the
%   test uses a NON-SQUARE map where a missed transpose is detectable.
%
%   UNITS. A Meadowlark-style file stores DN with 0..pixelmax spanning 2*pi. Reading
%   it as radians under-rotates by 2*pi/255 per count; reading radians as waves
%   over-rotates by 2*pi. The loader auto-detects the unambiguous cases and REFUSES
%   the ambiguous one, and that refusal is the behaviour under test -- a guess here
%   is a whole-correction scale error.
%
%   ZERNIKE INDEXING. Noll's single index puts |m| = 0 in ONE slot and every other
%   |m| in TWO. Treating them uniformly shifts every even-order term by one index,
%   which relabels defocus as astigmatism -- the exact class of mistake Phase 3 item
%   17 of the AO protocol warns costs a day. The table from that item is asserted
%   mode by mode.
%
%   See also: load_slm_correction, zernike_phase, phase_to_frame

    rng(11);   % deterministic; this test must not be flaky on the rig

    total = 0; failed = 0;
    tmp = tempname;
    mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's'));

    % ---- 1. DN image round trip, non-square so orientation is testable ---------
    nx = 8; ny = 5;                       % panel is (x, y) = 8 x 5
    want = zeros(nx, ny);
    for ix = 1:nx
        for iy = 1:ny
            want(ix, iy) = ix * 10 + iy;  % every pixel distinct
        end
    end
    dn = uint8(round(want));              % small integers -> unambiguous DN
    f = fullfile(tmp, 'corr.bmp');
    imwrite(dn.', f);                     % write as an IMAGE: (y, x) = 5 x 8

    [base, info] = load_slm_correction(f, 'Size', [nx ny], 'Verbose', false);
    [total, failed] = check(total, failed, 'image source comes back as (x, y)', ...
        isequal(size(base), [nx ny]));
    [total, failed] = check(total, failed, 'image source transposes by default', ...
        info.transpose);
    [total, failed] = check(total, failed, 'DN detected automatically', ...
        strcmp(info.units, 'dn'));

    % Piston is removed, so compare shapes rather than absolute values.
    expect = double(dn) / 255 * 2*pi;
    expect = expect - mean(expect(:));
    [total, failed] = check(total, failed, 'DN converts to radians correctly', ...
        max(abs(base(:) - expect(:))) < 1e-12);

    % The whole point of the non-square map: a missed transpose must not survive.
    [total, failed] = check(total, failed, 'Transpose false is rejected by the size check', ...
        threw(@() load_slm_correction(f, 'Size', [nx ny], 'Transpose', false, ...
                                      'Verbose', false), ...
              'load_slm_correction:sizeMismatch'));

    % ---- 2. sign flip is exactly a negation -----------------------------------
    neg = load_slm_correction(f, 'Size', [nx ny], 'Sign', -1, 'Verbose', false);
    [total, failed] = check(total, failed, 'Sign -1 negates the map', ...
        max(abs(neg(:) + base(:))) < 1e-12);

    % ---- 3. flips act on the right dimension ---------------------------------
    % Compared with a tolerance, not isequal: the loader removes the piston by
    % subtracting mean(phase(:)), and summing the same values in a different order
    % differs in the last bits. That is float addition, not a flip getting it wrong.
    fx = load_slm_correction(f, 'Size', [nx ny], 'FlipX', true, 'Verbose', false);
    [total, failed] = check(total, failed, 'FlipX flips SLM x (dimension 1)', ...
        max(abs(fx(:) - reshape(flip(base, 1), [], 1))) < 1e-12);
    fy = load_slm_correction(f, 'Size', [nx ny], 'FlipY', true, 'Verbose', false);
    [total, failed] = check(total, failed, 'FlipY flips SLM y (dimension 2)', ...
        max(abs(fy(:) - reshape(flip(base, 2), [], 1))) < 1e-12);

    % ---- 4. units: the ambiguous case must refuse to guess --------------------
    small = (rand(nx, ny) * 2 - 1) * 0.5;      % could be radians or waves
    [total, failed] = check(total, failed, 'ambiguous radians-vs-waves errors', ...
        threw(@() load_slm_correction(small, 'Size', [nx ny], 'Verbose', false), ...
              'load_slm_correction:ambiguousUnits'));

    big = (rand(nx, ny) * 2 - 1) * 5;          % beyond half a turn -> radians
    [~, ib] = load_slm_correction(big, 'Size', [nx ny], 'Verbose', false);
    [total, failed] = check(total, failed, 'large values auto-detect as radians', ...
        strcmp(ib.units, 'radians'));

    w = load_slm_correction(small, 'Size', [nx ny], 'Units', 'waves', 'Verbose', false);
    r = load_slm_correction(small, 'Size', [nx ny], 'Units', 'radians', 'Verbose', false);
    [total, failed] = check(total, failed, 'waves is 2*pi times radians', ...
        max(abs(w(:) - 2*pi*r(:))) < 1e-12);

    % ---- 5. bad inputs ------------------------------------------------------
    [total, failed] = check(total, failed, 'a missing file errors', ...
        threw(@() load_slm_correction(fullfile(tmp, 'nope.bmp'), 'Verbose', false), ...
              'load_slm_correction:noFile'));
    [total, failed] = check(total, failed, 'an unsupported extension errors', ...
        threw(@() local_load_ext(tmp), 'load_slm_correction:badExt'));
    [total, failed] = check(total, failed, 'a bad Sign errors', ...
        threw(@() load_slm_correction(big, 'Sign', 0, 'Verbose', false), ...
              'load_slm_correction:badSign'));

    % ---- 6. zernike_phase: defocus has the shape it should ---------------------
    % Z(2,0) is 2*rho^2 - 1 unnormalized: -1 at the centre, +1 at the pupil edge.
    N = 65;
    ph = zernike_phase(N, N, 1, 'Modes', [2 0], 'Mask', true);
    c = (N + 1) / 2;
    [total, failed] = check(total, failed, 'defocus is -1 at the pupil centre', ...
        abs(ph(c, c) + 1) < 1e-9);
    [total, failed] = check(total, failed, 'defocus is +1 at the pupil edge', ...
        abs(ph(N, c) - 1) < 1e-9);
    [total, failed] = check(total, failed, 'defocus is rotationally symmetric', ...
        abs(ph(c, N) - ph(N, c)) < 1e-9);
    [total, failed] = check(total, failed, 'outside the pupil is masked to zero', ...
        ph(1, 1) == 0);

    % Tilt must run along the axis it is named for: Z(1,1) is the cosine term, x.
    tx = zernike_phase(N, N, 1, 'Modes', [1 1], 'Mask', true);
    [total, failed] = check(total, failed, 'Z(1,+1) tilts along SLM x', ...
        tx(N, c) > 0.9 && abs(tx(c, N)) < 1e-9);
    ty = zernike_phase(N, N, 1, 'Modes', [1 -1], 'Mask', true);
    [total, failed] = check(total, failed, 'Z(1,-1) tilts along SLM y', ...
        ty(c, N) > 0.9 && abs(ty(N, c)) < 1e-9);

    % ---- 7. Noll and OSA single-index tables ---------------------------------
    % Noll j=1..15, from Phase 3 item 17's convention. j=1 is piston and is dropped
    % by default, so the table starts at 2.
    noll = [2 1 1; 3 1 -1; 4 2 0; 5 2 -2; 6 2 2; 7 3 -1; 8 3 1; 9 3 -3; ...
            10 3 3; 11 4 0; 12 4 2; 13 4 -2; 14 4 4; 15 4 -4];
    ok = true;
    for i = 1:size(noll, 1)
        j = noll(i,1); n = noll(i,2); m = noll(i,3);
        v = zeros(1, j); v(j) = 1;
        a = zernike_phase(17, 17, v, 'Convention', 'noll');
        b = zernike_phase(17, 17, 1, 'Modes', [n m]);
        if max(abs(a(:) - b(:))) > 1e-9
            fprintf('        Noll j=%d is not (n=%d, m=%+d)\n', j, n, m);
            ok = false;
        end
    end
    [total, failed] = check(total, failed, 'Noll j=2..15 matches the reference table', ok);

    % OSA/ANSI j=0..14. j=0 is piston.
    osa = [1 1 -1; 2 1 1; 3 2 -2; 4 2 0; 5 2 2; 6 3 -3; 7 3 -1; 8 3 1; ...
           9 3 3; 10 4 -4; 11 4 -2; 12 4 0; 13 4 2; 14 4 4];
    ok = true;
    for i = 1:size(osa, 1)
        j = osa(i,1); n = osa(i,2); m = osa(i,3);
        v = zeros(1, j + 1); v(j + 1) = 1;      % OSA starts at j = 0
        a = zernike_phase(17, 17, v, 'Convention', 'osa');
        b = zernike_phase(17, 17, 1, 'Modes', [n m]);
        if max(abs(a(:) - b(:))) > 1e-9
            fprintf('        OSA j=%d is not (n=%d, m=%+d)\n', j, n, m);
            ok = false;
        end
    end
    [total, failed] = check(total, failed, 'OSA j=0..14 matches the reference table', ok);

    % Normalisation is a real factor, not a no-op: Z(2,0) gains sqrt(3).
    un = zernike_phase(N, N, 1, 'Modes', [2 0], 'Normalized', false);
    nm = zernike_phase(N, N, 1, 'Modes', [2 0], 'Normalized', true);
    [total, failed] = check(total, failed, 'Normalized scales Z(2,0) by sqrt(3)', ...
        abs(nm(c, c) / un(c, c) - sqrt(3)) < 1e-9);

    % ---- 8. coefficients through the loader ----------------------------------
    s = struct('coefficients', 1.5, 'modes', [2 0], 'normalized', false);
    [zb, zi] = load_slm_correction(s, 'Size', [N N], 'Verbose', false);
    [total, failed] = check(total, failed, 'coefficient source needs no transpose', ...
        ~zi.transpose && isequal(size(zb), [N N]));
    [total, failed] = check(total, failed, 'coefficient source is reported as zernike', ...
        strcmp(zi.kind, 'zernike') && strcmp(zi.units, 'radians'));
    [total, failed] = check(total, failed, 'coefficients without Size error', ...
        threw(@() load_slm_correction(s, 'Verbose', false), ...
              'load_slm_correction:sizeRequired'));

    % ---- 9. the loaded map is usable by phase_to_frame -------------------------
    Setup.SLM.pixelmax = 255;
    Setup.SLM.correction = zb;
    frame = phase_to_frame(Setup, zeros(N, N));
    [total, failed] = check(total, failed, 'a loaded correction wraps to a uint8 frame', ...
        isa(frame, 'uint8') && isequal(size(frame), [N N]));

    fprintf('\ntest_slm_correction: %d/%d passed\n', total - failed, total);
    if failed > 0
        error('test_slm_correction:failed', '%d check(s) failed.', failed);
    end
end

% -------------------------------------------------------------------------
function local_load_ext(tmp)
    f = fullfile(tmp, 'corr.xyz');
    fid = fopen(f, 'w'); fwrite(fid, 'x'); fclose(fid);
    load_slm_correction(f, 'Verbose', false);
end

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
