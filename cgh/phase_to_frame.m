function Hologram = phase_to_frame(Setup, phase)
%PHASE_TO_FRAME Wrap a CGH phase map to the SLM's 8-bit frame, adding any static correction.
%   Hologram = PHASE_TO_FRAME(Setup, phase) is the ONE place phase (radians)
%   becomes the uint8 image the panel displays. It replaces the line
%
%       Hologram = uint8(floor((Setup.SLM.pixelmax*(Holo.phase+pi)/(2*pi))));
%
%   which was copy-pasted identically into four CGH functions
%   (function_Make_3D_SHOT_Holos, _disks_KCZ, _tuneZ, function_makeAB). Those four
%   are upstream of EVERY hologram this repo makes -- the listener's compile path,
%   new_alignment_tools, get_psf_no_power, calibrate_DE_powermeter, the AO bead
%   scripts and the tests -- so a base-phase layer added here reaches all of them
%   and cannot be applied to stimulation but forgotten for the alignment and DE
%   measurements that the power model is calibrated against.
%
%   THE STATIC CORRECTION. If Setup.SLM.correction is a non-empty Nx-by-Ny array of
%   RADIANS it is added to every hologram before wrapping. This is Phase 4 item 18
%   of the AO bead protocol ("add the correction phase to the CGH before wrapping,
%   as a static layer applied to every hologram"), and it is the hook
%   ao_bead_calibration/README.md said would be needed. Build the array with
%   load_slm_correction, which handles Meadowlark-style .bmp files and Zernike
%   coefficients and reports what it did.
%
%   WHY THE ADDITION IS IN RADIANS AND NOT IN DN. Two traps, both silent:
%
%     * uint8 arithmetic in MATLAB SATURATES, it does not wrap.
%       uint8(200) + uint8(100) is 255, not 44. Adding a correction to an
%       already-compiled uint8 stack therefore clips a few percent of the panel to
%       255 -- a flat region of maximum phase, which shows up as a diffuse haze in
%       the reconstruction rather than as an obvious failure.
%     * the DN modulus is pixelmax, NOT 256. The forward map below puts 2*pi
%       across pixelmax (255) counts, so DN and DN+255 are the same phase. Wrapping
%       a DN-domain sum with mod(...,256) mis-scales by one count in 255 for every
%       pixel that wraps.
%
%   Adding before the wrap sidesteps both: mod(...,2*pi) is exact in radians.
%
%   BYTE-IDENTICAL WITHOUT A CORRECTION. The no-correction branch is the original
%   expression, unchanged, rather than the mod() form. They differ at exactly one
%   input -- phase == +pi gives DN 255 through the legacy expression and DN 0
%   through mod(), which are the same physical phase but not the same byte. Keeping
%   the old line means enabling this file changes nothing you could measure until a
%   correction is actually set, so it can go on the rig ahead of the mask.
%
%   Reconstruction is deliberately NOT corrected: function_VolumeIntensity is fed
%   the raw Holo.phase by every caller, because it predicts the focus for the phase
%   you COMMANDED. Feeding it the correction would show a simulated aberration that
%   the real optics is supposed to be cancelling.
%
%   See also: load_slm_correction, zernike_phase, function_Make_3D_SHOT_Holos

    base = [];
    if isfield(Setup, 'SLM') && isfield(Setup.SLM, 'correction')
        base = Setup.SLM.correction;
    end

    if isempty(base)
        % The original line, preserved exactly. See BYTE-IDENTICAL above.
        Hologram = uint8(floor((Setup.SLM.pixelmax*(phase+pi)/(2*pi))));
        return
    end

    % A wrong orientation cannot be caught by a size check on a square panel (the
    % scope2k Meadowlark is 1024x1024), so this assert only catches the gross case.
    % load_slm_correction prints the transpose/flip it applied for that reason.
    assert(isequal(size(base), size(phase)), 'phase_to_frame:sizeMismatch', ...
        ['Setup.SLM.correction is %s but the phase map is %s.\nThe correction ' ...
         'must be one value per SLM pixel, laid out the same way as the hologram:\n' ...
         'dimension 1 is SLM x, dimension 2 is SLM y (function_Make_3D_SHOT_Holos\n' ...
         'builds Masks(x,y)), which is the TRANSPOSE of what imread returns for an\n' ...
         'image file. Load it with load_slm_correction.'], ...
        local_dims(base), local_dims(phase));

    assert(all(isfinite(base(:))), 'phase_to_frame:nonFiniteCorrection', ...
        ['Setup.SLM.correction contains %d non-finite value(s). NaN propagates ' ...
         'through the\nwrap into uint8(NaN) == 0, so a single bad pixel would ' ...
         'become a silent flat spot\nrather than an error.'], ...
        sum(~isfinite(base(:))));

    Hologram = uint8(floor(Setup.SLM.pixelmax * ...
                           mod(double(phase) + double(base) + pi, 2*pi) / (2*pi)));
end

% -------------------------------------------------------------------------
function s = local_dims(a)
    s = strjoin(arrayfun(@(n) sprintf('%d', n), size(a), 'UniformOutput', false), 'x');
end
