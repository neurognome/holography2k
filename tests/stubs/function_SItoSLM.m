function out = function_SItoSLM(targets, CoC)
%FUNCTION_SITOSLM  Offline test stub: CoC.att supplies the per-target attenuation.
%
%   The real function (repo root) maps ScanImage coordinates to SLM coordinates using a
%   measured CoC and returns [x y z attenuation], where column 4 is the diffraction
%   efficiency at that spot. Tests only need column 4 to be controllable, so this stub
%   takes it straight from CoC.att and leaves the coordinates alone.
%
%   Lives in tests/stubs/ so it is only on the path when a test puts it there -- it must
%   never shadow the real function during a run. Used by test_pattern_de.
    assert(isstruct(CoC) && isfield(CoC, 'att'), ...
        'function_SItoSLM stub needs CoC.att (per-target attenuation)');
    att = CoC.att(:);
    assert(numel(att) == size(targets, 1), ...
        'CoC.att has %d entries for %d targets', numel(att), size(targets, 1));
    out = [targets, att];
end
