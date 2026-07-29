function test_powerbias_delivery()
%TEST_POWERBIAS_DELIVERY  Invariant behind the generate_holograms_new powerbias fix.
%   Verifies that after the 1/DE inversion and the per-target powerbias multiply,
%   the power actually DELIVERED to each target (post-diffraction) is proportional
%   to powerbias, independent of each target's diffraction efficiency. This is the
%   property the graded-power calibration relies on. Run: >> test_powerbias_delivery
%
%   (Pure arithmetic mirror of generate_holograms_new.m's col-4 handling -- no SLM
%   hardware or CoC needed. A full integration check with a real CoC lives in the
%   end-to-end rig verification.)

    DE   = [0.8, 0.3];    % different diffraction efficiencies (position-dependent)
    bias = [1.0, 3.0];    % we want target 2 to receive 3x the power of target 1

    % --- reproduce generate_holograms_new.m col-4 handling ---
    col4 = 1 ./ DE(:);    % invert DE -> power boost   (generate_holograms_new.m:108)
    col4 = col4 .* bias(:);  % apply powerbias ALWAYS  (the fix at ~line 112)

    % Gerchberg-Saxton realizes requested intensity ~ col4 (normalised); diffraction
    % then attenuates target i by DE_i, so delivered_i ~ (col4_i / sum) * DE_i.
    requested = col4 / sum(col4);
    delivered = requested .* DE(:);

    ratio_delivered = delivered / sum(delivered);
    ratio_bias      = bias(:) / sum(bias);

    assert(max(abs(ratio_delivered - ratio_bias)) < 1e-12, ...
        'delivered power must be proportional to powerbias, independent of DE');

    % And with flat powerbias, equal delivered power to every target (DE-corrected).
    col4f = (1 ./ DE(:)) .* ones(2, 1);
    deliveredf = (col4f / sum(col4f)) .* DE(:);
    assert(max(abs(deliveredf - mean(deliveredf))) < 1e-12, ...
        'flat powerbias must deliver equal power to all targets');

    disp('PASS test_powerbias_delivery (delivered power ~ powerbias, DE-independent).');
end
