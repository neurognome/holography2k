function c = find_latest_calib(wavelength, calib_dir)
%FIND_LATEST_CALIB Load the newest dated calibration for a wavelength.
%   c = FIND_LATEST_CALIB(wavelength) picks the most recent
%   *_Calib_<wavelength>*.mat in the default calib folder (by file date) and
%   importdata's it. Falls back to the known-good hardcoded file for that
%   wavelength if none match. This replaces the hand-edited switch in
%   MsocketHolorequest2K.m so the holo prime never needs a code edit.
%
%   c = FIND_LATEST_CALIB(wavelength, calib_dir) overrides the folder.

    if nargin < 2 || isempty(calib_dir)
        calib_dir = 'C:\Users\holos\Documents\calibs';
    end

    files = dir(fullfile(calib_dir, sprintf('*_Calib_%d*.mat', wavelength)));
    if isempty(files)
        c = local_fallback(wavelength, calib_dir);
        return
    end
    [~, idx] = max([files.datenum]);   % newest by file modification date
    fname = files(idx).name;
    fprintf('  calib %dnm: %s\n', wavelength, fname);
    c = importdata(fullfile(files(idx).folder, fname));
    c = local_clean(c, wavelength);
end

function c = local_fallback(wavelength, calib_dir)
    switch wavelength
        case 900,  name = '20-Feb-2026_Calib_900_Nikon16x.mat';
        case 1100, name = '23-Feb-2026_Calib_1100_Nikon16x.mat';
        case 1030, name = '23-Jan-2025_Calib_1030.mat';
        case 607,  name = '12-Nov-2024_Calib_607.mat';
        otherwise, error('find_latest_calib: no calibration for %dnm', wavelength);
    end
    fprintf('  calib %dnm: (fallback) %s\n', wavelength, name);
    c = importdata(fullfile(calib_dir, name));
    c = local_clean(c, wavelength);
end

function c = local_clean(c, wavelength)
    % 1030 calib carries stray fit fields the pipeline doesn't want (see the
    % original MsocketHolorequest2K switch).
    if wavelength == 1030 && isstruct(c) && all(isfield(c, {'FitX', 'FitY', 'FitZ'}))
        c = rmfield(c, {'FitX', 'FitY', 'FitZ'});
    end
end
