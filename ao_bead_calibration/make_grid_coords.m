function gridCoords = make_grid_coords(source, wavelength, opts)
% MAKE_GRID_COORDS  Build/load the N-spot hologram grid for AO bead calibration.
%
%   gridCoords = make_grid_coords(source, wavelength, opts)
%
% Returns an N-by-4 array [x y z power] in NORMALIZED SLM coordinates
% (x,y in [0,1], z = SLM defocus units, power = relative weight), which is
% exactly what function_Make_3D_SHOT_Holos() consumes.
%
% Genesis provides the calibration pattern as a "target coordinate list, grid
% spacing, coordinate frame" (README item 2). This helper accepts that file in
% either coordinate frame, and falls back to a synthetic grid so the acquisition
% script is runnable BEFORE her file arrives.
%
% INPUTS
%   source     One of:
%              []                      -> synthetic ndgrid fallback (default)
%              path to .mat            -> variable 'coords' (N-by-3/4) or 'gridCoords'
%              path to .csv/.txt       -> numeric N-by-3/4, columns x y z [power]
%              N-by-3/4 numeric array  -> used directly
%   wavelength SLM wavelength (900 / 1030 / 1100). Only used when frame='SI',
%              to pick the correct calibration for function_SItoSLM.
%   opts       struct, optional:
%     .frame        'SLM' (default) coords already normalized-SLM, or
%                   'SI' coords are in ScanImage/FOV frame and must be converted.
%     .power        default per-target power weight if source has only 3 cols (default 1).
%     .n            number of points for the synthetic fallback (default 20).
%     .xrange       [lo hi] normalized-SLM x extent for fallback (default [0.25 0.75]).
%     .yrange       [lo hi] normalized-SLM y extent for fallback (default [0.25 0.75]).
%     .z            defocus plane for fallback (default 0).
%     .calib_dir    override calib dir for find_latest_calib (frame='SI' only).
%
% NOTE ON DROPPED POINTS: function_Make_3D_SHOT_Holos silently ignores any
% target with x or y outside [0,1] (it prints a "Bad Warning"). We validate here
% and error out instead, so N spots asked == N spots projected == N PSFs measured.

if nargin < 1, source = []; end
if nargin < 2 || isempty(wavelength), wavelength = 1030; end
if nargin < 3, opts = struct(); end

def = struct('frame','SLM', 'power',1, 'n',20, ...
             'xrange',[0.25 0.75], 'yrange',[0.25 0.75], 'z',0, 'calib_dir','');
opts = merge_defaults(opts, def);

% ---- 1. obtain raw coordinates -----------------------------------------------
if isempty(source)
    coords = synthetic_grid(opts);
    fprintf('make_grid_coords: no source given, using synthetic %d-point grid.\n', size(coords,1));
elseif isnumeric(source)
    coords = source;
elseif ischar(source) || isstring(source)
    coords = load_coords_file(char(source));
else
    error('make_grid_coords:badSource', 'source must be [], a numeric array, or a file path.');
end

% ---- 2. ensure a power column ------------------------------------------------
if size(coords,2) == 3
    coords(:,4) = opts.power;
elseif size(coords,2) ~= 4
    error('make_grid_coords:badShape', 'coords must be N-by-3 or N-by-4, got N-by-%d.', size(coords,2));
end

% ---- 3. convert SI -> SLM frame if requested ---------------------------------
if strcmpi(opts.frame, 'SI')
    if isempty(opts.calib_dir)
        try
            calib_dir = rig_remote_get('paths.calib_dir', '');
        catch
            calib_dir = '';
        end
    else
        calib_dir = opts.calib_dir;
    end
    CoC = find_latest_calib(wavelength, calib_dir);   % holorequest/find_latest_calib.m
    SLMXYZP = function_SItoSLM(coords(:,1:3), CoC);    % -> [x y z DE]
    coords  = [SLMXYZP(:,1:3), coords(:,4)];           % keep the requested power weights
    fprintf('make_grid_coords: converted %d targets SI->SLM via find_latest_calib(%d).\n', ...
            size(coords,1), wavelength);
end

% ---- 4. validate normalized-SLM range (fail loud, do NOT silently drop) ------
bad = coords(:,1) < 0 | coords(:,1) > 1 | coords(:,2) < 0 | coords(:,2) > 1;
if any(bad)
    error('make_grid_coords:outOfRange', ...
        ['%d of %d targets fall outside the normalized SLM range [0,1] in x/y ' ...
         '(rows: %s). function_Make_3D_SHOT_Holos would silently drop these, ' ...
         'changing the spot count. Fix the coordinate list or the frame setting.'], ...
        sum(bad), size(coords,1), mat2str(find(bad)'));
end

gridCoords = coords;
fprintf('make_grid_coords: returning %d targets, x[%.3f %.3f] y[%.3f %.3f] z[%.3g %.3g].\n', ...
    size(gridCoords,1), min(gridCoords(:,1)), max(gridCoords(:,1)), ...
    min(gridCoords(:,2)), max(gridCoords(:,2)), min(gridCoords(:,3)), max(gridCoords(:,3)));
end

% =============================================================================
function coords = synthetic_grid(opts)
% Square-ish grid of opts.n points over the requested normalized-SLM extent,
% at a single defocus plane. Kept away from the SLM edges and the zero order.
side = round(sqrt(opts.n));
xs = linspace(opts.xrange(1), opts.xrange(2), side);
ys = linspace(opts.yrange(1), opts.yrange(2), side);
[X, Y] = ndgrid(xs, ys);
coords = [X(:), Y(:), opts.z*ones(numel(X),1)];
coords = coords(1:min(opts.n, size(coords,1)), :);   % trim to exactly n
end

% =============================================================================
function coords = load_coords_file(fpath)
if ~isfile(fpath)
    error('make_grid_coords:noFile', 'Coordinate file not found: %s', fpath);
end
[~,~,ext] = fileparts(fpath);
switch lower(ext)
    case '.mat'
        S = load(fpath);
        if isfield(S,'gridCoords'), coords = S.gridCoords;
        elseif isfield(S,'coords'), coords = S.coords;
        else
            fn = fieldnames(S);
            num = fn(structfun(@(v) isnumeric(v) && size(v,2)>=3, S));
            if isscalar(num), coords = S.(num{1});
            else
                error('make_grid_coords:matVar', ...
                    'Could not find an N-by-3/4 numeric variable in %s (looked for ''gridCoords''/''coords'').', fpath);
            end
        end
    case {'.csv', '.txt', '.tsv'}
        coords = readmatrix(fpath);
    otherwise
        error('make_grid_coords:ext', 'Unsupported coordinate file type: %s', ext);
end
coords = double(coords);
end

% =============================================================================
function s = merge_defaults(s, def)
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(s, f{i}) || isempty(s.(f{i}))
        s.(f{i}) = def.(f{i});
    end
end
end
