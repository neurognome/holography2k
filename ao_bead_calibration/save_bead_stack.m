function outpaths = save_bead_stack(outdir, stem, dataUZ, bgd, refs, meta)
% SAVE_BEAD_STACK  Write the AO bead-calibration z-stack + references + metadata.
%
%   outpaths = save_bead_stack(outdir, stem, dataUZ, bgd, refs, meta)
%
% There is no TIFF writer anywhere in holography2k, so this is it. Everything is
% written for a downstream Python (tifffile) analysis pipeline.
%
% WHY float32, not uint8/uint16: the Basler is an 8-bit camera (bascam.camMax =
% 255), but each plane is the MEAN of several frames, so the useful data is
% fractional. We therefore save the raw averaged planes as 32-bit float, exactly
% as measured -- NO smoothing, NO background subtraction, NO clipping (all three
% corrupt phase retrieval / deconvolution; README item 11). Background and
% reference frames are saved separately so Genesis can subtract if she chooses.
%
% INPUTS
%   outdir   destination folder (created if needed).
%   stem     filename stem, e.g. '260814_beadcal_1030'.
%   dataUZ   H-by-W-by-nZ raw averaged z-stack (double). Axis order [y x z].
%   bgd      H-by-W beam-off background (double), or [] to skip.
%   refs     struct of extra 2-D/3-D reference images to write, e.g.
%              refs.reference_blank  (blank-SLM frame, README item 14a)
%              refs.reference_nolut  (existing-correction-disabled frame, item 14b)
%            Each field -> '<stem>_<fieldname>.tif'. Pass struct() for none.
%   meta     struct of metadata (README item 15). Written verbatim to JSON + MAT,
%            with acquisition/layout fields added here.
%
% OUTPUT
%   outpaths struct of the files written (stack, background, refs, json, mat).

if nargin < 5 || isempty(refs), refs = struct(); end
if nargin < 6 || isempty(meta), meta = struct(); end

if ~isfolder(outdir), mkdir(outdir); end
outpaths = struct();

% ---- main z-stack ------------------------------------------------------------
stack_path = fullfile(outdir, [stem '_grid_stack.tif']);
write_float_tiff(stack_path, single(dataUZ));
outpaths.stack = stack_path;
fprintf('save_bead_stack: wrote %s  (%d x %d x %d, float32)\n', ...
    stack_path, size(dataUZ,1), size(dataUZ,2), size(dataUZ,3));

% ---- background --------------------------------------------------------------
if ~isempty(bgd)
    bg_path = fullfile(outdir, [stem '_background.tif']);
    write_float_tiff(bg_path, single(bgd));
    outpaths.background = bg_path;
end

% ---- reference images --------------------------------------------------------
rf = fieldnames(refs);
for i = 1:numel(rf)
    img = refs.(rf{i});
    if isempty(img), continue; end
    rp = fullfile(outdir, [stem '_' rf{i} '.tif']);
    write_float_tiff(rp, single(img));
    outpaths.(rf{i}) = rp;
end

% ---- metadata (item 3 / item 15) ---------------------------------------------
% Fields added here describe the file itself; caller-supplied meta describes the
% acquisition. Anything the caller did not supply is left as a NEEDS_INPUT
% sentinel so it is obvious in the handoff what still has to be filled in.
meta.saved_by       = mfilename;
meta.stem           = stem;
meta.tiff_axis_order = 'YXZ (row, col, plane); planes = TIFF directories';
meta.tiff_dtype     = 'float32';
meta.camera_bit_depth = 8;      % Basler is 8-bit; planes are frame-averages
meta.camera_max_dn  = 255;
meta.n_planes       = size(dataUZ, 3);
meta.image_size_px  = [size(dataUZ,1), size(dataUZ,2)];
meta.files          = outpaths;
meta = ensure_fields(meta, { ...
    'z_step_um','z_range_um','pxPerMu','pixel_size_um','wavelength_nm', ...
    'power_mW','frames_averaged','camera_exposure','na','objective_mag', ...
    'pupil_fill_fraction','pattern_file','pupil_mapping','timestamp'});

json_path = fullfile(outdir, [stem '_metadata.json']);
fid = fopen(json_path, 'w');
fwrite(fid, jsonencode(meta, 'PrettyPrint', true));
fclose(fid);
outpaths.metadata_json = json_path;

mat_path = fullfile(outdir, [stem '_metadata.mat']);
save(mat_path, '-struct', 'meta');
outpaths.metadata_mat = mat_path;

fprintf('save_bead_stack: wrote metadata %s (+ .mat)\n', json_path);

% Warn loudly about anything still unfilled -- these block Genesis's analysis.
todo = fieldnames(meta);
todo = todo(structfun(@(v) ischar(v) && strcmp(v,'NEEDS_INPUT'), meta));
if ~isempty(todo)
    warning('save_bead_stack:incompleteMeta', ...
        'Metadata still NEEDS_INPUT for: %s', strjoin(todo, ', '));
end
end

% =============================================================================
function write_float_tiff(path, img)
% Multipage 32-bit float TIFF via the built-in Tiff class (no toolbox deps).
t = Tiff(path, 'w');
nZ = size(img, 3);
for k = 1:nZ
    ts.ImageLength   = size(img, 1);
    ts.ImageWidth    = size(img, 2);
    ts.Photometric   = Tiff.Photometric.MinIsBlack;
    ts.BitsPerSample = 32;
    ts.SampleFormat  = Tiff.SampleFormat.IEEEFP;
    ts.SamplesPerPixel = 1;
    ts.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    ts.Compression   = Tiff.Compression.None;
    setTag(t, ts);
    write(t, img(:,:,k));
    if k < nZ
        writeDirectory(t);
    end
end
close(t);
end

% =============================================================================
function s = ensure_fields(s, names)
for i = 1:numel(names)
    if ~isfield(s, names{i}) || isempty(s.(names{i}))
        s.(names{i}) = 'NEEDS_INPUT';
    end
end
end
