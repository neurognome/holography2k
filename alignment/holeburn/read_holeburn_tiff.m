function [img, n_used] = read_holeburn_tiff(file, chan_idx, n_chan)
%READ_HOLEBURN_TIFF Average one channel's frames out of a ScanImage tiff.
%   img = READ_HOLEBURN_TIFF(file) averages every page.
%   img = READ_HOLEBURN_TIFF(file, chan_idx, n_chan) averages only the pages
%   belonging to channel chan_idx, given that n_chan channels are interleaved.
%
%   ScanImage writes channels INTERLEAVED, one page per channel per frame, so
%   with two channels saved the green frames are pages 1,3,5,... Averaging every
%   page would silently mix the two PMTs into one image -- and a hole burn is
%   measured as a small local intensity DROP, which is exactly the kind of signal
%   that a wrong average buries.
%
%   Uses imread rather than bigread3 or the ScanImage tiff reader on purpose:
%   this runs on whichever machine is doing the fit, and imread is always there.
%   A calibration is a handful of frames, so the speed does not matter.
%
%   Returns double, not the native integer type: everything downstream subtracts
%   two images, and unsigned subtraction saturating at zero would throw away
%   exactly the negative half of the difference.
%
%   See also: detect_holeburns, fit_holeburn_offsets

    if nargin < 2 || isempty(chan_idx), chan_idx = 1; end
    if nargin < 3 || isempty(n_chan),   n_chan   = 1; end

    assert(exist(file, 'file') == 2, 'read_holeburn_tiff:noFile', ...
        'No such tiff:\n  %s', file);
    assert(n_chan >= 1 && chan_idx >= 1 && chan_idx <= n_chan, ...
        'read_holeburn_tiff:channel', ...
        'chan_idx %d is not within 1..%d.', chan_idx, n_chan);

    info = imfinfo(file);
    n_pages = numel(info);
    pages = chan_idx:n_chan:n_pages;
    assert(~isempty(pages), 'read_holeburn_tiff:noPages', ...
        ['%s has %d page(s), so there is no page for channel %d of %d.\n' ...
         'Was more than one channel actually saved?'], file, n_pages, chan_idx, n_chan);

    acc = [];
    for k = pages
        f = double(imread(file, k));
        if isempty(acc)
            acc = zeros(size(f));
        else
            assert(isequal(size(f), size(acc)), 'read_holeburn_tiff:ragged', ...
                'Page %d of %s is a different size from page %d.', k, file, pages(1));
        end
        acc = acc + f;
    end

    img = acc / numel(pages);
    n_used = numel(pages);
end
