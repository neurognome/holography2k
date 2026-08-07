function [opts, mode, gui] = parse_holo_listener_args(varargin)
%PARSE_HOLO_LISTENER_ARGS Split start_holo_listener's arguments into opts + mode + gui.
%   [opts, mode, gui] = PARSE_HOLO_LISTENER_ARGS(...)
%
%   The name/value options ('Wavelengths', 'CalibDir', 'Timeout', 'Offline') are
%   passed through untouched in OPTS, so every previously valid invocation keeps
%   meaning what it meant. 'blocking' / 'async' / 'gui' / 'nogui' are recognised as
%   standalone keywords in any position; matching is case-insensitive and tolerates
%   surrounding whitespace.
%
%       ()                            -> opts = {},                   async, gui
%       ('Wavelengths', [1100 900])   -> opts = {'Wavelengths', ...}, async, gui
%       ('nogui')                     -> opts = {},                   async, no gui
%       ('blocking')                  -> opts = {},                   blocking, gui
%       ('Timeout', 1700, 'nogui')    -> opts = {'Timeout', 1700},    async, no gui
%
%   gui only asks for the status window; whether one is possible is
%   start_holo_listener's call (blocking mode has no usable window, because the
%   loop owns the session and a Stop button could never be clicked).
%
%   Only EXACT keyword matches are stripped, so a name/value pair is never eaten:
%   the option names are all longer words, and their values are numbers or paths.
%
%   Split out of start_holo_listener so it can be tested -- that function calls
%   addpath(genpath(...)) on the whole repo and opens SLM hardware, so argument
%   handling was untestable in place. Deliberately the same shape as holodaq's
%   parse_si_listener_args, so both listeners take the same keywords.
%
%   See also: start_holo_listener, HoloListener, test_holo_listener_monitor,
%             parse_si_listener_args (holodaq)

    mode = 'async';
    gui  = true;
    args = varargin;
    keep = true(1, numel(args));
    for i = 1:numel(args)
        if ischar(args{i}) || isstring(args{i})
            switch lower(strtrim(char(args{i})))
                case 'blocking', mode = 'blocking'; keep(i) = false;
                case 'async',    mode = 'async';    keep(i) = false;
                case 'gui',      gui  = true;       keep(i) = false;
                case 'nogui',    gui  = false;      keep(i) = false;
            end
        end
    end
    opts = args(keep);

    % Catch the one way stripping a keyword can corrupt a name/value list: an odd
    % number of leftovers means a name lost its value (or vice versa). Better here,
    % with the original call visible, than as inputParser's generic complaint.
    assert(mod(numel(opts), 2) == 0, 'start_holo_listener:oddArgs', ...
        ['Expected name/value options plus ''blocking''/''async''/''gui''/''nogui''; ' ...
         'got %d leftover\nargument(s), which cannot be name/value pairs. Usage:\n' ...
         '    start_holo_listener(''Wavelengths'', [1100 900], ''nogui'')'], numel(opts));
end
