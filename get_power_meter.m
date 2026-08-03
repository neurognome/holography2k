function pm = get_power_meter(wavelength)
%GET_POWER_METER Connect to the Thorlabs power meter and set it for this wavelength.
%   pm = GET_POWER_METER(wavelength) returns a CONNECTED ThorlabsPowerMeter, with
%   its sensor wavelength already set. wavelength is the laser's nominal 900, 1100
%   or 1030 -- the same argument get_slm takes.
%
%   Wraps the connect handshake of the driver
%   (github.com/Tinyblack/Matlab-Driver-for-Thorlabs-power-meter), which is easy to
%   get wrong two ways, and was, in every script here:
%     1) the constructor alone gives you a DISCOVERY STUB, not a usable meter. You
%        must listdevices() then connect(), and connect returns a *copy* -- its
%        return value has to be assigned or you keep talking to the stub.
%     2) the setter names and units are not what you would guess:
%        setWaveLength (capital L), setAverageTime in SECONDS, setTimeout in
%        MILLISECONDS. Averaging and timeout stay with the caller, since those are
%        per-script knobs; only wavelength is a property of which laser you're on.
%
%   It is a .NET wrapper, not VISA: Thorlabs Optical Power Monitor must be installed
%   for Thorlabs.TLPM_64.Interop.dll. This errors early and by name rather than
%   letting an undefined-function surface in the middle of a half-finished sweep.
%
%   Reads are in WATTS via the meterPowerReading property. Always check
%   meterPowerUnit -- a head left in dBm from a previous session returns numbers
%   that look perfectly plausible and silently ruin a calibration.
%
%   See also: get_slm, AutoLaserPowerCalib_HWP, calibrate_DE_powermeter

    % Argument first, environment second, so a typo'd wavelength reports itself as
    % one even on a machine where the driver is missing too.
    switch wavelength
        case 900
            set_nm = 960;   % deliberately 960, as these calibrations always used
        case 1100
            set_nm = 1100;
        case 1030
            set_nm = 1030;
        otherwise
            error('get_power_meter:badWavelength', ...
                'Power meter wavelength not found (900, 1100, 1030).');
    end

    if isempty(which('ThorlabsPowerMeter'))
        error('get_power_meter:noDriver', ...
            ['ThorlabsPowerMeter is not on the MATLAB path.\n' ...
             'It is not part of this repo or holodaq -- it is the driver from\n' ...
             '  github.com/Tinyblack/Matlab-Driver-for-Thorlabs-power-meter\n' ...
             'Clone it anywhere inside this checkout (holo_paths genpaths the\n' ...
             'whole tree, so any subfolder works) or install the File Exchange\n' ...
             'add-on, and install Thorlabs Optical Power Monitor for its DLLs.']);
    end

    stub = ThorlabsPowerMeter();
    desc = stub.listdevices();      % also the first thing that touches the DLLs
    if isempty(desc)
        error('get_power_meter:noDevice', ...
            ['No Thorlabs power meter found. Check the head is plugged in and ' ...
             'that no other program (Optical Power Monitor itself) holds it.']);
    end

    pm = stub.connect(desc);        % MUST assign: connect returns a copy
    pm.setWaveLength(set_nm);

    fprintf('Loaded power meter at %dnm.\n', set_nm);
end
