function slm = get_slm(wavelength)
switch wavelength
    case 900
        slm = MeadowlarkOneK(1, 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6902_at900.lut');
    case 1100
        slm = MeadowlarkOneK(2, 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6490_at1100.lut');
    case 1030
        slm = MeadowlarkOneK(2, 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6490_at1030.lut');
    otherwise
        fprint('SLM wavelength not found (900, 1100, 1030).\n')
        return
end

fprintf('Loaded %dnm SLM.\n', wavelength);
end