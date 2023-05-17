classdef HoloeyePLUTO < SLM
    properties
        bit_depth = 8
        % max_transients = 10
        % external_pulse = 1;
        timeout_ms = 5000;
        Nx = 1920
        Ny = 1080
        psSLM = 8e-6;       % meters    SLM pixel dimensions
        wait_for_trigger = 0
        state = 0;
        pixelmax = 255;
        true_frames = 3;

        lut_file = '' % 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6490_at1064.LUT';
        lut
    end

    methods
        function obj = HoloeyePLUTO()
            obj = obj@SLM();
        end

        function start(obj)
            heds_init_slm;
            
            % get lut
            try 
                obj.lut = heds_load_data(importdata(obj.lut_file)); % should be a vector, uint8 of mappings
            catch
                obj.lut = [];
                fprintf('No LUT found at %s\n', obj.lut_file);
            end
        end

        function stop(obj)
            heds_close_slm;
        end

        function feed(obj, frame)
            data_handle = heds_load_data(frame);
            if ~isempty(obj.lut)
                % apply LUT if exists
                heds_datahandle_set_lookuptable(data_handle, obj.lut);
            end
            % show the hologram
            heds_show_datahandle(data_handle);
        end
    end
end