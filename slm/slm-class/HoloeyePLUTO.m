classdef HoloeyePLUTO < SLM
    properties
        bit_depth = 8
        % max_transients = 10
        % external_pulse = 1;
        timeout_ms = 5000;
        wait_for_trigger = 0
        state = 0;
        pixelmax = 255;
        true_frames = 3;

        lut_file = '' % 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6490_at1064.LUT';
        lut

        trigger = [];
    end

    methods
        function obj = HoloeyePLUTO()
            obj = obj@SLM();
            obj.Nx = 1080;
            obj.Ny = 1920;
            obj.psSLM = 8e-6;
        end

        function start(obj)
            % need an arduino
            obj.trigger = arduino('COM4'); % get out of hardcode later...
            obj.trigger.writeDigitalPin('D2', 0);
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

        function out = feed(obj, frame)
            data_handle = heds_load_data(frame);
            if ~isempty(obj.lut)
                % apply LUT if exists
                heds_datahandle_set_lookuptable(data_handle, obj.lut);
            end
            % show the hologram
            if obj.wait_for_trigger
                obj.wait();
            end
            out = heds_show_datahandle(data_handle.id);
            obj.send_flip();
        end

        function wait(obj)
            tic
            while ~obj.trigger.readDigitalPin('D13') 
                elapsed = toc;
                disp(elapsed)
            end
            disp('firing!');
        end

        function send_flip(obj)
            obj.trigger.writeDigitalPin('D2', 1);
            drawnow();
            obj.trigger.writeDigitalPin('D2', 0);
        end
    end
end