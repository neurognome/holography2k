classdef MeadowlarkOneK < SLM
    properties
        bit_depth = 12
        num_boards_found = libpointer('uint32Ptr', 0)
        constructed_okay = libpointer('int32Ptr', 0)
        is_nematic_type = 1
        RAM_write_enable = 1
        use_GPU = 0
        max_transients = 10
        external_pulse = 1;
        timeout_ms = 5000;
        Nx = 1024
        Ny = 1024
        wait_for_trigger = 0
        state = 0;
        pixelmax = 255;
        true_frames = 3;

        lut_file = 'C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\LUT Files\slm6490_at1064.LUT';
        reg_lu = libpointer('string');
    end

    methods
        function obj = MeadowlarkOneK()
            obj = obj@SLM();
        end

        function start(obj)
            if obj.state ~= 1
                if ~libisloaded('Blink_C_wrapper')
                    if SLM.is_onek
                        loadlibrary("C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK\Blink_C_wrapper.dll", "C:\Program Files\Meadowlark Optics\Blink OverDrive Plus\SDK\Blink_C_wrapper.h");
                    else
                        loadlibrary('Blink_C_wrapper.dll', 'Blink_C_wrapper.h');
                    end
                end

                calllib('Blink_C_wrapper', 'Create_SDK', obj.bit_depth, obj.num_boards_found, obj.constructed_okay, obj.is_nematic_type, obj.RAM_write_enable, obj.use_GPU, obj.max_transients, obj.reg_lut);
                if obj.constructed_okay.value ~= 1
                    disp('Blink SDK was not successfully constructed');
                    disp(calllib('Blink_C_wrapper', 'Get_last_error_message'));
                    calllib('Blink_C_wrapper', 'Delete_SDK');
                else
                    disp('Blink SDK was successfully constructed');
                    fprintf('Found %u SLM controller(s)\n', obj.num_boards_found.value);
                    % A linear LUT must be loaded to the controller for OverDrive Plus
                    calllib('Blink_C_wrapper', 'Load_LUT_file', 1, obj.lut_file);
                end

                if obj.wait_for_trigger == 1
                    disp('Triggering active')
                else
                    disp('Triggering not active')
                end
                disp('SLM ready to fire!')
                    obj.state = 1;
            else
                disp('SLM already loaded.')
            end
        end

        function stop(obj)
            try
                calllib('Blink_C_wrapper', 'Delete_SDK');
                disp('SLM has just been successfully turned OFF')
                if libisloaded('Blink_C_wrapper')
                    unloadlibrary('Blink_C_wrapper');
                    disp('SLM command library successfully unloaded.')
                end
            catch
                disp('Either SLM is already off, or warning for other issues...')
            end
            
            obj.state = 0;
        end

        function out = feed(obj, frame)
            calllib('Blink_C_wrapper', 'Write_image', 1, frame, obj.Nx*obj.Ny, obj.wait_For_trigger,0, 1, 0, obj.timeout_ms);
            out = calllib('Blink_C_wrapper', 'ImageWriteComplete', 1, obj.timeout_ms);
        end
    end
end