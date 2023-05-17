classdef HoloeyePLUTO < SLM
    properties
    end

    methods
        function obj = HoloeyePLUTO()
            obj = obj@SLM();
            add_heds_path; % importing the SDK
        end

        function start(obj)
            heds_init_slm;
        end

        function stop(obj)
            heds_close_slm;
        end

        function feed(obj, frame)
            hed_show_data(frame)
        end
    end
end