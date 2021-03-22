classdef tracer < water_particle
    %TRACER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        label_x = double(0);
        tails_x =[];
        tails_z =[];
        color_vector =[];
        no_frames = double(0);
        frame_freq = double(0);
        %frame_pointer = double(0);
        frame_length = double(0);
    end
    
    methods
        function obj = tracer(x,z,H,T,h,z_axis_lower,time_offset,label_x,no_frames,frame_length)
            
            if nargin == 0
                %fprintf('tracer (empty)\n');
                super_args = {};
            else
                %fprintf('tracer x:%.2f z:%.2f\n',x,z);
                super_args{1} = x;
                super_args{2} = z;
                super_args{3} = H;
                super_args{4} = T;
                super_args{5} = h;
                super_args{6} = z_axis_lower;
                super_args{7} = time_offset;
            end
            obj@water_particle(super_args{:});
            
            if nargin ~= 0
                obj.label_x = label_x;
                obj.no_frames = no_frames;
                obj.tails_x = zeros(no_frames,1);
                obj.tails_z = zeros(no_frames,1);
                obj.frame_length = frame_length;              
                for i = 1:no_frames
                    obj.tails_x(i) = -1000;
                    obj.tails_z(i) = -1000;
                    obj.color_vector(i) = 1;
                end
                %obj.frame_pointer = 0;
            end
        end
        
        function [do_plot,return_position_x,return_position_z,x_array,z_array] = plot_position(obj,t,time_per_meter)
            [do_plot,return_position_x,return_position_z] = plot_position@water_particle(obj,t,time_per_meter);
            
            %obj.frame_pointer = obj.frame_pointer + 1;
            %obj.tails_x(obj.frame_pointer) = return_position_x;
            %obj.tails_z(obj.frame_pointer) = return_position_z;
            
            for i = 1:obj.no_frames;
                [do_plot_tail,return_position_x_tail,return_position_z_tail] = plot_position@water_particle(obj,t - (i * obj.frame_length),time_per_meter);             
                obj.tails_x(i) = return_position_x_tail;
                obj.tails_z(i) = return_position_z_tail;                
            end         
            x_array = obj.tails_x;
            z_array = obj.tails_z;
            
        end
    end  
end

