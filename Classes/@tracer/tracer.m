classdef tracer < position
    %TRACER subclass of the position class
    %   Provides a way to draw tracers.
    
    properties
        label_x = double(0);
        tails_x =[];
        tails_z =[];
        color_vector =[];
        no_frames = double(0);
        frame_freq = double(0);
        frame_length = double(0);
    end
    
    methods
        function obj = tracer(x,y,z,H,T,h,z_axis_lower,label_x,no_frames,frame_length)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = x;
                super_args{2} = y;
                super_args{3} = z;
                super_args{4} = H;
                super_args{5} = T;
                super_args{6} = h;
                super_args{7} = z_axis_lower;
            end
            obj@position(super_args{:});
            
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
            end
        end
        
        function [do_plot,return_position_x,return_position_y,return_position_z,x_array,z_array] = plot_position(obj,t)
            [do_plot,return_position_x,return_position_y,return_position_z] = plot_position@water_particle(obj,t);                      
            
            for i = 1:obj.no_frames;
                [do_plot_tail,return_position_x_tail,return_position_y_tail,return_position_z_tail] = plot_position@water_particle(obj,t - (i * obj.frame_length));             
                obj.tails_x(i) = return_position_x_tail;
                obj.tails_z(i) = return_position_z_tail;                
            end
            
            x_array = obj.tails_x;
            z_array = obj.tails_z;
            
            %obj.t = t;
            %obj.particle_displacement_amplitudes;
        end
    end  
end

