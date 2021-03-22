classdef water_particle < coastal_object
    %WATER_PARTICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        position_x = double(0);             % Current x particle position.
        position_y = double(0);             % Current y particle position.   
        position_z = double(0);             % Current z particle position.
        position_x_original = double(0);    % Rotation point x co-ord.
        position_y_original = double(0);    % Rotation point y co-ord.
        position_z_original = double(0);    % Rotation point z co-ord.
        %sea_bed_z = double(0);              %     
        z_axis_lower = double(0);           % Lowest point on the seafloor.
    end
    
    methods
        function obj = water_particle(x,y,z,H,T,h,z_axis_lower)
            % Constructor for the object.
            if nargin == 0 % Object is created empty.
                %fprintf('water_particle (empty)\n');
                super_args = {};  
            else % Object is with values.
                %fprintf('water_particle x:%.2f z:%.2f H:%.2f T:%.2f h:%.2f\n',x,z,H,T,h);
                super_args{1} = x;
                super_args{2} = y;
                super_args{3} = z;
                super_args{4} = H;
                super_args{5} = T;
                super_args{6} = h;
            end
            obj@coastal_object(super_args{:});
            
            if nargin ~= 0
               obj.z_axis_lower = z_axis_lower;
               obj.position_x_original = x;
               obj.position_y_original = y;
               obj.position_z_original = z;             
            end
        end
        
        function [do_plot,return_position_x,return_position_y,return_position_z] = plot_position(obj,t)       
            obj.t = t;
            if obj.position_z_original > -obj.h
                do_plot = 1;
%                 obj.particle_velocites;  
%                 return_position_x = double(obj.position_x_original + obj.u);
%                 return_position_z = double(obj.position_z_original + obj.w);               
%                 obj.velocity_amplitudes;
%                 return_position_x = double(obj.position_x_original) + (obj.u_hat);
%                 return_position_z = double(obj.position_z_original) + (obj.w_hat);               
                obj.particle_displacement_amplitudes;
                return_position_x = double(obj.position_x_original) + (obj.ALPHAdisp);
                return_position_z = double(obj.position_z_original) + (obj.BETAdisp);
                return_position_y = obj.y;
            else
                do_plot = 0;
                return_position_x = obj.position_x_original;
                return_position_z = obj.position_z_original;
                return_position_y = obj.y;
            end          
        end             
    end
end

