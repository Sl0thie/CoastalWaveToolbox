classdef position < numeretical_functions
    %POSITION A subclass of NUMERETICAL_FUNCTIONS
    %   This class provides a set of numeretical functions relative to a
    %   poticular position.
    %
    % See also NUMERETICAL_FUNCTIONS TRACER
    
    properties (Access = public)
        index = double(0);                  % Index. Currently unused.
        position_x = double(0);             % Current x particle position.
        position_y = double(0);             % Current y particle position.
        position_z = double(0);             % Current z particle position.
        position_x_original = double(0);    % Rotation point x co-ord.
        position_y_original = double(0);    % Rotation point y co-ord.
        position_z_original = double(0);    % Rotation point z co-ord.
        z_axis_lower = double(0);           % Lowest point on the seafloor.
    end
    
    methods
        function obj = position(x,y,z,H,T,h,z_axis_lower)
        %POSITION constructor for the object.
        %   This function is used by MatLAB to construct this object when
        %   it is created.       
            persistent current_index          
            if nargin == 0 % Object is created empty.
                super_args = {};
            else % Object is with values.
                super_args{1} = x;
                super_args{2} = y;
                super_args{3} = z;
                super_args{4} = H;
                super_args{5} = T;
                super_args{6} = h;
            end
            obj@numeretical_functions(super_args{:});
            
            if nargin ~= 0
                obj.z_axis_lower = z_axis_lower;
                obj.position_x_original = x;
                obj.position_y_original = y;
                obj.position_z_original = z;
            end
            % incremented index for collections.
            obj.index = current_index;current_index = current_index + 1;
        end
        
        function obj = set.position_x_original(obj,value)
            obj.x = value;
            obj.position_x = value;
            obj.position_x_original = value;
        end
        
        function obj = set.position_z_original(obj,value)
            if value <= 0
                obj.z = value;
                obj.position_z = value;
                obj.position_z_original = value;
            else
                error('z value set to position above water level');
            end
        end
        
        function obj = set.position_y_original(obj,value)
            obj.y = value;
            obj.position_y = value;
            obj.position_y_original = value;
        end
        
        function [do_plot,return_position_x,return_position_y,return_position_z] = plot_position(obj,t)
            obj.t = t;
            if obj.position_z_original > -obj.h
                do_plot = 1;
                obj.particle_displacement_amplitudes;
                return_position_x = double(obj.position_x_original) + (obj.ALPHA_disp);
                return_position_z = double(obj.position_z_original) + (obj.BETA_disp);
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

