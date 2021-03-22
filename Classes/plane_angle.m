classdef plane_angle < matlab.mixin.SetGet
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        DEBUG = 0;
    end
    
    properties(Dependent)
        radian = double(0);
        degree = double(0);
        rad = double(0);
        deg = double(0);
    end
    
    properties(Access=private)
        baseradian = double(0);
    end
    
    methods
        function obj = plane_angle(parm1)
            if nargin > 0
                obj.baseradian = parm1;
                if obj.DEBUG
                    obj.show_object;
                end
            end
        end
        
        %<--- Metric ------------
        % rad
        function obj = set.radian(obj,value)
            obj.baseradian = value;
            obj.Validate
        end
        function value = get.radian(obj)
            value = obj.baseradian;
        end
        function obj = set.rad(obj,value)
            obj.baseradian = value;
            obj.Validate
        end
        function value = get.rad(obj)
            value = obj.baseradian;
        end
        
        % deg
        function obj = set.degree(obj,value)
            obj.baseradian = value / 57.295777937149167208732472433712;
            obj.Validate
        end
        function value = get.degree(obj)
            value = obj.baseradian * 57.295777937149167208732472433712;
        end
        function obj = set.deg(obj,value)
            obj.baseradian = value / 57.295777937149167208732472433712;
            obj.Validate
        end
        function value = get.deg(obj)
            value = obj.baseradian * 57.295777937149167208732472433712;
        end
        
        function show_object(obj)
            fprintf('plane_angle:constructor : %0.0f\n',obj.baseradian)
        end
        
        function Validate(obj)
            % Function to validate the data of the object. 
            
        end
    end
    
end

