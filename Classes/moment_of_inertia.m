classdef moment_of_inertia < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        DEBUG = 0;
    end
    
    properties(Dependent)
        meter = double(0);
        m = double(0);
        centimeter = double(0);
        cm = double(0);
        millimeter = double(0);
        mm = double(0);
        foot = double(0);
        inch = double(0);
    end
    
    properties(Access=private)
        basemeter = double(0);
    end
    
    methods
        function obj = moment_of_inertia(parm1)
            if nargin > 0
                obj.basemeter = parm1;
                if obj.DEBUG
                    obj.show_object;
                end
            end
        end
        
        function obj = set.meter(obj,value)
            obj.basemeter = value;
            obj.Validate
        end
        function value = get.meter(obj)
            value = obj.basemeter;
        end
        function obj = set.m(obj,value)
            obj.basemeter = value;
            obj.Validate
        end
        function value = get.m(obj)
            value = obj.basemeter;
        end
        
        function obj = set.centimeter(obj,value)
            obj.basemeter = value / 1E8;
            obj.Validate
        end
        function value = get.centimeter(obj)
            value = obj.basemeter * 1E8;
        end
        function obj = set.cm(obj,value)
            obj.basemeter = value / 1E8;
            obj.Validate
        end
        function value = get.cm(obj)
            value = obj.basemeter * 1E8;
        end
        
        function obj = set.millimeter(obj,value)
            obj.basemeter = value / 1E12;
            obj.Validate
        end
        function value = get.millimeter(obj)
            value = obj.basemeter * 1E12;
        end
        function obj = set.mm(obj,value)
            obj.basemeter = value / 1E12;
            obj.Validate
        end
        function value = get.mm(obj)
            value = obj.basemeter * 1E12;
        end
        
        % ft
        function obj = set.foot(obj,value)
            obj.basemeter = value / 115.8617675;
            obj.Validate
        end
        function value = get.foot(obj)
            value = obj.basemeter * 115.8617675;
        end
        
        % in
        function obj = set.inch(obj,value)
            obj.basemeter = value / 2402509.61;
            obj.Validate
        end
        function value = get.inch(obj)
            value = obj.basemeter * 2402509.61;
        end
        
        
        function show_object(obj)
            fprintf('moment_of_inertia:constructor : %0.0f\n',obj.basemeter)
        end
        
        function Validate(obj)
            % Function to validate the data of the object. 
        end
        
    end
    
    
end

