classdef (ConstructOnLoad = true) torque < matlab.mixin.SetGet
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        DEBUG = 0;
    end
    
    properties(Dependent)
        newton_meter = double(0);
        Nm = double(0);
        kilonewton_meter = double(0);
        
        inch_poundforce = double(0);
        foot_poundforce = double(0);
    end
    
    properties(Access=private)
        basenewton_meter = double(0);
    end
    
    
    methods
        function obj = torque(parm1)
            
            if nargin == 0
                obj.basenewton_meter = 0;
            else
                obj.basenewton_meter = parm1;
            end
            if obj.DEBUG;fprintf('torque:constructor = %0.0f  %0.0f\n',obj.basenewton_meter,nargin);end;
        end
        
        %<--- Metric ------------
        % Nm
        function obj = set.newton_meter(obj,value)
            obj.basenewton_meter = value;
            obj.Validate;
        end
        function value = get.newton_meter(obj)
            value = obj.basenewton_meter;
        end
        function obj = set.Nm(obj,value)
            obj.basenewton_meter = value;
            obj.Validate;
        end
        function value = get.Nm(obj)
            value = obj.basenewton_meter;
        end
        
        function obj = set.kilonewton_meter(obj,value)
            obj.basenewton_meter = value / 0.001;
            obj.Validate;
        end
        function value = get.kilonewton_meter(obj)
            value = obj.basenewton_meter * 0.001;
        end
        
        %<--- US (Imperial) Measure ------------
        % in.lbf
        function obj = set.inch_poundforce(obj,value)
            obj.basenewton_meter = value / 8.850745791327184;
            obj.Validate;
        end
        function value = get.inch_poundforce(obj)
            value = obj.basenewton_meter * 8.850745791327184;
        end
        
        function obj = set.foot_poundforce(obj,value)
            obj.basenewton_meter = value / 0.7375621492772656;
            obj.Validate;
        end
        function value = get.foot_poundforce(obj)
            value = obj.basenewton_meter * 0.7375621492772656;
        end
        
        function Validate(obj)
            
        end
    end
end

