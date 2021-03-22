classdef (ConstructOnLoad = true) force < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties(Constant)
        DEBUG = 0;
    end
    
    properties(Dependent)
        newton = double(0);
        N = double(0);
        kilonewton = double(0);
        
        kilopoundforce = double(0);
        poundforce = double(0);
        ounceforce = double(0);
    end
    
    properties(Access=private)
        basenewton = double(0);
    end
    
    methods
        function obj = force(parm1)
            
            if nargin == 0
                obj.basenewton = 0;
            else
                obj.basenewton = parm1;
            end
            if obj.DEBUG;fprintf('force:constructor = %0.0f  %0.0f\n',obj.basenewton,nargin);end;
        end
        
        %<--- Metric ------------
        % kN
        function obj = set.kilonewton(obj,value)
            obj.basenewton = value / 0.001;
            obj.Validate
        end
        function value = get.kilonewton(obj)
            value = obj.basenewton * 0.001;
        end
        
        % N
        function obj = set.newton(obj,value)
            obj.basenewton = value;
            obj.Validate
        end
        function value = get.newton(obj)
            value = obj.basenewton;
        end
        function obj = set.N(obj,value)
            obj.basenewton = value;
            obj.Validate
        end
        function value = get.N(obj)
            value = obj.basenewton;
        end
        
        %<--- US (Imperial) Measure ------------
        % kip, kipf, klbf
        function obj = set.kilopoundforce(obj,value)
            obj.basenewton = value / 2.248089430997105e-4;
            obj.Validate
        end
        function value = get.kilopoundforce(obj)
            value = obj.basenewton * 2.248089430997105e-4;
        end
        
        % lbf
        function obj = set.poundforce(obj,value)
            obj.basenewton = value / 0.2248089430997105;
            obj.Validate
        end
        function value = get.poundforce(obj)
            value = obj.basenewton * 0.2248089430997105;
        end
        
        % ozf
        function obj = set.ounceforce(obj,value)
            obj.basenewton = value / 3.596943089595368;
            obj.Validate
        end
        function value = get.ounceforce(obj)
            value = obj.basenewton * 3.596943089595368;
        end
        
        function Validate(obj)
            % Function to validate the data of the object. 
            
        end
    end
end

