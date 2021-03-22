classdef (ConstructOnLoad = true) pressure < matlab.mixin.SetGet
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent)
        pascal = double(0);
        Pa = double(0);
        kilopascal = double(0);
        megapascal = double(0);
        gigapascal = double(0);
        
        poundforce_squarefoot = double(0);
        poundforce_squareinch = double(0);
        kilopoundforce_squareinch = double(0);
        ksi = double(0);
        
        bar = double(0);
        standard_atmosphere = double(0);
    end
    
    properties(Access=private)
        basepascal = double(0);
    end
    
    methods
        function obj = pressure(parm1)
            
            switch nargin
                case 0
                    % Create empty.
                    obj.basepascal = 0;
                case 1
                    obj.basepascal = parm1;
                otherwise
                   warning('Unhandled ammount of parameters')     
            end
        end
        
        %<--- Metric ------------
        function obj = set.gigapascal(obj,value)
            obj.basepascal = value / 0.000000001;
            obj.Validate
        end
        function value = get.gigapascal(obj)
            value = obj.basepascal * 0.000000001;
        end
        
        function obj = set.megapascal(obj,value)
            obj.basepascal = value / 0.000001;
            obj.Validate
        end
        function value = get.megapascal(obj)
            value = obj.basepascal * 0.000001;
        end
        
        % bar
        function obj = set.bar(obj,value)
            obj.basepascal = value / 0.00001;
            obj.Validate
        end
        function value = get.bar(obj)
            value = obj.basepascal * 0.00001;
        end
        
        function obj = set.kilopascal(obj,value)
            obj.basepascal = value / 0.001;
            obj.Validate
        end
        function value = get.kilopascal(obj)
            value = obj.basepascal * 0.001;
        end
        
        function obj = set.pascal(obj,value)
            obj.basepascal = value;
            obj.Validate
        end
        function value = get.pascal(obj)
            value = obj.basepascal;
        end
        function obj = set.Pa(obj,value)
            obj.basepascal = value;
            obj.Validate
        end
        function value = get.Pa(obj)
            value = obj.basepascal;
        end
        
        %<--- US (Imperial) Measure ------------
        % psf
        function obj = set.poundforce_squarefoot(obj,value)
            obj.basepascal = value / 0.0208854337883712;
            obj.Validate
        end
        function value = get.poundforce_squarefoot(obj)
            value = obj.basepascal * 0.0208854337883712;
        end
        
        % ksi
        function obj = set.kilopoundforce_squareinch(obj,value)
            obj.basepascal = value / 1.450377438972831e-7;
            obj.Validate
        end
        function value = get.kilopoundforce_squareinch(obj)
            value = obj.basepascal * 1.450377438972831e-7;
        end
        function obj = set.ksi(obj,value)
            obj.basepascal = value / 1.450377438972831e-7;
            obj.Validate
        end
        function value = get.ksi(obj)
            value = obj.basepascal * 1.450377438972831e-7;
        end
        
        % psi
        function obj = set.poundforce_squareinch(obj,value)
            obj.basepascal = value / 1.450377438972831e-4;
            obj.Validate
        end
        function value = get.poundforce_squareinch(obj)
            value = obj.basepascal * 1.450377438972831e-4;
        end
        
        % atm
        function obj = set.standard_atmosphere(obj,value)
            obj.basepascal = value / 9.869232667160128e-6;
            obj.Validate
        end
        function value = get.standard_atmosphere(obj)
            value = obj.basepascal * 9.869232667160128e-6;
        end
        
        function Validate(obj)
            
        end
    end
end

