classdef square_area < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        square_nanometer = double(0);
        square_micrometer = double(0);
        square_millimeter = double(0);
        square_centimeter = double(0);
        square_meter = double(0);
        sm = double(0);
        hectare = double(0);
        square_kilometer = double(0);
        
        square_mile = double(0);
        square_yard = double(0);
        square_foot = double(0);
        square_inch = double(0);
    end
    
    properties(Access=private)
        square_basemetre = double(0)
    end
   
    methods
        function obj = square_area(parm1)
            switch nargin
                case 0
                    % Create empty.
                    obj.square_basemetre = 0;
                case 1
                    obj.square_basemetre = parm1;
                otherwise
                   warning('Unhandled ammount of parameters')     
            end
        end
        
        %<--- Metric ------------
        % nm^2
        function obj = set.square_nanometer(obj,value)
            obj.square_basemetre = value / 1E18;
            obj.Validate
        end
        function value = get.square_nanometer(obj)
            value = obj.square_basemetre * 1E18;
        end
        
        % um^2
        function obj = set.square_micrometer(obj,value)
            obj.square_basemetre = value / 1E12;
            obj.Validate
        end
        function value = get.square_micrometer(obj)
            value = obj.square_basemetre * 1E12;
        end
        
        % mm^2
        function obj = set.square_millimeter(obj,value)
            obj.square_basemetre = value / 1E6;
            obj.Validate
        end
        function value = get.square_millimeter(obj)
            value = obj.basemetre * 1E6;
        end
        
        % cm^2
        function obj = set.square_centimeter(obj,value)
            obj.basemetre = value / 10000;
            obj.Validate
        end
        function value = get.square_centimeter(obj)
            value = obj.square_basemetre * 10000;
        end
        
        % m^2
        function obj = set.square_meter(obj,value)
            obj.square_basemetre = value;
            obj.Validate
        end
        function value = get.square_meter(obj)
            value = obj.square_basemetre;
        end
        function obj = set.sm(obj,value)
            obj.square_basemetre = value;
            obj.Validate
        end
        function value = get.sm(obj)
            value = obj.square_basemetre;
        end
        
        % ha
        function obj = set.hectare(obj,value)
            obj.square_basemetre = value / 0.0001;
            obj.Validate
        end
        function value = get.hectare(obj)
            value = obj.square_basemetre * 0.0001;
        end
        
        % km^2
        function obj = set.square_kilometer(obj,value)
            obj.square_basemetre = value / 0.000001;
            obj.Validate
        end
        function value = get.square_kilometer(obj)
            value = obj.basemetre * 0.000001;
        end
        
        %<--- US (Imperial) Measure ------------
        % sq mi
        function obj = set.square_mile(obj,value)
            obj.basemetre = value / 3.861021585424458e-7;
            obj.Validate
        end
        function value = get.square_mile(obj)
            value = obj.square_basemetre * 3.861021585424458e-7;
        end
        
        % sq yd
        function obj = set.square_yard(obj,value)
            obj.square_basemetre = value / 1.19599004630108;
            obj.Validate
        end
        function value = get.square_yard(obj)
            value = obj.square_basemetre * 1.19599004630108;
        end
        
        % sq ft
        function obj = set.square_foot(obj,value)
            obj.square_basemetre = value / 10.76391041670972;
            obj.Validate
        end
        function value = get.square_foot(obj)
            value = obj.square_basemetre * 10.76391041670972;
        end
        
        % sq in
        function obj = set.square_inch(obj,value)
            obj.square_basemetre = value / 1550.0031000062;
            obj.Validate
        end
        function value = get.square_inch(obj)
            value = obj.square_basemetre * 1550.0031000062;
        end
        
        
        
        
        function Validate(obj)
            % Function to validate the data of the object. A length cannot be a
            % negative number. Raise warning if the object contains a negative
            % number.
            if obj.square_basemetre < 0
                warning('Area is less than zero!')
            end
        end
        
    end
    
end

