classdef capacity_volume < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        cubic_nanometer = double(0);
        cubic_micrometer = double(0);
        cubic_millimeter = double(0);
        cubic_centimeter = double(0);
        cubic_meter = double(0);
        hectare = double(0);
        cubic_kilometer = double(0);
        
        cubic_mile = double(0);
        cubic_yard = double(0);
        cubic_foot = double(0);
        cubic_inch = double(0);
    end
    
    properties(Access=private)
        basecubic_metre = double(0)
    end
   
    methods
        function obj = capacity_volume(parm1)
            switch nargin
                case 0
                    % Create empty.
                    obj.basecubic_metre = 0;
                case 1
                    obj.basecubic_metre = parm1;
                otherwise
                   warning('Unhandled ammount of parameters')     
            end
        end
        
        %<--- Metric ------------
        
        function obj = set.cubic_nanometer(obj,value)
            obj.basecubic_metre = value / 1E27;
            obj.Validate
        end
        function value = get.cubic_nanometer(obj)
            value = obj.basecubic_metre * 1E27;
        end
        
        function obj = set.cubic_micrometer(obj,value)
            obj.basecubic_metre = value / 1E18;
            obj.Validate
        end
        function value = get.cubic_micrometer(obj)
            value = obj.basecubic_metre * 1E18;
        end
        
        function obj = set.cubic_millimeter(obj,value)
            obj.basecubic_metre = value / 1E9;
            obj.Validate
        end
        function value = get.cubic_millimeter(obj)
            value = obj.basecubic_metre * 1E9;
        end
        
        function obj = set.cubic_centimeter(obj,value)
            obj.basecubic_metre = value / 1000000;
            obj.Validate
        end
        function value = get.cubic_centimeter(obj)
            value = obj.basecubic_metre * 1000000;
        end
        
        function obj = set.cubic_meter(obj,value)
            obj.basecubic_metre = value;
            obj.Validate
        end
        function value = get.cubic_meter(obj)
            value = obj.basecubic_metre;
        end
        
        function obj = set.cubic_kilometer(obj,value)
            obj.basecubic_metre = value / 0.000000001;
            obj.Validate
        end
        function value = get.cubic_kilometer(obj)
            value = obj.basecubic_metre * 0.000000001;
        end
        
        %<--- US (Imperial) Measure ------------
        
        function obj = set.cubic_mile(obj,value)
            obj.basecubic_metre = value / 2.399127586E-10;
            obj.Validate
        end
        function value = get.cubic_mile(obj)
            value = obj.basecubic_metre * 2.399127586E-10;
        end
        
        function obj = set.cubic_yard(obj,value)
            obj.basecubic_metre = value / 1.307950619;
            obj.Validate
        end
        function value = get.cubic_yard(obj)
            value = obj.basecubic_metre * 1.307950619;
        end
        
        function obj = set.cubic_foot(obj,value)
            obj.basecubic_metre = value / 35.314666721;
            obj.Validate
        end
        function value = get.cubic_foot(obj)
            value = obj.basecubic_metre * 35.314666721;
        end
        
        function obj = set.cubic_inch(obj,value)
            obj.basecubic_metre = value / 61023.744094732;
            obj.Validate
        end
        function value = get.cubic_inch(obj)
            value = obj.basecubic_metre * 61023.744094732;
        end
        
        function Validate(obj)
            % Function to validate the data of the object. A length cannot be a
            % negative number. Raise warning if the object contains a negative
            % number.
            if obj.basecubic_metre < 0
                warning('Capacity Volume is less than zero!')
            end
        end
        
    end
    
end

