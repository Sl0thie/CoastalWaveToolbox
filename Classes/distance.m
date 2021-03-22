classdef distance < matlab.mixin.SetGet
    %…DISTANCE CLASS…
    %   This class defines a Length measurement object. A Length must
    %   contain a positive value.
    properties(Constant)
        DEBUG = 0;
    end
    
    properties (Dependent)
        nanometer = double(0);
        micron = double(0);
        micrometer = double(0);
        millimeter = double(0);
        centimeter = double(0);
        meter = double(0);
        m = double(0);
        kilometer = double(0);
        
        mile = double(0);
        yard = double(0);
        foot = double(0);
        inch = double(0);
        milliinch = double(0);
        
        lightyear = double(0);
        nauticalmile = double(0);
    end
    
    properties(Access=private)
        basemetre = double(0)
    end
   
    methods
        function obj = distance(parm1)
            if nargin > 0
                obj.basemetre = parm1;
                if obj.DEBUG
                    obj.show_object;
                end
            end
        end
        
        %<--- Metric ------------
        
        function obj = set.nanometer(obj,value)
            obj.basemetre = value / 1000000000;
            obj.Validate
        end
        function value = get.nanometer(obj)
            value = obj.basemetre * 1000000000;
        end
        
        % µ
        function obj = set.micron(obj,value)
            obj.basemetre = value / 1000000;
            obj.Validate
        end
        function value = get.micron(obj)
            value = obj.basemetre * 1000000;
        end
        
        % µm
        function obj = set.micrometer(obj,value)
            obj.basemetre = value / 1000000;
            obj.Validate
        end
        function value = get.micrometer(obj)
            value = obj.basemetre * 1000000;
        end
        
        % mm
        function obj = set.millimeter(obj,value)
            obj.basemetre = value / 1000;
            obj.Validate
        end
        function value = get.millimeter(obj)
            value = obj.basemetre * 1000;
        end
        
        % cm
        function obj = set.centimeter(obj,value)
            obj.basemetre = value / 100;
            obj.Validate
        end
        function value = get.centimeter(obj)
            value = obj.basemetre * 100;
        end
        
        % m
        function obj = set.meter(obj,value)
            obj.basemetre = value;
            obj.Validate
        end
        function value = get.meter(obj)
            value = obj.basemetre;
        end
        function obj = set.m(obj,value)
            obj.basemetre = value;
            obj.Validate
        end
        function value = get.m(obj)
            value = obj.basemetre;
        end
        
        % km
        function obj = set.kilometer(obj,value)
            obj.basemetre = value / 0.001;
            obj.Validate
        end
        function value = get.kilometer(obj)
            value = obj.basemetre * 0.001;
        end
        
        %<--- US (Imperial) Measure ------------
        %mi
        function obj = set.mile(obj,value)
            obj.basemetre = value / 0.00062137119;
            obj.Validate
        end
        function value = get.mile(obj)
            value = obj.basemetre * 0.00062137119;
        end
        
        % yd
        function obj = set.yard(obj,value)
            obj.basemetre = value / 1.093613298337708;
            obj.Validate
        end
        function value = get.yard(obj)
            value = obj.basemetre * 1.093613298337708;
        end
        
        % ft
        function obj = set.foot(obj,value)
            obj.basemetre = value / 3.280839895013123;
            obj.Validate
        end
        function value = get.foot(obj)
            value = obj.basemetre * 3.280839895013123;
        end
        
        % in
        function obj = set.inch(obj,value)
            obj.basemetre = value / 39.37007874015748;
            obj.Validate
        end
        function value = get.inch(obj)
            value = obj.basemetre * 39.37007874015748;
        end
        
        % mil, thou
        function obj = set.milliinch(obj,value)
            obj.basemetre = value / 39370.07874015748;
            obj.Validate
        end
        function value = get.milliinch(obj)
            value = obj.basemetre * 39370.07874015748;
        end
        
        %<--- Other Measure ------------
        % ly
        function obj = set.lightyear(obj,value)
            obj.basemetre = value / 1.057000834024615e-16;
            obj.Validate
        end
        function value = get.lightyear(obj)
            value = obj.basemetre * 1.057000834024615e-16;
        end
        
        % NM, nmi
        function obj = set.nauticalmile(obj,value)
            obj.basemetre = value / 5.399568034557235e-4;
            obj.Validate
        end
        function value = get.nauticalmile(obj)
            value = obj.basemetre * 5.399568034557235e-4;
        end
        
        function show_object(obj)
            fprintf('distance:basemetre : %0.0f\n',obj.basemetre);
        end
        
        function Validate(obj)
            
        end
    end
end

