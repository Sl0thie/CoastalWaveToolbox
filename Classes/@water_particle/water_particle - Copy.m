classdef water_particle < handle
    %WATER_PARTICLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant,Access = public)
        % Constants.
        g = double(9.81);
        RHO = double(1026);
        RHOair = double(1.2);
        MU = double(0.001);
        NU = double(0.000001);
        
    end
    
    properties (Access = public)
        % Variables.
        Cr = double(0);
        c = double(0);      % wave celerity (speed at which the wave form moves) = L/T
        cg = double(0);     % Wave Group Velosity.        
        Ef = double(0);     % Energy Flux.
        E = double(0);      % Energy Density.        
        f = double(0);      % average frequency.       
        F = double(0);      % Fetch length.
        ff = double(0);     % Coriolis parameter.
        Fstar = double(0);  %        
        Fstar_eff = double(0); %             
        H = double(0);      % Wave Height (Vertical distance between creats and trough)
        Hs = double(0);     %  
        Ho = double(0);     % ? Originasl H height.
        Hb = double(0);     % the breaker height.
        Hr = double(0);     %
        Hi = double(0);     %
        H2 = double(0);     %
        Hsb = double(0);    %
        hb = double(0);     %
        Hsig = double(0)    % significant wave height. (Average height of the highest 1/3 of waves)
        Havg = double(0);   % Average wave height.
        Hrms = double(0);   % Root mean square wave height = sqrt(Havg^2)
        Hmax = double(0);   % Maximum wave height.
        Hn = double(0);     % number exceeding.
        Hbar = double(0);   % 
        Hper = double(0);   %
        h = double(0);      % water depth (vertical distance between seabed and the MWS)
        Hmo = double(0);    % Significant wave height derived from the zeroth moment.
        Hstar_mo = double(0); %                
        k = double(0);      % Wave Number.
        Kp = double(0);     %
        Kr = double(0);     % is and increasing function of the surf similarity parameter.
        Ks = double(0);     %
        L = double(0);      % wave length (horizontal distance between two consecutive crests or troughs)
        Lo = double(0);     % Original wavge length.      
        m = double(0);      % number of events per year.
        MWS = double(0);    % Mean Water Surface is the time averaged water surface level      
        N = double(0);      % Number of waves.
        NeHu = double(0)    % Number of events exceeding threshold.        
        n = double(0);      % Count of something.       
        PH = double(0);     % probability of a wave height.      
        ps = double(0);     % Hydrostatic Pressure (ps = -RHO * g * z)
        pd = double(0);     % Dynamic Pressure (pd = RHO * g * ETA * Kp)
        pmax = double(0);   % Max pessure (Crest)
        pmin = double(0);   % Min pressure (Trough)      
        record_length = double(0);  % timespan of data collection.        
        T = double(0);      % wave period (time taken for two wave crests to pass a fixed point in space)
        Tstarp = double(0); %
        Tpstar = double(0); %      
        Tp = double(0);     % The period corosponding to the peak of the spectrum.
        t = double(0);      % time. in seconds.      
        td = double(0);     % Duration the wind blows for.
        tstar = double(0);  %       
        Tz = double(0);     % average zero crossing period, length of record/number of waves.
        U10 = double(0);      % Wind speed 10m above surface.       
        Ug = double(0);     %
        u = double(0);      % Threshold. %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       
        w = double(0);      % Offset location? <<<<<<<<<<<<<<<<<<<<<<<<<<<      
        x = double(0);      % horizontal coordinate positive in the direction of wave propagation
        xtv = double(0);    % x threshold value (Percentage?).
        XT = double(0);     %
        yT = double(0);     % 
        z = double(0);      % vertical coordinate positive upwards from the mean water surface (ie. z = 0 at the MWS and z = -h at the seabed)       
        
        ALPHA = double(0);  % Offset x from original point due to 
        BETA = double(0);   % Offset Z from original point due to      
        ETA = double(0);    % “eta” is the instantaneous water surface elevation relative to the MWS        
        OMEGA = double(0);  % angular frequency.       
        NUmax = double(0);  %       
        ZETAu = double(0)   % Is the probability of event exceeding u.
        SIGMA = double(0);  % Scale parameter.
        XI = double(0);     % Shape parameter.
        XIo = double(0);    % Surf Similarity Parameter.
        THETAb = double(0); %
        MUmax = double(0);  %
        
        % owned properties.
        position_x = double(0); %
        position_z = double(0); %        
        position_x_original = double(0); %
        position_z_original = double(0); %        
        sea_bed_z = double(0);  %
        u_hat = double(0);      %
        w_hat = double(0);      %
        z_axis_lower = double(0); %
        time_offset = double(0); %
    end
    
    methods
        function obj = water_particle(x,z,H,T,h,z_axis_lower,time_offset)
            % Constructor for the object.
            if nargin == 0 % Object is created empty.
                % Set object properties.
                obj.position_x = 0;
                obj.position_z = double(0);
                obj.position_x_original = 0; %
                obj.position_z_original = double(0); %
                obj.z = double(0);
                obj.H = 0;
                obj.T = 0;
                obj.h = 1;
                obj.ETA = 1/2;
                obj.z_axis_lower = 1;
                obj.time_offset = 1;              
            else % Object is with values.
                obj.position_x = x;
                obj.position_z = double(z);
                obj.position_x_original = x; %
                obj.position_z_original = double(z); %
                obj.z = double(z);
                obj.H = H;
                obj.T = T;
                obj.h = h;
                obj.ETA = h/2; %%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                obj.z_axis_lower = z_axis_lower;
                obj.time_offset = time_offset;
                % Calculate dependant properties.
                obj.wave_length;
                obj.wave_number;
                obj.angular_frequency;
                obj.wave_celerity;
                obj.velocity_amplitudes;
            end
            
        end
        
        function [do_plot,return_position_x,return_position_z] = plot_position(obj,t,time_per_meter)
            
            if obj.position_z_original > -obj.h
                do_plot = 1;
                obj.t = t - (time_per_meter * obj.time_offset);
                obj.sinusoidal_wave_shape;
                
                obj.particle_velocites;
                return_position_x = double(obj.position_x_original) + (obj.u);
                return_position_z = double(obj.position_z_original) + (obj.w);
                
%                 obj.velocity_amplitudes;
%                 return_position_x = double(obj.position_x_original) + (obj.u_hat);
%                 return_position_z = double(obj.position_z_original) + (obj.w_hat);
                
%                 obj.particle_displacement_amplitudes;
%                 return_position_x = double(obj.position_x_original) + (obj.ALPHA);
%                 return_position_z = double(obj.position_z_original) + (obj.BETA);
            else
                do_plot = 0;
                return_position_x = obj.position_x_original;
                return_position_z = obj.position_z_original;
            end
            
        end      
        
        function [L] = wave_length(obj)         
            obj.L = ((obj.g * obj.T^2)/(2 * pi)) * (tanh( (((2 * pi)/obj.T) * sqrt(obj.h/obj.g))^(3/2) ))^(2/3);
            L = obj.L;     
        end
        
        function [k] = wave_number(obj)        
            obj.k = (pi * 2) / obj.L;
            k = obj.k;
        end
        
        function [OMEGA] = angular_frequency(obj)           
            obj.OMEGA = (pi * 2) / obj.T;
            OMEGA = obj.OMEGA;
        end
        
        function [c] = wave_celerity(obj)      
            obj.c = obj.L/obj.T;
            c = obj.c;
        end
        
        function [u_hat,w_hat] = velocity_amplitudes(obj)           
            obj.u_hat = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            obj.w_hat = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            u_hat = obj.u_hat;
            w_hat = obj.w_hat;
        end
        
        function [u,w] = particle_velocites(obj)         
            obj.u = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.w = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));        
            u = obj.u;
            w = obj.w;
        end
        
        function [ALPHA,BETA] = particle_displacement_amplitudes(obj)           
            obj.ALPHA =  (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.BETA =  (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));          
            ALPHA = obj.ALPHA;
            BETA = obj.BETA;
        end
        
        function [ETA] = sinusoidal_wave_shape(obj)           
            obj.ETA = (obj.H/2) * cos((obj.k * obj.x) - (obj.OMEGA * obj.t));        
            ETA = obj.ETA;
        end
        
    end
    
end

