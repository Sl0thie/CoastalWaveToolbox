classdef coastal_object < handle 
    %COASTAL_OBJECT class to manage coastal functions.
    %   This class provides a simple wrapper for the coast functions used
    %   within the formula sheet used for 6110ENG - Coastal Engineering and
    %   modeling.
    
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
        
        
        
    end
    
    methods
        function obj = coastal_object
            % Constructor for the object.
            obj.reset_data;
        end
        
        
        %% A.4 Surf Zone
        function [MUmax] = longshore_current(obj)
            
           obj.MUmax = some_const * sqrt(obj.g * obj.Hb) * tan(obj.BETA) * sin(2 * obj.THETAb);
           
           MUmax = obj.MUmax;
        end
        
        function [XIo] = surf_similarity_parameter(obj)
            
            obj.XIo = tan(obj.BETA/ sqrt(obj.Ho/obj.Lo));
            if obj.XIo >= 4
                disp('Surfzone is considered reflective');
            else
                disp('Surfzone is considered dissipative');
            end
            XIo = obj.XIo;
        end
        
        % A.3.3 Reflection coefficient.
        function [Cr] = find_Cr(obj)          
            obj.Cr = obj.Hr / obj.Hi;          
            Cr = obj.Cr;
        end
        
        function [Kr] = find_Kr(obj)           
            if obj.H > obj.H2
               obj.Kr = sqrt(obj.cg1/obj.cg2);           
                Kr = obj.Kr;
            else
                warning('H1 must be larger than H2');
            end           
        end
        
        
        
        
        %A.3.1 Shoaling and Refraction.
        function [H2] = find_H2(obj)      
            obj.H2 = obj.H * obj.Kr * obj.Ks;          
        end
        
        % A.2.3 Wave Hindcasting.
        function [Hmo,Tp] = JONSWAP(obj)       
            % Requires F,t,U10
            % upper limits H*mo 0.243, T*p 8.13
            
            obj.tstar = (obj.g * obj.td)/(obj.U10);
            
            
            obj.Fstar = (obj.g * obj.F)/(obj.U10^2);  
            
            
            obj.Fstar_eff = (obj.tstar/68.8)^(3/2);
            fprintf('Data = %f\n' , obj.Fstar_eff);
            
            F_eff = 0;
            if obj.Fstar > obj.Fstar_eff
                F_eff = obj.Fstar_eff;
            else
                F_eff = obj.Fstar;
            end
            
            obj.Hstar_mo = 0.0016 * ((F_eff)^(1/2));
            
            obj.Hmo = (obj.Hstar_mo * (obj.U10^2))/(obj.g);
            
            obj.Tpstar = 0.286 * (F_eff^(1/3));
            
            obj.Tp = (obj.Tpstar * obj.U10) / obj.g;
            
            
            Tp = obj.Tp
            Hmo = obj.Hmo;
            
        end
        
        function [XT] = generalised_pareto_distribution(obj)
            
            obj.ZETAu = obj.NeHu/obj.N;
            
            if obj.XI == 0
                obj.XT = obj.u + (obj.SIGMA * log(obj.m * obj.T * obj.ZETAu));
                
            else
                obj.XT = obj.u + ((obj.SIGMA / obj.XI) * (((obj.m * obj.T * obj.ZETAu)^obj.XI) - 1));             
            end
            
            XT = obj.XT;
        end
               
        function [XT] = generalised_value_distribution(obj)
            
            obj.yT = -log(1 - (1/obj.T))
                      
            if obj.XI == 0            
                obj.XT = obj.Mu - (obj.SIGMA * log(obj.yT));          
            else               
                obj.XT = obj.MU - ((obj.SIGMA / obj.XI) * (1 - obj.yT^-obj.XI));               
            end
            
            XT = obj.XT;
        end
        
        function [Hper] = find_Hper(obj)
            
            obj.Hper = obj.Hrms * ((-log(obj.xtv/100))^(1/2));
            Hper = obj.Hper;
            
        end
        
        function [Hbar] = find_Hbar(obj)
           
            obj.Hbar = 0.886 * obj.Hrms;
            Hbar = obj.Hbar;
            
        end
        
        function [Hsig] = find_Hsig(obj)
           
            obj.Hsig = (2)^(1/2) * obj.Hrms;
            Hsig = obj.Hsig;
                   
        end
        
        function [Hmax] = find_Hmax(obj)
            
           obj.Hmax = obj.Hrms * (log(obj.N));         
           Hmax = obj.Hmax;
            
        end
        
        function [Hn] = wave_height_exceedance(obj)
            
            obj.Hn = obj.Hrms * (-log(obj.n/obj.N))^(1/2);
            
            Hn = obj.Hn;
        end
        
        function [PH] = probability_of_a_wave_height(obj)
            
           obj.PH = exp(-(obj.xtv/obj.Hrms)^2);
            
           PH = obj.PH;
        end
        
        function [H] = find_wave_height_from_Ap(obj,Ap)
            
            H = (2 * Ap * sinh(obj.k * obj.h))/(cosh(obj.k *(obj.z + obj.h)))
            
        end
        
        function [u,w] = particle_velocites(obj)
            
            obj.u = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.w = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            
            u = obj.u;
            w = obj.w;
        end
        
        function [u_hat,w_hat] = velocity_amplitudes(obj)
            
            obj.L = ((9.81 * obj.T^2)/(2 * pi)) * (tanh( (((2 * pi)/obj.T) * sqrt(obj.h/9.81))^(3/2) ))^(2/3);
            obj.k = (pi * 2) / obj.L;
            obj.OMEGA = (pi * 2) / obj.T;
            
            u_hat = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            w_hat = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            
        end
        
        function [ETA] = find_wave_ETA(obj)
            % Requires h,z,T,pmax
            
            obj.L = ((9.81 * obj.T^2)/(2 * pi)) * (tanh( (((2 * pi)/obj.T) * sqrt(obj.h/9.81))^(3/2) ))^(2/3);      
            obj.k = (pi * 2) / obj.L;     
            obj.Kp = (cosh(obj.k * (obj.z + obj.h)))/(cosh(obj.k * obj.h));
            
            obj.ETA =  (obj.pmax + (obj.RHO * obj.g * obj.z))/(obj.RHO * obj.g * obj.Kp);
            
            
            ETA = obj.ETA;
            
            
        end
        
        function [z] = find_depth_from_hydrostatic_pressure(obj)
            
            obj.z = obj.ps / -obj.RHO / obj.g;
            
            
            z = obj.z;
        end
        
        function [pd] = dynamic_pressure(obj)
            
            obj.pd = obj.RHO * obj.g * obj.ETA * obj.Kp;
                       
            pd = obj.pd;
            
        end
        
        function [ps] = hydrostatic_pressure(obj)
            
           obj.ps = -obj.RHO * obj.g * obj.z;
           
           ps = obj.ps; 
        end        
        
        function [Hn] = height_probibility(obj,n)
           % n is the percentage. eg, 5% = 5.
           
           obj.Hn = obj.Hrms * sqrt(-ln(n/obj.N)); 
           
        end
                
        function [Hs] = estimate_significant_wave_height(obj,LineB)
           % LineA is from the left (n).
           % LineB is from the bottom (Havg).
           
           obj.Hs = LineB;
           
           obj.Hrms = obj.Havg / LineB;
           
           
           obj.Hs = obj.Hs * obj.Hrms;
           
           
           
           
            
        end
        
        function [E] = energy_density(obj)
            
            obj.E = (obj.RHO * obj.g * obj.H^2)/8;
            
            
            E = obj.E;
        end
        
        function [Kp] = pressure_coefficient(obj)
            % Requires k,z,h
            
            
            
            
            obj.Kp = (cosh(obj.k * (obj.z + obj.h)))/(cosh(obj.k * obj.h));
            
            
            Kp = obj.Kp;
        end     
        
        function [cg] = group_velosity(obj)
            
            obj.cg = (obj.c/2) * (1 + ((2 * obj.k * obj.h)/(sinh(2 * obj.k * obj.h))));
            
            cg = obj.cg;
        end
        
        function [c] = wave_celerity(obj)
            % Requires L (wave length) and T (wave period).

            if(isnan(obj.L))
                warning('L is required');
                return;
            end
            
            if(isnan(obj.T))
                warning('T is required');
                return;
            end
            
            obj.c = obj.L/obj.T;
            
            c = obj.c;
        end
        
        function [k] = wave_number(obj)
            
            if(isnan(obj.L))
                warning('L is required');
                return;
            end
            
            obj.k = (pi * 2) / obj.L;
            
            k = obj.k;
        end
        
        function [f] = average_frequency(obj)
            
            if(isnan(obj.T))
                warning('T is required');
                return;
            end
                      
            obj.f = 1 / obj.T;
            
            f = obj.f;
        end
        
        function [OMEGA] = angular_frequency(obj)
            
            if(isnan(obj.T))
                warning('T is required');
                return;
            end
            
            obj.OMEGA = (pi * 2) / obj.T;
            
            OMEGA = obj.OMEGA;
        end
        
        function [ALPHA,BETA] = particle_displacement_amplitudes(obj)
            
            obj.ALPHA =  (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.BETA =  (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            
            ALPHA = obj.ALPHA;
            BETA = obj.BETA;
        end
        
        function [ETA] = sinusoidal_wave_shape(obj)
            % Requires L,T
            obj.k = Wave_Number(obj.L);
            obj.OMEGA = Angular_Frequency(obj.T);
            obj.ETA = (obj.H/2) * cos((obj.k * obj.x) - (obj.OMEGA * obj.t));           
            ETA = obj.ETA;
        end
        
        function [PSI] = sine_wave_velocity_potential(obj)
            % Requires L,T          
            obj.k = Wave_Number(obj.L);
            obj.OMEGA = Angular_Frequency(obj.T);
            obj.PSI = (-9.81/obj.OMEGA) * (obj.H/2) * (cosh(obj.k * (obj.z + obj.h))/cosh(obj.k * obj.h)) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));           
            PSI = obj.PSI;
        end
        
        function [L] = wave_length(obj)
            % Requires T,h      
            obj.L = ((obj.g * obj.T^2)/(2 * pi)) * (tanh( (((2 * pi)/obj.T) * sqrt(obj.h/obj.g))^(3/2) ))^(2/3);           
            L = obj.L;
        end
        
        function reset_data(obj)
            obj.H = NaN;
            obj.Hs = NaN;
            obj.Hsig = NaN;
            obj.Havg = NaN;
            obj.Hrms = NaN;
            obj.Hmax = NaN;
            obj.Hn = NaN;
            obj.Hbar = NaN;
            obj.Hper = NaN;
            obj.OMEGA = NaN;
            obj.MWS = NaN;
            obj.L = NaN;
            obj.ETA = NaN;
            obj.h = NaN;
            obj.z = NaN;
            obj.T = NaN;
            obj.Tz = NaN;
            obj.N = NaN;
            obj.c = NaN;
            obj.x = NaN;            
            obj.f = NaN;
            obj.k = NaN;
            obj.Kp = NaN;
            obj.Ef = NaN;
            obj.E = NaN;            
            obj.F = NaN;
            obj.td = NaN;
            obj.u = NaN;
            obj.SIGMA = NaN;
            obj.XI = NaN;
            obj.ps = NaN;
            obj.pmax = NaN;
            obj.pmin = NaN;            
            obj.U10 = NaN;
            obj.n = NaN;
            obj.cg = NaN;
            obj.xtv = NaN;
            obj.PH = NaN;
            obj.pd = NaN;
            obj.yT = NaN;
            obj.XT = NaN;
            obj.w = NaN;
            obj.NeHu = NaN;
            obj.record_length = NaN;
            obj.m = NaN;
            obj.Hmo = NaN;
            obj.Hstar_mo = NaN;
            obj.Ug = NaN;
            obj.ff = NaN;
            obj.Fstar = NaN;
            obj.tstar = NaN;
            obj.Fstar_eff = NaN;
            obj.Tstarp = NaN;
            obj.Tpstar = NaN;
            obj.Tp = NaN;
            obj.t = NaN;
            obj.td = NaN;
            obj.ZETAu = NaN;
            obj.ALPHA = NaN;
            obj.BETA = NaN;
        end  
    end    
end

