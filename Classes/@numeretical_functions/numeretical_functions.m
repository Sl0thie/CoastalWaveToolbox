classdef numeretical_functions < handle
    %numeretical_functions class to manage coastal numeretical functions.
    %This class provides a simple wrapper for the coast functions used
    %within the formula sheet used for 6110ENG - Coastal Engineering and
    %modeling.
    %
    % numeretical_functions Properties:
    %   a - Owen formula slope parameter a.
    %   b - Owen formula slope parameter b.
    %
    % numeretical_functions Methods:
    %   find_SLR - Finds the SLR (shoreline retreat)
    %   find_Q_grad - Finds the Q gradient
    %
    % See also POSITION TRACER
    
    properties (Access = public, Hidden = false)
        a = double(NaN);      % Owen formula slope parameter
        b = double(NaN);      % Owen formula slope parameter
        B = double(NaN);      % Bearm height.
        C_r = double(NaN);
        C_10 = double(NaN)    % Drag co-efficent.
        c = double(NaN);      % wave celerity (speed at which the wave form moves) = L/T
        c_g = double(NaN);     % Wave Group Velosity.
        D_n50 = double(NaN);  % Length of the side of a hypothetical cude of the same weight.
        E_f = double(NaN);     % Energy Flux.
        E = double(NaN);      % Energy Density.
        f = double(NaN);      % average frequency.
        F = double(NaN);      % Fetch length.
        ff = double(NaN);     % Coriolis parameter.
        Fstar = double(NaN);  %
        Fstar_eff = double(NaN);%
        GSLR = double(NaN);   % <<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>>>>>
        H = double(NaN);      % Wave Height (Vertical distance between creats and trough)
        H_s = double(NaN);     % Design significant wave height.
        H_stoe = double(NaN); % Design significant wave height at toe.
        H_o = double(NaN);     % Original H height.
        H_b = double(NaN);     % the breaker height.
        H_d = double(NaN);     % Diffracted wave height.
        H_r = double(NaN);     %
        H_i = double(NaN);     % Incident wave height.
        H2 = double(NaN);     %
        H_sb = double(NaN);    % <<<<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>
        h_b = double(NaN);     %
        H_BW = double(NaN);   % Height at the break water.
        H_sig = double(NaN);   % significant wave height. (Average height of the highest 1/3 of waves)
        H_avg = double(NaN);   % Average wave height.
        H_rms = double(NaN);   % Root mean square wave height = sqrt(H_avg^2)
        H_max = double(NaN);   % Maximum wave height.
        h_cl = double(NaN);   % depth at offshore point for bruun rule.
        H_n = double(NaN);     % number exceeding.
        H_bar = double(NaN);   %
        H_per = double(NaN);   %
        h = double(NaN);      % water depth (vertical distance between seabed and the MWS)
        h_toe = double(NaN);  % Design depth at structure toe at high water.
        H_mo = double(NaN);    % Significant wave height derived from the zeroth moment.
        Hstar_mo = double(NaN); %
        K_D = double(NaN);    % Is an empirical stability coefficient.
        k = double(NaN);      % Wave Number.
        k_o = double(NaN);     % Orignial Wave Number. <<<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>>>>>
        K_p = double(NaN);     %
        K_prime = double(NaN); % Diffraction coefficient.
        K_r = double(NaN);     % is and increasing function of the surf similarity parameter.
        K_s = double(NaN);     %
        L = double(NaN);      % wave length (horizontal distance between two consecutive crests or troughs) <<<<<<<< Double >>>>>>>> overloaded in last lecture
        L_basin = double(NaN);% Length of closed basin.
        L_bay = double(NaN);  % Length of open bay.
        Lo = double(NaN);     % Original wavge length.
        L_om = double(NaN);   %
        L_r = double(NaN);     % Vertical Length Scale.
        m = double(NaN);      % number of events per year.
        MWS = double(NaN);    % Mean Water Surface is the time averaged water surface level
        N = double(NaN);      % Number of waves.
        N_s = double(NaN);    % Sability ratio.
        NeHu = double(NaN);   % Number of events exceeding threshold.
        n = double(NaN);      % Count of something.
        P = double(NaN);      % Probability. ***************
        p = double(NaN);      % sediment porosity
        PH = double(NaN);     % probability of a wave height.
        ps = double(NaN);     % Hydrostatic Pressure (ps = -RHO * g * z)
        pd = double(NaN);     % Dynamic Pressure (pd = RHO * g * ETA * K_p)
        p_max = double(NaN);   % Max pessure (Crest)
        p_min = double(NaN);   % Min pressure (Trough)
        Qstar = double(NaN);  % Non-diamensional discharge
        Q_in = double(NaN);   % Sediment inflow.
        Q_out = double(NaN);  % Sediment outflow.
        Q_sinks = double(NaN);% Are any losses of sediment such as dredging.
        Q_grad = double(NaN); % is the loneshore gradient in alongshore transport rate. partial Q_y/ partial y.
        Q_x = double(NaN);    % is the shore-normal sediment transport rate per unit length of coastline.
        Q_y = double(NaN);    % is the alone shore sediment transport rate.
        Q_sources = double(NaN);% are additional sources of sediment such as beach nourishment or riverine discharge.
        q = double(NaN);      % Maxium allowable average overtopping discharge;
        R = double(NaN);      % Wave Runup.
        R_c = double(NaN);    % Creast elevation.
        R_star = double(NaN); %
        R_n = double(NaN);     %
        r = double(NaN);      % Distance from corner of breakwater (barrier) and point of interest.
        record_length = double(NaN);  % timespan of data collection.
        s = double(NaN)       % is the specific weight RHO_s/RHO.
        S_om = double(NaN);   % wave steepness
        SWL = double(NaN);    % still water level.
        SLR = double(NaN);    % Shore line retreat.
        T = double(NaN);      % wave period (time taken for two wave crests to pass a fixed point in space)
        Tstar_p = double(NaN); %
        T_p = double(NaN);     % The period corosponding to the peak of the spectrum.
        T_om = double(NaN);   %
        T_m = double(NaN);    % Mean period 0.76 * T_p.
        t = double(NaN);      % time. in seconds.
        t_d = double(NaN);     % Duration the wind blows for.
        tstar = double(NaN);  %
        T_z = double(NaN);     % average zero crossing period, length of record/number of waves. <<<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>
        T_n = double(NaN);    %
        U_10 = double(NaN);      % Wind speed 10m above surface.
        U_g = double(NaN);     % <<<<<<<<<<<<<<<<<<<< unused >>>>>>>>>>>>>>>>>>>>>>>>>>>
        u = double(NaN);      % Threshold. %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        u_hat = double(NaN);  %
        W = double(NaN)       % Width of the shelf between water depths.
        W_n50 = double(NaN)   % Median Weight.
        w = double(NaN);      % Offset location? <<<<<<<<<<<<<<<<<<<<<<<<<<<
        w_hat = double(NaN);  %
        x = double(NaN);      % horizontal coordinate positive in the direction of wave propagation.
        y = double(NaN);      %
        xtv = double(NaN);    % x threshold value (Percentage?).
        X_T = double(NaN);     %
        y_T = double(NaN);     %
        z = double(NaN);      % vertical coordinate positive upwards from the mean water surface (ie. z = 0 at the MWS and z = -h at the seabed)
        z_100 = double(NaN);   % highest point transgressed by 100% of the waves.
        zwm = double(NaN);    % max runup elevation <<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>>>>>>>>
        
        
        ALPHA = double(NaN);  % is the structure slope and (ie tanALPHA is the slope)
        ALPHA_disp = double(NaN);  % Offset x from original point due to
        ALPHA_NU = double(NaN);   % 1 <= ALPHA_NU <= 1.5
        BETA = double(NaN);   % Angle between line of breakwater (barrier) and wave ray.
        BETA_disp = double(NaN);   % Offset Z from original point due to
        tanBETA = double(NaN);% Bottom Slope. (Value only not function tan)
        DELTA = double(NaN);  % (RHO_s/RHO_W) - 1 is the relative buoyant density of the armour.
        DELTA_x = double(NaN);% messure of shoreline retreat.
        DELTA_z = double(NaN);%
        ETA = double(NaN);    % “eta” is the instantaneous water surface elevation relative to the MWS
        ETAbar = double(NaN)  %
        ETAbar_wind = double(NaN)%
        GAMMA_b = double(NaN); % (?b) Breaker Index H_b/h_b.
        GAMMA_r = double(NaN);% surface roughness redution factor.
        OMEGA = double(NaN);  % angular frequency.
        MU_max = double(NaN);  %
        NU_max = double(NaN);  %  <<<<<<<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        PHI = double(NaN)     %
        PSI = double(NaN);    %
        RHO_s = double(NaN);  % Density of armour unit material.
        SIGMA = double(NaN);  % Scale parameter.
        TAO_w = double(NaN);   % Wind shear stress.
        TAO_b = double(NaN);   % Bottom shear stress.
        THETA = double(NaN);  % Angle that the wave approches the shore.
        THETA_b = double(NaN); % Wave breaker angle.
        THETA_o = double(NaN); %
        XI = double(NaN);     % Shape parameter.
        XI_o = double(NaN);    % Surf Similarity Parameter. (Iribarren Number)
        XI_b = double(NaN);    % Surf Similarity Parameter. for breaker type. <<<<<<<<<<<<<<<< Unused >>>>>>>>>>>>>>>>>>>>>>>>>>
        ZETA_u = double(NaN)   % Is the probability of event exceeding u.
        
        
        deep_water = double(NaN); % Flag for functions.
        wave_is_broken = double(NaN);%
    end
    
    properties (Constant,Access = public, Hidden = true)
        % Constants.
        g = double(9.81);
        RHO = double(1026);     % density of water.
        RHO_air = double(1.2);  % density of air.
        MU = double(0.001);
        NU = double(0.000001);
        RCP_2_5_SLR_min = double(0.28);%
        RCP_2_5_SLR_max = double(0.61);%
        RCP_8_5_SLR_min = double(0.53);%
        RCP_8_5_SLR_max = double(0.98);%
        K = double(0.77); % is an empirical coefficient witch has a weak dependance on grain size.
    end
    
    methods (Access = public, Hidden = true)
        %% Object Methods
        function obj = coastal.numeretical_functions(x,y,z,H,T,h)
            % Constructor for the object.            
            if nargin == 0
                obj.x = NaN;
                obj.z = NaN;
                obj.H = NaN;
                obj.T = NaN;
                obj.h = NaN;
                obj.y = NaN;
            else
                obj.x = x;
                obj.z = z;
                obj.H = H;
                obj.T = T;
                obj.h = h;
                obj.y = y;
                obj.basic_wave_properties;
            end
        end
        
        function basic_wave_properties(obj)
            obj.wave_length;
            obj.wave_number;
            obj.angular_frequency;
            obj.wave_celerity;
            obj.group_velosity;
            obj.velocity_amplitudes;
            obj.energy_density;
        end
        
        function reset_data(obj)
            obj.a = NaN;
            obj.b = NaN;
            obj.B = NaN;
            obj.C_r = NaN;
            obj.C_10 = NaN;
            obj.c = NaN;
            obj.c_g = NaN;
            obj.D_n50 = NaN;
            obj.E_f = NaN;
            obj.E = NaN;
            obj.f = NaN;
            obj.F = NaN;
            obj.ff = NaN;
            obj.Fstar = NaN;
            obj.Fstar_eff = NaN;
            obj.GSLR = NaN;
            obj.H = NaN;
            obj.H_s = NaN;
            obj.H_stoe = NaN;
            obj.H_o = NaN;
            obj.H_b = NaN;
            obj.H_d = NaN;
            obj.H_r = NaN;
            obj.H_i = NaN;
            obj.H2 = NaN;
            obj.H_sb = NaN;
            obj.h_b = NaN;
            obj.H_BW = NaN;
            obj.H_sig = NaN;
            obj.H_avg = NaN;
            obj.H_rms = NaN;
            obj.H_max = NaN;
            obj.h_cl = NaN;
            obj.H_n = NaN;
            obj.H_bar = NaN;
            obj.H_per = NaN;
            obj.h = NaN;
            obj.h_toe = NaN;
            obj.H_mo = NaN;
            obj.Hstar_mo = NaN;
            obj.K_D = NaN;
            obj.k = NaN;
            obj.k_o = NaN;
            obj.K_p = NaN;
            obj.K_prime = NaN;
            obj.K_r = NaN;
            obj.K_s = NaN;
            obj.L = NaN;
            obj.L_basin = NaN;
            obj.L_bay = NaN;
            obj.Lo = NaN;
            obj.L_om = NaN;
            obj.L_r = NaN;
            obj.m = NaN;
            obj.MWS = NaN;
            obj.N = NaN;
            obj.N_s = NaN;
            obj.NeHu = NaN;
            obj.n = NaN;
            obj.P = NaN;
            obj.p = NaN;
            obj.PH = NaN;
            obj.ps = NaN;
            obj.pd = NaN;
            obj.p_max = NaN;
            obj.p_min = NaN;
            obj.Qstar = NaN;
            obj.Q_in = NaN;
            obj.Q_out = NaN;
            obj.Q_sinks = NaN;
            obj.Q_grad = NaN;
            obj.Q_x = NaN;
            obj.Q_y = NaN;
            obj.Q_sources = NaN;
            obj.q = NaN;
            obj.R = NaN;
            obj.R_c = NaN;
            obj.R_star = NaN;
            obj.R_n = NaN;
            obj.r = NaN;
            obj.record_length = NaN;
            obj.s = NaN;
            obj.S_om = NaN;
            obj.SWL = NaN;
            obj.SLR = NaN;
            obj.T = NaN;
            obj.Tstar_p = NaN;
            obj.T_p = NaN;
            obj.T_om = NaN;
            obj.T_m = NaN;
            obj.t = NaN;
            obj.t_d = NaN;
            obj.tstar = NaN;
            obj.T_z = NaN;
            obj.T_n = NaN;
            obj.U_10 = NaN;
            obj.U_g = NaN;
            obj.u = NaN;
            obj.u_hat = NaN;
            obj.W = NaN;
            obj.W_n50 = NaN;
            obj.w = NaN;
            obj.w_hat = NaN;
            obj.x = NaN;
            obj.y = NaN;
            obj.xtv = NaN;
            obj.X_T = NaN;
            obj.y_T = NaN;
            obj.z = NaN;
            obj.z_100 = NaN;
            obj.zwm = NaN;
            
            
            obj.ALPHA = NaN;
            obj.ALPHA_disp = NaN;
            obj.ALPHA_NU = NaN;
            obj.BETA = NaN;
            obj.BETA_disp = NaN;
            obj.tanBETA = NaN;
            obj.DELTA = NaN;
            obj.DELTA_x = NaN;
            obj.DELTA_z = NaN;
            obj.ETA = NaN;
            obj.ETAbar = NaN;
            obj.ETAbar_wind = NaN;
            obj.GAMMA_b = NaN;
            obj.GAMMA_r = NaN;
            obj.OMEGA = NaN;
            obj.MU_max = NaN;
            obj.NU_max = NaN;
            obj.PHI = NaN;
            obj.PSI = NaN;
            obj.RHO_s = NaN;
            obj.SIGMA = NaN;
            obj.TAO_w = NaN;
            obj.TAO_b = NaN;
            obj.THETA = NaN;
            obj.THETA_b = NaN;
            obj.THETA_o = NaN;
            obj.XI = NaN;
            obj.XI_o = NaN;
            obj.XI_b = NaN;
            obj.ZETA_u = NaN;
            
            
            obj.deep_water = NaN;
            obj.wave_is_broken = NaN;
        end
        
%         function disp(obj)
%             
%             if ~isnan(obj.H)
%                 disp(['H wave height = ',num2str(obj.H)])               
%             end
%         end
    end
    
    methods (Access = public, Hidden = false)
        %%
        function [SLR] = find_SLR(obj)
        %find_SLR class to manage coastal functions.
        %   This class provides a simple wrapper for the coast functions used
        %   within the formula sheet used for 6110ENG - Coastal Engineering and
        %   modeling.
        %
        % See also POSITION TRACER
            if(isnan(obj.Q_grad))
                error('Q_grad is required');
                return;
            end
            if(isnan(obj.p))
                error('p is required');
                return;
            end
            if(isnan(obj.h_cl))
                error('h_cl is required');
                return;
            end
            if(isnan(obj.B))
                error('B is required');
                return;
            end
            if(isnan(obj.Q_x))
                error('Q_x is required');
                return;
            end
            grad_per_year = obj.Q_grad * 365 * 24 * 60 * 60;
            obj.SLR = (obj.Q_x + grad_per_year)/((1-obj.p)*(obj.h_cl + obj.B));
            SLR = obj.SLR;
            
        end
        
        function [Q_grad] = find_Q_grad(obj,refQ_y,refy)
        %find_Q_grad class to manage coastal functions.
        %   This class provides a simple wrapper for the coast functions used
        %   within the formula sheet used for 6110ENG - Coastal Engineering and
        %   modeling.
        %
        % See also POSITION TRACER
            if(isnan(obj.Q_y))
                error('Q_y is required');
                return;
            end
            if(isnan(obj.y))
                error('y is required');
                return;
            end
            obj.Q_grad = (obj.Q_y - refQ_y)/(obj.y - refy);
            Q_grad = obj.Q_grad;
            
        end
              
        function [Q_x] = find_Q_x(obj)
            if(isnan(obj.Q_grad))
                error('Q_grad is required');
                return;
            end
            if(isnan(obj.p))
                error('p is required');
                return;
            end
            if(isnan(obj.h_cl))
                error('h_cl is required');
                return;
            end
            if(isnan(obj.B))
                error('B is required');
                return;
            end
            if(isnan(obj.SLR))
                error('SLR is required');
                return;
            end
            grad_per_year = obj.Q_grad * 365 * 24 * 60 * 60;
            obj.Q_x = (grad_per_year)-(1-obj.p)*(obj.h_cl + obj.B)*(obj.SLR);
            Q_x = obj.Q_x;
            
        end
        
        function [Q_y] = find_Q_y(obj)
            if(isnan(obj.K))
                error('K is required');
                return;
            end
            if(isnan(obj.s))
                error('s is required');
                return;
            end
            if(isnan(obj.H_s))
                error('H_s is required');
                return;
            end
            if(isnan(obj.GAMMA_b))
                error('GAMMA_b is required');
                return;
            end
            if(isnan(obj.THETA_b))
                error('THETA_b is required');
                return;
            end
            obj.Q_y = -(obj.K/(16 * (obj.s-1) * sqrt(obj.GAMMA_b))) * (sqrt(obj.g) * (obj.H_b^(5/2))) * (sind(2 * obj.THETA_b));
            Q_y = obj.Q_y;
            
        end
        
        %%
        
        function[DELTA_z] = find_DELTA_z_high_emissions(obj)
            
            obj.DELTA_z = (obj.RCP_8_5_SLR_min + obj.RCP_8_5_SLR_max)/2;
            DELTA_z = obj.DELTA_z;
            
        end
        
        function[DELTA_z] = find_DELTA_z_low_emissions(obj)
            
            obj.DELTA_z = (obj.RCP_2_5_SLR_min + obj.RCP_2_5_SLR_max)/2;
            DELTA_z = obj.DELTA_z;
            
        end
        
        function [DELTA_x]  =bruuns_rule(obj)
            if(isnan(obj.DELTA_z))
                error('DELTA_z is required');
                return;
            end
            if(isnan(obj.W))
                error('W is required');
                return;
            end
            if(isnan(obj.B))
                error('B is required');
                return;
            end
            if(isnan(obj.h_cl))
                error('h_cl is required');
                return;
            end
            obj.DELTA_x =  obj.DELTA_z * (obj.W/(obj.B + obj.h_cl));
            DELTA_x = obj.DELTA_x;
            
        end
        
        %% A.7.1 Rubble Mound Structures
        
        function [R_c] = find_R_c(obj)
            if(isnan(obj.R_star))
                error('R_star is required');
                return;
            end
            if(isnan(obj.H_s))
                error('H_s is required');
                return;
            end
            if(isnan(obj.GAMMA_r))
                error('GAMMA_r is required');
                return;
            end
            if(isnan(obj.S_om))
                error('S_om is required');
                return;
            end
            obj.R_c = (obj.R_star * obj.H_s * obj.GAMMA_r)/sqrt(obj.S_om/(2 * pi)) ;
            R_c= obj.R_c;
            
        end
        
        function[L_om] = find_L_om(obj)
            if(isnan(obj.T_om))
                error('T_om is required');
                return;
            end
            obj.L_om = (obj.g * obj.T_om^2)/(2 * pi) ;
            L_om = obj.L_om;
            
        end
        
        function[S_om] = find_wave_steepness(obj,refL_om)
            if(isnan(obj.H_stoe))
                error('H_stoe is required');
                return;
            end
            obj.S_om = (obj.H_stoe/refL_om) ;
            S_om = obj.S_om;
            
        end
        
        function[Qstar] = owens_formula_Q_star(obj)
            if(isnan(obj.H_s))
                error('H_s is required');
                return;
            end
            if(isnan(obj.T_om))
                error('T_om is required');
                return;
            end
            obj.Qstar = obj.q/(obj.g * obj.H_s * obj.T_om);
            Qstar = obj.Qstar;
            
        end
        
        function[R_star] = owens_formula_R_star(obj)
            if(isnan(obj.a))
                error('a is required');
                return;
            end
            if(isnan(obj.b))
                error('b is required');
                return;
            end
            if(isnan(obj.Qstar))
                error('Qstar is required');
                return;
            end
            obj.R_star = (-1/obj.b) * log(obj.Qstar/obj.a);
            %obj.R_star = (obj.R_c/obj.H_s) * sqrt(obj.S_om/(2 * pi)) * (1/obj.GAMMA_y);
            R_star = obj.R_star;
            
        end
        
        function [H_b] = breaker_criteria(obj)
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            if(isnan(obj.GAMMA_b))
                error('GAMMA_b is required');
                return;
            end
            
            obj.H_b = obj.h * obj.GAMMA_b;
            if obj.H_BW <= obj.H_b
                obj.wave_is_broken = 0;
            else
                
                obj.wave_is_broken = 1;
                obj.H_b = obj.H_BW;
            end
            H_b = obj.H_b;
            
        end
        
        function [DELTA] = find_DELTA(obj)
            if(isnan(obj.RHO_s))
                error('RHO_s is required');
                return;
            end
            obj.DELTA = (obj.RHO_s/obj.RHO) - 1;
            DELTA = obj.DELTA;
            
        end
        
        function [W_n50] = weight_of_armour_unit(obj)
            if(isnan(obj.PHO_s))
                error('RHO_s is required');
                return;
            end
            if(isnan(obj.D_n50))
                error('D_n50 is required');
                return;
            end
            obj.W_n50 = obj.RHO_s * obj.g * obj.D_n50^3;
            W_n50 = obj.W_n50;
            
        end
        
        function [D_n50] = hudsons_formula_for_D_n50(obj)
            if(isnan(obj.DELTA))
                error('DELTA is required');
                return;
            end
            if(isnan(obj.K_D))
                error('K_D is required');
                return;
            end
            if(isnan(obj.ALPHA))
                error('ALPHA is required');
                return;
            end
            obj.D_n50 = (obj.H)/(obj.DELTA * (obj.K_D * obj.ALPHA)^(1/3));
            D_n50 = obj.D_n50;
            
        end
        
        %% A.6.4 Seiches
        function [T_n] = natural_period_oscillation_closed_basin(obj)
            if(isnan(obj.L_basin))
                error('L_basin is required');
                return;
            end
            if(isnan(obj.n))
                error('n is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            obj.T_n = (2 * obj.L_basin)/(obj.n * sqrt(obj.g * obj.h));
            T_n = obj.T_n;
            
        end
        
        function [T_n] = natural_period_oscillation_open_bay(obj)
            if(isnan(obj.L_basin))
                error('L_basin is required');
                return;
            end
            if(isnan(obj.n))
                error('n is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            obj.T_n = (4 * obj.L_basin)/((2 * obj.n + 1) * sqrt(obj.g * obj.h));
            T_n = obj.T_n;
            
        end
        
        %% A.6.3 Storm Surges
        function [ETAbar_wind] = wind_setup(obj,refh)
            if(isnan(obj.ALPHA_NU))
                error('ALPHA_NU is required');
                return;
            end
            if(isnan(obj.TAO_w))
                error('TAO_w is required');
                return;
            end
            if(isnan(obj.ALPHA_NU))
                error('ALPHA_NU is required');
                return;
            end
            if(isnan(refh))
                error('refh is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            obj.ETAbar_wind = obj.ALPHA_NU *(obj.TAO_w/(obj.RHO * obj.g)) * (obj.W/refh) * log(refh/obj.h);
            ETAbar_wind = obj.ETAbar_wind;
            
        end
        
        function [TAO_w] = find_TAO_w(obj)
            if(isnan(obj.C_10))
                error('C_10 is required');
                return;
            end
            if(isnan(obj.U_10))
                error('U_10 is required');
                return;
            end
            if(isnan(obj.PHI))
                error('PHI is required');
                return;
            end
            obj.TAO_w = 0.5 * obj.RHO_air * obj.C_10 * (obj.U_10 * cosd(obj.PHI))^2;
            TAO_w = obj.TAO_w;
            
        end
        
        %% A.5.3 Wave Runup
        function [P] = natural_waves_rayleigh_distribution_probability(obj)
            if(isnan(obj.R_n))
                error('R_n is required');
                return;
            end
            if(isnan(obj.L_r))
                error('L_r is required');
                return;
            end
            obj.P = exp(-(obj.R_n/obj.L_r)^2);
            P = obj.P;
        end
        
        function [R_n] = natural_waves_rayleigh_distribution_number(obj)
            if(isnan(obj.L_r))
                error('L_r is required');
                return;
            end
            if(isnan(obj.n))
                error('n is required');
                return;
            end
            if(isnan(obj.N))
                error('N is required');
                return;
            end
            if(isnan(obj.z_100))
                error('z_100 is required');
                return;
            end
            obj.R_n =  obj.L_r * sqrt(-log(obj.n/obj.N)) + obj.z_100;
            R_n = obj.R_n;
            
        end
        
        function [L_r] = find_vertical_length_scale(obj)
            if(isnan(obj.tanBETA))
                error('tanBETA is required');
                return;
            end
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            if(isnan(obj.Lo))
                error('Lo is required');
                return;
            end
            if obj.tanBETA > 0.1
                obj.L_r = 0.6 * obj.tanBETA * sqrt(obj.H_rms * obj.Lo);
                %obj.L_r = 0.06 * sqrt(obj.H_rms * obj.Lo);
            else
                %obj.L_r = 0.6 * obj.tanBETA * sqrt(obj.H_rms * obj.Lo);
                obj.L_r = 0.06 * sqrt(obj.H_rms * obj.Lo);
            end
            
            L_r = obj.L_r;
        end
        
        %% A.5 Surf Zone
        function [H_b] = find_H_b(obj)
            if(isnan(obj.H_o))
                error('H_o is required');
                return;
            end
            if(isnan(obj.THETA_o))
                error('THETA_o is required');
                return;
            end
            if(isnan(obj.GAMMA_b))
                error('GAMMA_b is required');
                return;
            end
            %obj.H_b = obj.H_o * sqrt(cosd(obj.THETA_o)) * (1/(nthroot((4 * obj.k * (obj.H_o / obj.GAMMA_b)),4)));
            obj.H_b = ((obj.H_o^4 * cosd(obj.THETA_o) * cosd(obj.THETA_o))/(4 * (obj.k / obj.GAMMA_b )))^(1/5);
            H_b = obj.H_b;
        end
        
        function [K_s] = nielsens_explicit_approximations_K_s(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            obj.K_s = 1/(nthroot((4 * obj.k * obj.h),4));
            K_s = obj.K_s;
            
        end
        
        function [K_r] = nielsens_explicit_approximations_K_r(obj)
            if(isnan(obj.THETA_o))
                error('THETA_o is required');
                return;
            end
            obj.K_r = sqrt(cosd(obj.THETA_o));
            K_r = obj.K_r;
            
        end
        
        function [MU_max] = longshore_current(obj)
            if(isnan(obj.H_b))
                error('H_b is required');
                return;
            end
            if(isnan(obj.BETA))
                error('BETA is required');
                return;
            end
            if(isnan(obj.THETA_b))
                error('THETA_b is required');
                return;
            end
            obj.MU_max = some_const * sqrt(obj.g * obj.H_b) * tan(obj.BETA) * sin(2 * obj.THETA_b);
            MU_max = obj.MU_max;
        end
        
        function [XI_o] = surf_similarity_parameter(obj)
            if(isnan(obj.tanBETA))
                error('tanBETA is required');
                return;
            end
            if(isnan(obj.H_o))
                error('H_o is required');
                return;
            end
            if(isnan(obj.Lo))
                error('Lo is required');
                return;
            end
            obj.XI_o = obj.tanBETA / sqrt(obj.H_o/obj.Lo);
            if obj.XI_o >= 4
                disp('Surfzone is considered reflective');
            else
                disp('Surfzone is considered dissipative');
            end
            XI_o = obj.XI_o;
        end
        
        %% A.4.3 Reflection coefficient.
        function [C_r] = refraction_coefficient(obj)
            if(isnan(obj.H_r))
                error('H_r is required');
                return;
            end
            if(isnan(obj.H_i))
                error('H_i is required');
                return;
            end
            obj.C_r = obj.H_r / obj.H_i;
            C_r = obj.C_r;
        end
        
        %% A.4.2 Diffraction coefficient.
        function [K_prime] = diffration_coefficient(obj)
            if(isnan(obj.H_d))
                error('H_d is required');
                return;
            end
            if(isnan(obj.H_i))
                error('H_i is required');
                return;
            end
            obj.K_prime = obj.H_d / obj.H_i;
            K_prime = obj.K_prime;
        end
        
        %% A.4.1 Shoaling and Refraction.
        function [H] = find_H(obj,ref_H)
            if(isnan(obj.K_r))
                error('K_r is required');
                return;
            end
            if(isnan(obj.K_s))
                error('K_s is required');
                return;
            end
            obj.H = ref_H * obj.K_r * obj.K_s;
            H = obj.H;
            % Now need to set new values for c,k.
            %obj.basic_wave_properties;
        end
        
        function [K_s] = shoaling_coefficent(obj,ref_c_g)
            if(isnan(obj.c_g))
                error('c_g is required');
                return;
            end
            obj.K_s = sqrt(ref_c_g/obj.c_g);
            K_s = obj.K_s;
        end
        
        function [K_r] = refraction_coefficent(obj,ref_THETA)
            if(isnan(obj.THETA))
                error('THETA is required');
                return;
            end
            obj.K_r = sqrt(cosd(ref_THETA)/cosd(obj.THETA));
            K_r = obj.K_r;
        end
        
        function [THETA] = snells_law(obj,ref_c,ref_THETA)
            if(isnan(obj.c))
                error('c is required');
                return;
            end
            obj.THETA = asind(sind(ref_THETA) * (obj.c/ref_c));
            THETA = obj.THETA;
        end
        
        function [THETA_1] = snells_law_for_THETA_1(obj,ref_c)
            if(isnan(obj.c))
                error('c is required');
                return;
            end
            if(isnan(obj.THETA))
                error('THETA is required');
                return;
            end
            %THETA_1 = 1/asind(obj.c/ref_c)/sind(obj.THETA);
            THETA_1 = asind(1/((obj.c/ref_c)/sind(obj.THETA)));
        end
        
        %% A.3.3 Wave Hindcasting.
        function [H_mo,T_p] = JONSWAP(obj)
            if(isnan(obj.t_d))
                error('t_d is required');
                return;
            end
            if(isnan(obj.F))
                error('F is required');
                return;
            end
            if(isnan(obj.U_10))
                error('U_10 is required');
                return;
            end
            
            % Requires F,t,U_10
            % upper limits H*mo 0.243, T*p 8.13
            obj.tstar = (obj.g * obj.t_d)/(obj.U_10);
            obj.Fstar = (obj.g * obj.F)/(obj.U_10^2);
            obj.Fstar_eff = (obj.tstar/68.8)^(3/2);
            fprintf('Data = %f\n' , obj.Fstar_eff);
            F_eff = 0;
            if obj.Fstar > obj.Fstar_eff
                F_eff = obj.Fstar_eff;
            else
                F_eff = obj.Fstar;
            end
            obj.Hstar_mo = 0.0016 * ((F_eff)^(1/2));
            obj.H_mo = (obj.Hstar_mo * (obj.U_10^2))/(obj.g);
            obj.Tstar_p = 0.286 * (F_eff^(1/3));
            obj.T_p = (obj.Tstar_p * obj.U_10) / obj.g;
            T_p = obj.T_p;
            H_mo = obj.H_mo;
        end
        
        %% A.3.2 Long-term record.
        function [X_T] = generalised_pareto_distribution(obj)
            if(isnan(obj.NeHu))
                error('NeHu is required');
                return;
            end
            if(isnan(obj.N))
                error('N is required');
                return;
            end
            if(isnan(obj.SIGMA))
                error('SIGMA is required');
                return;
            end
            if(isnan(obj.m))
                error('m is required');
                return;
            end
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            if(isnan(obj.u))
                error('u is required');
                return;
            end
            obj.ZETA_u = obj.NeHu/obj.N;
            if obj.XI == 0
                obj.X_T = obj.u + (obj.SIGMA * log(obj.m * obj.T * obj.ZETA_u));
            else
                obj.X_T = obj.u + ((obj.SIGMA / obj.XI) * (((obj.m * obj.T * obj.ZETA_u)^obj.XI) - 1));
            end
            X_T = obj.X_T;
        end
        
        function [X_T] = generalised_value_distribution(obj)
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            if(isnan(obj.SIMGA))
                error('SIGMA is required');
                return;
            end
            if(isnan(obj.XI))
                error('XI is required');
                return;
            end
            
            obj.y_T = -log(1 - (1/obj.T))
            if obj.XI == 0
                obj.X_T = obj.Mu - (obj.SIGMA * log(obj.y_T));
            else
                obj.X_T = obj.MU - ((obj.SIGMA / obj.XI) * (1 - obj.y_T^-obj.XI));
            end
            X_T = obj.X_T;
        end
        
        %% A.3.1 Short-term record - Rayleigh Distribution.
        function [H_rms] = find_H_rms(obj)
            if(isnan(obj.H_s))
                error('H_s is required');
                return;
            end
            obj.H_rms = obj.H_s / 1.42;
            H_rms = obj.H_rms;
            
        end
        
        function [H_per] = find_H_per(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            if(isnan(obj.xtv))
                error('xtv is required');
                return;
            end
            obj.H_per = obj.H_rms * ((-log(obj.xtv/100))^(1/2));
            H_per = obj.H_per;
        end
        
        function [H_bar] = find_H_bar(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            obj.H_bar = 0.886 * obj.H_rms;
            H_bar = obj.H_bar;
        end
        
        function [H_sig] = find_H_sig(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            obj.H_sig = (2)^(1/2) * obj.H_rms;
            H_sig = obj.H_sig;
        end
        
        function [H_max] = find_H_max(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            if(isnan(obj.N))
                error('N is required');
                return;
            end
            obj.H_max = obj.H_rms * (log(obj.N));
            H_max = obj.H_max;
        end
        
        function [H_n] = wave_height_exceedance(obj)
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            if(isnan(obj.n))
                error('n is required');
                return;
            end
            if(isnan(obj.N))
                error('N is required');
                return;
            end
            obj.H_n = obj.H_rms * (-log(obj.n/obj.N))^(1/2);
            H_n = obj.H_n;
        end
        
        function [PH] = probability_of_a_wave_height(obj)
            if(isnan(obj.xtv))
                error('xtv is required');
                return;
            end
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            obj.PH = exp(-(obj.xtv/obj.H_rms)^2);
            PH = obj.PH;
        end
        
        %% A.2 Small Amplittude (Linear) wave theory.
        function [L] = wave_length(obj)
            if obj.deep_water == 1
                if(isnan(obj.T))
                    error('T is required');
                    return;
                end
                obj.L = ((obj.g * obj.T^2)/(2 * pi));
                L = obj.L;
            else
                if(isnan(obj.T))
                    error('T is required');
                    return;
                end
                if(isnan(obj.h))
                    error('h is required');
                    return;
                end
                obj.L = ((obj.g * obj.T^2)/(2 * pi)) * (tanh( (((2 * pi)/obj.T) * sqrt(obj.h/obj.g))^(3/2) ))^(2/3);
                L = obj.L;
            end
        end
        
        function [c] = wave_celerity(obj)
            if obj.deep_water == 1
                if(isnan(obj.T))
                    error('T is required');
                    return;
                end
                %obj.c = obj.L/obj.T;
                obj.c = (obj.g * obj.T)/(2 * pi);
                c = obj.c;
            else
                if(isnan(obj.L))
                    error('L is required');
                    return;
                end
                if(isnan(obj.T))
                    error('T is required');
                    return;
                end
                obj.c = obj.L/obj.T;
                c = obj.c;
            end
        end
        
        function [c_g] = group_velosity(obj)
            if obj.deep_water == 1
                if(isnan(obj.c))
                    error('c is required');
                    return;
                end
                obj.c_g = (obj.c/2);
                c_g = obj.c_g;
            else
                if(isnan(obj.c))
                    error('c is required');
                    return;
                end
                if(isnan(obj.k))
                    error('k is required');
                    return;
                end
                if(isnan(obj.h))
                    error('h is required');
                    return;
                end
                obj.c_g = (obj.c/2) * (1 + ((2 * obj.k * obj.h)/(sinh(2 * obj.k * obj.h))));
                c_g = obj.c_g;
            end
            
        end
        
        function [ALPHA_disp,BETA_disp] = particle_displacement_amplitudes(obj)
            if(isnan(obj.H))
                error('H is required');
                return;
            end
            if(isnan(obj.k))
                error('k is required');
                return;
            end
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            if(isnan(obj.x))
                error('x is required');
                return;
            end
            if(isnan(obj.t))
                error('t is required');
                return;
            end
            if(isnan(obj.OMEGA))
                error('OMEGA is required');
                return;
            end
            obj.ALPHA_disp =  (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.BETA_disp =  (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            ALPHA_disp = obj.ALPHA_disp;
            BETA_disp = obj.BETA_disp;
        end
        
        function [u,w] = particle_velocites(obj)
            if(isnan(obj.OMEGA))
                error('OMEGA is required');
                return;
            end
            if(isnan(obj.H))
                error('H is required');
                return;
            end
            if(isnan(obj.k))
                error('k is required');
                return;
            end
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            if(isnan(obj.x))
                error('x is required');
                return;
            end
            if(isnan(obj.t))
                error('t is required');
                return;
            end
            obj.u = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            obj.w = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h))) * (cos((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            u = obj.u;
            w = obj.w;
            %t = obj.t
        end
        
        function [u_hat,w_hat] = velocity_amplitudes(obj)
            if(isnan(obj.L))
                error('L is required');
                return;
            end
            if(isnan(obj.k))
                error('k is required');
                return;
            end
            if(isnan(obj.OMEGA))
                error('OMEGA is required');
                return;
            end
            obj.u_hat = obj.OMEGA * (obj.H/2) * ((cosh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            obj.w_hat = obj.OMEGA * (obj.H/2) * ((sinh(obj.k * (obj.z + obj.h)))/(sinh(obj.k * obj.h)));
            u_hat = obj.u_hat;
            w_hat = obj.w_hat;
        end
        
        function [E] = energy_density(obj)
            if(isnan(obj.H))
                error('H is required');
                return;
            end
            obj.E = (obj.RHO * obj.g * obj.H^2)/8;
            E = obj.E;
        end
        
        function [H] = find_wave_height_from_Ap(obj,Ap)
            if(isnan(obj.k))
                error('k is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            obj.H = (2 * Ap * sinh(obj.k * obj.h))/(cosh(obj.k *(obj.z + obj.h)));
            H = obj.H;
        end
        
        function [ETA] = find_wave_ETA(obj)
            if(isnan(obj.K_p))
                error('K_p is required');
                return;
            end
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            if(isnan(obj.p_max))
                error('p_max is required');
                return;
            end
            obj.ETA =  (obj.p_max + (obj.RHO * obj.g * obj.z))/(obj.RHO * obj.g * obj.K_p);
            ETA = obj.ETA;
        end
        
        function [z] = find_depth_from_hydrostatic_pressure(obj)
            if(isnan(obj.ps))
                error('ps is required');
                return;
            end
            obj.z = obj.ps / -obj.RHO / obj.g;
            z = obj.z;
        end
        
        function [pd] = dynamic_pressure(obj)
            if(isnan(obj.ETA))
                error('ETA is required');
                return;
            end
            if(isnan(obj.K_p))
                error('K_p is required');
                return;
            end
            obj.pd = obj.RHO * obj.g * obj.ETA * obj.K_p;
            pd = obj.pd;
        end
        
        function [ps] = hydrostatic_pressure(obj)
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            obj.ps = -obj.RHO * obj.g * obj.z;
            ps = obj.ps;
        end
        
        function [H_n] = height_probibility(obj,n)
            % n is the percentage. eg, 5% = 5.
            if(isnan(obj.H_rms))
                error('H_rms is required');
                return;
            end
            if(isnan(obj.N))
                error('N is required');
                return;
            end
            obj.H_n = obj.H_rms * sqrt(-ln(n/obj.N));
            H_n = obj.H_n;
        end
        
        function [H_s] = estimate_significant_wave_height(obj,LineB)
            if(isnan(obj.H_avg))
                error('H_avg is required');
                return;
            end
            if(isnan(obj.H_s))
                error('H_s is required');
                return;
            end
            % LineA is from the left (n).
            % LineB is from the bottom (H_avg).
            obj.H_s = LineB; %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            obj.H_rms = obj.H_avg / LineB;
            obj.H_s = obj.H_s * obj.H_rms;
            H_s = obj.H_s;
        end
        
        function [K_p] = pressure_coefficient(obj)
            % Requires k,z,h
            if(isnan(obj.k))
                error('k is required');
                return;
            end
            if(isnan(obj.z))
                error('z is required');
                return;
            end
            if(isnan(obj.h))
                error('h is required');
                return;
            end
            obj.K_p = (cosh(obj.k * (obj.z + obj.h)))/(cosh(obj.k * obj.h));
            K_p = obj.K_p;
        end
        
        function [ETA] = sinusoidal_wave_shape(obj)
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            if(isnan(obj.L))
                error('L is required');
                return;
            end
            if(isnan(obj.OMEGA))
                error('OMEGA is required');
                return;
            end
            obj.ETA = (obj.H/2) * cos((obj.k * obj.x) - (obj.OMEGA * obj.t));
            ETA = obj.ETA;
        end
        
        function [PSI] = sine_wave_velocity_potential(obj)
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            if(isnan(obj.L))
                error('L is required');
                return;
            end
            if(isnan(obj.OMEGA))
                error('OMEGA is required');
                return;
            end
            obj.PSI = (-obj.g/obj.OMEGA) * (obj.H/2) * (cosh(obj.k * (obj.z + obj.h))/cosh(obj.k * obj.h)) * (sin((obj.k * obj.x)-(obj.OMEGA * obj.t)));
            PSI = obj.PSI;
        end
        
        %% A.1 Basic Wave Properties.
        function [k] = wave_number(obj)
            if(isnan(obj.L))
                error('L is required');
                return;
            end
            obj.k = (pi * 2) / obj.L;
            k = obj.k;
        end
        
        function [OMEGA] = angular_frequency(obj)
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            obj.OMEGA = (pi * 2) / obj.T;
            OMEGA = obj.OMEGA;
        end
        
        function [f] = average_frequency(obj)
            if(isnan(obj.T))
                error('T is required');
                return;
            end
            obj.f = 1 / obj.T;
            f = obj.f;
        end
           
    end
end

