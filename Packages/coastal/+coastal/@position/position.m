classdef position < coastal.numeretical_functions & matlab.mixin.CustomDisplay
    %POSITION A subclass of NUMERETICAL_FUNCTIONS
    %   This class provides a set of numeretical functions relative to a
    %   poticular position.
    %
    % See also NUMERETICAL_FUNCTIONS TRACER
    
    properties (Access = public)
        index = double(0);                  % Index. Currently unused.
        position_x = double(0);             % Current x particle position.
        position_y = double(0);             % Current y particle position.
        position_z = double(0);             % Current z particle position.
        position_x_original = double(0);    % Rotation point x co-ord.
        position_y_original = double(0);    % Rotation point y co-ord.
        position_z_original = double(0);    % Rotation point z co-ord.
        z_axis_lower = double(0);           % Lowest point on the seafloor.
    end
    
    methods
        function obj = position(x,y,z,H,T,h,z_axis_lower)
            %POSITION constructor for the object.
            %   This function is used by MatLAB to construct this object when
            %   it is created.
            persistent current_index
            if nargin == 0 % Object is created empty.
                super_args = {};
            else % Object is with values.
                super_args{1} = x;
                super_args{2} = y;
                super_args{3} = z;
                super_args{4} = H;
                super_args{5} = T;
                super_args{6} = h;
            end
            obj@coastal.numeretical_functions(super_args{:});
            
            if nargin ~= 0
                obj.z_axis_lower = z_axis_lower;
                obj.position_x_original = x;
                obj.position_y_original = y;
                obj.position_z_original = z;
            end
            % incremented index for collections.
            obj.index = current_index;current_index = current_index + 1;
        end
        
    end
    
    methods (Access = protected)
        
        function displayScalarObject(obj)
            
            header = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            header = sprintf('%s\n',header);
            disp(header);
            
            
            if ~isnan(obj.a);disp(['a = ',num2str(obj.a)]);end
            if ~isnan(obj.b);disp(['b = ',num2str(obj.b)]);end
            if ~isnan(obj.B);disp(['B = ',num2str(obj.B)]);end
            if ~isnan(obj.C_r);disp(['C_r = ',num2str(obj.C_r)]);end
            if ~isnan(obj.C_10);disp(['C_10 = ',num2str(obj.C_10)]);end
            if ~isnan(obj.c);disp(['c = ',num2str(obj.c)]);end
            if ~isnan(obj.c_g);disp(['c_g = ',num2str(obj.c_g)]);end
            if ~isnan(obj.D_n50);disp(['D_n50 = ',num2str(obj.D_n50)]);end
            if ~isnan(obj.E_f);disp(['E_f = ',num2str(obj.E_f)]);end
            if ~isnan(obj.E);disp(['E = ',num2str(obj.E)]);end
            if ~isnan(obj.f);disp(['f = ',num2str(obj.f)]);end
            if ~isnan(obj.F);disp(['F = ',num2str(obj.F)]);end
            if ~isnan(obj.ff);disp(['ff = ',num2str(obj.ff)]);end
            if ~isnan(obj.Fstar);disp(['Fstar = ',num2str(obj.Fstar)]);end
            if ~isnan(obj.Fstar_eff);disp(['Fstar_eff = ',num2str(obj.Fstar_eff)]);end
            if ~isnan(obj.GSLR);disp(['GSLR = ',num2str(obj.GSLR)]);end
            if ~isnan(obj.H);disp(['H = ',num2str(obj.H)]);end
            if ~isnan(obj.H_s);disp(['H_s = ',num2str(obj.H_s)]);end
            if ~isnan(obj.H_stoe);disp(['H_stoe = ',num2str(obj.H_stoe)]);end
            if ~isnan(obj.H_o);disp(['H_o = ',num2str(obj.H_o)]);end
            if ~isnan(obj.H_b);disp(['H_b = ',num2str(obj.H_b)]);end
            if ~isnan(obj.H_d);disp(['H_d = ',num2str(obj.H_d)]);end
            if ~isnan(obj.H_r);disp(['H_r = ',num2str(obj.H_r)]);end
            if ~isnan(obj.H_i);disp(['H_i = ',num2str(obj.H_i)]);end
            if ~isnan(obj.H2);disp(['H2 = ',num2str(obj.H2)]);end
            if ~isnan(obj.H_sb);disp(['H_sb = ',num2str(obj.H_sb)]);end
            if ~isnan(obj.h_b);disp(['h_b = ',num2str(obj.h_b)]);end
            if ~isnan(obj.H_BW);disp(['H_BW = ',num2str(obj.H_BW)]);end
            if ~isnan(obj.H_sig);disp(['H_sig = ',num2str(obj.H_sig)]);end
            if ~isnan(obj.H_avg);disp(['H_avg = ',num2str(obj.H_avg)]);end
            if ~isnan(obj.H_rms);disp(['H_rms = ',num2str(obj.H_rms)]);end
            if ~isnan(obj.H_max);disp(['H_max = ',num2str(obj.H_max)]);end
            if ~isnan(obj.h_cl);disp(['h_cl = ',num2str(obj.h_cl)]);end
            if ~isnan(obj.H_n);disp(['H_n = ',num2str(obj.H_n)]);end
            if ~isnan(obj.H_bar);disp(['H_bar = ',num2str(obj.H_bar)]);end
            if ~isnan(obj.H_per);disp(['H_per = ',num2str(obj.H_per)]);end
            if ~isnan(obj.h);disp(['h = ',num2str(obj.h)]);end
            if ~isnan(obj.h_toe);disp(['h_toe = ',num2str(obj.h_toe)]);end
            if ~isnan(obj.H_mo);disp(['H_mo = ',num2str(obj.H_mo)]);end
            if ~isnan(obj.Hstar_mo);disp(['Hstar_mo = ',num2str(obj.Hstar_mo)]);end
            if ~isnan(obj.K_D);disp(['K_D = ',num2str(obj.K_D)]);end
            if ~isnan(obj.k);disp(['k = ',num2str(obj.k)]);end
            if ~isnan(obj.k_o);disp(['k_o = ',num2str(obj.k_o)]);end
            if ~isnan(obj.K_p);disp(['K_p = ',num2str(obj.K_p)]);end
            if ~isnan(obj.K_prime);disp(['K_prime = ',num2str(obj.K_prime)]);end
            if ~isnan(obj.K_r);disp(['K_r = ',num2str(obj.K_r)]);end
            if ~isnan(obj.K_s);disp(['K_s = ',num2str(obj.K_s)]);end
            if ~isnan(obj.L);disp(['L = ',num2str(obj.L)]);end
            if ~isnan(obj.L_basin);disp(['L_basin = ',num2str(obj.L_basin)]);end
            if ~isnan(obj.L_bay);disp(['L_bay = ',num2str(obj.L_bay)]);end
            if ~isnan(obj.Lo);disp(['Lo = ',num2str(obj.Lo)]);end
            if ~isnan(obj.L_om);disp(['L_om = ',num2str(obj.L_om)]);end
            if ~isnan(obj.L_r);disp(['L_r = ',num2str(obj.L_r)]);end
            if ~isnan(obj.m);disp(['m = ',num2str(obj.m)]);end
            if ~isnan(obj.MWS);disp(['MWS = ',num2str(obj.MWS)]);end
            if ~isnan(obj.N);disp(['N = ',num2str(obj.N)]);end
            if ~isnan(obj.N_s);disp(['N_s = ',num2str(obj.N_s)]);end
            if ~isnan(obj.NeHu);disp(['NeHu = ',num2str(obj.NeHu)]);end
            if ~isnan(obj.n);disp(['n = ',num2str(obj.n)]);end
            if ~isnan(obj.P);disp(['P = ',num2str(obj.P)]);end
            if ~isnan(obj.p);disp(['p = ',num2str(obj.p)]);end
            if ~isnan(obj.PH);disp(['PH = ',num2str(obj.PH)]);end
            if ~isnan(obj.ps);disp(['ps = ',num2str(obj.ps)]);end
            if ~isnan(obj.pd);disp(['pd = ',num2str(obj.pd)]);end
            if ~isnan(obj.p_max);disp(['p_max = ',num2str(obj.p_max)]);end
            if ~isnan(obj.p_min);disp(['p_min = ',num2str(obj.p_min)]);end
            if ~isnan(obj.Qstar);disp(['Qstar = ',num2str(obj.Qstar)]);end
            if ~isnan(obj.Q_in);disp(['Q_in = ',num2str(obj.Q_in)]);end
            if ~isnan(obj.Q_out);disp(['Q_out = ',num2str(obj.Q_out)]);end
            if ~isnan(obj.Q_sinks);disp(['Q_sinks = ',num2str(obj.Q_sinks)]);end
            if ~isnan(obj.Q_grad);disp(['Q_grad = ',num2str(obj.Q_grad)]);end
            if ~isnan(obj.Q_x);disp(['Q_x = ',num2str(obj.Q_x)]);end
            if ~isnan(obj.Q_y);disp(['Q_y = ',num2str(obj.Q_y)]);end
            if ~isnan(obj.Q_sources);disp(['Q_sources = ',num2str(obj.Q_sources)]);end
            if ~isnan(obj.q);disp(['q = ',num2str(obj.q)]);end
            if ~isnan(obj.R);disp(['R = ',num2str(obj.R)]);end
            if ~isnan(obj.R_c);disp(['R_c = ',num2str(obj.R_c)]);end
            if ~isnan(obj.R_star);disp(['R_star = ',num2str(obj.R_star)]);end
            if ~isnan(obj.R_n);disp(['R_n = ',num2str(obj.R_n)]);end
            if ~isnan(obj.r);disp(['r = ',num2str(obj.r)]);end
            if ~isnan(obj.record_length);disp(['record_length = ',num2str(obj.record_length)]);end
            if ~isnan(obj.s);disp(['s = ',num2str(obj.s)]);end
            if ~isnan(obj.S_om);disp(['S_om = ',num2str(obj.S_om)]);end
            if ~isnan(obj.SWL);disp(['SWL = ',num2str(obj.SWL)]);end
            if ~isnan(obj.SLR);disp(['SLR = ',num2str(obj.SLR)]);end
            if ~isnan(obj.T);disp(['T = ',num2str(obj.T)]);end
            if ~isnan(obj.Tstar_p);disp(['Tstar_p = ',num2str(obj.Tstar_p)]);end
            if ~isnan(obj.T_p);disp(['T_p = ',num2str(obj.T_p)]);end
            if ~isnan(obj.T_om);disp(['T_om = ',num2str(obj.T_om)]);end
            if ~isnan(obj.T_m);disp(['T_m = ',num2str(obj.T_m)]);end
            if ~isnan(obj.t);disp(['t = ',num2str(obj.t)]);end
            if ~isnan(obj.t_d);disp(['t_d = ',num2str(obj.t_d)]);end
            if ~isnan(obj.tstar);disp(['tstar = ',num2str(obj.tstar)]);end
            if ~isnan(obj.T_z);disp(['T_z = ',num2str(obj.T_z)]);end
            if ~isnan(obj.T_n);disp(['T_n = ',num2str(obj.T_n)]);end
            if ~isnan(obj.U_10);disp(['U_10 = ',num2str(obj.U_10)]);end
            if ~isnan(obj.U_g);disp(['U_g = ',num2str(obj.U_g)]);end
            if ~isnan(obj.u);disp(['u = ',num2str(obj.u)]);end
            if ~isnan(obj.u_hat);disp(['u_hat = ',num2str(obj.u_hat)]);end
            if ~isnan(obj.W);disp(['W = ',num2str(obj.W)]);end
            if ~isnan(obj.W_n50);disp(['W_n50 = ',num2str(obj.W_n50)]);end
            if ~isnan(obj.w);disp(['w = ',num2str(obj.w)]);end
            if ~isnan(obj.w_hat);disp(['w_hat = ',num2str(obj.w_hat)]);end
            if ~isnan(obj.x);disp(['x = ',num2str(obj.x)]);end
            if ~isnan(obj.y);disp(['y = ',num2str(obj.y)]);end
            if ~isnan(obj.xtv);disp(['xtv = ',num2str(obj.xtv)]);end
            if ~isnan(obj.X_T);disp(['X_T = ',num2str(obj.X_T)]);end
            if ~isnan(obj.y_T);disp(['y_T = ',num2str(obj.y_T)]);end
            if ~isnan(obj.z);disp(['z = ',num2str(obj.z)]);end
            if ~isnan(obj.z_100);disp(['z_100 = ',num2str(obj.z_100)]);end
            if ~isnan(obj.zwm);disp(['zwm = ',num2str(obj.zwm)]);end
            
            
            if ~isnan(obj.ALPHA);disp(['ALPHA = ',num2str(obj.ALPHA)]);end
            if ~isnan(obj.ALPHA_disp);disp(['ALPHA_disp = ',num2str(obj.ALPHA_disp)]);end
            if ~isnan(obj.ALPHA_NU);disp(['ALPHA_NU = ',num2str(obj.ALPHA_NU)]);end
            if ~isnan(obj.BETA);disp(['BETA = ',num2str(obj.BETA)]);end
            if ~isnan(obj.BETA_disp);disp(['BETA_disp = ',num2str(obj.BETA_disp)]);end
            if ~isnan(obj.tanBETA);disp(['tanBETA = ',num2str(obj.tanBETA)]);end
            if ~isnan(obj.DELTA);disp(['DELTA = ',num2str(obj.DELTA)]);end
            if ~isnan(obj.DELTA_x);disp(['DELTA_x = ',num2str(obj.DELTA_x)]);end
            if ~isnan(obj.DELTA_z);disp(['DELTA_z = ',num2str(obj.DELTA_z)]);end
            if ~isnan(obj.ETA);disp(['ETA = ',num2str(obj.ETA)]);end
            if ~isnan(obj.ETAbar);disp(['ETAbar = ',num2str(obj.ETAbar)]);end
            if ~isnan(obj.ETAbar_wind);disp(['ETAbar_wind = ',num2str(obj.ETAbar_wind)]);end
            if ~isnan(obj.GAMMA_b);disp(['GAMMA_b = ',num2str(obj.GAMMA_b)]);end
            if ~isnan(obj.GAMMA_r);disp(['GAMMA_r = ',num2str(obj.GAMMA_r)]);end
            if ~isnan(obj.OMEGA);disp(['OMEGA = ',num2str(obj.OMEGA)]);end
            if ~isnan(obj.MU_max);disp(['MU_max = ',num2str(obj.MU_max)]);end
            if ~isnan(obj.NU_max);disp(['NU_max = ',num2str(obj.NU_max)]);end
            if ~isnan(obj.PHI);disp(['PHI = ',num2str(obj.PHI)]);end
            if ~isnan(obj.PSI);disp(['PSI = ',num2str(obj.PSI)]);end
            if ~isnan(obj.RHO_s);disp(['RHO_s = ',num2str(obj.RHO_s)]);end
            if ~isnan(obj.SIGMA);disp(['SIGMA = ',num2str(obj.SIGMA)]);end
            if ~isnan(obj.TAO_w);disp(['TAO_w = ',num2str(obj.TAO_w)]);end
            if ~isnan(obj.TAO_b);disp(['TAO_b = ',num2str(obj.TAO_b)]);end
            if ~isnan(obj.THETA);disp(['THETA = ',num2str(obj.THETA)]);end
            if ~isnan(obj.THETA_b);disp(['THETA_b = ',num2str(obj.THETA_b)]);end
            if ~isnan(obj.THETA_o);disp(['THETA_o = ',num2str(obj.THETA_o)]);end
            if ~isnan(obj.XI);disp(['XI = ',num2str(obj.XI)]);end
            if ~isnan(obj.XI_o);disp(['XI_o = ',num2str(obj.XI_o)]);end
            if ~isnan(obj.XI_b);disp(['XI_b = ',num2str(obj.XI_b)]);end
            if ~isnan(obj.ZETA_u);disp(['ZETA_u = ',num2str(obj.ZETA_u)]);end
            
            
            
            if ~isnan(obj.deep_water);disp(['deep_water = ',num2str(obj.deep_water)]);end
            if ~isnan(obj.wave_is_broken);disp(['wave_is_broken = ',num2str(obj.wave_is_broken)]);end
            
            
            disp('');
            
        end
        
    end
    
    methods
        function obj = set.position_x_original(obj,value)
            obj.x = value;
            obj.position_x = value;
            obj.position_x_original = value;
        end
        
        function obj = set.position_z_original(obj,value)
            if value <= 0
                obj.z = value;
                obj.position_z = value;
                obj.position_z_original = value;
            else
                error('z value set to position above water level');
            end
        end
        
        function obj = set.position_y_original(obj,value)
            obj.y = value;
            obj.position_y = value;
            obj.position_y_original = value;
        end
        
        function [do_plot,return_position_x,return_position_y,return_position_z] = plot_position(obj,t)
            obj.t = t;
            if obj.position_z_original > -obj.h
                do_plot = 1;
                obj.particle_velocites;
                obj.particle_displacement_amplitudes;
                return_position_x = double(obj.position_x_original) + (obj.ALPHA_disp);
                return_position_z = double(obj.position_z_original) + (obj.BETA_disp);
                return_position_y = obj.y;
            else
                do_plot = 0;
                return_position_x = obj.position_x_original;
                return_position_z = obj.position_z_original;
                return_position_y = obj.y;
            end
        end
    end
end


