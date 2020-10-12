classdef CartPendulum < handle
    %{
    Implemenation of inverted pendulum on a cart.
    
                     O m
                  l /
            ______ /__
           |      /   | M
      f -->|     /    |
           |____|_____|
            0   |     0
                |
                |<--\theta = 0
                O
       |-----> X
    
    The state is position of the cart X, angle of the pendulum \theta
    the velocity of the cart dX and the angular velocity d\theta.
    The input is the force f acting on the cart.
    
               x:=[ X, \theta, dX, d\theta ]^\top
               u:= f
    %}
    
    properties
        
        n          = 4     % state dimension
        p          = 1     % input dimension
        samplerate = 0.05  % sample rate
        
    end
    
    methods
        
        % === Constructor ===
        function obj = CartPendulum(the_state)
            obj.state = the_state;
        end
        
        % === Setters ===
        function set_samplerate(obj,the_samplerate)
            obj.samplerate = the_samplerate;
        end
        
        % === Getters ===
        function y = measure(obj)
            %{
            TODO: add output here
            e.g. y = [X, \theta]^\top
            and add noise
            %}
            y = obj.state;
        end
        
        % === methods ===
        function x_plus = integrate(obj, u)
            t  = 0;             % init time
            dt = obj.intg_step; % length of ruku integration
            x  = obj.state;
            
            while t < obj.samplerate
                k1 = obj.f(x,         u);
                k2 = obj.f(x+dt/2*k1, u);
                k3 = obj.f(x+dt/2*k2, u);
                k4 = obj.f(x+dt*k3,   u);
                x_plus = x + dt/6*(k1+2*k2+2*k3+k4);
                t = t + dt;
                x = x_plus;
            end
            
            obj.state = x;
            
        end % integrate
        
        function [A, B] = linearize(obj, x_bar, u_bar)
            % declare indices
            Xind        = obj.X_ind;
            tind        = obj.theta_ind;
            dXind       = obj.dX_ind;
            dtind       = obj.dtheta_ind;
            % declare parameter for readability
            m           = obj.m;
            M           = obj.M;
            l           = obj.l;
            g           = obj.g;
            d           = obj.d;
            % precompute
            dtheta2     = x_bar(dtind)^2;
            sin_th      = sin( x_bar(tind) );
            cos_th      = cos( x_bar(tind) );
            mll         = m*l*l;
            tmp_1       = mll*(M+m*(1-cos_th^2));
            tmp_2       = m*l*mll*sin_th;
            tmp_3       = m*mll*sin_th^2;
            tmp_4       = m*mll*cos_th^2;
            tmp_5       = m*mll*cos_th*sin_th;
            % compute derivatives
            g0          = 1/tmp_1;
            dg0_dtheta  = -1/(tmp_1^2)*(2*m*mll*cos_th*sin_th);
            g1          = m*mll*g*cos_th*sin_th + tmp_2*dtheta2 - mll*d*x_bar(dXind);
            dg1_dtheta  = -tmp_3*g + tmp_4*g + m*l*mll*cos_th*dtheta2;
            g2          = -(m+M)*m*g*l*sin_th - tmp_5*dtheta2 + m*l*cos_th*d*x_bar(dXind) ;
            dg2_dtheta  = -(m+M)*m*g*l*cos_th + tmp_3*dtheta2 - tmp_4*dtheta2 - m*l*sin_th*d*x_bar(dXind);
            
            
            df3_dtheta  = (g1 + mll*u_bar)*dg0_dtheta + g0*dg1_dtheta;
            df3_ddtheta = 2*tmp_2*g0*x_bar(dtind);
            df3_ddX     = -mll*d*g0;
            df4_dtheta  = (g2 - m*l*cos_th*u_bar)*dg0_dtheta + (dg2_dtheta + m*l*sin_th*u_bar)*g0;
            df4_ddtheta = -2*tmp_5*g0*x_bar(dtind);
            df4_ddX     = m*l*cos_th*d*g0;
            
            df3_du      = mll*g0;
            df4_du      = -m*l*cos_th*g0;
            
            % get jacobians (hard coded for now)
            A = [ 
              0,          0,       1,           0;
              0,          0,       0,           1;
              0, df3_dtheta, df3_ddX, df3_ddtheta;
              0, df4_dtheta, df4_ddX, df4_ddtheta
            ];
        
            B = [
              0;
              0;
              df3_du;
              df4_du
            ];
                
        end % linearize
        
        function visualize(obj)
            x           = obj.state;
            Xind        = obj.X_ind;
            tind        = obj.theta_ind;
            % define some dimensions
            y_line      = [ 0, 0];
            x_line      = [-2, 2];
            car_width   = 0.3;
            car_hight   = 0.2;
            ball_radius = 0.05;
            % define some collors
            my_red      = obj.car_color;
            my_green    = obj.ball_color;
            
            plot(x_line, y_line, 'k')
            xlabel('X[m]','Interpreter','Latex')
            axis equal
            hold on
            % plot car
            rectangle('Position', [x(Xind) - car_width/2, 0, car_width, car_hight] ,'FaceColor', my_red);
            % plot rod
            x_rod = [ x(Xind)    , x(Xind) + obj.l*sin(x(tind))];
            y_rod = [ car_hight/2, car_hight/2-obj.l*cos(x(tind))];
            plot( x_rod, y_rod, 'k' ,'LineWidth',2);
            % plot ball
            r = rectangle('Position', [ x_rod(2) - ball_radius, y_rod(2) - ball_radius, 2*ball_radius, 2*ball_radius]);
            r.Curvature = [1, 1];
            r.FaceColor = my_green;

            xlim([-2,2])
            ylim([-1,1])
            drawnow
            hold off
            
        end % visualizze
    
    end % public methods
    
    properties(Access = private)
        m = 0.3   % mass of pendulum
        M = 0.8   % mass of cart
        l = 0.6   % length of rod
        d = 2.0   % damping on the cart
        g = 9.81  % gravity
        
        intg_step = 0.001; % step size of the integrator
        
        X_ind      = 1;
        theta_ind  = 2;
        dX_ind     = 3;
        dtheta_ind = 4;
        
        state              % current nominal state
        
        car_color     = [153, 0, 0]/255;
        ball_color   = [0, 153, 0]/255;
        
    end
    
    methods(Access = private)
        
        function dx = f(obj, x, u)
            % declare indices
            Xind         = obj.X_ind;
            tind         = obj.theta_ind;
            dXind        = obj.dX_ind;
            dtind        = obj.dtheta_ind;
            % declare parameter for readability
            m            = obj.m;
            M            = obj.M;
            l            = obj.l;
            g            = obj.g;
            d            = obj.d;
            % precalculate terms
            sin_th       = sin( x(tind) );
            cos_th       = cos( x(tind) );
            mll          = m*l*l;
            tmp          = mll*( M+m*(1-cos_th^2) );
            
            dx           = zeros(obj.n, 1);
            
            dx(Xind, 1)  = x(dXind);
            dx(tind, 1)  = x(dtind);
            dx(dXind,1)  = (1/tmp)*( m*mll*g*cos_th*sin_th + mll*(m*l*x(dtind)^2*sin_th - d*x(dXind)) ) + mll/tmp*u;
            dx(dtind, 1) = (1/tmp)*( -(m+M)*m*g*l*sin_th - m*l*cos_th*(m*l*x(dtind)^2*sin_th - d*x(dXind)) ) - m*l*cos_th/tmp*u;
            
        end % f
                
    end % private methods
end

