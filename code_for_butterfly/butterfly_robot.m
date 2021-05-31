classdef butterfly_robot
    %Differential equation:
    %M(q)dot_dot_q + C(q, dot_q)dot_q + G(q) = [u;0]
    properties
        a = 0.1148; %Constant for describing the butterfly frame
        b = 0.04022;%Constant for describing the butterfly frame

        %c = 0.49; %Constant for virtual holonomic constraint.
        g = 9.81; %Gravity
        %J_s = 4.44556451612904e-07;%5.48e-7;%5.48e-7; %Mass moment of inertia of the ball 
        %J_f = 1.58117e-3; %Mass moment of inertia of the frame, from article.
        J_f = 8.8015e-04;
        %J_f = 1.8e-3/2; %From lab in modrob.
        m_b = 3e-3; %Mass of the ball
        %m_f = 0.15; %Mass of the frame
        r_b = 16.6;
        %R_b = sqrt(16.55e-3^2-12.5e-3^2); From paper.
        R_b = sqrt(16.6e-3^2-13.5e-3^2); %Measured on the robot
        J_s = 2/3*3e-3*16.6e-3^2;
        delta
        tau_delta
        rho
        length_rho
        bad_rho
        kappa_frame
        %alpha
        curve_for_phi        
        ds
        dds
        dsf
        ddsf
        tau
        kappa
        R
        diff_R
        tau_diff 
        theta
        d_theta
        dd_theta
        all_theta
        K
        function_for_X
        function_for_dphi
        spline_dphi
        %% Change these to tune the LQR
        Gamma =@(t) 1
        Q_integral = @(t)[40 0 0 0;0 40 0 0;0 0 20 0;0 0 0 0.1];
        Q = @(t) [10 0 0;0 10 0;0 0 10];
        %%
        X
        phi_test
        A
        B
        Phi
        dPhi
        ddPhi
        dddPhi
        normal_vector
        g_fun
        dg_fun
        ddg_fun
        drho_test 
        ys
        test
        spline_delta
        P
        function_for_P
    end
    methods
        function obj = butterfly_robot(calculate_riccati,with_integral)
            %% phi = f(varphi), varphi = g(varphi) 
            e_phi =@(phi) [sin(phi);cos(phi);0];
            d_e_phi =@(phi) [cos(phi);-sin(phi);0];

            %% Delta
            obj.spline_delta = load('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\shape.mat');
%             delta = @(x) ppval(obj.spline_delta.spline_delta,x);
%             d_delta = @(x) ppval(obj.spline_delta.spline_d_delta,x);
            %delta =@(phi) (obj.a-obj.b*cos(2*phi));
            a =      0.1148;
            b =    -0.04023;
            c =   -0.001535;
            d =    0.001152;
            e =   0.0008215;
            delta =  @(phi) a+b*cos(2*phi)+c*cos(4*phi)+d*cos(6*phi)+e*cos(8*phi);
            d_delta =  @(phi) -b*2*sin(2*phi)-c*4*sin(4*phi)-d*6*sin(6*phi)-e*8*sin(8*phi);

            delta_v = @(phi) delta(phi).*e_phi(phi);
            %d_delta = @(phi) obj.b*2*sin(2*phi);
            d_delta_v = @(phi) d_delta(phi).*e_phi(phi)+delta(phi).*d_e_phi(phi);
            norm_d_delta_v = @(phi) norm(d_delta_v(phi));



            %% Tau
            tau =@(phi) d_delta_v(phi)./norm_d_delta_v(phi);
            range_for_functions = linspace(0,2*pi,500);
            function_g = @(x) atan2([1 0 0]*delta_v(x)-obj.R_b*[0 1 0]*tau(x),...
                    [0 1 0]*delta_v(x)+obj.R_b*[1 0 0]*tau(x));
            res_fun_g = zeros(length(range_for_functions),1);
            res_fun_dg = zeros(length(range_for_functions),1);
            for i = 1:length(range_for_functions)
                res_fun_g(i) = function_g(range_for_functions(i));
            end
            res_fun_g = unwrap(res_fun_g);
            figure

            plot(res_fun_g,range_for_functions)
            hold on;
            grid on;
            k = spline(res_fun_g,range_for_functions);
            
            result_spline = @(x)ppval(k,mod(x,2*pi));
            subplot(1,2,1)
            plot(range_for_functions,result_spline(range_for_functions));
            grid on;
            obj.Phi = result_spline;
                   
            
            syms phi delta d_delta dd_delta ddd_delta;
            %% e_phi 
            e_phi = [sin(phi);cos(phi);0];
            d_e_phi = [cos(phi);-sin(phi);0];
            dd_e_phi = [-sin(phi);-cos(phi);0];
            ddd_e_phi = [-cos(phi);sin(phi);0];

            %% Delta
            delta_v =  delta*e_phi;
            d_delta_v =  d_delta*e_phi+delta*d_e_phi;
            dd_delta_v =  dd_delta*e_phi+2*d_delta*d_e_phi+delta*dd_e_phi;
            ddd_delta_v =  ddd_delta*e_phi+ 3*dd_delta*d_e_phi + ...
                3*d_delta*dd_e_phi + delta*ddd_e_phi;
            norm_d_delta_v = sqrt(d_delta^2+delta^2);
            obj.delta = matlabFunction(delta_v);

            %% Tau
            tau =  d_delta_v/norm_d_delta_v;
            d_tau =  dd_delta_v/norm_d_delta_v - d_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3;
            dd_tau =  ddd_delta_v/norm_d_delta_v - dd_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3 - ...
                (dd_delta_v*d_delta_v'*dd_delta_v+d_delta_v*(dd_delta_v'*dd_delta_v)+d_delta_v*d_delta_v'*ddd_delta_v)/...
                norm_d_delta_v^3 + d_delta_v*d_delta_v'*dd_delta_v*3*d_delta_v'*dd_delta_v/norm_d_delta_v^5;
            obj.tau = matlabFunction(tau);

            %% Normal vector:
            normal =  [0 -1 0;1 0 0;0 0 0]*tau;
            d_normal =  [0 -1 0;1 0 0;0 0 0]*d_tau;
            dd_normal =  [0 -1 0;1 0 0;0 0 0]*dd_tau;
            obj.normal_vector = matlabFunction(normal);                

            %% Rho :
            rho =  delta_v + obj.R_b*normal;
            length_rho =  norm(rho);
            d_rho =  d_delta_v + obj.R_b*d_normal;
            dd_rho =  dd_delta_v + obj.R_b*dd_normal;

            obj.rho = matlabFunction(rho);
            obj.length_rho = matlabFunction(length_rho);

            %% Need to find varphi = g(phi)
            gamma =  delta_v'*[1;0;0]-obj.R_b*tau'*[0;1;0];
            psi =  delta_v'*[0;1;0]+obj.R_b*tau'*[1;0;0];
            d_gamma =  d_delta_v'*[1;0;0]-obj.R_b*d_tau'*[0;1;0];
            d_psi =  d_delta_v'*[0;1;0]+obj.R_b*d_tau'*[1;0;0];
            dd_gamma =  dd_delta_v'*[1;0;0]-obj.R_b*dd_tau'*[0;1;0];
            dd_psi =  dd_delta_v'*[0;1;0]+obj.R_b*dd_tau'*[1;0;0];
            %func_g = atan2(gamma, psi);
            d_func_g = (d_gamma*psi-d_psi*gamma)/(gamma^2+psi^2);
            dd_func_g =  (dd_gamma*psi-dd_psi*gamma)/(psi^2+gamma^2)-(d_gamma*psi-d_psi*gamma)*(2*gamma*d_gamma+2*psi*d_psi)/(gamma^2+psi^2)^2;

            %% s
            d_s =  norm(d_rho)/d_func_g;
            dd_s =  d_rho'*dd_rho/norm(d_rho)/d_func_g^2 - norm(d_rho)*dd_func_g/d_func_g^3;
            obj.ds = matlabFunction(d_s);
            obj.dds = matlabFunction(dd_s);


            %% s_f:
            d_sf =  norm_d_delta_v/d_func_g;
            dd_sf =  d_delta_v'*dd_delta_v/norm_d_delta_v/d_func_g - norm_d_delta_v*dd_func_g/d_func_g^3;
            obj.dsf = matlabFunction(d_sf);
            obj.ddsf = matlabFunction(dd_sf);

            %% Kappa:
            kappa =  dd_rho/d_func_g^2/d_s^2;
            obj.kappa = matlabFunction(kappa);

            %% Rotations: 
            syms theta
            obj.R = matlabFunction([cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1]);
            obj.diff_R = matlabFunction([-sin(theta) -cos(theta) 0;cos(theta) -sin(theta) 0;0 0 0]);

            c1= -1.207115e+01; c2= 4.626067e-02; c3= 9.333719e-03; c4= -4.626870e-04; c5= -2.898753e-03;
            c1= -4.954183e+01; c2= 1.170221e-02; c3= 2.955707e-03; c4= 4.730799e-04; c5= -5.471393e-05;
            %%c1= -6.926763e+00; c2= 8.305011e-02; c3= 1.438304e-02; c4= -2.502948e-03; c5 = -4.531121e-03; c6= 8.731885e-04;

            %c1= 4.982470e+01; c2= -1.206215e-02; c3= -2.596353e-03; c4= 1.293764e-04; c5= 4.283631e-04;
%             c1= -6.227222e-01; c2= 1.474808e-02; c3= 1.555738e-01; c4= 1.927990e-01; c5= -1.102467e-04, c6= 6.171664e-02;
%             c1 =     -0.4189;
%             c2 =        2.67;
%             c3 =      0.9033;
%             c4 =     -0.1537;
%             c5 =     -0.2513;
%             c6 =     -0.0803;
            
            a= -2.739474e+01; b= 2.156690e-02; c= -1.323878e-03; d= -2.646793e-03; e= -1.770084e-03;
            syms varphi;
            Theta = c1*atan(c2*sin(2*varphi)+c3*sin(4*varphi)+c4*sin(6*varphi)+c5*sin(8*varphi))+varphi;
            dTheta = c1*(2*c2*cos(2*varphi) + 4*c3*cos(4*varphi) + 6*c4*cos(6*varphi) + 8*c5*cos(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1) + 1;
            ddTheta = c1*(-4*c2*sin(2*varphi) - 16*c3*sin(4*varphi) - 36*c4*sin(6*varphi) - 64*c5*sin(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1) - 2*c1*(2*c2*cos(2*varphi) + 4*c3*cos(4*varphi) + 6*c4*cos(6*varphi) + 8*c5*cos(8*varphi))^2*(c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1)^2;
            obj.theta = matlabFunction(Theta);
            obj.d_theta = matlabFunction(dTheta);
            obj.dd_theta = matlabFunction(ddTheta);
            
            %% Plots phase plane of alpha betta gamma function
            figure
            a = @(x) obj.alpha_beta_gamma(x);
            k = linspace(0,pi, 400);
            linear_term = @(x) -obj.get_dgamma(x(1))/([1 0 0]*a(x(1)));
            res = zeros(400,1);
            for i = 1:400
               res(i) = linear_term(k(i));%[0 0 1]*a(k(i))/([1 0 0]*a(k(i))); 
            end
            plot(k,res);
            f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];%f_plane = @(t,x) [x(2);-obj.get_dgamma(0)/([1 0 0]*a(0))*x(1)];
            %% Finding solution to periodic Riccati equation
            if calculate_riccati
                options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                [ts,ys] = ode45(f_plane,[0,10],[0;0.1455], options);
                i = 1;
                phi_dot = [0 0];
                length(ys)
                while ys(i,1) <= 2*pi
                    phi_dot(i,:) = ys(i,:);
                    i = i+1;
                    if i >= length(ys)
                        break
                    end
                end
                phi_dot(i,:) = ys(i,:);
                %% Interpolated 
                spline_dphi = spline(phi_dot(:,1), phi_dot(:,2));
                obj.spline_dphi = spline_dphi;
                %interpolation_for_dphi = @(x) interp1(phi_dot(:,1),phi_dot(:,2), x);
                obj.function_for_dphi = @(phi) ppval(spline_dphi,(mod(phi,pi)));
                figure();
                plot(ys(:,1),ys(:,2))
                if(with_integral)
                    [A, B] = obj.get_linearization(phi, obj.function_for_dphi,true,with_integral);          
                    obj.A = @(phi) A(phi);
                    obj.B = @(phi) B(phi);
                    [X,phi] = sdp_riccati(obj.A,obj.B,obj.Q_integral,obj.Gamma,0,pi,700,20,4);
                    obj.X = X;    
                    interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x) interp1(phi,reshape(X(1,4,:),1,length(X)),x);
                          interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x) interp1(phi,reshape(X(2,4,:),1,length(X)),x);
                          interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x) interp1(phi,reshape(X(3,4,:),1,length(X)),x);
                          interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x) interp1(phi,reshape(X(4,4,:),1,length(X)),x)];
                    obj.function_for_X = @(x) interpolation_for_X(mod(x,pi));
                else
                    [A, B, C] = obj.get_linearization(phi, obj.function_for_dphi,true,with_integral);          
                    obj.A = @(phi) A(phi);
                    obj.B = @(phi) B(phi);
                    [X,phi] = sdp_riccati(obj.A,obj.B,obj.Q,obj.Gamma,0,pi,400,20,3);
                    obj.X = X;
                    interpolation_for_X = @(x)[interp1(phi,reshape(X(1,1,:),1,length(X)),x) interp1(phi,reshape(X(1,2,:),1,length(X)),x) interp1(phi,reshape(X(1,3,:),1,length(X)),x);
                                               interp1(phi,reshape(X(2,1,:),1,length(X)),x) interp1(phi,reshape(X(2,2,:),1,length(X)),x) interp1(phi,reshape(X(2,3,:),1,length(X)),x);
                                               interp1(phi,reshape(X(3,1,:),1,length(X)),x) interp1(phi,reshape(X(3,2,:),1,length(X)),x) interp1(phi,reshape(X(3,3,:),1,length(X)),x)];
                    obj.function_for_X = @(x) interpolation_for_X(mod(x,pi));
                    A_kalman = @(t) obj.A(t)-obj.B(t)/obj.Gamma(t)*obj.B(t)'*obj.function_for_X(t);
                    E = @(t)[1 0 0;0 1 0;0 0 1];
                    R = @(t) [0.001^2 0 0;0 0.01^2 0;0 0 0.01^2]*obj.function_for_dphi(t);
                    G = @(t) eye(3)/obj.function_for_dphi(t);
                    ERE = @(t) E(t)*R(t)*E(t)';
                    Q = @(t)[0.0001^2 0 0;0 0.001^2 0;0 0 0.01^2]/obj.function_for_dphi(t);
                    GQG = @(t) G(t)*Q(t)*G(t)';
                    [P,phi] = sdp_riccati_kalman(A,C,GQG,ERE,0,pi,400,20,3);
                    obj.P = P;    
                    interpolation_for_P = @(x)[interp1(phi,reshape(P(1,1,:),1,length(P)),x) interp1(phi,reshape(P(1,2,:),1,length(X)),x) interp1(phi,reshape(P(1,3,:),1,length(P)),x);
                                               interp1(phi,reshape(P(2,1,:),1,length(P)),x) interp1(phi,reshape(P(2,2,:),1,length(X)),x) interp1(phi,reshape(P(2,3,:),1,length(P)),x);
                                               interp1(phi,reshape(P(3,1,:),1,length(P)),x) interp1(phi,reshape(P(3,2,:),1,length(X)),x) interp1(phi,reshape(P(3,3,:),1,length(P)),x)];
                    obj.function_for_P = @(x) interpolation_for_P(mod(x,pi));
                end
            end
        end
        
        function delta = get_delta(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,~,~,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            delta = obj.delta(delta,obj.Phi(varphi));
        end
        function tau = get_tau(obj,varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,~,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
            tau = obj.tau(d_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function normal_vector = get_normal_vector(obj,varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,~,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
            normal_vector = obj.normal_vector(d_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function rho = get_rho(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,~,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
           rho = obj.rho(d_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function l_rho = get_l_rho(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,~,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
            l_rho = obj.length_rho(d_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function kappa = get_kappa(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,dd_delta,ddd_delta] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
            %dd_delta = ppval(obj.spline_delta.spline_dd_delta,varphi);
            %ddd_delta = ppval(obj.spline_delta.spline_ddd_delta,varphi);
            kappa = obj.kappa(d_delta, ...
                dd_delta, ...
                ddd_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function ds = get_ds(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,dd_delta,~] = delta_5_cosine(varphi);
            %delta = ppval(obj.spline_delta.spline_delta,varphi);
            %d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
            %dd_delta = ppval(obj.spline_delta.spline_dd_delta,varphi);
            ds = obj.ds(d_delta, ...
                dd_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function dds = get_dds(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,dd_delta,ddd_delta] = delta_5_cosine(varphi);
%             delta = ppval(obj.spline_delta.spline_delta,varphi);
%             d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
%             dd_delta = ppval(obj.spline_delta.spline_dd_delta,varphi);
%             ddd_delta = ppval(obj.spline_delta.spline_ddd_delta,varphi);
            dds = obj.dds(d_delta, ...
                dd_delta, ...
                ddd_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function dsf = get_dsf(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,dd_delta,~] = delta_5_cosine(varphi);
%             delta = ppval(obj.spline_delta.spline_delta,varphi);
%             d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
%             dd_delta = ppval(obj.spline_delta.spline_dd_delta,varphi);
            dsf = obj.dsf(d_delta, ...
                dd_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function ddsf = get_ddsf(obj, varphi)
            varphi = mod(varphi,2*pi);
            [delta,d_delta,dd_delta,ddd_delta] = delta_5_cosine(varphi);
%             delta = ppval(obj.spline_delta.spline_delta,varphi);
%             d_delta = ppval(obj.spline_delta.spline_d_delta,varphi);
%             dd_delta = ppval(obj.spline_delta.spline_dd_delta,varphi);
%             ddd_delta = ppval(obj.spline_delta.spline_ddd_delta,varphi);
            ddsf = obj.ddsf(d_delta, ...
                dd_delta, ...
                ddd_delta, ...
                delta, ...
                obj.Phi(varphi));
        end
        function theta = get_theta(obj, varphi)
            theta = obj.theta(varphi);
            %[theta,~,~] = theta_paper(varphi,false);
        end
        function d_theta = get_dtheta(obj, varphi)
            d_theta = obj.d_theta(varphi);
            %[~,d_theta, ~] = theta_paper(varphi,false);
        end
        function dd_theta = get_ddtheta(obj, varphi)
            dd_theta = obj.dd_theta(varphi);
            %[~,~,dd_theta] = theta_paper(varphi,false);
        end
        function d_alpha = get_dalpha(obj,varphi)
            rxt = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            rxk = [0 0 1]*cross(obj.get_rho(varphi),obj.get_kappa(varphi));
            d_s = obj.get_ds(varphi);d_sf = obj.get_dsf(varphi);
            dd_s = obj.get_dds(varphi);dd_sf = obj.get_ddsf(varphi);
            dth = obj.get_dtheta(varphi);
            ddth = obj.get_ddtheta(varphi);
            JRm = obj.J_s/obj.m_b/obj.R_b;
            d_alpha = obj.m_b*((dd_s*rxt+rxk*d_s^2-dd_sf*JRm)*dth + ...
                (d_s*rxt-d_sf*JRm)*ddth + 2*dd_s + 2*dd_sf*JRm/obj.R_b);
        end
        function energy = get_energy(obj, q, dq)
            v_b = obj.diff_R(q(1))*dq(1)*obj.get_rho(q(2))...
                +obj.R(q(1))*obj.get_tau(q(2))*obj.get_ds(q(2))*dq(2);
            omega_b = [0;0;dq(1)-obj.get_dsf(q(2))/obj.R_b*dq(2)];
            omega_frame = dq(1);
            kinetic = 1/2*obj.m_b*(v_b')*v_b+1/2*obj.J_s*(omega_b')*omega_b+obj.J_s*omega_frame^2;
            potential = obj.m_b*[0 obj.g 0]*obj.R(q(1))*obj.get_rho(q(2));
            energy = [potential;kinetic];
        end
        
        function gamma_diff = get_dgamma(obj,varphi)
            gamma_diff = obj.m_b*obj.g*[0 1 0]*(...
                obj.diff_R(obj.get_theta(varphi))*obj.get_dtheta(varphi)*obj.get_tau(varphi)*obj.get_ds(varphi)+...
                obj.R(obj.get_theta(varphi))*obj.get_kappa(varphi)*obj.get_ds(varphi)^2+...
                obj.R(obj.get_theta(varphi))*obj.get_tau(varphi)*obj.get_dds(varphi));
        end
        
        function M = get_M(obj, q)
            varphi = q(2);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            m11 = obj.get_l_rho(varphi)^2+obj.J_f/obj.m_b+obj.J_s/obj.m_b;
            m12 = obj.get_ds(varphi)*rhoxtau-obj.get_dsf(varphi)*obj.J_s/obj.m_b/obj.R_b;
            m22 = obj.get_ds(varphi)^2+obj.get_dsf(varphi)^2*obj.J_s/obj.m_b/obj.R_b^2;
            M =   obj.m_b*[m11 m12;
                           m12 m22];
        end
        
        function C = get_C(obj, q, dq)
            varphi = q(2);
            ds = obj.get_ds(varphi);
            dsf = obj.get_dsf(varphi);
            dds = obj.get_dds(varphi);
            
            ddsf = obj.get_ddsf(varphi);
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi), obj.get_tau(varphi));
            taudotrho = obj.get_tau(varphi)'*obj.get_rho(varphi);
            rhoxkappa = [0 0 1]*cross(obj.get_rho(varphi),obj.get_kappa(varphi));
            c11 = obj.m_b*taudotrho*ds*dq(2);
            c12 = obj.m_b*(taudotrho*ds*dq(1) ...
                +(dds*rhoxtau-ddsf*obj.J_s/obj.m_b/obj.R_b ...
                +ds^2*rhoxkappa)*dq(2));
            c21 = -obj.m_b*ds*taudotrho*dq(1);
            c22 = obj.m_b*(ds*dds+dsf*ddsf*obj.J_s/obj.R_b^2/obj.m_b)*dq(2);
            
            C = [c11 c12;
                 c21 c22];
        end
        
        function G = get_G(obj, q)
            varphi = q(2);
            G = [obj.m_b*[0;obj.g;0]'*obj.diff_R(q(1))*obj.get_rho(varphi);
                 obj.m_b*[0;obj.g;0]'*obj.R(q(1))*obj.get_tau(varphi)*obj.get_ds(varphi)];
        end
        
        function ddq = calculate_ddq(obj, q, dq, u)
            C = obj.get_C(q,dq);
            G = obj.get_G(q);
            M = obj.get_M(q);
            ddq = (M^-1*(-C*dq-G+[u;0]));
        end
  
        function plot_constraints(obj)
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            figure 
            k = linspace(0,pi,500);
            %constraint_eq_1 = obj.get_ds(0)*(1+obj.J_s/obj.m_b*obj.get_p(0)^2)/([0 0 1]*cross(obj.get_rho(0),obj.get_tau(0))+obj.J_s/obj.m_b*obj.get_p(0));
            %constraint_eq_2 = obj.get_ds(pi/2)*(1+obj.J_s/obj.m_b*obj.get_p(pi/2)^2)/([0 0 1]*cross(obj.get_rho(pi/2),obj.get_tau(pi/2))+obj.J_s/obj.m_b*obj.get_p(pi/2));
            JRm = obj.J_s/obj.R_b/obj.m_b;
            not_eq_constraint = @(x) -(obj.get_ds(x)^2+obj.get_dsf(x)^2*JRm/obj.R_b)/ ...
                (obj.get_ds(x)*([0 0 1]*cross(obj.get_rho(x),obj.get_tau(x))) - obj.get_dsf(x)*JRm);
            hold on;
            constraint_gamma =@(x) -obj.get_tau(x)'*obj.diff_R(x)'*(obj.R(x)*obj.get_kappa(x)*obj.get_ds(x)-obj.R(x)*obj.get_tau(x)*obj.get_dds(x)/obj.get_ds(x));
            diff_gamma = @(x) [0 1 0]*(obj.diff_R(obj.get_theta(x))*obj.get_dtheta(x)*obj.get_tau(x)*obj.get_ds(x)+...
                obj.R(obj.get_theta(x))*obj.get_kappa(x)*obj.get_ds(x)^2+obj.R(obj.get_theta(x))*obj.get_tau(x)*obj.get_dds(x));
            results = zeros(length(k),1);
            theta = zeros(length(k),1);
            diff_theta = zeros(length(k),1);
            dd_theta = zeros(length(k),1);
            find_eq_points_atan_tau = zeros(length(k),1);
            constraint = zeros(length(k),1);
            for i = 1:length(k)
                find_eq_points_atan_tau(i,1) = atan2(-[0 1 0]*obj.get_tau(k(i)),[1 0 0]*obj.get_tau(k(i)));
                results(i,1) = not_eq_constraint(k(i));
                diff_theta(i,1) = obj.get_dtheta(k(i));
                theta(i,1) = obj.get_theta(k(i));
                dd_theta(i,1) = obj.get_ddtheta(k(i));
                constraint(i) = constraint_gamma(k(i));
            end
            plot(k,theta);
            plot(k,dd_theta);
            %plot(k,obj.get_theta(k));
            %plot(k, unwrap(find_eq_points_atan_tau));
            eq_points = [0;pi/16;pi/4;pi/2;pi/2+pi/4;pi-pi/16];
            scatter(eq_points,arrayfun(constraint_gamma, eq_points));
            not_eq_constraint(pi/2)
            constraint_gamma(0)
            constraint_gamma(pi/2)
            plot(k,diff_theta);
            plot(k,constraint);
            grid on;
            legend("Constraint all $\varphi$", "Constraint eq. point  $\varphi = n\pi$", "$\Theta'(\varphi)$");
            title("Constraints for $\Theta'(\varphi)$")
            
            set(gca,'XTick',0:pi/4:pi) 
            set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$'})
        end

        function normal_force = get_normal_force(obj, q)
            %Using a simplified model. 
            % Assuming the normal vector is 90 deg anti-clockwise to tau.
            n = [0 -1 0;1 0 0;0 0 1]*obj.get_tau(q(2));
            normal_force = [0;9.81;0]'*obj.R(q(1))*n;
        end
        function position_ball = get_position_ball(obj, q)
            ball_position = obj.get_rho(q(2));
            position_ball = obj.R(q(1))*ball_position;
        end
        function velocity_ball = get_velocity_ball(obj, q, dq)
            txRrho = cross([0;0;dq(1)],obj.R(q(1))*obj.get_rho(q(2)));
            velocity_ball = txRrho + obj.R(q(1))*obj.get_tau(q(2))*obj.get_ds(q(2))*dq(2);
        end
        function M = make_movie_butterfly_robot(obj, results, fps)
            h = figure;
            time = results.tout(end);
            number_of_frames = round(time*fps);
            x_frame = zeros(1,length(results.tout));
            y_frame = zeros(1,length(results.tout));
            x_ball = obj.R_b*cos(linspace(0, 2*pi,20));
            y_ball = obj.R_b*sin(linspace(0, 2*pi,20));
            xlim([-0.2 0.2])
            ylim([-0.2 0.2])
            axis equal;
            grid on;
            M(number_of_frames) = struct('cdata',[],'colormap',[]);
            %h.Visible = 'off';
            j = 0;
            for i = linspace(0,2*pi)
                j = j + 1;
                position = obj.delta(i);
                x_frame(j) = position(1);
                y_frame(j)= position(2);
            end
            frame = hgtransform;
            ball = hgtransform;
            patch('XData',x_frame,'YData',y_frame,'FaceColor','yellow','Parent',frame) 
            patch('XData',x_ball,'YData',y_ball,'FaceColor','red','Parent',ball) 
            
            %start_falling = length(results.after_fall(:,1))-length(find(results.after_fall(:,1)));
            acumulated_time = 0;
            current_frame = 0;
            for i = 2:length(results.tout);
                acumulated_time = acumulated_time + results.tout(i) - results.tout(i-1);
                if acumulated_time >= 1/fps
                    ball_position = obj.rho(results.q(i,2));
                    ball_position_inertial = obj.R(results.q(i,1))*ball_position;
                    frame.Matrix = makehgtform('zrotate',results.q(i,1));
                    ball.Matrix = makehgtform('translate', ball_position_inertial);
                    drawnow
                    current_frame = current_frame + 1
                    M(current_frame) = getframe;
                    acumulated_time = acumulated_time - 1/fps;
                end
            end
%             for i = start_falling+1:length(results.after_fall(:,1))
%                acumulated_time = acumulated_time + results.tout(i)-results.tout(i-1);
%                if acumulated_time >= 1/fps
%                     frame.Matrix = makehgtform('zrotate',results.after_fall(i,1));
%                     ball.Matrix = makehgtform('translate', [results.after_fall(i,2);results.after_fall(i,3);0]);
%                     drawnow
%                     current_frame = current_frame + 1
%                     M(current_frame) = getframe;
%                     acumulated_time = acumulated_time - 1/fps;
%                 end
%             end
%             h.Visible = 'on';
            movie(M,1,fps);
        end
        
        
        
        %% Controller Stuff
        
        function a = alpha_beta_gamma(obj, varphi)
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            rhodottau = obj.get_rho(varphi)'*obj.get_tau(varphi);
            ds_ = obj.get_ds(varphi); dsf_ = obj.get_dsf(varphi);
            dds_ = obj.get_dds(varphi); ddsf_ = obj.get_ddsf(varphi);
            JRm = obj.J_s/obj.R_b/obj.m_b;
            alpha = obj.m_b*((ds_*rhoxtau-dsf_*JRm)*obj.get_dtheta(varphi)+ds_^2+dsf_^2*JRm/obj.R_b);
            beta = obj.m_b*((ds_*rhoxtau-dsf_*JRm)*obj.get_ddtheta(varphi)-rhodottau*ds_*obj.get_dtheta(varphi)^2 ...
                + ds_*dds_ + dsf_*ddsf_*JRm/obj.R_b);
            gamma = obj.m_b*[0 obj.g 0]*obj.R(obj.get_theta(varphi))*obj.get_tau(varphi)*ds_;
            
            a = [alpha;beta;gamma];
        end
        function g = get_g_w_y_on_trajectory(obj, varphi, d_varphi)
            rhoxtau = [0 0 1]*cross(obj.get_rho(varphi),obj.get_tau(varphi));
            g_w = -obj.m_b*(obj.get_ds(varphi)*rhoxtau-obj.get_dsf(varphi)*obj.J_s/obj.m_b/obj.R_b);
            g_y_dot = obj.m_b*obj.get_ds(varphi)*obj.get_tau(varphi)'*obj.get_rho(varphi)*(2*obj.get_dtheta(varphi)*d_varphi);
            g_y = -obj.get_ds(varphi)*obj.m_b*obj.g*[cos(obj.get_theta(varphi)) -sin(obj.get_theta(varphi)) 0]*obj.get_tau(varphi);
            g = [g_w;g_y_dot;g_y];
        end
        function [A, B,C] = get_linearization(obj, phi, phi_dot, function_handles,with_integral)
            if function_handles == true
                abg = @(x)obj.alpha_beta_gamma(x);
                gwy = @(x)obj.get_g_w_y_on_trajectory(x, phi_dot(x));
                if with_integral
                    A = @(x)[0 1 0 0;
                             0 0 0 0;
                            [0 0 1]*gwy(x)/([1 0 0]*abg(x)) [0 1 0]*gwy(x)/([1 0 0]*abg(x)) ([0 0 1]*abg(x)-[0 1 0]*abg(x)*phi_dot(x)^2)/([1 0 0]*abg(x))/phi_dot(x) 0;
                             0 0 1 0];
                    B = @(x) [0;1;[1 0 0]*gwy(x)/([1 0 0]*abg(x));0];
                else
                    A = @(x)[0 1 0;
                             0 0 0;
                            [0 0 1]*gwy(x)/([1 0 0]*abg(x)) [0 1 0]*gwy(x)/([1 0 0]*abg(x)) ([0 0 1]*abg(x)-[0 1 0]*abg(x)*phi_dot(x)^2)/([1 0 0]*abg(x))/phi_dot(x)]/phi_dot(x);
                    B = @(x) [0;1;[1 0 0]*gwy(x)/([1 0 0]*abg(x))]/phi_dot(x);
                    C = @(x) [1 0 0;0 1 0;0 0 1];
%                     A = @(x)[0 1/phi_dot(x) 0;
%                              0 0 0;
%                             2*[0 0 1]*gwy(x)/([1 0 0]*abg(x)) 2*[0 1 0]*gwy(x)/([1 0 0]*abg(x)) (2*[0 1 0]*abg(x))/([1 0 0]*abg(x))];
%                     B = @(x) [0;1/phi_dot(x);[1 0 0]*gwy(x)/([1 0 0]*abg(x))];
%                     C = @(x) [1 0 0;0 1 0;0 0 1];
                end
            else
                abg = obj.alpha_beta_gamma(phi);
                gwy = obj.get_g_w_y_on_trajectory(phi, phi_dot);
                A = [0 1 0;
                     0 0 0;
                     gwy(3)/abg(1) gwy(2)/abg(1) (abg(3)-abg(2)*phi_dot^2)/abg(1)/phi_dot];
                B = [0;1/phi_dot(x);2*gwy(1)/abg(1)];
            end
        end
        function w = get_w(obj, riccati_sol, q, epsilon)
            w = -(1/obj.Gamma)*obj.B(q(2))'*riccati_sol*epsilon;
        end
        function out = riccati_times_epsilon(obj, q, epsilon)
            out = obj.function_for_X(mod(q(2),2*pi))*epsilon;
        end
        function u = get_u(obj, q, dq, epsilon)
            L = [1 obj.get_dtheta(q(2));
                 0 1];
            N = [obj.get_ddtheta(q(2))*dq(2)^2;
                 0];
            M = obj.get_M([obj.get_theta(q(2));q(2)]);
            C = obj.get_C([obj.get_theta(q(2));q(2)],[obj.get_dtheta(q(2))*dq(2);dq(2)]);
            G = obj.get_G([obj.get_theta(q(2));q(2)]);
            denom = (L^-1)*(M^-1);
            numerator = L^-1*(N+(M^-1)*C*L*[epsilon(2);dq(2)]+(M^-1)*G)*1;
            numerator = numerator(1)+obj.get_w(obj.function_for_X(q(2)),q, epsilon);
            u = numerator/denom(1,1);
        end
        function u = get_my_u(obj, q, dq, epsilon)
            M = obj.get_M([obj.get_theta(q(2));q(2)]);
            C = obj.get_C([obj.get_theta(q(2));q(2)],[obj.get_dtheta(q(2))*obj.function_for_dphi(q(2));obj.function_for_dphi(q(2))]);
            G = obj.get_G([obj.get_theta(q(2));q(2)]);
            k = [1 0]*M*[obj.get_dtheta(q(2));1]/([obj.get_dtheta(q(2)) 1]*M*[0;1]);
            u = k*([k^-1 -1]*M*[obj.get_w(obj.function_for_X(q(2)),q, epsilon)+obj.get_ddtheta(q(2))*obj.function_for_dphi(q(2))^2;0]+[k^-1 -1]*C*[(obj.get_dtheta(q(2)))*obj.function_for_dphi(q(2));dq(2)]+[k^-1 -1]*G);
        end
        function epsilon = get_epsilon(obj, q, dq)
            epsilon = [q(1)-obj.get_theta(q(2));dq(1)-obj.get_dtheta(q(2))*dq(2);dq(2)-obj.function_for_dphi(q(2))];
        end
        function plot_phase_plane(obj, x_0,x_end,y_0,y_end,title,time)
            a = @(x) obj.alpha_beta_gamma(x);
            f_plane = @(t,x) [x(2);-1/(a(x(1))'*[1;0;0])*((a(x(1))'*[0;1;0])*x(2)^2+a(x(1))'*[0;0;1])];
            phase_plot_2_interactive(f_plane,[x_0 x_end;y_0 y_end],time,title,[20,20],0.1);
        end
        function plot_abg(obj)
            set(groot, 'defaultAxesTickLabelInterpreter','latex');
            set(groot, 'defaultLegendInterpreter','latex');
            phi = linspace(0,2*pi,1000);
            abg = zeros(3,1,length(phi));
            for i = 1:length(phi)
                abg(:,1,i) = obj.alpha_beta_gamma(phi(i));
            end
            figure
            k = 1
            y_label = ["$\alpha$" "$\beta$" "$\gamma$"]
            for i = 1:3
                subplot(3,1,i)
                hold on;
                plot(phi,reshape(abg(i,1,:),1,length(phi)));
                set(gca,'XTick',0:pi/4:2*pi) 
                set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$' ...
                    ,'$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
                k = k+1;
                grid on;
                xlabel('$\varphi$', 'interpreter', 'latex')
                ylabel(y_label(i), 'interpreter', 'latex')
            end
            grid on
            sgtitle('$\alpha, \beta, \gamma$ of $\phi$', 'Interpreter','latex');
        end  
    end
end