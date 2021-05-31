function [u,current_epsilon,epsilon_measured] = controller(theta,varphi,dtheta,dvarphi,t,with_luenberger,luenberger) %#codegen
    q = [theta;varphi];
    dq = [dtheta;dvarphi];
    persistent spline_phi spline_riccati spline_dphi s_delta s_d_delta s_dd_delta s_ddd_delta prev_time epsilon spline_kalman phi_prev
    if isempty(spline_phi)
        %phi_prev = q(2);
        prev_time = t;
        %integral_term = 0;
        solutions = coder.load('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\solutions_riccati_and_dphi.mat');
        shape = coder.load('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\shape.mat');
        [breaks,coefs,~,~,dim] = unmkpp(solutions.varphi_to_phi);
        spline_phi = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(solutions.spline_riccati_sol);
        spline_riccati = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(solutions.spline_dphi);
        spline_dphi = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(shape.spline_delta);
        s_delta = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(shape.spline_d_delta);
        s_d_delta = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(shape.spline_dd_delta);
        s_dd_delta = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(shape.spline_ddd_delta);
        s_ddd_delta = mkpp(breaks,coefs,dim);
        [breaks,coefs,~,~,dim] = unmkpp(solutions.spline_kalman);
        spline_kalman = mkpp(breaks,coefs,dim);
    end
%     a = 0.1095; %Constant for describing the butterfly frame
%     b = 0.0405; %Constant for describing the butterfly frame
%     a = 0.1150;
%     b = 0.04;
%     a = 0.1148; %Constant for describing the butterfly frame
%     b = 0.04022;
    g = 9.81; %Gravity
    %J_s = 4.44556451612904e-07;%5.48e-7;%5.48e-7; %Mass moment of inertia of the ball 
    J_f = 8.8015e-04; %Mass moment of inertia of the frame, from article.
    m_b = 3e-3; %Mass of the ball
    %R_b = sqrt(16.55e-3^2-12.5e-3^2); From paper.
    r_b = 16.6e-3;
    J_s = 2/3*m_b*r_b^2;
    R_b = sqrt(r_b^2-13.5e-3^2); %Measured on the robot


    varphi = q(2);

    %c1= -6.012854e-01; c2= 1.415366e+00; c3= 3.964811e-01; c4= 5.645448e-02; c5= 2.783807e-03;
    %c1= -1.207115e+01; c2= 4.626067e-02; c3= 9.333719e-03; c4= -4.626870e-04; c5= -2.898753e-03;
    c1= -6.926763e+00; c2= 8.305011e-02; c3= 1.438304e-02; c4= -2.502948e-03; c5 = -4.531121e-03;
    Theta = c1*atan(c2*sin(2*varphi)+c3*sin(4*varphi)+c4*sin(6*varphi)+c5*sin(8*varphi))+varphi;
    d_Theta = c1*(2*c2*cos(2*varphi) + 4*c3*cos(4*varphi) + 6*c4*cos(6*varphi) + 8*c5*cos(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1) + 1;
    dd_Theta = c1*(-4*c2*sin(2*varphi) - 16*c3*sin(4*varphi) - 36*c4*sin(6*varphi) - 64*c5*sin(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1) - 2*c1*(2*c2*cos(2*varphi) + 4*c3*cos(4*varphi) + 6*c4*cos(6*varphi) + 8*c5*cos(8*varphi))^2*(c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))/((c2*sin(2*varphi) + c3*sin(4*varphi) + c4*sin(6*varphi) + c5*sin(8*varphi))^2 + 1)^2;
    
    phi = ppval(spline_phi,mod(varphi,2*pi));
    phi = phi(1,1); %Need to do this because of simulink coder...
    d_phi_star = ppval(spline_dphi, mod(varphi,pi));
    %d_phi_star = d_phi_star(1,1);
    riccati_sol = ppval(spline_riccati, mod(varphi,pi));
    sol_riccati = reshape(riccati_sol,3,3);
    %riccati_sol = reshape(riccati_sol,3,3);
    %% e_phi 
    e_phi = [sin(phi);cos(phi);0];
    d_e_phi = [cos(phi);-sin(phi);0];
    dd_e_phi = [-sin(phi);-cos(phi);0];
    ddd_e_phi = [-cos(phi);sin(phi);0];
    
    %% Delta
    delta = ppval(s_delta,mod(varphi,2*pi));
    d_delta = ppval(s_d_delta,mod(varphi,2*pi));
    dd_delta = ppval(s_dd_delta,mod(varphi,2*pi));
    ddd_delta = ppval(s_ddd_delta,mod(varphi,2*pi));
    %[delta,d_delta,dd_delta,ddd_delta] = delta_5_cosine(varphi);
    %delta = (a-b*cos(2*phi));
    delta_v = delta*e_phi;
%     d_delta = b*2*sin(2*phi);
%     dd_delta = b*4*cos(2*phi);
%     ddd_delta = -b*8*sin(2*phi);
    d_delta_v = d_delta*e_phi+delta*d_e_phi;
    dd_delta_v = dd_delta*e_phi+2*d_delta*d_e_phi+delta*dd_e_phi;
    ddd_delta_v = ddd_delta*e_phi+ 3*dd_delta*d_e_phi + ...
        3*d_delta*dd_e_phi + delta*ddd_e_phi;
    norm_d_delta_v = sqrt(d_delta^2+delta^2);
    
    %% Tau
    tau = d_delta_v/norm_d_delta_v;
    d_tau = dd_delta_v/norm_d_delta_v - d_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3;
    dd_tau = ddd_delta_v/norm_d_delta_v - dd_delta_v*d_delta_v'*dd_delta_v/norm_d_delta_v^3 - ...
        (dd_delta_v*d_delta_v'*dd_delta_v+d_delta_v*(dd_delta_v'*dd_delta_v)+d_delta_v*d_delta_v'*ddd_delta_v)/...
        norm_d_delta_v^3 + d_delta_v*d_delta_v'*dd_delta_v*3*d_delta_v'*dd_delta_v/norm_d_delta_v^5;
    %% Normal vector:
    normal = [0 -1 0;1 0 0;0 0 0]*tau;
    d_normal = [0 -1 0;1 0 0;0 0 0]*d_tau;
    dd_normal = [0 -1 0;1 0 0;0 0 0]*dd_tau;
    

    %% Rho :
    rho = delta_v + R_b*normal;
    d_rho = d_delta_v + R_b*d_normal;
    dd_rho = dd_delta_v + R_b*dd_normal;
    l_rho = norm(rho);
    
    %% Need to find varphi = g(phi)
    gamma = delta_v'*[1;0;0]-R_b*tau'*[0;1;0];
    psi = delta_v'*[0;1;0]+R_b*tau'*[1;0;0];
    d_gamma = d_delta_v'*[1;0;0]-R_b*d_tau'*[0;1;0];
    d_psi = d_delta_v'*[0;1;0]+R_b*d_tau'*[1;0;0];
    dd_gamma = dd_delta_v'*[1;0;0]-R_b*dd_tau'*[0;1;0];
    dd_psi = dd_delta_v'*[0;1;0]+R_b*dd_tau'*[1;0;0];
    %func_g = atan2(gamma, psi);
    d_func_g = (d_gamma*psi-d_psi*gamma)/(gamma^2+psi^2);
    dd_func_g = (dd_gamma*psi-dd_psi*gamma)/(psi^2+gamma^2)-(d_gamma*psi-d_psi*gamma)*(2*gamma*d_gamma+2*psi*d_psi)/(gamma^2+psi^2)^2;
    
    %% s
    d_s = norm(d_rho)/d_func_g;
    dd_s = d_rho'*dd_rho/norm(d_rho)/d_func_g^2 - norm(d_rho)*dd_func_g/d_func_g^3;
    
    %% s_f:
    d_sf = norm_d_delta_v/d_func_g;
    dd_sf = d_delta_v'*dd_delta_v/norm_d_delta_v/d_func_g - norm_d_delta_v*dd_func_g/d_func_g^3;
    
    %% Kappa:
    kappa = dd_rho/d_func_g^2/d_s^2;
%     
%     %% Rotations:
%     R = [cos(varphi) -sin(varphi) 0;sin(varphi) cos(varphi) 0;0 0 1];
%     dR = [-sin(varphi) -cos(varphi) 0;cos(varphi) -sin(varphi) 0;0 0 0];
%     ddR = [-cos(varphi) sin(varphi) 0;-sin(varphi) -cos(varphi) 0;0 0 0];
%     %% Theta = varphi + atan(gamma/psi)
%     c = -0.03;
%     
%     gamma = c*sin(2*varphi-pi)*[1 0 0]*R*tau-[0 1 0]*R*tau;
%     psi = c*sin(2*varphi-pi)*[0 1 0]*R*tau+[1 0 0]*R*tau;
%     d_gamma = c*2*cos(2*varphi-pi)*[1 0 0]*R*tau + c*sin(2*varphi-pi)*[1 0 0]*dR*tau + c*sin(2*varphi-pi)*[1 0 0]*R*d_tau - [0 1 0]*dR*tau-[0 1 0]*R*d_tau;
%     d_psi = c*2*cos(2*varphi-pi)*[0 1 0]*R*tau + c*sin(2*varphi-pi)*[0 1 0]*dR*tau + c*sin(2*varphi-pi)*[0 1 0]*R*d_tau + [1 0 0]*dR*tau + [1 0 0]*R*d_tau;
%     
%     dd_gamma = -c*4*sin(2*varphi-pi)*[1 0 0]*R*tau + c*2*cos(2*varphi-pi)*[1 0 0]*dR*tau + c*2*cos(2*varphi-pi)*[1 0 0]*R*d_tau ... 
%         +c*2*cos(2*varphi-pi)*[1 0 0]*dR*tau + c*sin(2*varphi-pi)*[1 0 0]*ddR*tau + c*sin(2*varphi-pi)*[1 0 0]*dR*d_tau ...
%         + c*2*cos(2*varphi-pi)*[1 0 0]*R*d_tau + c*sin(2*varphi-pi)*[1 0 0]*dR*d_tau + c*sin(2*varphi-pi)*[1 0 0]*R*dd_tau ...
%         - [0 1 0]*ddR*tau - 2*[0 1 0]*dR*d_tau -[0 1 0]*R*dd_tau;
%     
%     dd_psi = -c*4*sin(2*varphi-pi)*[0 1 0]*R*tau + c*2*cos(2*varphi-pi)*[0 1 0]*dR*tau + c*2*cos(2*varphi-pi)*[0 1 0]*R*d_tau ...
%         + c*2*cos(2*varphi-pi)*[0 1 0]*dR*tau + c*sin(2*varphi-pi)*[0 1 0]*ddR*tau + c*sin(2*varphi-pi)*[0 1 0]*dR*d_tau ...
%         + c*2*cos(2*varphi-pi)*[0 1 0]*R*d_tau + c*sin(2*varphi-pi)*[0 1 0]*dR*d_tau + c*sin(2*varphi-pi)*[0 1 0]*R*dd_tau ...
%         + [1 0 0]*ddR*tau + 2*[1 0 0]*dR*d_tau + [1 0 0]*R*dd_tau;
%     Theta = varphi + atan(gamma/psi);
%     d_Theta = (-gamma*d_psi + psi^2 + d_gamma*psi + gamma^2)*1/(psi^2 + gamma^2);
%     dd_Theta = ((psi^3 + psi*gamma^2)*dd_gamma + (-psi^2*gamma - gamma^3)*dd_psi - 2*(psi*d_psi + gamma*d_gamma)*(d_gamma*psi - gamma*d_psi))*1/(psi^2 + gamma^2)^2;
%     

    
    rhoxtau = cross(rho,tau);
    rhoxkappa = cross(rho,kappa);
    
     
    %% Rotations:
    R = [cos(Theta) -sin(Theta) 0;sin(Theta) cos(Theta) 0;0 0 1];
    dR = [-sin(Theta) -cos(Theta) 0;cos(Theta) -sin(Theta) 0;0 0 0];
    
    
    

    alpha = m_b*((d_s*rhoxtau(3)-d_sf*J_s/R_b/m_b)*d_Theta+d_s^2+d_sf^2*J_s/R_b/m_b/R_b);
    beta = m_b*((d_s*rhoxtau(3)-d_sf*J_s/R_b/m_b)*dd_Theta-rho'*tau*d_s*d_Theta^2 ...
         + d_s*dd_s + d_sf*dd_sf*J_s/R_b/m_b/R_b);
    gamma = m_b*[0 g 0]*R*tau*d_s;
    
    g_w = -m_b*(d_s*rhoxtau(3)-d_sf*J_s/m_b/R_b);
%     g_y_dot = m_b*d_s*tau'*rho*(2*d_Theta*d_varphi);
%     g_y = -get_ds(varphi)*m_b*g*[cos(get_theta(varphi)) -sin(get_theta(varphi)) 0]*get_tau(varphi);
%     
    %B = [0;1;g_w/alpha];%0];
    B = [0;1;g_w/alpha]/d_phi_star;
    Gamma = 1;
    %integral_term = integral_term + d_phi_star*0.001;
    %epsilon = [q(1)-Theta;dq(1)-d_Theta*dq(2);dq(2)^2-ppval(spline_I,mod(q(2),pi))];%q(2)-integral_term];

    
    m_11 = m_b*(rho'*rho+J_s/m_b+J_f/m_b);
    m_12 = m_b*(d_s*rhoxtau(3)-d_sf*J_s/R_b/m_b);
    m_22 = m_b*(d_s^2+d_sf^2*J_s/R_b^2/m_b);
    M = [m_11 m_12;m_12 m_22];
    
    c_11 = m_b*rho'*tau*d_s*d_phi_star;
    c_12 = m_b*(rho'*tau*d_s*d_Theta*(d_phi_star)+(dd_s*rhoxtau(3)-dd_sf*J_s/R_b/m_b+d_s^2*rhoxkappa(3))*d_phi_star);
    c_21 = -m_b*rho'*tau*d_s*d_Theta*(d_phi_star);
    c_22 = m_b*(d_s*dd_s+d_sf*dd_sf*J_s/R_b^2/m_b)*d_phi_star;
    C = [c_11 c_12;c_21 c_22];
    
    G = m_b*[[0 g 0]*dR*rho;[0 g 0]*R*tau*d_s];
    
    epsilon_measured = [q(1)-Theta;dq(1)-d_Theta*dq(2);dq(2)-d_phi_star];

    if isempty(epsilon)
        epsilon = epsilon_measured;
    end
    if with_luenberger
        w = -(1/Gamma)*B'*sol_riccati*epsilon;
    else
        w = -(1/Gamma)*B'*sol_riccati*epsilon_measured;
    end
    k = [1 0]*M*[d_Theta;1]/([d_Theta 1]*M*[0;1]);
    u = k*([k^-1 -1]*M*[w+dd_Theta*d_phi_star^2;0]+[k^-1 -1]*C*[d_Theta*d_phi_star+epsilon(2);d_phi_star]+[k^-1 -1]*G);
 
    g_w = -m_b*(d_s*rhoxtau(3)-d_sf*J_s/m_b/R_b);
    g_y_dot = m_b*d_s*tau'*rho*(2*d_Theta*d_phi_star);
    g_y = -d_s*m_b*g*[cos(Theta) -sin(Theta) 0]*tau;
   
    %A = [0 1 0;0 0 0;g_y/alpha g_y_dot/alpha (-beta*d_phi_star^2+gamma)/alpha]/d_phi_star;
    A = [0 1 0;
         0 0 0;
         g_y/alpha g_y_dot/alpha (gamma-beta*d_phi_star^2)/alpha/d_phi_star];
    B = [0;1;g_w/alpha];

    %B = [0;1;g_w/alpha]/d_phi_star;%0];
    C = [1 0 0;0 1 0;0 0 1];
    L = [luenberger(1) 0 0;0 luenberger(2) 0;0 0 luenberger(3)];
    %kalman = ppval(spline_kalman,mod(phi,pi));
    %K = reshape(kalman,3,3);
    %K(1,3) = 0; K(3,1) = 0;
    %R = [1e-4  0 0;0 1e-4 0;0 0 1e-5];
    %L = R\K*C';
    phi_prev = q(2);
    dt = t-prev_time;
    prev_time = t;
    epsilon = epsilon + (A*epsilon+B*u+L*(C*epsilon_measured-C*epsilon))*dt; 
    current_epsilon = epsilon;
    if dq(1) > 0
        u = u + 0.002;
    elseif dq(1) <= 0
        u = u - 0.003;
    end     

        