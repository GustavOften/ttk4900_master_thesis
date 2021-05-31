syms Jf Js m_b rho rhoxtau rhoxn rhoxkappa d_s dd_s d_sf dd_sf norm_k tau n g R R_dot R_b real
syms u q1 q2 q3 q4 dq1 dq2 dq3 dq4 real
syms ddq1 ddq2 ddq3 ddq4 real
M = [Jf+Js+m_b*rho*rho m_b*d_s*rhoxtau Js m_b*rhoxn;
    m_b*d_s*rhoxtau m_b*d_s^2 0 0;
    Js 0 Js 0;
    m_b*rhoxn 0 0 m_b];
C = m_b*[rho*tau*d_s*dq1+rho*n*dq4 rho*tau*d_s*dq2+(dd_s*rhoxtau+rhoxkappa)*dq2+1/2*(d_s-norm_k*rhoxtau*d_s)*dq4 0 rho*n*dq1+1/2*(d_s-norm_k*rhoxtau*d_s)*dq2;
    rho*tau*d_s*dq1 dd_s*d_s*dq2 0 0;
    0 0 0 0;
    -rho*tau*d_s*dq1+1/2*(d_s-norm_k*rhoxtau*d_s)*dq2 1/2*(d_s-norm_k*rhoxtau*d_s)*dq1 0 0];
G =[R_dot*g*m_b*rho;
  R*d_s*g*m_b*tau;
             0;
       R*g*m_b*n];
 
B = [1;0;0;0];
f = [0 0;0 d_sf;0 R_b;1 0];
f_L = inv(f'*f)*f';
dq = [dq1;dq2;dq3;dq4];
ddq = [ddq1;ddq2;ddq3;ddq4];
dynamics = M*ddq+C*dq+G-B*u;
lambda = f_L*dynamics;
lambda_with_constraints = simplify(subs(lambda,[dq3 ddq3 dq4 ddq4],[-d_sf/R_b*dq2 -dd_sf/R_b*dq2^2-d_sf/R_b*ddq2 0 0]));


% lambda_1 = m*(ddq1*rhoxn + R*g*n + dq1*dq2*d_s + dq1*dq2*rhoxkappa - dq1^2*d_s*rho*tau - dq1*dq2*d_s*rhoxn);
% lambda_2 = ((Js*R_b + d_s*d_sf*m*rhoxtau)/(R_b^2 + d_sf^2))*ddq1 + (-(- d_sf*m*d_s^2 + Js*d_sf)/(R_b^2 + d_sf^2))*ddq2 + (d_sf*(d_s*m*rho*tau*dq1^2 + dd_s*d_s*m*dq2^2 + R*d_s*g*m*tau) - Js*dd_sf*dq2^2)/(R_b^2 + d_sf^2);
