a = 0.1095; %Constant for describing the butterfly frame
b = 0.0405; %Constant for describing the butterfly frame
g = 9.81; %Gravity
J_s = 5.48e-7; %Mass moment of inertia of the ball
J_f = 1.581e-2; %Mass moment of inertia of the frame, from article.
m_b = 3.0e-3; %Mass of the ball
R_b = sqrt(16.55e-3^2-12.5e-3^2); %Radius of the ball
%% phi = f(varphi), varphi = g(varphi) 
delta_curve_fit = @(phi)(a - b*cos(2*phi))*[sin(phi);cos(phi);0];
tau_curve_fit = @(phi)(2*b*sin(2*phi)*[sin(phi);cos(phi);0]+(a - b*cos(2*phi))*[cos(phi);-sin(phi);0]) ...
    /sqrt(sum((2*b*sin(2*phi)*[sin(phi);cos(phi);0]+(a - b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));

range_for_functions = linspace(0,2*pi,500);
function_g = @(x) atan2([1 0 0]*delta_curve_fit(x) - R_b*[0 1 0]*tau_curve_fit(x), ...
                    [0 1 0]*delta_curve_fit(x)+R_b*[1 0 0]*tau_curve_fit(x));
res_fun_g = zeros(length(range_for_functions),1);
res_fun_dg = zeros(length(range_for_functions),1);
for i = 1:length(range_for_functions)
    res_fun_g(i) = function_g(range_for_functions(i));
end
res_fun_g = unwrap(res_fun_g);

k = spline(res_fun_g,range_for_functions);
result_spline = @(x)ppval(k,x);   


n_varphi = 5000;
n_riccati = 5000;
range_varphi_to_phi = linspace(0,2*pi,n_varphi);
range_sol_riccati = linspace(0,pi,n_riccati);
d_phi_star = zeros(1,n_riccati);
riccati_sol = zeros(3,3,n_riccati);
varphi_to_phi = zeros(1,n_varphi);
for i = 1:n_varphi
    varphi_to_phi(1,i) = result_spline(range_varphi_to_phi(i));
end
for i = 1:n_riccati
    d_phi_star(1,i) = bf.function_for_dphi(range_sol_riccati(i));
    riccati_sol(:,:,i) = bf.function_for_X(range_sol_riccati(i));
end

save('solutions_riccati_and_dphi.mat','d_phi_star','riccati_sol','varphi_to_phi');

