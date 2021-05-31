delete C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\solutions_riccati_and_dphi.mat;
%spline_dphi = load('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\shape.mat');
shape = load('shape.mat');

% 
a = 0.1095; %Constant for describing the butterfly frame
b = 0.0405; %Constant for describing the butterfly frame
a0 =    0.1148;
a1 =    -0.04023;
b1 =   0.0001534;
a2 =   -0.001519;
b2 =  -4.019e-05;
a3 =    0.001161;
b3 =  -6.347e-06;
a4 =   0.0008238;
b4 =  -1.364e-06;
a5 =   0.0002184;
b5 =   6.756e-07;
a6 =  -0.0001084;
b6 =   1.232e-06;
w =       1.999;
% a = 0.1148; %Constant for describing the butterfly frame, measured
% b = 0.0390; %Constant for describing the butterfly frame, measured
% a = 0.1150;
% b = 0.04;
a = 0.1148; %Constant for describing the butterfly frame
b = 0.04022;
a =      0.1148;
b =    -0.04023;
c =   -0.001535;
d =    0.001152;
e =   0.0008215;
g = 9.81; %Gravity
%J_s = 4.44556451612904e-07;%5.48e-7;%5.48e-7; %Mass moment of inertia of the ball 
%J_f = 1.58117e-3; %Mass moment of inertia of the frame, from article.
J_f = 8.8015e-04;
m_b = 3.0e-3; %Mass of the ball
%R_b = sqrt(16.55e-3^2-12.5e-3^2); %From paper.
r_b = 16.6e-3;
R_b = sqrt(r_b^2-13.5e-3^2); %Measured on the robot
J_s = 2/3*m_b*r_b^2;
n_riccati = 400;
%% phi = f(varphi), varphi = g(varphi) 
%delta_curve_fit = @(phi)(a - b*cos(2*phi))*[sin(phi);cos(phi);0];
%tau_curve_fit = @(phi)(2*b*sin(2*phi)*[sin(phi);cos(phi);0]+(a - b*cos(2*phi))*[cos(phi);-sin(phi);0]) ...
%    /sqrt(sum((2*b*sin(2*phi)*[sin(phi);cos(phi);0]+(a - b*cos(2*phi))*[cos(phi);-sin(phi);0]).^2));

e_phi =@(phi) [sin(phi);cos(phi);0];
d_e_phi =@(phi) [cos(phi);-sin(phi);0];
dd_e_phi = @(phi) [-sin(phi);-cos(phi);0];

%% Delta
% delta = @(x) (a0 + a1*cos(x*w) + b1*sin(x*w) + ...
%          a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ... 
%          a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
%          a6*cos(6*x*w) + b6*sin(6*x*w));
% 
% d_delta = @(x) (-w*a1*sin(x*w) + w*b1*cos(x*w) + ...
%          -2*w*a2*sin(2*x*w) + 2*b2*w*cos(2*x*w) - 3*w*a3*sin(3*x*w) + 3*w*b3*cos(3*x*w) + ... 
%          -4*w*a4*sin(4*x*w) + 4*w*b4*cos(4*x*w) -5*w*a5*sin(5*x*w) + 5*w*b5*cos(5*x*w) + ...
%          -6*w*a6*sin(6*x*w) + 6*w*b6*cos(6*x*w));
delta =  @(phi)a+b*cos(2*phi)+c*cos(4*phi)+d*cos(6*phi)+e*cos(8*phi);
d_delta =  @(phi)-b*2*sin(2*phi)-c*4*sin(4*phi)-d*6*sin(6*phi)-e*8*sin(8*phi);


[breaks,coefs,~,~,dim] = unmkpp(shape.spline_delta);
spline_delta = mkpp(breaks,coefs,dim);
[breaks,coefs,~,~,dim] = unmkpp(shape.spline_d_delta);
spline_d_delta = mkpp(breaks,coefs,dim);
[breaks,coefs,~,~,dim] = unmkpp(shape.spline_dd_delta);
spline_dd_delta = mkpp(breaks,coefs,dim);
[breaks,coefs,~,~,dim] = unmkpp(shape.spline_ddd_delta);
spline_ddd_delta = mkpp(breaks,coefs,dim);
% spline_delta = spline(X,len);
% spline_d_delta = spline(X_dlen,dlen);
% spline_dd_delta = spline(X_ddlen,ddlen);
% spline_ddd_delta = spline(X_dddlen,dddlen);

delta = @(x) ppval(spline_delta,x);
d_delta = @(x) ppval(spline_d_delta,x);
a =      0.1148;
b =    -0.04023;
c =   -0.001535;
d =    0.001152;
e =   0.0008215;
% delta =  @(phi) a+b*cos(2*phi)+c*cos(4*phi)+d*cos(6*phi)+e*cos(8*phi);
% d_delta =  @(phi) -b*2*sin(2*phi)-c*4*sin(4*phi)-d*6*sin(6*phi)-e*8*sin(8*phi);

%delta =@(phi) (a-b*cos(2*phi));
delta_v = @(phi) delta(phi)*e_phi(phi);
%d_delta = @(phi) b*2*sin(2*phi);
d_delta_v = @(phi) d_delta(phi)*e_phi(phi)+delta(phi)*d_e_phi(phi);
norm_d_delta_v = @(phi) norm(d_delta_v(phi));
dd_delta = @(phi) b*4*cos(2*phi);
curvature = @(phi) ((dd_delta(phi)*e_phi(phi) + 2*d_delta(phi)*d_e_phi(phi)+delta(phi)*dd_e_phi(phi))/norm(d_delta_v(phi))^2 - ...
    d_delta_v(phi)*d_delta_v(phi)'*(dd_delta(phi)*e_phi(phi) + 2*d_delta(phi)*d_e_phi(phi)+delta(phi)*dd_e_phi(phi))/norm(d_delta_v(phi))^4)/norm((dd_delta(phi)*e_phi(phi) + 2*d_delta(phi)*d_e_phi(phi)+delta(phi)*dd_e_phi(phi))/norm(d_delta_v(phi))^2 - ...
    d_delta_v(phi)*d_delta_v(phi)'*(dd_delta(phi)*e_phi(phi) + 2*d_delta(phi)*d_e_phi(phi)+delta(phi)*dd_e_phi(phi))/norm(d_delta_v(phi))^4);
s_curvature = @(phi) sum(curvature(phi));


%% Tau
tau =@(phi) d_delta_v(phi)/norm_d_delta_v(phi);
range_for_functions = linspace(0,2*pi,500);
function_g_one = @(x) atan2([1 0 0]*delta_v(x)-sign([0 1 0]*tau(x)*sin(x))*R_b*[0 1 0]*tau(x),...
                    [0 1 0]*delta_v(x)+R_b*[1 0 0]*tau(x));
function_g_two = @(x) atan2([1 0 0]*delta_v(x)+R_b*[0 1 0]*tau(x),...
                    [0 1 0]*delta_v(x)+R_b*[1 0 0]*tau(x));
function_g = @(x) atan2([1 0 0]*delta_v(x)-R_b*[0 1 0]*tau(x),...
                    [0 1 0]*delta_v(x)+R_b*[1 0 0]*tau(x));
res_fun_g = zeros(length(range_for_functions),1);
res_fun_g_one = zeros(length(range_for_functions),1);
res_fun_dg = zeros(length(range_for_functions),1);
res_fun_g_two = zeros(length(range_for_functions),1);
sign_of_kappa = zeros(length(range_for_functions),1);
tau_frame = zeros(length(range_for_functions),3);
normal_frame = zeros(length(range_for_functions),3);
kapp = zeros(length(range_for_functions),3);
normal =@(x) [-[0 1 0]*tau(x);[1 0 0]*tau(x);0];
dot_product_kappa_normal = zeros(length(range_for_functions),1);
func =@(x) sign([0 1 0]*tau(x)*sin(x));
test = zeros(length(range_for_functions),3);
for i = 1:length(range_for_functions)
    res_fun_g(i) = function_g(range_for_functions(i));
    res_fun_g_one(i) = function_g_one(range_for_functions(i));
    res_fun_g_two(i) = function_g_two(range_for_functions(i));
    sign_of_kappa(i) = sign(s_curvature(range_for_functions(i)));
    normal_frame(i,:) = normal(range_for_functions(i));
    tau_frame(i,:) = tau(range_for_functions(i)); 
    kapp(i,:) = curvature(range_for_functions(i));
    dot_product_kappa_normal(i) = kapp(i,:)*normal_frame(i,:)';
    test(i) = func(range_for_functions(i));
end
res_fun_g = unwrap(res_fun_g);
res_fun_g_two = unwrap(res_fun_g_two);
res_fun_g_one = unwrap(res_fun_g_one);
figure
plot(range_for_functions,res_fun_g);
hold on;
%plot(range_for_functions,res_fun_g_one);
%plot(range_for_functions,res_fun_g_two);
plot(range_for_functions,res_fun_g_one);
plot(range_for_functions,range_for_functions);
legend('g(\phi)','f(varphi','line')
varphi_to_phi = spline(res_fun_g,range_for_functions);
%[breaks,coefs,L,order,dim] = unmkpp(varphi_to_phi);
%varphi_to_phi = mkpp(breaks,coefs,dim);
%%result_spline = @(x)ppval(k,x);   
%range_varphi_to_phi = linspace(0,2*pi,n_varphi);
range_sol_riccati = linspace(0,pi,n_riccati);
%d_phi_star = zeros(1,n_riccati);
riccati_sol = zeros(3,3,n_riccati);
% varphi_to_phi = zeros(1,n_varphi);
% for i = 1:n_varphi
%     varphi_to_phi(1,i) = result_spline(range_varphi_to_phi(i));
% end
delta_frame = zeros(2,500);
rho_frame = zeros(2,500);
rho = zeros(1,500);
for i = 1:length(range_for_functions)
    delta_frame(:,i) = [sin(range_for_functions(i))*delta(range_for_functions(i));
                        cos(range_for_functions(i))*delta(range_for_functions(i))];
    rho_frame(:,i) = [sin(range_for_functions(i))*delta(range_for_functions(i))-R_b*[0 1 0]*tau(range_for_functions(i));
                        cos(range_for_functions(i))*delta(range_for_functions(i))+R_b*[1 0 0]*tau(range_for_functions(i))];
end
figure
hold on;
xlim([-0.2 0.2])
ylim([-0.2 0.2])
plot(delta_frame(1,:),delta_frame(2,:))
plot(rho_frame(1,:),rho_frame(2,:))
%plot(out.simout2.Data(:,1),out.simout2.Data(:,2))

%% Finding interpolation for the transverse coordinate I
abg =@(x) bf.alpha_beta_gamma(x);
%Integral beta/alpha:

f = @(t,x) ([0 1 0]*abg(t))/([1 0 0]*abg(t));
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
ba_integral = ode45(f,[0,pi],0,options);
ba_spline = spline(ba_integral.x,ba_integral.y);
psi_minus = @(phi) exp(-2*ppval(ba_spline,phi));
psi_pluss = @(phi) exp(2*ppval(ba_spline,phi));
f = @(t,x) psi_pluss(t)*2*[0 0 1]*abg(t)/([1 0 0]*abg(t));
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
big_integral = ode45(f,[0, pi],0,options);
big_spline = spline(big_integral.x,big_integral.y);
I = @(phi) psi_minus(phi)*(ppval(bf.spline_dphi,0)^2-ppval(big_spline,phi));

k = linspace(0,pi,500);
res = zeros(1,500);
for i = 1:500
    res(1,i) = I(k(i));
end
spline_I = spline(k,res);
spline_dphi = bf.spline_dphi;
spline_kalman = spline(range_sol_riccati,reshape(bf.P,9,n_riccati));
spline_riccati_sol = spline(range_sol_riccati,reshape(bf.X,9,n_riccati));
save('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\solutions_riccati_and_dphi.mat','spline_dphi','spline_riccati_sol','varphi_to_phi','spline_delta','spline_d_delta','spline_dd_delta','spline_ddd_delta','spline_kalman','spline_I');




