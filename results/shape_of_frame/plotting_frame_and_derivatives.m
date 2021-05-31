close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

fname = 'shape.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
x = zeros(length(val.shape)/2,1);
y = zeros(length(val.shape)/2,1);
j = 1;
for i = 1:length(val.shape)
    if mod(i,2) == 0
        y(j) = val.shape(i);
    else 
        x(j) = val.shape(i); 
        j = j + 1;
    end
end
a = 0.1148; %Constant for describing the butterfly frame
b = -0.04022;
func = @(x) (a+b*cos(2*x));
d_func = @(x) -b*2*sin(2*x);
dd_func = @(x) -b*4*cos(2*x);
ddd_func = @(x) b*8*cos(2*x);
a0 =      0.1148;  
a1 =    -0.04024;
b1 =    0.000155;
a2 =   -0.001519;
b2 =  -4.007e-05;
a3 =    0.001161;
b3 =  -6.491e-06;
a4 =   0.0008238;
b4 =    -1.5e-06;
w =       2;

a =      0.1148;
b =    -0.04023;
c =   -0.001535;
d =    0.001152;
e =   0.0008215;

delta = @(x) a+b*cos(2*x)+c*cos(4*x)+d*cos(6*x)+e*cos(8*x);
ddelta = @(x) -b*2*sin(2*x)-c*4*sin(4*x)-d*6*sin(6*x)-e*8*sin(8*x);
dddelta = @(x) -b*4*cos(2*x)-c*16*cos(4*x)-d*36*cos(6*x)-e*64*cos(8*x);
ddddelta = @(x) b*8*sin(2*x)+c*16*4*sin(4*x)+d*36*6*sin(6*x)+e*64*8*sin(8*x);
fourier = @(x) (a0+a1*cos(x*w) + b1*sin(x*w) + ...
               a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ... 
               a4*cos(4*x*w) + b4*sin(4*x*w));
a0 =      0.1148;
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

len = zeros(length(val.shape)/2-2,1);
for i =1:length(val.shape)/2-2
    len(i) = sqrt(x(i+1)^2+y(i+1)^2);
end
X = linspace(0,2*pi,length(len));
x = x(2:end-1);
y = y(2:end-2);
phi = unwrap(atan2(x,y));
figure
k = linspace(0,2*pi,length(x));
          
f = @(x) (a0 + a1*cos(x*w) + b1*sin(x*w) + ...
         a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ... 
         a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
         a6*cos(6*x*w) + b6*sin(6*x*w));

df = @(x) (-w*a1*sin(x*w) + w*b1*cos(x*w) + ...
         -2*w*a2*sin(2*x*w) + 2*b2*w*cos(2*x*w) - 3*w*a3*sin(3*x*w) + 3*w*b3*cos(3*x*w) + ... 
         -4*w*a4*sin(4*x*w) + 4*w*b4*cos(4*x*w) -5*w*a5*sin(5*x*w) + 5*w*b5*cos(5*x*w) + ...
         -6*w*a6*sin(6*x*w) + 6*w*b6*cos(6*x*w));
ddf = @(x) (-w^2*a1*cos(x*w) - w^2*b1*sin(x*w) + ...
         -4*w^2*a2*cos(2*x*w) - w^2*4*b2*sin(2*x*w) - 9*w^2*a3*cos(3*x*w) - 9*w^2*b3*sin(3*x*w) + ... 
         -16*w^2*a4*cos(4*x*w) - 16*w^2*b4*sin(4*x*w) -25*w^2*a5*cos(5*x*w) - 25*w^2*b5*sin(5*x*w) + ...
         -36*w^2*a6*cos(6*x*w) - 36*w^2*b6*sin(6*x*w));
dddf = @(x) (w^3*a1*sin(x*w) - w^3*b1*cos(x*w) + ...
         4*2*w^3*a2*sin(2*x*w) - w^3*4*2*b2*cos(2*x*w) + 9*3*w^3*a3*sin(3*x*w) - 9*3*w^3*b3*cos(3*x*w) + ... 
         16*4*w^3*a4*sin(4*x*w) - 16*4*w^3*b4*cos(4*x*w) + 25*5*w^3*a5*sin(5*x*w) - 25*5*w^3*b5*cos(5*x*w) + ...
         36*6*w^3*a6*sin(6*x*w) - 36*6*w^3*b6*cos(6*x*w));
           
figure
plot(linspace(0,2*pi,length(x)),len);
hold on;
k = func(linspace(0,2*pi,4000));
g = f(linspace(0,2*pi,4000));
dlen = zeros(length(len)-1,1);
X_dlen = zeros(length(len)-1,1);
for i = 2:length(len)
    dlen(i-1) = (len(i)-len(i-1))/(phi(i)-phi(i-1));
    X_dlen(i-1) = (phi(i)+phi(i-1))/2;
end
spline_delta = spline(X,len);
spline_d_delta = spline(X_dlen,dlen);

figure 
hold on
k = d_func(linspace(0,2*pi,4000));
g = df(linspace(0,2*pi,4000));
plot(linspace(0,2*pi,4000),k(1,:));
plot(linspace(0,2*pi,4000),g(1,:));
plot(linspace(0,2*pi,length(dlen)),dlen);
plot(linspace(0,2*pi,4000),ddelta(linspace(0,2*pi,4000)));
legend('cos','fourier','numerical derivation')

ddlen = zeros(length(len)-2,1);

for i = 2:length(len)-1
    ddlen(i-1) = (len(i+1)-2*len(i)+len(i-1))/((phi(i)-phi(i-1))*(phi(i+1)-phi(i)));
    
end
ddlen = zeros(length(dlen)-1,1);
X_ddlen = zeros(length(dlen)-2,1);
for i = 2:length(dlen)
    ddlen(i-1) = (dlen(i)-dlen(i-1))/(phi(i)-phi(i-1));
    X_ddlen(i-1) = (phi(i)+phi(i-1))/2;
end
spline_dd_delta = spline(X_ddlen,ddlen);
figure 
hold on
k = dd_func(linspace(0,2*pi,4000));
g = ddf(linspace(0,2*pi,4000));
plot(linspace(0,2*pi,4000),k(1,:));
plot(linspace(0,2*pi,4000),g(1,:));
plot(linspace(0,2*pi,length(ddlen)),ddlen);
plot(linspace(0,2*pi,4000),dddelta(linspace(0,2*pi,4000)));
plot(linspace(0,2*pi,4000),ppval(spline_dd_delta,linspace(0,2*pi,4000)))
legend('cos','fourier','numerical_derivation','cos sum','spline');

dddlen = zeros(length(ddlen)-1,1);
X_dddlen = zeros(length(ddlen)-2,1);
for i = 2:length(ddlen)
    dddlen(i-1) = (ddlen(i)-ddlen(i-1))/(phi(i)-phi(i-1));
    X_dddlen(i-1) = (phi(i)+phi(i-1))/2;
end
spline_ddd_delta = spline(X_dddlen,dddlen);


figure 
hold on
k = ddd_func(linspace(0,2*pi,4000));
g = dddf(linspace(0,2*pi,4000));


plot(linspace(0,2*pi,4000),k(1,:));
plot(linspace(0,2*pi,4000),g(1,:));
plot(linspace(0,2*pi,length(dddlen)),dddlen);
plot(linspace(0,2*pi,4000),ddddelta(linspace(0,2*pi,4000)));
plot(linspace(0,2*pi,4000),ppval(spline_ddd_delta,linspace(0,2*pi,4000)))
legend('cos','fourier','numerical_derivation')

save('C:\Users\g-oft\OneDrive\Dokumenter\master\ttk4900_master_thesis\shape.mat','spline_delta','spline_d_delta','spline_dd_delta','spline_ddd_delta');

figure
subplot(4,1,1)
hold on;
grid on;
plot(phi,len)
plot(phi,func(phi))
plot(phi,delta(phi))
set(gca,'XTick',0:pi/4:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
legend('Frame','1 cosine','4 cosine')
xlim([0,2*pi])
%xlabel('$\phi$','fontsize',14,'interpreter','latex')
ylabel('$\delta(\phi)$','fontsize',10,'interpreter','latex')
%ylabel('','fontsize',14,'interpreter','latex')


subplot(4,1,2)
hold on;
grid on;
plot(X_dlen,dlen)
plot(X_dlen,d_func(X_dlen))
plot(X_dlen,ddelta(X_dlen))
set(gca,'XTick',0:pi/4:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
legend('Finite difference','1 cosine','4 cosine')
xlim([0,2*pi])
%xlabel('$\phi$','fontsize',14,'interpreter','latex')
%ylabel('$\frac{\partial \delta(\phi)}{\partial \phi}$','fontsize',14,'interpreter','latex')
ylabel('First derivative','fontsize',10,'interpreter','latex')


subplot(4,1,3)
hold on;
grid on;
plot(X_ddlen,ddlen)
plot(X_ddlen,dd_func(X_ddlen))
plot(X_ddlen,dddelta(X_ddlen))
set(gca,'XTick',0:pi/4:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
legend('Finite difference','1 cosine','4 cosine')
xlim([0,2*pi])
%xlabel('$\phi$','fontsize',14,'interpreter','latex')
%ylabel('$\frac{\partial^2 \delta(\phi)}{\partial \phi^2}$','fontsize',14,'interpreter','latex')
ylabel('Second derivative','fontsize',10,'interpreter','latex')

subplot(4,1,4)
hold on;
grid on;
plot(X_dddlen,dddlen)
plot(X_dddlen,ddd_func(X_dddlen))
plot(X_dddlen,ddddelta(X_dddlen))
set(gca,'XTick',0:pi/4:2*pi) 
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
legend('Finite difference','1 cosine','4 cosine')
xlabel('$\phi$','fontsize',14,'interpreter','latex')
%ylabel('$\frac{\partial^3 \delta(\phi)}{\partial \phi^3}$','fontsize',14,'interpreter','latex')
ylabel('Third derivative','fontsize',10,'interpreter','latex')
xlim([0,2*pi])
sgtitle('Functions for frame and their derivatives', 'Interpreter','latex');




