close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
epsilon = load('final_luenberger_no_noise_epsilon.txt');
states = load('final_luenberger_no_noise.txt');
eps_measured = load('final_luenberger_no_noise_eps_measured.txt');

phi = mod(states(:,4),2*pi);
jumps = [];
k = 1;
for i = 1:length(phi)-1
    if phi(i+1)-phi(i) < -3
        jumps(k) = i;
        k = k+1;
    end
end

%% varphi
figure
grid on;
hold on;
plot(states(:,1),states(:,4))
xlabel('Time[s]','fontsize',14,'interpreter','latex')
xlim([5,30])
ylabel('$\varphi$[rad]','fontsize',14,'interpreter','latex')
sgtitle('$\varphi$ as a function of time')

%% Theta and vartheta
figure
grid on;
hold on;
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),3))
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),3)-eps_measured(jumps(12)+1:jumps(13),1))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
legend('$\vartheta$','$\Theta(\varphi)$');
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
ylabel('[rad]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');

%% dot_theta and dot Theta
figure
grid on;
hold on;
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),5))
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),5)-eps_measured(jumps(12)+1:jumps(13),2))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
legend('$\dot{\vartheta}$','$\Theta''(\varphi)\dot{\varphi}$');
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
ylabel('[$\frac{rad}{s}$]','fontsize',14,'interpreter','latex');
sgtitle('$\dot{\vartheta}$ and the derivative of the VHC $\Theta''(\varphi)\dot{\varphi}$','Interpreter','latex');

%% dot_theta and dot Theta
figure
grid on;
hold on;
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),6))
plot(phi(jumps(12)+1:jumps(13)),states(jumps(12)+1:jumps(13),6)-eps_measured(jumps(12)+1:jumps(13),3))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
legend('$\dot{\varphi}$','$\dot{\varphi}_\star$');
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex')
ylabel('[$\frac{rad}{s}$]','fontsize',14,'interpreter','latex');
sgtitle('$\dot{\varphi}$ and $\dot{\varphi}_\star$','Interpreter','latex');


%% Epsilon_measured
figure
subplot(3,1,1)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),eps_measured(jumps(12)+1:jumps(13),1))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$y_1$[rad]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,2)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),eps_measured(jumps(12)+1:jumps(13),2))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$y_2$[$\frac{{rad}}{{s}}$]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,3)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),eps_measured(jumps(12)+1:jumps(13),3))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
ylabel('$z [\frac{{rad}}{{s}}]$','fontsize',14,'interpreter','latex');
sgtitle('Measured transverse coordinates','Interpreter','latex');

%% Epsilon_observed
figure
subplot(3,1,1)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),epsilon(jumps(12)+1:jumps(13),1))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\hat{y}_1$[rad]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,2)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),epsilon(jumps(12)+1:jumps(13),2))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\hat{y}_2$[$\frac{{rad}}{{s}}$]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,3)
grid on; hold on;
plot(phi(jumps(12)+1:jumps(13)),epsilon(jumps(12)+1:jumps(13),3))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
ylabel('$\hat{z} [\frac{{rad}}{{s}}]$','fontsize',14,'interpreter','latex');
sgtitle('Estimated transverse coordinates','Interpreter','latex');

%% Epsilon_observed_many
figure
subplot(3,1,1)
grid on; hold on;
scatter(phi(:),epsilon(:,1),1)
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\hat{y}_1$[rad]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,2)
grid on; hold on;
scatter(phi(:),epsilon(:,2),1)
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\hat{y}_2$[$\frac{{rad}}{{s}}$]','fontsize',14,'interpreter','latex');
sgtitle('$\vartheta$ and the VHC $\Theta(\varphi)$','Interpreter','latex');
subplot(3,1,3)
grid on; hold on;
scatter(phi(:),epsilon(:,3),1)
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
xlabel('mod($\varphi$,2$\pi$)[rad]','fontsize',14,'interpreter','latex');
ylabel('$\hat{z} [\frac{{rad}}{{s}}]$','fontsize',14,'interpreter','latex');
sgtitle('Estimated transverse coordinates for entire run','Interpreter','latex');

%% Input to the system
figure
grid on;
hold on;
plot(states(:,1),states(:,2))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
ylabel('u[Nm]','fontsize',14,'interpreter','latex');
sgtitle('Input to the Butterfly Robot','Interpreter','latex');

ts = timeseries([states(:,3) states(:,4) states(:,5) states(:,6)], states(:,1));
torque = timeseries(states(:,2), states(:,1));
xy = timeseries( [states(:,7) states(:,8)], states(:,1));

phi_2 = mod(out.lambda_and_varphi.data(:,1),2*pi);
l1 = out.lambda_and_varphi.data(:,2);
l2 = out.lambda_and_varphi.data(:,3);
jumps = [];
k = 1;
for i = 1:length(phi_2)-1
    if phi_2(i+1)-phi_2(i) < -3
        jumps(k) = i;
        k = k+1;
    end
end
figure
subplot(2,1,1)
grid on; hold on;
plot(phi_2(jumps(12)+1:jumps(13)),l1(jumps(12)+1:jumps(13)))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\lambda_1$[N]','fontsize',14,'interpreter','latex');
sgtitle('$\lambda_1$ as a function of $\varphi$','Interpreter','latex');
subplot(2,1,2)
grid on; hold on;
plot(phi_2(jumps(12)+1:jumps(13)),l2(jumps(12)+1:jumps(13)))
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\lambda_2$[N]','fontsize',14,'interpreter','latex');
sgtitle('Constraint forces as functions of $\varphi$','Interpreter','latex');

figure
grid on; hold on;
plot(phi_2(jumps(12)+1:jumps(13)),l1(jumps(12)+1:jumps(13)));
plot(phi_2(jumps(12)+1:jumps(13)),l2(jumps(12)+1:jumps(13)));
set(gca,'XTick',0:pi/4:2*pi)
set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
xlim([0 2*pi])
ylabel('$\lambda_2$[N]','fontsize',14,'interpreter','latex');
sgtitle('Constraint forces as functions of $\varphi$','Interpreter','latex');

ts = timeseries([states(:,3) states(:,4) states(:,5) states(:,6)], states(:,1));
torque = timeseries(states(:,2), states(:,1));
xy = timeseries( [states(:,7) states(:,8)], states(:,1));

