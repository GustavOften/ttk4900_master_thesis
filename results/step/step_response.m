step_001 = load('step_response_tau_001_matlab.txt');
step_005 = load('step_response_tau_005_matlab.txt');
steps_001 = [1;find(step_001(:,1) > 9.999)];
steps_005 = [1;687;1373;2064;2750];
figure
hold on;
a001 = 0;
for i = 1:4
    plot(step_001(steps_001(i)+1:steps_001(i+1),1), step_001(steps_001(i)+1:steps_001(i+1),5))
    a001 = a001 + (step_001(steps_001(i+1),5)-step_001(steps_001(i)+1,5))/(step_001(steps_001(i+1),1)-step_001(steps_001(i)+1,1));
end
a001 = a001/4;
figure
hold on;
a = 0;
for i = 1:4
    plot(step_005(steps_005(i)+1:steps_005(i+1),1), step_005(steps_005(i)+1:steps_005(i+1),5))
    a_new = (step_005(steps_005(i+1),5)-step_005(steps_005(i)+1,5))/(step_005(steps_005(i+1),1)-step_005(steps_005(i)+1,1));
    a = a + a_new;
    a_new
end
grid on;
xlabel('Time[s]','fontsize',14,'Interpreter','latex')
ylabel('$\dot{\vartheta}$ [rad/s]','fontsize',14,'Interpreter','latex')
title('Step response $\tau$ = 0.05[Nm]','fontsize',14,'Interpreter','latex')
a = a/4;
J = (0.05-0.0022)/a
J001 = (0.01-0.0022)/a001
