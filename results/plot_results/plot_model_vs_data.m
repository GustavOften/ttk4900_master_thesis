figure
time = 9;
l = 1;
for i = 1:4
    time = time + 1.2;
    time_n = find(ts.Time > time,1);
    run_i = sim('butterfly_matlab_to_c');
    phi = mod(run_i.simout(:,2),2*pi);
    dtheta = run_i.simout(:,3);
    dphi = run_i.simout(:,4);
    mod_dtheta = run_i.simout(:,5);
    mod_dphi = run_i.simout(:,6);
        jumps = [];
    k = 1;
    for j = 1:length(phi)-1
    if phi(j+1)-phi(j) < -3
    jumps(k) = j;
    k = k+1;
    end
    end
    subplot(4,2,l)
    l = l+1;
    grid on;
    hold on;
    plot(phi(1:jumps(1)),dtheta(1:jumps(1)),'LineWidth',1)
    plot(phi(1:jumps(1)),mod_dtheta(1:jumps(1)),'LineWidth',1)
    set(gca,'XTick',0:pi/4:2*pi)
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
    xlim([0 2*pi])
    if i == 4
        xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
    end
    ylabel('$\dot{\vartheta}$','fontsize',14,'interpreter','latex')
    subplot(4,2,l)
    l = l+1;
    grid on;
    hold on;
    plot(phi(1:jumps(1)),dphi(1:jumps(1)),'LineWidth',1)
    plot(phi(1:jumps(1)),mod_dphi(1:jumps(1)),'LineWidth',1)
    set(gca,'XTick',0:pi/4:2*pi)
    set(gca,'XTickLabel',{'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'})
    xlim([0 2*pi])
    ylabel('$\dot{\varphi}$','fontsize',14,'interpreter','latex')
    if i == 4
        xlabel('$\varphi$[rad]','fontsize',14,'interpreter','latex');
    end
end
sgtitle('Data from measurements insearted in model','Interpreter','latex');

