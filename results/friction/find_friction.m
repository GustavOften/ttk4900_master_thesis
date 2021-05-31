fric_res = load('friction_rising_10_to_matlab.txt');
fric_v_1 = load('friction_v_ref_1_to_matlab.txt');
fric_v_n_1 = load('friction_v_ref_negative_1_to_matlab.txt');
J = 8.8500e-04;
[~,time_10_index] = min(abs(fric_res(:,1)-10));
[~,time_30_index] = min(abs(fric_res(:,1)-30));
[~,time_31_index] = min(abs(fric_res(:,1)-31));
time_50_index = length(fric_res);
a = 0.4;
u_negative = sum(fric_res(time_10_index:time_30_index,2)) ...
    /length(fric_res(time_10_index:time_30_index,2));
u_positive = sum(fric_res(time_31_index:time_50_index,2)) ...
    /length(fric_res(time_31_index:time_50_index,2));
fric_negative = a*J - u_negative;
fric_positive = a*J - u_positive;
figure
hold on;
plot(fric_res(time_10_index:end,1),fric_res(time_10_index:end,2));
plot(fric_res(time_10_index:end,1),ones(length(fric_res(time_10_index:end,1)),1)*fric_positive);
plot(fric_res(time_10_index:end,1),ones(length(fric_res(time_10_index:end,1)),1)*fric_negative);
figure
plot(fric_res(time_10_index:end,1),fric_res(time_10_index:end,5))
figure
plot(fric_res(time_10_index:end,1),fric_res(time_10_index:end,3));

t_10 = find(fric_v_1 > 10, 1);
figure
plot(fric_v_1(t_10:end,1),fric_v_1(t_10:end,2))
fric_positive = sum(fric_v_1(t_10:end,2))/(length(fric_v_1(t_10:end,2)));
figure
plot(fric_v_1(t_10:end,1),fric_v_1(t_10:end,5))


t_10 = find(fric_v_n_1 > 10, 1);
figure
plot(fric_v_n_1(t_10:end,1),fric_v_n_1(t_10:end,2))
fric_negative = sum(fric_v_n_1(t_10:end,2))/(length(fric_v_n_1(t_10:end,2)));
figure
plot(fric_v_n_1(t_10:end,1),fric_v_n_1(t_10:end,5))

figure
subplot(1,2,1)
plot(fric_v_n_1(t_10:end,1),fric_v_n_1(t_10:end,2))
subplot(1,2,2)
plot(fric_v_1(t_10:end,1),fric_v_1(t_10:end,2))

highpass(fric_res(:,5),5,300)
