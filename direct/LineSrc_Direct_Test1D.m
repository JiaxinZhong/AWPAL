% Test 2D field
x = 0;
y = linspace(0, 1.2, 4e2);
a = 0.05;
c0 = 343;
f = 1e3;
k = 2*pi*f/c0 + 1i*AbsorpAttenCoef(f);

tic
[prs, vel_x, vel_y, lag] = LineSrc_Direct(k, a, x, y, ...
    'src_norm', 'surf_vel',...
    'is_cal_vel', true,...
    'is_cal_lag', true);
toc

spl = PrsToSpl(prs*1e-1);

figure;
plot(y, spl);
ylim([110,150])
title('SPL');

figure
plot(y, real(prs));
hold on
plot(y, imag(prs), '--')
title('Pressure')

figure 
plot(y, real(vel_x));
hold on
plot(y, imag(vel_x), '--');
title('Vel x')

figure 
plot(y, real(vel_y));
hold on
plot(y, imag(vel_y), '--');
title('Vel y')

figure;
plot(y, lag);
title('Lagnrangian');
