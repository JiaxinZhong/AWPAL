% Test 1D field
phi = pi/2;
rho = linspace(0, 1.2, 4e2);
a = 0.05;
c0 = 343;
f = 1e3;
k = 2*pi*f/c0 + 1i*AbsorpAttenCoef(f);

tic
[prs, vel_rho, vel_phi, lag] = LineSrc_CWE(k, a, rho, phi, ...
    'src_norm', 'surf_vel',...
    'is_cal_vel', true,...
    'is_cal_lag', true);
toc

[vel_x, vel_y] = VectorFieldPolar2Cart(...
        vel_rho, vel_phi, phi);
spl = PrsToSpl(prs*1e-1);

figure;
plot(rho, spl);
ylim([110,150])
title('SPL');

figure
plot(rho, real(prs));
hold on
plot(rho, imag(prs), '--')
title('Pressure')

figure 
plot(rho, real(vel_x));
hold on
plot(rho, imag(vel_x), '--');
title('Vel x')
% 
figure 
plot(rho, real(vel_y));
hold on
plot(rho, imag(vel_y), '--');
title('Vel y')

figure;
plot(y, lag);
title('Lagnrangian');
