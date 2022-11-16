% Test conv
clear all
c0 = 343;
fu = 40e3;
fa = 2e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.1;

%% The field points 
rho = [linspace(0,0.2,8e1).'; linspace(0.2, 1.2, 1e2).'];
phi = linspace(0, pi, 5e1);
[x, y] = Polar2Cart(rho, phi);
m1_max = 8e1;
ma_max = 50;
src_norm = 'surf_vel';
eqn = 'Westervelt';
profile = 'steerable';
steer_angle = pi/3;

tic
prs = PalLineSrc_CWE(...
    k1, k2, ka, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'src_norm', src_norm,...
    'eqn', eqn, ...
    'profile', profile,...
    'steer_angle', steer_angle);
toc

spl = PrsToSpl(prs*(1e-2)^2);

fig = Figure;
pcolor(y, x+0.5*exp(1)/pi, spl);
fig.Init
xlim([0,1]);
ylim([0,1]*exp(1)/pi)
fprintf("ka*a = %-g\n", ka*a)
