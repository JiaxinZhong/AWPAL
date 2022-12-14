% Calculate the 2D field generated by a source with a continuous steerable profile
clear all

trunc_term = 2e2;

c0 = 343;
f = 40e3;
k = 2*pi*f/c0 + 1i*AbsorpAttenCoef(f);

rho = linspace(0, 2.5, 2e2).';
phi = linspace(0, pi, 5e1);
% phi = pi/2;
[x, y] = Polar2Cart(rho, phi);
% 
% x = linspace(-1, 1, 2e2).';
% y = linspace(0, 2, 3e2);
% [rho, phi] = Cart2Polar(x, y);

a = .1;

profile = 'steerable';
src_norm = 'surf_vel';
steer_angle = [pi/4;pi/2;2*pi/3];
tic
prs = LineSrc_CWE(k, a, rho, phi, ...
    'trunc_term', trunc_term,...
    'profile', profile, ...
    'src_norm', src_norm, ...
    'profile', profile, ...
    'steer_angle', steer_angle);
toc
spl = prs2spl(prs*.1e-2);

fig = Figure;
pcolor(y.', x.', spl.')
fig.Init;
ylim([-1,1])
xlim([0,2]);
% caxis([90,140]);
xlabel('y(m)');
ylabel('x(m)');
fig.SetColorbarTitle('SPL (dB)')
% fig.ExportTikz('filename', 'CWE/fig/LineSrc_CWE_Test2DUniform_220625B_.tex')
