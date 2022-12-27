% Test conv
% variable: ma_max
clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.05;

%% The field points 
rho = [linspace(0,0.2,8e1).'; linspace(0.2, 1.2,   1e2).'];
phi = linspace(0, pi, 5e1);
[x, y] = Polar2Cart(rho, phi);
m1_max = 8e1;
ma_max = 20;
src_norm = 'surf_vel';

tic
[prs, prs_cum] = PalLineSrc_CWE(...
    k1, k2, ka, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'src_norm', src_norm,...
    'eqn', 'Westervelt');
toc

spl = PrsToSpl(prs*(1e-2)^2);
spl_cum = PrsToSpl(prs_cum*(1e-2)^2);

fig_cum = Figure;
pcolor(y, x+0.5*exp(1)/pi, spl_cum);
fig_cum.Init
xlim([0,1]);
ylim([0,1]*exp(1)/pi)
caxis([35,65])
% fig_cum.ExportTikz('filename', 'CWE/fig/PalLineSrc_CWE_Test2D_220627A_Cum_.tex');

fig = Figure;
pcolor(y, x+0.5*exp(1)/pi, spl);
fig.Init
xlim([0,1]);
ylim([0,1]*exp(1)/pi)
caxis([35,65])
% fig.ExportTikz('filename', 'CWE/fig/PalLineSrc_CWE_Test2D_220627A_Local_.tex');

fig_diff = Figure;
pcolor(y, x+0.5*exp(1)/pi, spl-spl_cum);
fig_diff.Init
xlim([0,1]);
ylim([0,1]*exp(1)/pi)
caxis([-1,1])
% fig_diff.ExportTikz('filename', 'CWE/fig/PalLineSrc_CWE_Test2D_220627A_Diff_.tex');

% save('CWE/data/PalLineSrc_CWE_Test2D_220627A_.mat')
fprintf("ka*a = %-g\n", ka*a)
