clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.1;

%% The field points 
rho = logspace(-2,2,1e2).';
phi = pi/2;
[x, y] = Polar2Cart(rho, phi);
m1_max = ceil(1.2*real(k1*a));
ma_max = 20;
src_norm = 'surf_vel';
eqn = 'Westervelt';

tic
prs_Westervelt = PalLineSrc_CWE(...
    k1, k2, ka, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'src_norm', src_norm,...
    'eqn', eqn, ...
    'is_farfield', false);
toc
tic
prs_far = PalLineSrc_CWE(...
    k1, k2, ka, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'src_norm', src_norm,...
    'eqn', eqn, ...
    'is_farfield', true);
toc

spl_Westervelt = PrsToSpl(prs_Westervelt*(1e-2)^2);
spl_far = PrsToSpl(prs_far*(1e-2)^2);

figure 
subplot(211)
semilogx(rho, spl_Westervelt)
hold on
plot(rho, spl_far)
subplot(212)
semilogx(rho, spl_far-spl_Westervelt)
ylim([-3,3])