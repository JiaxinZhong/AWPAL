clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;% + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1)*1;
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2)*1;

a = .1;
b = inf;
m1 = 0;
% m1 = (0:150).';
ma = 10;
% ma = 10;

tic
R = PalLineSrc_CWE_Radial_IntJ_Matlab(...
    k1, k2, ka, ...
    a, a, b, ...
    'm1', m1, 'ma', ma);
toc