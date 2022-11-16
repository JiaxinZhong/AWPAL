clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;% + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1)*1;
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2)*1;

a = 0.1;
b = inf;
m1 = 100;
m1 = (0:150).';
ma = 0;
% ma = 10;

tic
radial = PalLineSrc_CWE_RadialJ(...
    k1, k2, ka, ...
    a, a, b, m1, ma);
toc

% fig = Figure;
% pcolor(ma, m1, log10(abs(real(radial))));
% fig.Init;

figure
plot(m1, real(radial));
