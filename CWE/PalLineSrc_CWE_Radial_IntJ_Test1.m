clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);
a = 0.2;

m1 = nan;
ma = nan;
m1_max = 80;
ma_max = 30;
radial = PalLineSrc_CWE_Radial_IntJ(...
    k1, k2, ka, a, m1, ma, ...
    'm1_max', m1_max, 'ma_max', ma_max, ...
    'lower_limit', a, 'upper_limit', inf);
