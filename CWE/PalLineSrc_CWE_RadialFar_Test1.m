% Calculate a single point

c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;% + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.1;

%% The field points 
rho = [0.2; 1];
% phi = [0, pi/3];

tic
radial = PalLineSrc_CWE_RadialFar(k1, k2, ka, a, rho);
toc