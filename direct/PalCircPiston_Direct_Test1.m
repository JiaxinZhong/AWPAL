c0 = 343;
fa = 1e3;
fu = 40e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

x = 0;
y = 0;
z = 1;
a = 0.01;

tic
prs = PalCircPiston_Direct(ka, k1, k2, a, x, y, z)
spl = prs2spl(prs)
toc
