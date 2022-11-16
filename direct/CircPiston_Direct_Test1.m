% Plot 2D fields
clear all
c0 = 343;
f = 40e3;
k = 2*pi*f/c0 + 1i*AbsorpAttenCoef(f);
a = 0.04;
x = linspace(-1,1, 3e2).';
y = 0;
z = linspace(0, 2, 5e2);

tic
prs = CircPiston_Direct(k, a, x, y, z);
toc

spl = prs2spl(prs);

fig = Figure;
pcolor(z, x, spl);
fig.Init;
