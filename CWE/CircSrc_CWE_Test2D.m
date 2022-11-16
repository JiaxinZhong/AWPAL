clear all

f = 40e3;
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);
a = 0.025;

x = linspace(-0.25, .25, 1e2).'*1;
% x = 0;
rho = abs(x);
z = linspace(0, .7, 1e3)*1;

tic
prs = CircSrc_CWE(k, a, rho, z, 'int_num', 4e2);
toc
spl = prs2spl(prs*1.21*343/(pi*a.^2)*1e-4);

if numel(x) == 1
    fig = Figure;
    plot(z, spl+12)
else
    fig = Figure;
    pcolor(z, x, spl)
    fig.Init;
    caxis([80,130])
end