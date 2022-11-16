% Test convergence
clear all
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

gauss_num = (5:5:200).';

prs = 0 * gauss_num;
for i = 1:length(gauss_num)
    fprintf('Processing %d of %d...\n', i, length(gauss_num));
    tic
    prs(i) = PalCircPiston_Direct(ka, k1, k2, a, x, y, z, ...
        'gauss_num', gauss_num(i));
    toc
end

spl = prs2spl(prs);

fig = Figure;
plot(gauss_num, spl)
