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
b = logspace(log10(a*2), 3, 1e2).';
m1 = 0;
% m1 = (0:150).';
ma = 100;
% ma = 10;

R = 0*b;
figure;
hold on
for i = 1:length(b)
    tic
    R(i) = PalLineSrc_CWE_Radial_IntJ_Matlab(...
        k1, k2, ka, ...
        a, a, b(i), ...
        'm1', m1, 'ma', ma);
    toc
    plot(b(i), real(R(i)), 'k.', 'markersize', 10)
    drawnow
end

