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

m1 = 0;
ma = 100;

int_num = (2:200).';

int_method = 'gauss';

int = int_num * 0;
for i = 1:length(int_num)
    tic
    fprintf("Processing %d of %d ...\n", i, length(int_num));
    int(i) = PalLineSrc_CWE_RadialH2(k1, k2, ka, a, m1, ma, ...
        'int_num', int_num(i),...
        'int_method', int_method, ...
        'is_log', true);
    int(i) = exp(int(i));
    toc
end
 int0 = PalLineSrc_CWE_RadialH2(k1, k2, ka, a, m1, ma, ...
        'int_method', 'matlab');
rho_vsrc = logspace(-1, 3, 1e2).';

figure;
plot(int_num, real(int));
hold on
plot(int_num, imag(int), '--');

fprintf("Exact: %g %+gi\n", real(int0), imag(int0));
fprintf("Calculated: %g %+gi\n", real(int(end)), imag(int(end)));

