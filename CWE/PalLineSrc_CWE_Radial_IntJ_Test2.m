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
m1 = 0;
% m1 = (0:150).';
ma = 100;
% ma = 10;

tic
radial_matlab = PalLineSrc_CWE_Radial_IntJ_Matlab(...
    k1, k2, ka, ...
    a, a, b, ...
    m1, ma);
toc

tic
radial = ...
    exp(PalLineSrc_CWE_Radial_IntJ(...
    k1, k2, ka, a,  ...
    'lower_limit', a,...
    'upper_limit', inf,...
    'm1', m1,...
    'ma', ma, ...
    'is_log', true));
toc

[~, rel_err] = Error(radial_matlab, radial);
fprintf("Calculated by matlab: %g %+gi.\n", real(radial_matlab), imag(radial_matlab));
fprintf("Calculated by me: %g %+gi.\n", real(radial), imag(radial));
fprintf("Relative error: %g %+gi.\n", real(rel_err), imag(rel_err));


% fig = Figure;
% pcolor(ma, m1, log10(abs(real(radial))));
% fig.Init;

% figure
% plot(m1, real(radial));
