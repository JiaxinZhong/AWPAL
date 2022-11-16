% Test conv
% vary rho
clear all
c0 = 343;
fu = 40e3;
fa = 1e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);


a = 0.1;

%% The field points 
rho = logspace(0, 3, 1e2).';
% phi = [0, pi/3];
m1_max = 50;
ma_max = 20;
ma = (-ma_max:ma_max).';

% radial = zeros(1, 1, 2*ma_max+1, length(m1_max));
tic
radial_near = PalLineSrc_CWE_Radial(k1, k2, (ka), a, rho, ...
    'm1_max', m1_max, ...
    'ma_max', ma_max, ...
    'is_farfield', false);
toc
radial_near = squeeze(radial_near);

tic
radial_far = PalLineSrc_CWE_Radial(k1, k2, (ka), a, rho, ...
    'm1_max', m1_max, ...
    'ma_max', ma_max, ...
    'is_farfield', true);
toc
radial_far = squeeze(radial_far);


idx = ma_max+1;

figure
subplot(211)
semilogx(rho, real(radial_near(:,idx)));
hold on
plot(rho, real(radial_far(:,idx)), '--');

subplot(212)
semilogx(rho, imag(radial_near(:,idx)));
hold on
plot(rho, imag(radial_far(:,idx)), '--');

% fig2 = Figure;
% pcolor(m1_max, ma, real(radial));
% fig2.Init;
% xlabel('m1 max')
% ylabel('ma')
% 
% fprintf("k1*a = %-g\n", k1*a)
