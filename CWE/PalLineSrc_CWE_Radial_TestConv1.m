% Test conv
% vary m1_max
clear all
c0 = 343;
fu = 40e3;
fa = 2e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

is_farfield = true;

a = 0.1;

%% The field points 
rho = 1;
% phi = [0, pi/3];
m1_max = (10:5:180);
ma_max = 20;
ma = (-ma_max:ma_max).';

radial = zeros(1, 1, 2*ma_max+1, length(m1_max));
for i = 1:length(m1_max)
    fprintf("Processing %d of %d.\n", i, length(m1_max));
    tic
    radial(1,1,:,i) = PalLineSrc_CWE_Radial(k1, k2, (ka), a, rho, ...
        'm1_max', m1_max(i), ...
        'ma_max', ma_max, ...
        'is_farfield', is_farfield);
    toc
end

radial = squeeze(radial);
idx = ma_max+1;
fig = Figure;
plot(m1_max, real(radial(idx,:)));
hold on
plot(m1_max, imag(radial(idx,:)), '--');

fig2 = Figure;
pcolor(m1_max, ma, real(radial));
fig2.Init;
xlabel('m1 max')
ylabel('ma')

fprintf("k1*a = %-g\n", k1*a)
