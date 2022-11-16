% Calculate a single point
% Variable: gauss_num
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
rho = [0.2];
% phi = [0, pi/3];
ma_max = 30;
ma = (-ma_max:ma_max);

int_num = (2:1:5e1).';
radial = zeros(length(int_num), 1, 2*ma_max+1);
for i = 1:length(int_num)
    fprintf("Processing %d of %d..\n", i, length(int_num));
    tic
    radial(i,1,:) = PalLineSrc_CWE_Radial(k1, k2, ka, a, rho, ...
        'ma_max', ma_max, ...
        'int_num', int_num(i));
    toc
end
radial = squeeze(radial);

idx_ma = ma_max+3;

figure;
plot(int_num, real(radial(:,idx_ma)));
hold on
plot(int_num, imag(radial(:,idx_ma)))

fig = Figure;
pcolor(ma, int_num, real(radial));
fig.Init;