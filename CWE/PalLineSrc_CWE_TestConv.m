% Test conv
% variable: ma_max
clear all
c0 = 343;
fu = 40e3;
fa = 4e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.1;

%% The field points 
rho = 1;
phi = linspace(0, pi, 1e2);
m1_max = 8e1;
ma_max = 200;
ma = (0:ma_max).';

is_farfield = true;
tic
[prs, ~, prs_m] = PalLineSrc_CWE(...
    k1, k2, ka, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'is_farfield', is_farfield);
toc
prs_cum = zeros(ma_max+1,length(phi));
for i = 0:ma_max
    prs_cum(i+1,:) = sum(prs_m(:,:,(-i:i)+ma_max+1), 3);
end

spl_cum = PrsToSpl(prs_cum*1e-4);

fig = Figure;
pcolor(phi/pi*180, ma, spl_cum)
fig.Init;

fprintf("ka*a = %-g\n", ka*a)
