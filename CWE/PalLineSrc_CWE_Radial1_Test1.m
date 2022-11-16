% Calculate a single point

c0 = 343;
fu = 40e3;
fa = 4e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0 + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1);
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2);

a = 0.1;
ma_max = 2e2;
ma = (-ma_max:ma_max).';

% The field points 
rho = [0.2];
% phi = [0, pi/3];

profile = 'steerable';
steer_angle = pi/4;
is_farfield = true;

int_num = (2:1e2);
radial = 0 * ma .* int_num;
for i = 1:length(int_num)
    fprintf('Processing %d of %d...\n', i, length(int_num));
    tic
    radial(:,i) = PalLineSrc_CWE_Radial1(k1, k2, ka, a, rho, ...
        'int_num', int_num(i), ...
        'ma_max', ma_max, ...
        'profile', profile, 'steer_angle', steer_angle, ...
        'is_farfield', is_farfield);
    toc
end
radial_diff = radial(:,end) - radial(:,end-1);
radial_diff = real(radial_diff)./real(radial(:,end)) ...
    + 1i * imag(radial_diff)./imag(radial(:,end));

% fig = Figure;
% subplot(211)
% plot(int_num, real(radial));
% title('Real');
% subplot(212)
% plot(int_num, imag(radial));
% title('Imag');

fig = Figure;
plot(ma, real(radial_diff));
hold on
plot(ma, imag(radial_diff));

idx_ma = ma_max + 1 + 11;
fig = Figure;
plot(int_num, real(radial(idx_ma,:)));
hold on
plot(int_num, imag(radial(idx_ma,:)));
