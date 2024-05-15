%
clear all

%% parameters
% source
src.wav.freq = 40e3;
src.wav.num = 2*pi*src.wav.freq/343 + 1i*AbsorpAttenCoef(src.wav.freq, 'temperature', 20, 'humidity', 70);
src.r = 0.1;

% uniform profile
src.prf.name = 'uniform';
src.prf.val = @(rhos) 1;

% focusing profile
% src.prf.name = 'focus';
% src.prf.focal_dist = 0.2; % focal distance
% src.prf.val = @(rs) exp(-1i * real(src.wav.num) * sqrt(rs.^2 + src.prf.focal_dist^2));


% field points
fp.x = linspace(-1, 1, 1e2).';
fp.y = 0;
fp.z = linspace(0, 3, 1.1e2);

%% main function
% profile on -memory

[prs, vel] = DIM3D(src, fp, 'int_num', 200, 'int_coord', 'polar');

% profile viewer

spl = 20*log10(abs(prs)/20e-6/sqrt(2));

%% 
figure;
% pcolor(fp.z, fp.x, abs(prs));
pcolor(fp.z, fp.x, spl);
colormap(MyColor('vik'));
shading interp
colorbar