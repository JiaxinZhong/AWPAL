% INTRO
%   - Plot 2D sound fields
% clear all

%% parameters settings
src.wav.freq = 40e3; % frequency
src.wav.num = 2*pi*src.wav.freq/343 + 1i*AbsorpAttenCoef(src.wav.freq, 'temperature', 20, 'humidity', 70); % wavenumber
src.r = 0.1; % radius of the source
src.prf.phi_m = 0; % azimuthal order of the profile
% uniform profile
src.prf.name = 'uniform';
src.prf.val = @(rs) 1;
% focusing profile
% src.prf.name = 'focus';
% src.prf.focal_dist = 0.2; % focal distance
% src.prf.val = @(rs) exp(-1i * real(src.wav.num) * sqrt(rs.^2 + src.prf.focal_dist^2));

% field points
fp.r = logspace(-3, 1, 2e2).';
fp.theta = linspace(0, pi/2, 1e2);
fp.phi = 0;
[fp.x, fp.y, fp.z] = Sph2Cart(fp.r, fp.theta, fp.phi);

%% main function
prs = SWE3D(src, fp);

%%
z = [flipud(fp.z.'); fp.z.'];
x = [-flipud(fp.x.'); fp.x.'];
prs_show = [flipud(prs.'); prs.'];
spl_show = 20*log10(abs(prs_show) / sqrt(2) /20e-6);

%% 2D field
figure;
% pcolor(z, x, abs(prs_show))
pcolor(z, x, spl_show)
colormap(MyColor('vik'))
shading interp
colorbar
xlim([0,5])
ylim([-3.09,3.09]/2)
% clim([40,120])
xlabel('z (m)')
ylabel('x (m)')