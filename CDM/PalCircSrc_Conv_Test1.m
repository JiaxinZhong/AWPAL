clear all
src.shape = 'circle';
src.radius = 0.05;
src.prf.name = 'steerable';
src.prf.theta = 10/180*pi;
src.prf.phi = 0;
pal = PalWave('audio_freq', 4e3, 'ultra_freq', 60e3, 'src', src);

%% The field points 
theta = linspace(0, pi/2, 2e2).';
phi = [0, pi];

tic
dir_convention = PalPlanarSrc_Conv(...
    pal, theta, phi,  ...
    'type', 'conventional');
dir_modified =  PalPlanarSrc_Conv(...
     pal, theta, phi,  ...
    'type', 'modified');
toc
dir_convention = [flip(dir_convention(:,2)); dir_convention(:,1)];
dir_modified = [flip(dir_modified(:,2)); dir_modified(:,1)];
angle = [-flip(theta); theta ]/pi*180;

k1 = pal.ultra_low.num;
k2 = pal.ultra_high.num;
ka = pal.audio.num;
ampl = 1.2*ka*src.radius./(imag(k1)+imag(k2)) /pi *sqrt(k1*k2);
dir_lvl_covention = PrsToSpl(dir_convention/ampl);
dir_lvl_covention = dir_lvl_covention - max(dir_lvl_covention);
dir_lvl_modified = PrsToSpl(dir_modified/ampl);
dir_lvl_modified = dir_lvl_modified - max(dir_lvl_modified);

figure 
plot(angle, dir_lvl_covention)
hold on
plot(angle, dir_lvl_modified, '--')
% ylim([-90,0])