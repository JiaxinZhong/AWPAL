clear all
src.shape = 'line';
src.radius = 0.1;
src.prf.name = 'cosine';
src.prf.phi = (90+15)/180*pi;
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
phi = linspace(0, pi, 2e2).';

tic
dir_convention = PalLineSrc_Conv(pal, phi, ...
    'type', 'conventional');
dir_improved = PalLineSrc_Conv(pal, phi, ...
    'type', 'improved');
toc

dir_lvl_covention = PrsToSpl(dir_convention);
dir_lvl_covention = dir_lvl_covention - max(dir_lvl_covention);
dir_lvl_improved = PrsToSpl(dir_improved);
dir_lvl_improved = dir_lvl_improved - max(dir_lvl_improved);

figure 
plot(phi/pi*180, dir_lvl_covention)
hold on
plot(phi/pi*180, dir_lvl_improved, '--')
