clear all

%% wave info
prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'cosine', 'order', 1);
% prf = SrcProfile('name', 'steerable', 'phi', (90+15)/180*pi);
% prf = SrcProfile('name', 'cosine_steerable', 'order', 1, 'phi', (90+15)/180*pi);
src = LineSrc('radius', .05, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
fp = Point2D('rho', 1, 'phi', linspace(0, pi, 2e2));
fp.Polar2Cart();
ma_max = 250;

tic
dir_exact = PalLineSrc_CWE( ...
    pal, fp, ...
    'ma_max', ma_max,...
    'is_farfield', true);
toc
dir_exact_lvl = PrsToSpl(dir_exact);
dir_exact_lvl = dir_exact_lvl - max(dir_exact_lvl);

%% Convolution model
dir_direct = PalLineSrc_Conv(pal, fp.phi, 'type', 'direct', 'is_norm_dB', true);
[dir_improved, dir_aperture] = PalLineSrc_Conv(pal, fp.phi, 'type', 'improved', 'is_norm_dB', true);
dir_westervelt = PalLineSrc_Conv(pal, fp.phi, 'type', 'Westervelt', 'is_norm_dB', true);

fig = Figure ;
hold on
plot(fp.phi/pi*180, dir_exact_lvl)
plot(fp.phi/pi*180, dir_direct, '-.')
plot(fp.phi/pi*180, dir_improved, '--')
% plot(fp.phi/pi*180, dir_westervelt, ':')
legend('Exact solution', 'Direct conv.', 'Improved conv.', 'Westervelt')

% fn = sprintf('PalLineSrc_CWE_TestDir_220715A_%dkHz_%dHz_%dcm_%s_', ...
%     pal.ultra.freq/1e3, pal.audio.freq, pal.src_ultra.radius*1e2, ...
%     src.prf.name);
% fig.ExportTikz('filename', ['CWE/fig/', fn, '.tex']);

fig_diff = Figure;
hold on
plot(fp.phi/pi*180, abs(dir_direct - dir_exact_lvl), '-.')
plot(fp.phi/pi*180, abs(dir_improved - dir_exact_lvl), '--')
% plot(fp.phi/pi*180, abs(dir_westervelt - dir_exact_lvl), ':')
title(sprintf('f = %d kHz, a = %dcm', pal.audio.freq/1e3, pal.src_ultra.radius*1e2))
legend('Direct conv.', 'Improved conv.', 'Westervelt')


% save(['CWE/data/', fn, '.mat'])