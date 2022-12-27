%% calculate the directivity
clear all

%% wave info
prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'cosine', 'order', 1);
% prf = SrcProfile('name', 'steerable', 'phi', (90+15)/180*pi);
% prf = SrcProfile('name', 'cosine_steerable', 'order', 1, 'phi', (90+15)/180*pi);
src = LineSrc('radius', .2, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
fp = Point2D('rho', 50, 'phi', linspace(-pi/2, pi/2, 3e2));
fp.Polar2Cart();
ma_max = 120;

tic
dir = PalLineSrc_CWE(pal, fp, 'ma_max', ma_max, 'is_farfield', false);
toc
dir_dB = PrsToSpl(dir);
dir_dB = dir_dB - max(dir_dB(:));

%% plot results
fig = Figure;
plot(fp.phi/pi*180, dir_dB);
xlim([-90, 90]);
xlabel('phi')
ylabel('Directivity (dB)')
fig.Init;