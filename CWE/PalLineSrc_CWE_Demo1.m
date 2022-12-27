clear all

%% wave info
prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'cosine', 'order', 1);
% prf = SrcProfile('name', 'steerable', 'phi', (90+15)/180*pi);
% prf = SrcProfile('name', 'cosine_steerable', 'order', 1, 'phi', (90+15)/180*pi);
src = LineSrc('radius', .05, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
fp = Point2D('rho', linspace(0, 3.4, 5.1e1).', 'phi', linspace(-pi/2, pi/2, 6e1));
fp.Polar2Cart();
ma_max = 120;

tic
prs = PalLineSrc_CWE(pal, fp, 'ma_max', ma_max);
toc
spl = PrsToSpl(prs);

%% plot results
fig = Figure;
pcolor(fp.x, fp.y, spl);
xlim([0,3]);
ylim([-1.5,1.5])
fig.Init;