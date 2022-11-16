clear all
s = 17.15/1e3;
xc = (-3.5:1:3.5).' * s;
a = 0.005;

wav = SoundWave('freq', 39e3);
p.name = 'steerable';
p.steer_angle = 70/180*pi;
prf = LineSrcArrayProfile('xc', xc, 'a', a);
prf.CreateProfile(wav, p);

trunc_term = 2e2;

rho = linspace(0, 2.5, 2e2).';
phi = linspace(0, pi, 5e1);
% phi = pi/2;
[x, y] = Polar2Cart(rho, phi);

tic
prs = LineSrc_CWE(wav.num, a, rho, phi, ...
    'trunc_term', trunc_term,...
    'array', prf);
toc
spl = prs2spl(prs*.1e-2);

%% calculate the directivity
dist = 3;
phi_dir = linspace(0, pi, 5e2);
dir = LineSrc_CWE(wav.num, a, dist, phi_dir, ...
    'trunc_term', trunc_term, ...
    'array', prf);
dir = PrsToSpl(dir);
dir = dir-max(dir(:));

fig = Figure;
pcolor(y.', x.', spl.')
fig.Init;
ylim([-1,1])
xlim([0,2]);
% caxis([-20,30]);
xlabel('y(m)');
ylabel('x(m)');
fig.SetColorbarTitle('SPL (dB)')
% fig.ExportTikz('filename', 'CWE/fig/LineSrc_CWE_Test2DUniform_220625B_.tex')

fig_dir = Figure;
plot(phi_dir/pi*180, dir)