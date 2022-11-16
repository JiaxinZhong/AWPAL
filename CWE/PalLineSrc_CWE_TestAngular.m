% Test conv
clear all
pal = PalWave('audio_freq', 2e3, 'ultra_freq', 40e3);

a = 0.08;

%% The field points 
rho = 3;
phi = linspace(0, pi, 2e2);
[x, y] = Polar2Cart(rho, phi);
m1_max = 1e2;
ma_max = 80;
src_norm = 'surf_vel';
eqn = 'Westervelt';
profile = 'steerable';
steer_angle = 70/180*pi;

s = 17.15/1e3;
xc = (-3.5:1:3.5).' * s;
a = 0.005;
wav = SoundWave('freq', 39e3);
p.name = 'steerable';
p.steer_angle = 70/180*pi;
prf = LineSrcArrayProfile('xc', xc, 'a', a);
prf.CreateProfile(wav, p);

tic
prs = PalLineSrc_CWE(...
    pal, a, rho, phi, ...
    'm1_max', m1_max,...
    'ma_max', ma_max,...
    'src_norm', src_norm,...
    'eqn', eqn, ...
    'profile', profile,...
    'steer_angle', steer_angle,...
    'array', prf);
toc

spl = PrsToSpl(prs*(1e-2)^2);

fig = Figure;
plot(phi/pi*180, spl)