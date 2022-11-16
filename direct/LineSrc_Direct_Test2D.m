% Test 2D field
x = linspace(-1,1, 2e2).';
y = linspace(0, 2.2, 4e2);
a = 0.05;
wav = SoundWave('freq', 40e3);
src = LineSrc('wav', wav, 'radius', a);

tic
prs = LineSrc_Direct(src, x, y);
toc

spl = PrsToSpl(prs*1e-1);

fig = Figure;
pcolor(y, x, spl);
fig.Init;