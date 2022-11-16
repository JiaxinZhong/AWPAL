% Calculate the 2d sound field 
clear all
a = 0.03;
pal = PalWave('audio_freq', 0, 'ultra_freq', 40e3);

z = linspace(0, 1.6, 1e2);
x = linspace(-.8, .8, 1.4e2).';

p0 = 0.12*1.21*343;

tic
q = PalCircSrc_GBE_SrcDensity(pal, a, abs(x), z);
toc
Q = 20*log10(abs(q));

figure;
pcolor(z, x, Q);
colormap jet
shading interp
colorbar
caxis([-40,0]+max(Q(:)))