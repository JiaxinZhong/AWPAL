% Calculate the 2d sound field 
clear all
a = 0.03;
pal = PalWave('audio_freq', 500, 'ultra_freq', 40e3);

z = linspace(0, 1.6, 1e2);
x = linspace(-.8, .8, 1.4e2).';

p0 = 0.12*1.21*343;

tic
prs = PalCircSrc_GBE(pal, a, abs(x), z, 'int_num', 1e3);
toc
spl = prs2spl(prs*p0^2);

figure;
pcolor(z, x, spl);
colormap jet
shading interp
colorbar
caxis([20,45])
