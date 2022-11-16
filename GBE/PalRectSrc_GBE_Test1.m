% Calculate the 2d sound field 
clear all
ax = 0.04;
ay = 0.02;
pal = PalWave('audio_freq', 500, 'ultra_freq', 40e3);

x = linspace(-.8, .8, .9e2).';
y = 0;
z = linspace(0, 1.6, 1e2);

p0 = 0.12*1.21*343;

tic
prs = PalRectSrc_GBE(pal, ax, ay, x, y, z, 'int_num', 1e3);
toc
spl = prs2spl(prs * p0^2);

figure;
pcolor(z, x, spl);
colormap jet
shading interp
colorbar
caxis([20,45]-3)
