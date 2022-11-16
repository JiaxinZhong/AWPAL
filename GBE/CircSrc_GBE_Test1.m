% Calculate the 2d sound field 
a = 0.1;
wav = SoundWave('freq', 40e3);

z = linspace(0, 5, 2e2);
rho = linspace(0, 2, 2.2e2).';

prs = CircSrc_GBE(wav.num, a, rho, z);
spl = prs2spl(prs);

figure;
pcolor(z, rho, spl);
colormap jet
shading interp
colorbar
caxis([0,100])