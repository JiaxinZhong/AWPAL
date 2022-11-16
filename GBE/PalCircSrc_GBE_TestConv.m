% Calculate the 2d sound field 
clear all
a = 0.04;
pal = PalWave('audio_freq', 500, 'ultra_freq', 40e3);

rho = .3;
z = 0.3;

p0 = 0.12*1.21*343;

int_num = (2e2:2e3).';
prs = int_num * 0;
for i = 1:length(int_num)
    tic
    fprintf("Processing %d of %d...\n", i, length(int_num));
    prs(i) = PalCircSrc_GBE(pal, a, rho, z, 'int_num', int_num(i));
    toc
end
spl = prs2spl(prs * p0^2);

figure;
plot(int_num, spl);
