% Calculate the 2d sound field 
clear all
a = 0.1;
wav = SoundWave('freq', 5);

z = logspace(-3, 1, 2e2);
rho = 0;

prs = CircSrc_GBE(wav.num, a, rho, z);
prs_exact = CircPiston_Axis_Exact(wav.num, a, z);
spl = prs2spl(prs*1.21*343);
spl_exact = PrsToSpl(prs_exact);

spl_delta = 1;
for i = length(z):-1:1
    if abs(spl_exact(i) - spl(i)) > spl_delta
        fprintf("Error at z = %gm.\n", z(i));
        break;
    end
end

figure;
semilogx(z, spl_exact)
hold on
plot(z, spl, '--');