clear all

pal = PalWave('audio_freq', 10, 'ultra_freq', 40e3);
a = 0.03;

x = linspace(-0.25, .25, 1e2).'*1;
% x = 0;
rho = abs(x);
z = linspace(0, .7, 1e3)*1;

tic
prs1 = CircSrc_CWE(pal.ultra_low.num, a, rho, z, 'int_num', 4e2);
prs2 = CircSrc_CWE(pal.ultra_high.num, a, rho, z, 'int_num', 4e2);
toc
q = conj(prs1) .* prs2;
q_phase = angle(q);
Q = 20*log10(abs(q));

fig = Figure;
pcolor(z, x, Q)
fig.Init;
caxis([-70,0]+max(Q(:)))

fig2 = Figure;
pcolor(z, x, q_phase)
fig2.Init;