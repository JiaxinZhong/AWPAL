wav = SoundWave('freq', 40e3);
a = 0.005;
k = wav.num;

theta = linspace(0, pi/2, 2e2).';
jinc = Jinc(real(k) .* a .* sin(theta));

fig = Figure;
subplot(211)
plot(theta/pi*180, jinc);
subplot(212)
plot(theta/pi*180, 20*log10(jinc));

