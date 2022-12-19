clear all

src = LineSrc('radius', 0.1);
pal = PalSrc('audio_freq', 1e3, 'ultra_freq', 40e3, 'src', src);

kax = linspace(0, 10*real(pal.audio.num), 2e2);

k1x = linspace(-1.5*real(pal.ultra_low.num), 1.5*real(pal.ultra_low.num), 1e3).';

k2x = k1x + kax;
k1y = sqrt(pal.ultra_low.num^2 - k1x.^2);
k2y = sqrt(pal.ultra_high.num^2 - k2x.^2);
kay = sqrt(pal.audio.num^2 - kax.^2);

u1 = sinc(k1x .* src.radius / pi);
u2 = sinc(k2x .* src.radius / pi);

y = 0.01;

I = (exp(1i*(k2y - conj(k1y)).*y) + exp(1i*kay.*y)) ...
    ./ 1i ./ (conj(k1y) - k2y - kay) ...
    + (exp(1i*(k2y - conj(k1y)).*y) - exp(1i*kay.*y)) ...
    ./ 1i ./ (k2y - conj(k1y) - kay);

% int = I .* conj(u1) .* u2 ./ conj(k1y) ./ k2y;
int = conj(u1) ./ conj(k1y) + 0*kax;

fig = Figure;
pcolor(kax/real(pal.audio.num), k1x/real(pal.ultra.num), (abs(int)));
fig.Init
xlabel('kax')
ylabel('k1x')