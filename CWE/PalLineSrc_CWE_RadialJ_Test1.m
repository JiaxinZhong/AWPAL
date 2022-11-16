clear all;

src.radius = 0.1;
src.shape = 'line';
src.prf.name = 'uniform';
pal = PalSrc('audio_freq', 1e3, 'ultra_freq', 40e3, 'src', src);

rho = 1;

ma_max = 100;
m1_max = 200;

R = 0;
% rho_part = [0; 1; 10; 100; 1e3; inf]*src.radius;
rho_part = [1e3; inf]*src.radius;
for i = 1:length(rho_part)-1
    R = R + PalLineSrc_CWE_RadialJ(...
        pal, ma_max, m1_max, rho_part(i), rho_part(i+1), rho, ...
        'is_farfield', true, 'int_num', 2e3);
end

R = R(:);
IsInvalid(R);
