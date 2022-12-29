%% Test the convergence of the radial components

prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'quadratic', 'order', 1);
src = CircSrc('radius', .2, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

n = (0:2:200).';
r = 100;

int_num = (2:10:4e2);

int = 0 * n .* int_num;

for i = 1:length(int_num)
    int(:,i) = CircSrc_SWE_RadialInt(pal.src_high, n, 'j',  0, src.radius, r, 'int_num', int_num(i));
end

figure;
subplot(211)
plot(int_num, real(int));
subplot(212);
plot(int_num, imag(int));