% clear all
c0 = 343;
fu = 40e3;
fa = 4e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;% + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1)*1;
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2)*1;

a = 0.2;

m1 = 50;
m1 = (-2e2:2e2).';
ma = 0;

int_num = (5e2:5.01e2);

int_method = 'deepest';
% int_method = 'direct_gauss';

int = ma .* m1 .* int_num * 0;
for i = 1:length(int_num)
    tic
    fprintf("Processing %d of %d ...\n", i, length(int_num));
    int(:,i) = PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma, 'int_num', int_num(i),...
        'int_method', int_method, 'is_log', true);
    toc
end

tic
int0 = log(PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma, ...
    'int_method', 'direct'));
toc

err = int(:,end) - int0;
rel_err = real(err)./real(int0) + 1i*imag(err)./imag(int0);
if numel(m1) * numel(ma) == 1
    fprintf("Calculated using MATLAB: %g%+gi\n", real(int0), imag(int0));
    fprintf("Calculated using %s: %g%+gi\n", int_method, real(int(end)), imag(int(end)));
    fprintf("Error: %g%+gi; Relative error: %g%+gi\n", real(err), imag(err), real(rel_err), imag(rel_err))
    % rho_vsrc = logspace(-1, 3, 1e2).';
end

% radial0 = conj(besselh(m1, 1, k1*rho_vsrc, 1)) ...
%     .* besselh(m2, 1, k2*rho_vsrc, 1) ...
%     .* besselh(ma, 1, ka*rho_vsrc, 1) ...
%     .* exp(-(imag(k1)+imag(k2)).*rho_vsrc) ...
%     .* ka^2 .* rho_vsrc;
% figure;
% semilogx(rho_vsrc, real(radial0));
% hold on
% plot(rho_vsrc, imag(radial0));

if numel(m1) * numel(ma) == 1
    figure;
    plot(int_num, real(int));
    hold on
    plot(int_num, imag(int), '--');
    title(sprintf('Method %s', int_method))
elseif numel(m1) > 1
    fig = Figure;
    pcolor(int_num, m1, log10(abs(real(int))));
    fig.Init;
    xlabel('int num')
    ylabel('m1')

    fig = Figure;
    pcolor(int_num, m1, log10(abs(imag(int))));
    fig.Init;
    xlabel('int num')
    ylabel('m1')

    figure;
    plot(m1, rel_err)
elseif numel(ma) > 1
    fig = Figure;
    pcolor(int_num, ma, log10(abs(real(int))));
    fig.Init;
    xlabel('int num')
    ylabel('ma')

    fig = Figure;
    pcolor(int_num, ma, log10(abs(imag(int))));
    fig.Init;
    xlabel('int num')
    ylabel('ma')

    figure;
    subplot(311);
    plot(ma, real(int(:,end)));
    hold on
    plot(ma, imag(int(:,end)));

    subplot(312)
    plot(ma, real(int0));
    hold on
    plot(ma, imag(int0));

    subplot(313)
    plot(ma, real(rel_err))
    hold on
    plot(ma, imag(rel_err));
    legend('real', 'imag')
    title('Relative error');
end
