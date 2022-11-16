% variable: m1 and ma
% clear all
c0 = 343;
fu = 40e3;
fa = 4e3;
f1 = fu - fa/2;
f2 = fu + fa/2;
ka = 2*pi*fa/c0;% + 1i*AbsorpAttenCoef(fa);
k1 = 2*pi*f1/c0 + 1i*AbsorpAttenCoef(f1)*1;
k2 = 2*pi*f2/c0 + 1i*AbsorpAttenCoef(f2)*1;

a = 0.1;

ma = permute((-2e2:2e2).', [3,2,1]);
m1 = permute((-1.1e2:1.1e2).', [4,2,3,1]);
% ma = (-30:30).';

int_num = (2e2:2.01e2).';

int_deepest = int_num .* m1 .* ma .*0;
int_gauss = int_num .* m1 .* ma .*0;
int_hybrid = int_num .* m1 .* ma .*0;
for i = 1:length(int_num)
    tic
    fprintf("Processing %d of %d ...\n", i, length(int_num));
    int_deepest(i,:,:,:) = PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma, 'int_num', int_num(i),...
        'int_method', 'deepest', 'is_log', true);
    int_gauss(i,:,:,:) = PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma, 'int_num', int_num(i),...
        'int_method', 'direct_gauss', 'is_log', true);
    int_hybrid(i,:,:,:) = PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma, 'int_num', int_num(i),...
        'int_method', 'hybrid', 'is_log', true);
    toc
end

int0 = ma .* m1 .* 0;
for i = 1:length(ma)
    fprintf("Processing %d of %d...\n", i, length(ma));
    tic
    int0(:,:,i,:) = log(PalLineSrc_CWE_RadialH(k1, k2, ka, a, m1, ma(i), ...
        'int_method', 'matlab'));
    toc
end

ma = permute(ma, [3,2,1]);
m1 = permute(m1, [1,4,3,2]);
int0 = squeeze(int0);
int_deepest = permute(permute(int_deepest, [3,2,1,4]), [1,4,3,2]);
int_gauss = permute(permute(int_gauss, [3,2,1,4]), [1,4,3,2]);
int_hybrid = permute(permute(int_hybrid, [3,2,1,4]), [1,4,3,2]);

[err_deepest, rel_err_deepest] = Error(int0, int_deepest(:,:,end), 'input_is_log', true);
[err_gauss, rel_err_gauss] = Error(int0, int_gauss(:,:,end), 'input_is_log', true);
[err_hybrid, rel_err_hybrid] = Error(int0, int_hybrid(:,:,end), 'input_is_log', true);

caxis_range = [-12,0];

%% convergence of integral
conv_deepest = int_deepest(:,:,end) - int_deepest(:,:,end-1);
conv_deepest = real(conv_deepest)./real(int_deepest(:,:,end)) ...
    + 1i * imag(conv_deepest) ./  imag(int_deepest(:,:,end));
conv_gauss = int_gauss(:,:,end) - int_gauss(:,:,end-1);
conv_gauss = real(conv_gauss)./real(int_gauss(:,:,end)) ...
    + 1i * imag(conv_gauss) ./  imag(int_gauss(:,:,end));

fig = Figure;
pcolor(m1, ma, log10(abs(real(conv_gauss))));
fig.Init;
title(sprintf('Convergence test, log10, real, %s', 'Gauss'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(imag(conv_gauss))));
fig.Init;
title(sprintf('Convergence test, log10, imag, %s', 'Gauss'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(real(conv_deepest))));
fig.Init;
title(sprintf('Convergence test, log10, real, %s', 'Deepest'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(imag(conv_deepest))));
fig.Init;
title(sprintf('Convergence test, log10, imag, %s', 'Deepest'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

%%
fig = Figure;
pcolor(m1, ma, log10(abs(real(rel_err_gauss))));
fig.Init;
title(sprintf('Real relative error, log10, %s', 'Gauss'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(imag(rel_err_gauss))));
fig.Init;
title(sprintf('Imag, relative error, log10, %s\n213', 'Gauss'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(real(rel_err_deepest))));
fig.Init;
title(sprintf('Real relative error, log10, %s', 'Deepest'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(imag(rel_err_deepest))));
fig.Init;
title(sprintf('Imag, relative error, log10, %s\n213', 'Deepest'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(real(rel_err_hybrid))));
fig.Init;
title(sprintf('Real relative error, log10, %s', 'Hybrid'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

fig = Figure;
pcolor(m1, ma, log10(abs(imag(rel_err_hybrid))));
fig.Init;
title(sprintf('Imag, relative error, log10, %s\n213', 'Hybrid'))
xlabel('m1');
ylabel('ma');
caxis(caxis_range);

% save('CWE/data/PalLineSrc_CWE_RadialH_Test2_.mat');
