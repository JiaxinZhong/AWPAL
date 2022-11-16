clear all
f = 40e3;
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);

a = 0.3;
lower_limit = 0;
upper_limit = a;

M = 200;
m = (-M:M).';

is_log = false; 

tic
int_matlab = LineSrc_CWE_IntJ(k, a, lower_limit, upper_limit, M, ...
    'profile', 'steerable', 'int_method', 'matlab', ...
    'is_log', is_log);
toc

int_num = (2:200);
int = int_num .* m .* 0;
for i = 1:length(int_num)
    tic 
    fprintf("Processing %d of %d...\n", i, length(int_num));
    int(:,i) = LineSrc_CWE_IntJ(k, a, lower_limit, upper_limit, M, ...
        'profile', 'steerable', 'int_method', 'gauss', ...
        'int_num', int_num(i), ...
        'is_log',  is_log);
    toc
end

[~, rel_err] = Error(int_matlab, int(:,end), 'input_is_log', is_log);

fig = Figure;
plot(int_num, real(int));
hold on
plot(int_num, imag(int))
legend('Real', 'Imag')
xlabel('int num')

% fig = Figure;
% pcolor(int_num, m, real(int));
% fig.Init;
% title('Real part');
% xlabel('int num');
% ylabel('m');

% fig = Figure;
% pcolor(int_num, m, imag(int));
% fig.Init;
% title('Imaginary part');
% xlabel('int num');
% ylabel('m');

% compare to the exact value
figure;
subplot(211)
plot(m, real(int(:,end)))
hold on
plot(m, real(int_matlab), '--')
title('Real')
subplot(212)
plot(m, imag(int(:, end)))
hold on
plot(m, imag(int_matlab),'--')

figure
title('Relative error');
plot(m, real(rel_err));
hold on
plot(m, imag(rel_err));

