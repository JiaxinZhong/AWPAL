% one point
clear all
f = 40e3;
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);

a = 0.1;
lower_limit = 0;
upper_limit = a;

m = 6;

is_log = true; 

tic
int_matlab = exp(LineSrc_CWE_IntJ(k, a, lower_limit, upper_limit, ...
    'order', m, ...
    'profile', 'uniform', 'int_method', 'gauss', ...
    'is_log', is_log))
toc
