clear all
f = 40e3;
max_order = 100;
rho = [0.05,0.01;0.2,0.5];
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);
a = 0.1;

[R, R_prime] = LineSrc_CWE_Radial(max_order, k, a, rho, ...
    'is_cal_prime', true);
