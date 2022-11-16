% Calculate the single value

clear all

a = 0.3;
prf.name = 'uniform';
src = LineSrc('freq', 40e3, 'radius', a, 'prf', prf);

rho1 = 0.1;
rho2 = a;
rho = 1;

m = (-100:100).';

int = LineSrc_CWE_Int(src, m, 'H', rho1, rho2, rho)