clear all

f = 40e3;
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);
a = 0.025;

rho = 0;
z = linspace(0, .7, 1e3).';

prs_cwe = CircSrc_CWE((k), a, rho, z, 'int_num', 4e2) *1.21*343/(pi*a.^2)*2e-4;
spl_cwe = PrsToSpl(prs_cwe);

prs_exact = CircPiston_Axis_Exact(k, a, z)/(pi*a.^2)*2e-4;
spl_exact = PrsToSpl(prs_exact);

fig = Figure;
plot(z, spl_cwe);
hold on
plot(z, spl_exact, '--')