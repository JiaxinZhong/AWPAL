clear all

f = 40e3;
k = 2*pi*f/343 + 1i*AbsorpAttenCoef(f);
a = 0.025;

rho = 1;
z = 1;
int_num = (2:5e2).';

prs = 0*int_num;
tic
for i = 1:length(int_num)
%     fprintf("Processing %d of %d...\n", i, length(int_num));
%     tic
    prs(i) = CircSrc_CWE((k), a, rho, z, 'int_num', int_num(i));
%     toc
end
toc
spl = prs2spl(prs);

fig = Figure;
plot(int_num, spl);
