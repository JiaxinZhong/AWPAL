clear all
j1 = 2*54;
j2 = 2*83;
j3 = 2*85;
m1 = 0;
m2 = 0;
ma = m2-m1;

w = Wigner3j(j1, j2, j3, -m1, m2, -ma)

% hypergeom([-j1-j2-j3-1, -j1+m1, -j3-m3], [-j1-j2-m3, -j2-j3+m1], 1)
if (m1 == 0) && (m2 == 0) && (m3 == 0)
    w0 = Wigner3j000(j1, j2, j3)
end