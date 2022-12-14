clear all

M = 10;
z = 2e1+1i*1;
% z = 0;

% m = (0:M).';
nu0 = 0;
tic
[H, H_prime] = HankelHLog_Rothwell(M, z, 'nu0', nu0, ...
    'is_cal_derivative', true);
H = exp(H);
H_prime = exp(H_prime);
H_sel = H(end);
% J = real(J) + wrapToPi(imag(J));
toc
H_exact = (besselh(M+nu0, z));

err = H_sel - H_exact;
[~, rel_err] = Error(H_exact, H_sel);

fprintf("H_exact = %g %+gi\n", real(H_exact), imag(H_exact))
fprintf("H = %g %+gi\n", real(H_sel), imag(H_sel))
fprintf("Error = %g %+gi\n", real(err), imag(err))
fprintf("Relative error = %g %+gi\n", real(rel_err), imag(rel_err))

