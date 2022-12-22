% Calculate the directivity
clear all

prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'steerable', 'theta', pi/3);
src = CircSrc('radius', .1, 'prf', prf, 'freq', 40e3);

fp = Point3D('r', 1, ...
    'theta', linspace(0, pi/2, 4e2), ...
    'phi', permute([0; pi], [3, 2, 1]));
fp.Sph2Cart();

dir_analytic = src.CalDirectivity(fp.theta, fp.phi);
dir_analytic_dB = PrsToSpl(dir_analytic);
dir_analytic_dB = dir_analytic_dB - max(dir_analytic_dB(:));

dir_swe = CircSrc_SWE(src, fp, 'is_farfield', true);
dir_swe_dB = PrsToSpl(dir_swe);
dir_swe_dB = dir_swe_dB- max(dir_swe_dB(:));

angle = [-flip(fp.theta), fp.theta];
dir_analytic_dB_show = [flip(dir_analytic_dB(1,:,2)), dir_analytic_dB(1,:,1)];
dir_swe_dB_show = [flip(dir_swe_dB(1,:,2)), dir_swe_dB(1,:,1)];

fig = Figure;
plot(angle/pi*180, dir_analytic_dB_show);
hold on
plot(angle/pi*180, dir_swe_dB_show, '--');
fig.Init;
legend({'Analytical results', 'SWE'})