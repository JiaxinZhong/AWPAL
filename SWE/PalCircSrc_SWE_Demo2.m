%% close all
clear all

prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'quadratic', 'order', 1);
src = CircSrc('radius', .12, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

fp = Point3D('r', 1, ...
    'theta', linspace(0, pi/2, 6e1), ...
    'phi', permute([0; pi], [3, 2, 1]));
fp.Sph2Cart();
 
prs = PalCircSrc_SWE(pal, fp, 'la_max', 70, 'is_farfield', true);
spl = PrsToSpl(prs);

angle = [-flip(fp.theta(:)); fp.theta(:)];
spl_show = [flip(spl(:,:,2)), spl(:,:,1)];
% normalization
spl_show = spl_show(:) - max(spl_show(:));

%% Directivity
fig = Figure;
plot(angle/pi*180, spl_show);
xlabel('z (m)')
ylabel('SPL (dB)')

