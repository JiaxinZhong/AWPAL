% close all
clear all

prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'quadratic', 'order', 1);
src = CircSrc('radius', .05, 'prf', prf);
pal = PalSrc('audio_freq', 1e3, 'ultra_freq', 40e3, 'src', src);

fp = Point3D('r', linspace(0, 3.4, 5e1).', ...
    'theta', linspace(0, pi/2, 6e1), ...
    'phi', permute([0; pi], [3, 2, 1]));
fp.Sph2Cart();
 
prs = PalCircSrc_SWE(pal, fp, 'la_max', 40);
spl = PrsToSpl(prs);

fp_z_show = [flipud(fp.z(2:end,:)); fp.z];
fp_x_show = [flipud(fp.x(2:end,:,2)); fp.x(:,:,1)];
spl_show = [flipud(spl(2:end,:,2)); spl(:,:,1)];

%% 2D field
fig = Figure;
% pcolor(fp.z(:,:,1), fp.x(:,:,1), spl(:,:,1));
pcolor(fp_z_show, fp_x_show, spl_show);
fig.Init;
ylim([-1.5,1.5]);
xlim([0,3])

%% Axial field
fig_axis = Figure;
plot(fp.r, spl(:,1,1));
xlabel('z (m)')
ylabel('SPL (dB)')