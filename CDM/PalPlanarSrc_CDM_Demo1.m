clear all

%% settings
prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'quadratic', 'order', 1);
prf = SrcProfile('name', 'steerable', 'phi', 0, 'theta', 15/180*pi);

% src = CircSrc('radius', .1, 'prf', prf);
src = RectSrc('prf', prf, 'ax', 0.1, 'ay', 0.05);

pal = PalSrc('audio_freq', 1e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
fp = Point3D('theta', linspace(0, pi/2, 2e2).', 'phi', [0,pi]);

tic
dir_Westervelt = PalPlanarSrc_CDM(pal, fp, 'type', 'Westervelt', 'is_norm_dB', true);
dir_direct = PalPlanarSrc_CDM(pal, fp, 'type', 'direct', 'is_norm_dB', true);
dir_modified =  PalPlanarSrc_CDM(pal, fp, 'type', 'modified', 'is_norm_dB', true);
toc

dir_Westervelt = [flip(dir_Westervelt(:,2)); dir_Westervelt(:,1)];
dir_direct = [flip(dir_direct(:,2)); dir_direct(:,1)];
dir_modified = [flip(dir_modified(:,2)); dir_modified(:,1)];
angle = [-flip(fp.theta); fp.theta]/pi*180;

dir_Westervelt = dir_Westervelt - max(dir_Westervelt(:));
dir_direct = dir_direct - max(dir_direct(:));
dir_modified = dir_modified - max(dir_modified(:));

fig = Figure ;
plot(angle, dir_Westervelt)
hold on
plot(angle, dir_direct)
plot(angle, dir_modified);
legend({'Westervelt directivity', 'Direct CDM', 'Modified CDM'}, 'location', 'south')
xlabel("Angle (\circ)")
ylabel("Directivity (dB)")
xlim([-90,90]);
fig.Init;
