clear all
% prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'cosine', 'order', 1);
% prf = SrcProfile('name', 'steerable', 'phi', 15/180*pi);
src = LineSrc('radius', 0.1, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

%% The field points 
phi = linspace(-pi/2, pi/2, 2e2).';

%% main procedure
tic
dir_Westervelt = PalLineSrc_CDM(pal, phi, 'type', 'Westervelt', 'is_norm_dB', true);
dir_direct = PalLineSrc_CDM(pal, phi, 'type', 'direct', 'is_norm_dB', true);
dir_modified = PalLineSrc_CDM(pal, phi, 'type', 'modified', 'is_norm_dB', true);
toc

%% plot results
fig = Figure; 
plot(phi/pi*180, dir_Westervelt)
hold on
plot(phi/pi*180, dir_direct)
plot(phi/pi*180, dir_modified)
legend({'Westervelt directivity', 'Direct CDM', 'Modified CDM'}, 'location', 'south')
xlabel("Angle (\circ)")
ylabel("Directivity (dB)")
xlim([-90,90])
fig.Init;