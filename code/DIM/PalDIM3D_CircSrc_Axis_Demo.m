clear all

%% parameters
pal.src1.wav.freq = 39.75e3;
pal.src1.wav.num = 2*pi*pal.src1.wav.freq / 343 + 1i*AbsorpAttenCoef(pal.src1.wav.freq, 'temperature', 20, 'humidity', 70);
pal.src1.r = 0.1;

pal.src2.wav.freq = 40.25e3;
pal.src2.wav.num = 2*pi*pal.src2.wav.freq / 343 + 1i*AbsorpAttenCoef(pal.src2.wav.freq, 'temperature', 20, 'humidity', 70);
pal.src2.r = 0.1;

% uniform profile
pal.src1.prf.val = @(rs) 1;
pal.src2.prf.val = @(rs) 1;

% focusing profile
% pal.src1.prf.focal_dist = 0.2; % focal distance
% pal.src1.prf.val = @(rs) exp(-1i * real(pal.src1.wav.num) * sqrt(rs.^2 + pal.src1.prf.focal_dist^2));
% pal.src2.prf.focal_dist = 0.2; % focal distance
% pal.src2.prf.val = @(rs) exp(-1i * real(pal.src2.wav.num) * sqrt(rs.^2 + pal.src2.prf.focal_dist^2));

pal.srca.wav.freq = pal.src2.wav.freq - pal.src1.wav.freq;

% field points
% x and y coordinates must be zero
fp.x = 0; fp.y = 0; 
% fp.z must lie in the first dimension, i.e., column vector
fp.z = logspace(-4, 1, 7e2).';

%% main function
v0 = 0.12;
p0 = 1.21 * 343 * v0;
% include local effects
is_incl_local = true;

profile on -memory
[prs, prs_tot] = PalDIM3D_CircSrc_Axis(pal, fp, ...
    'z_int_num', 300, 'rho_int_num', 200, ...
    'ultra_int_num', 200, ...
    'z_vsrc', [-15, -10,-5,-3,-1, 0, 0.1, 0.3, 0.5, 1, 3,5,10, 15], ... % integration ranges for z coordinate of virtual sources
    'is_incl_local', is_incl_local, ...
    'int_method', 'Gauss');
profile viewer

prs = prs * p0^2;
prs_tot = prs_tot * p0^2;
spl = 20*log10(abs(prs)/20e-6/sqrt(2));
spl_tot = 20*log10(abs(prs_tot)/20e-6/sqrt(2));

%% save data
save('DIM/data/PalDIM3D_CircSrc_Axis_Demo_.mat')

%%
C = MyColor('BrewerSet1');

figure;
semilogx(fp.z, spl, 'linewidth', 2, 'color', C(1,:));
hold on
if is_incl_local
    plot(fp.z, spl_tot, 'linewidth', 2, 'color', C(2,:));
end

xlabel('Axial distance (m)')
ylabel('SPL (dB)');
xlim([1e-4, 1e1])
ylim([20, 80])
leg = legend({'W/o local', 'W/ local'});
leg.Position = [0.547, 0.2306, 0.2423, 0.1281];
set(gca, 'linewidth', 1.5)
set(gca, 'fontsize', 20);
set(gca, 'xtick', 10.^(-4:1))
set(gca, 'fontname', 'times new roman')

% Export the figure
% exportgraphics(gcf, 'DIM/fig/PalDIM3D_CircSrc_Axis_Demo_.pdf', 'ContentType','vector') 
% exportgraphics(gcf, 'DIM/fig/PalDIM3D_CircSrc_Axis_Demo_.png', 'resolution', 300) 