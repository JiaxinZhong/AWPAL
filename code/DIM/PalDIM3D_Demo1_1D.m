% clear all
fn0 = mfilename;

%% parameters
pal.ultra1.wav.freq = 39.75e3;
pal.ultra1.wav.num = 2*pi*pal.ultra1.wav.freq / 343 + 1i*AbsorpAttenCoef(pal.ultra1.wav.freq, 'temperature', 20, 'humidity', 70);
pal.ultra1.r = 0.1;

pal.ultra2.wav.freq = 40.25e3;
pal.ultra2.wav.num = 2*pi*pal.ultra2.wav.freq / 343 + 1i*AbsorpAttenCoef(pal.ultra2.wav.freq, 'temperature', 20, 'humidity', 70);
pal.ultra2.r = 0.1;

% uniform profile
% pal.ultra1.prf.val = @(rs) 1;
% pal.ultra2.prf.val = @(rs) 1;

% focusing profile
pal.ultra1.prf.focal_dist = 0.2; % focal distance
pal.ultra1.prf.val = @(rs) exp(-1i * real(pal.ultra1.wav.num) * sqrt(rs.^2 + pal.ultra1.prf.focal_dist^2));
pal.ultra2.prf.focal_dist = 0.2; % focal distance
pal.ultra2.prf.val = @(rs) exp(-1i * real(pal.ultra2.wav.num) * sqrt(rs.^2 + pal.ultra2.prf.focal_dist^2));

pal.audio.wav.freq = pal.ultra2.wav.freq - pal.ultra1.wav.freq;

% field points
% x and y coordinates must be zero
fp.x = 0; fp.y = 0; 
% fp.z must lie in the first dimension, i.e., column vector
fp.z = logspace(-4, 1, 30).';

%% main function
v0 = 0.12;
p0 = 1.21 * 343 * v0;
% include local effects
is_incl_local = false;

% profile on -memory
prs_tot = fp.z * 0;
for i = 1:length(fp.z)
    fp_now = fp;
    fp_now.z = fp.z(i);
    fprintf('Processing %d of %d.\n', i, length(fp.z));
%     [prs_tot(i), ~, ~] = PalDIM3D(pal, fp_now, ...
%     'ultra_int_num', 100, ...
%     'x_vsrc_int_num', 50, ...
%     'y_vsrc_int_num', 50, ...
%     'z_vsrc_int_num', 100, ...
%     'x_vsrc', [-4, 4], ...
%     'y_vsrc', [-4, 4], ...
%     'z_vsrc', [-10, -5,-3,-1, 0, 1, 3, 5, -10], ... % integration ranges for z coordinate of virtual sources
%     'ultra_int_coord', 'polar', ...
%     'is_incl_local', is_incl_local);
    [prs_tot(i), ~, ~] = PalDIM3D(pal, fp_now, ...
    'ultra_int_num', 150, ...
    'rho_vsrc_int_num', 50, ...
    'phi_vsrc_int_num', 60, ...
    'z_vsrc_int_num', 100, ...
    'rho_vsrc', [0, 1, 2, 4], ...
    'phi_vsrc', [0, 2*pi], ...
    'z_vsrc', [-15, -10, -5,-3,-2,-1, 0, 1, 2,3, 5, 10, 15], ... % integration ranges for z coordinate of virtual sources
    'ultra_int_coord', 'polar', ...
    'int_coord', 'cylindrical', ...
    'is_incl_local', is_incl_local);
end
% profile viewer
prs_tot = prs_tot * p0^2;
spl_tot = 20*log10(abs(prs_tot)/sqrt(2)/20e-6);

%% save data
save(sprintf('DIM/data/%s_.mat', fn0));

%% plot results
figure;
semilogx(fp.z, spl_tot, 'linewidth', 2)

xlabel('Axial distance (m)')
ylabel('SPL (dB)');
xlim([1e-4, 1e1])
ylim([20, 80])
% leg = legend({'W/o local', 'W/ local'});
leg.Position = [0.547, 0.2306, 0.2423, 0.1281];
set(gca, 'linewidth', 1.5)
set(gca, 'fontsize', 20);
set(gca, 'xtick', 10.^(-4:1))
set(gca, 'fontname', 'times new roman');