% INTRO
%   - Plot 2D sound fields
clear all

%% parameters settings
src.wav.freq = 40e3; % frequency
src.wav.num = 2*pi*src.wav.freq/343 + 1i*AbsorpAttenCoef(src.wav.freq, 'temperature', 20, 'humidity', 70); % wavenumber
src.r = 0.1; % radius of the source
src.prf.phi_m = 0; % azimuthal order of the profile
% uniform profile
src.prf.name = 'uniform';
src.prf.val = @(rs) 1;

% field points
fp.r = logspace(-4, 2, 2e2).';
fp.theta = 0;
fp.phi = 0;
[fp.x, fp.y, fp.z] = Sph2Cart(fp.r, fp.theta, fp.phi);

%% main function
prs_SWE = SWE3D(src, fp);
[prs_ClosedForm, vel_ClosedForm] = ClosedForm_CircPiston(src, fp);
[prs_DIM, vel_DIM] = DIM3D(src, fp, 'int_num', 100, 'int_coord', 'polar', 'is_cal_vel', true);

prs_err_SWE = log10(abs((prs_SWE - prs_ClosedForm) ./ prs_ClosedForm));
prs_err_DIM = log10(abs((prs_DIM - prs_ClosedForm) ./ prs_ClosedForm));

vel_err_DIM.x = log10(abs((vel_DIM.x - vel_ClosedForm.x)));
vel_err_DIM.y = log10(abs((vel_DIM.y - vel_ClosedForm.y)));
vel_err_DIM.z = log10(abs((vel_DIM.z - vel_ClosedForm.z) ./ vel_ClosedForm.z));

%% plot results for sound pressure
figure;
subplot(211);
semilogx(fp.r, abs(prs_ClosedForm), 'linewidth', 2);
hold on
semilogx(fp.r, abs(prs_SWE), '--', 'linewidth', 2)
semilogx(fp.r, abs(prs_DIM), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|p/\rho_0c_0v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')

subplot(212)
semilogx(fp.r, prs_err_SWE, 'linewidth', 2)
hold on
plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')



%% plot results for velocity x
figure;
subplot(221);
semilogx(fp.r, real(vel_ClosedForm.x), 'linewidth', 2);
hold on
semilogx(fp.r, real(vel_DIM.x), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')

subplot(222);
semilogx(fp.r, abs(vel_ClosedForm.x), 'linewidth', 2);
hold on
semilogx(fp.r, abs(vel_DIM.x), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')


subplot(223)
semilogx(fp.r, vel_err_DIM.x , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')


subplot(224)
semilogx(fp.r, vel_err_DIM.x , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')



%% plot results for velocity y
figure;
subplot(221);
semilogx(fp.r, real(vel_ClosedForm.y), 'linewidth', 2);
hold on
semilogx(fp.r, real(vel_DIM.y), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')

subplot(222);
semilogx(fp.r, abs(vel_ClosedForm.x), 'linewidth', 2);
hold on
semilogx(fp.r, abs(vel_DIM.y), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')


subplot(223)
semilogx(fp.r, vel_err_DIM.y , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')


subplot(224)
semilogx(fp.r, vel_err_DIM.y , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')



%% plot results for velocity z
figure;
subplot(221);
semilogx(fp.r, real(vel_ClosedForm.z), 'linewidth', 2);
hold on
semilogx(fp.r, real(vel_DIM.z), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')

subplot(222);
semilogx(fp.r, abs(vel_ClosedForm.z), 'linewidth', 2);
hold on
semilogx(fp.r, abs(vel_DIM.z), '-.', 'linewidth', 2)

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
ylabel('$|vel_z/v_0|$', 'Interpreter','latex')
xlabel('Axial distance (m)')
set(gca, 'fontsize', 18)
set(gca, 'FontName', 'times new roman')


subplot(223)
semilogx(fp.r, vel_err_DIM.z , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')


subplot(224)
semilogx(fp.r, vel_err_DIM.z , 'linewidth', 2)
% hold on
% plot(fp.r, prs_err_DIM, 'linewidth', 2);

xlim([1e-4, 1e2])
set(gca, 'xtick', 10.^(-4:2))
% legend({'SWE', 'DIM'})
ylabel('log10(Rel. error)')
set(gca, 'fontsize', 18)
xlabel('Axial distance (m)')
set(gca, 'FontName', 'times new roman')
