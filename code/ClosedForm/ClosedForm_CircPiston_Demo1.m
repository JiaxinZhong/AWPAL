%
clear all

%% parameters
src.r = 0.1;
src.wav.freq = 40e3;
src.wav.num = 2*pi*src.wav.freq/343 + 1i*AbsorpAttenCoef(src.wav.freq, 'temperature', 20, 'humidity', 70);
fp.z = logspace(-4, 1, 1e3).';

% main process
[prs, vel] = ClosedForm_CircPiston(src, fp);

%%
figure;
subplot(211)
semilogx(fp.z, abs(prs), 'linewidth', 2);

title('Pressure')
ylabel('$|p_0/\rho_0c_0v_0$', 'Interpreter','latex')

xlabel('Axial distance (m)');
xlim([1e-4, 1e1])
set(gca, 'xtick', 10.^(-4:1));
set(gca, 'FontName', 'times new roman');
set(gca, 'fontsize', 18)

subplot(212)
semilogx(fp.z, abs(vel.z), 'linewidth', 2)

title('Velocity in z direction')
ylabel('$|v_z/v_0|$', 'Interpreter','latex')

xlabel('Axial distance (m)');
xlim([1e-4, 1e1])
set(gca, 'xtick', 10.^(-4:1));
set(gca, 'FontName', 'times new roman');
set(gca, 'fontsize', 18)
