% Test convergence
clear all
c0 = 343;
f = 40e3;
k = 2*pi*f/c0 + 1i*AbsorpAttenCoef(f);
a = 0.01;
% field points
x = 2;
y = 2;
z = linspace(0, 2, 2e2).';

gauss_num = (2:1:100);
prs = 0 * z .* gauss_num;
for i = 1:length(gauss_num)
    fprintf("Processing %d of %d..\n", i, length(gauss_num));
    tic
    prs(:,i) = CircPiston_Direct(k, a, x, y, z, 'gauss_num', gauss_num(i),...
        'profile_type', 'uniform', 'steering_angle', 75/180*pi);
    toc
end
spl = PrsToSpl(prs);

% tolerrance
spl_tol = 0.01;
% 
% % the required Gauss number to meet the SPL tolerrance
num = 0 * z;
spl_diff = diff(spl, 1, 2);
for i = 1:length(z)
    for j = length(gauss_num)-1:-1:1
        if (abs(spl_diff(i,j)) >= spl_tol)
            num(i) = gauss_num(j);
            break
        end
    end
end


%
fig_spl = Figure;
pcolor(gauss_num, z, spl);
fig_spl.Init;
hold on
plot(num, z, 'k')
xlabel('Gauss degree')
ylabel('z (m)')

% save('direct/data/CircPiston_Direct_TestConverg_.mat')