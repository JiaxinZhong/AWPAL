% INTRO
%   - Calculate the on-axial sound field
% INTPUT
% OUTPUT
%   - prs: audio sound pressure without local effects; 
%   - prs_tot: audio sound pressure with local effects
%   - all pressure normalized to p0^2, where p0 = rho0v0c0
% DIMENSION
%   - fp.z => rho_vsrc => z_vsrc

function [prs, prs_tot] = PalDIM3D_CircSrc_Axis(pal, fp, varargin)

    ip = inputParser();
    % integration range in rho direction
    ip.addParameter('rho_vsrc_max', 5);
    % integration range in z direction
    ip.addParameter('z_vsrc', [-4, 5]);
    % integration number for the virtual source
    ip.addParameter('rho_int_num', 1.9e2);
    ip.addParameter('z_int_num', 2e2);
    % integration number for the ultrasound field
    ip.addParameter('ultra_int_num', 2e2);
    % include local effects
    ip.addParameter('is_incl_local', false);
    % integration method
    ip.addParameter('int_method', 'Gauss');
    ip.parse(varargin{:});
    ip = ip.Results;

    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    pal.srca.angfreq = 2*pi*pal.srca.wav.freq;
    pal.srca.wav.num = pal.srca.angfreq / c0;

    switch ip.int_method
        case 'Gauss'
            prs = 0;

            for i = 1:length(ip.z_vsrc) - 1
                t0 = tic;
                fprintf('Calculating over z_vsrc from %g m to %g m...', ip.z_vsrc(i), ip.z_vsrc(i+1));

                % virtual source points
                [vsrc.rho, weight_rho_vsrc] = GaussLegendreQuadParam(ip.rho_int_num, 0, ip.rho_vsrc_max, 'dim', 2);
                [vsrc.z, weight_z_vsrc] = GaussLegendreQuadParam(ip.z_int_num, ip.z_vsrc(i), ip.z_vsrc(i+1), 'dim', 3);
                vsrc.x = vsrc.rho;
                vsrc.y = 0;
        
            % calculate ultrasound pressure field
%             prs1 = DIM3D_CircSrc(pal.src1, vsrc, 'int_num', ip.ultra_int_num);
%             prs2 = DIM3D_CircSrc(pal.src2, vsrc, 'int_num', ip.ultra_int_num);
        
            % audio sound prssure
%             dist = sqrt(vsrc.rho.^2 + (fp.z - vsrc.z).^2);
%             prs = - beta .* pal.srca.angfreq^2 / (2 * rho0 * c0^4) ...
%                 .* sum(sum(conj(prs1) .* prs2 .* exp(1i * pal.srca.wav.num .* dist) ./ dist ...
%                 .* vsrc.rho .* weight_z_vsrc, 3) .* weight_rho_vsrc, 2);
                prs = prs - beta .* pal.srca.angfreq^2 / (2 * rho0 * c0^4) ...
                    .* sum(sum(Integrand(pal, fp, vsrc.rho, vsrc.z, ip.ultra_int_num) ...
                    .* weight_z_vsrc, 3) .* weight_rho_vsrc, 2);

                fprintf(' Finished. Elpased time: %g s.\n', toc(t0));
            end

        case 'builtin'
            prs = fp.z * 0;
            for i = 1:length(fp.z)
                fprintf("Processing %d of %d, fp.z = %g m...\n", i, length(fp.z), fp.z(i));
                fp_now.x = 0;
                fp_now.y = 0;
                fp_now.z = fp.z(i);
                prs(i) = - beta .* pal.srca.angfreq^2 / (2 * rho0 * c0^4) ...
                    .* quad2d(@(rho_vsrc, z_vsrc) Integrand(pal, fp_now, rho_vsrc, z_vsrc, ip.ultra_int_num), ...
                    0, ip.rho_vsrc_max, ip.z_vsrc_min, ip.z_vsrc_max);
            end
    end
    
    %% 
    if ~ip.is_incl_local
        prs_tot = nan;
        return
    end
    [prs1, vel1] = DIM3D_CircSrc(pal.src1, fp, 'int_num', ip.ultra_int_num, 'is_cal_vel', true);
    [prs2, vel2] = DIM3D_CircSrc(pal.src2, fp, 'int_num', ip.ultra_int_num, 'is_cal_vel', true);

    vel1.z = vel1.z / (rho0 * c0);
    vel2.z = vel2.z / (rho0 * c0);
    prs_tot = prs - (rho0/2 * conj(vel1.z) .* vel2.z ...
        - (pal.src1.wav.freq / pal.src2.wav.freq + pal.src2.wav.freq / pal.src1.wav.freq - 1) ...
        .* conj(prs1) .* prs2 ./ (2 * rho0 * c0^2));

end


function int = Integrand(pal, fp, rho_vsrc, z_vsrc, int_num)
    vsrc.x = rho_vsrc;
    vsrc.y = 0;
    vsrc.z = z_vsrc;
    vsrc.rho = rho_vsrc;
    % calculate ultrasound pressure field
    prs1 = DIM3D_CircSrc(pal.src1, vsrc, 'int_num', int_num);
    prs2 = DIM3D_CircSrc(pal.src2, vsrc, 'int_num', int_num);

    % audio sound prssure
    dist = sqrt(vsrc.rho.^2 + (fp.z - vsrc.z).^2);
    int = conj(prs1) .* prs2 .* exp(1i * pal.srca.wav.num .* dist) ./ dist .* vsrc.rho;
end