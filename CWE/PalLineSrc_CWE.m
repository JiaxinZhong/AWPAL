% =========================================================================
% INTRO
%   - Calcualte the radiation generated by a PAL field using the CWE method 
% -------------------------------------------------------------------------
% INPUT
%   - pal: infor about the pal
%   - fp: field points
% DIMENSION
%   - 1: fp.rho
%   - 2: fp.phi
%   - 3: ma
% =========================================================================
function prs = PalLineSrc_CWE(pal, fp, varargin)

    ip = inputParser;
    ip.addParameter('ma_max', 3e1, @(x)validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
    ip.addParameter('m1_max', ceil(1.2*real(pal.ultra_low.num * pal.src_low.radius)));
    ip.addParameter('src_norm', 'vel', @(x)any(validatestring(x, {'vel', 'surf_vel', 'prs'})));
    ip.addParameter('eqn', 'Westervelt', @(x)any(validatestring(x, {'Westervelt', 'WesterveltCorrection'})));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    ma = permute((-ip.ma_max:ip.ma_max).', [3,2,1]);

    %% calcualte the radial component
    R = PalLineSrc_CWE_Radial(...
        pal, fp.rho, ip.ma_max, ip.m1_max, ...
        'is_farfield', ip.is_farfield);
    
    %% calculate the results 
    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    prs = beta * pi / (2*1i*rho0*c0^2) * sum(R .* exp(1i*ma.*fp.phi), 3);

    %% include the local effects
    switch ip.eqn
        case 'WesterveltCorrection'
            % [prs1, vel_rho1, vel_phi1] = LineSrc_CWE(k1, a, ...
                % fp.rho, fp.phi, ...
                % 'is_cal_vel', true, ...
                % 'src_norm', ip.src_norm);
            % [prs2, vel_rho2, vel_phi2] = LineSrc_CWE(k2, a, ...
                % fp.rho, fp.phi, ...
                % 'is_cal_vel', true, ...
                % 'src_norm', ip.src_norm);
            % lag = rho0/2.*(conj(vel_rho1) .* vel_rho2 + conj(vel_phi1) .* vel_phi2) ...
                % - (real(k1)/real(k2) + real(k2)/real(k1) - 1) ...
                % .* conj(prs1) .* prs2 / (2 * rho0 * c0^2);
            % prs = prs - lag;
    end
end
