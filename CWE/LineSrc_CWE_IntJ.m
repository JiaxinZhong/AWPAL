% =========================================================================
% INTRO
%   - Calculate the radial part for a baffled line source in a 2D problem
%   - Evaluati int_0^a bar{u}(rho_s) J_n(k rho_<) H_n(k rho_>) k drho_src
% -------------------------------------------------------------------------
% INPUT
%   - k: wavenumber
%   - a: half-width (radius)
%   - M: maximum order
% -------------------------------------------------------------------------
% OUTPUT
%   - dimension 1: order (0:M).'
%   - dimension 2 to N: lower_limit .* upper_limit
% =========================================================================
function int = LineSrc_CWE_IntJ(k, a, lower_limit, upper_limit, varargin)

    ip = inputParser;
    ip.addParameter('int_num', nan);
    ip.addParameter('profile', 'uniform', @(x)any(validatestring(x, {'uniform', 'steerable'})));
    % ip.addParameter('focus_dist', 0.2);
    ip.addParameter('steer_angle', pi/4);
    % ip.addParameter('is_farfield', 0);
    % ip.addParameter('is_discrete', 0);
    % ip.addParameter('discrete_size', 0.01);
    % ip.addParameter('is_effective', 0);
    % ip.addParameter('effective_size', 0.01);
    ip.addParameter('is_cal_prime', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('int_method', 'gauss');
    ip.addParameter('is_log', true);
    ip.addParameter('order', nan);
    ip.addParameter('max_order', nan);
    % normalization 
    parse(ip, varargin{:});
    ip = ip.Results;

    % to ensure the convergence at 40 kHz
    if isnan(ip.int_num)
        if a < 0.2
            ip.int_num = 1e2;
        elseif a < 0.3
            ip.int_num = 1.5e2;
        else
            ip.int_num = 2e2;
        end
    end
    m = ip.order;
    if ~isnan(ip.max_order)
        m = (-ip.max_order:ip.max_order).';
    end

    switch ip.int_method
        case 'matlab'
            int = integral(@(rho_src) Integrand(m, k, a, rho_src, ip.profile, ip.steer_angle), ...
                lower_limit, upper_limit, 'ArrayValue', true);
            if ip.is_log
                int = log(int);
            end
        case 'gauss'
            % [rho_src, weight] = GaussLegendreQuadParam(ip.int_num, 0, a, 'dim', 2);
            % u = LineSrc_CWE_Profile(m, k, a, rho_src, ip.profile, ...
                % 'steer_angle', ip.steer_angle);
            % J_log = BesselJLog(m, k.*rho_src);
            % int = sum(u .* exp(J_log) .* k .* weight, 2);
            int = GaussLegendreQuad(@(rho_src) ...
                IntegrandLog(m, k, a, rho_src, ip.profile, ip.steer_angle), ...
                0, a, 'int_num', ip.int_num, 'is_log', true);
            if ~ip.is_log
                int = exp(int);
            end
        otherwise
            error('Wrong integration method!');
    end
end

function int = Integrand(m, k, a, rho_src, profile, steer_angle)
    u = LineSrc_CWE_Profile(m, k, a, rho_src, profile, ...
        'steer_angle', steer_angle); 
    int = u .* BesselJ_Matlab(m, k*rho_src) .* k;
end

function int = IntegrandLog(m, k, a, rho_src, profile, steer_angle)
    u = LineSrc_CWE_Profile(m, k, a, rho_src, profile, 'steer_angle', steer_angle); 
    int = log(u) +  BesselJLog(m, k*rho_src) + log(k);
end
