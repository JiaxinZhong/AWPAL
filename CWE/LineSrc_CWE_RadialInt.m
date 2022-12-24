% =========================================================================
% INTRO
%   - Calculate the radial part for a baffled line source in a 2D problem
%   - Evaluati int_0^a bar{u}(rho_s) J_n(k rho_<) H_n(k rho_>) k drho_src
% -------------------------------------------------------------------------
% INPUT
%   - k: wavenumber
%   - a: half-width (radius)
%   - cyl: cylinder function type 'J' or 'H'
%   - rho1: the lower limit of the integral
%   - rho2: the upeer limit of the integral
%   - rho: the argument of the function outside the integral
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: m
%   - 2~N: rho1 .* rho2 .* rho
% =========================================================================
function int = LineSrc_CWE_RadialInt(src, m, cyl, rho1, rho2, rho, varargin)

    validateattributes(m, {'numeric'}, {'column'});
    validatestring(cyl, {'J', 'H'});
    validateattributes(rho1, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    validateattributes(rho2, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    validateattributes(rho, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    
    ip = inputParser;
    ip.addParameter('int_num', []);
    ip.addParameter('is_log', false);
    ip.addParameter('is_farfield', false)
    parse(ip, varargin{:});
    ip = ip.Results;

    k = src.wav.num;
    a = src.radius;

    % to ensure the convergence at 40 kHz
    if isempty(ip.int_num)
        if a < 0.2
            ip.int_num = 1e2;
        elseif a < 0.3
            ip.int_num = 1.5e2;
        else
            ip.int_num = 2e2;
        end
    end

    m_abs = abs(m);
    [m_abs_unique, ~, idx] = unique(m_abs);

    [rho_src, weight] = GaussLegendreQuadParam(ip.int_num, rho1, rho2, 'dim', 7);

%     u = (src.CalProfile(rho_src) + (-1).^m_abs_unique .* src.CalProfile(-rho_src))/2; 
    u = (1i.^(-m_abs_unique) .* src.CalProfile(rho_src) + 1i.^(m_abs_unique) .* src.CalProfile(-rho_src))/2; 
     
    switch cyl
        case 'J'
            if ip.is_farfield
                int = log(u) + log(k) ...
                    + BesselJ(m_abs_unique, k*rho_src, 'is_log', true) ...
                    + HankelH(m_abs_unique, k*rho, 'is_log', true, 'arg_is_large', ip.is_farfield);
                    % + 1i*k*rho + log(sqrt(2./(k*rho.*1i*pi))./1i.^m_abs_unique);
            else
                int = log(u) + log(k) ...
                    + BesselJ(m_abs_unique, k*rho_src, 'is_log', true) ...
                    + HankelH(m_abs_unique, k*rho, 'is_log', true);
            end
        case 'H'
            int = log(u) + log(k) ...
                + BesselJ(m_abs_unique, k*rho, 'is_log', true) ...
                + HankelH(m_abs_unique, k*rho_src, 'is_log', true);
        otherwise
            error('Wrong cylinder functions!')
    end
    int = sum(exp(int) .* weight, 7);
%     int = GaussLegendreQuad(@(rho_src) ...
%         Integrand(src, m_abs_unique, k, cyl, rho_src, rho), ...
%         rho1, rho2, 'int_num', ip.int_num, 'is_log', true, 'dim', 10);
    int = int(idx, :, :, :, :, :, :);
    int(m<0, :, :, :, :, :, :) = int(m<0, :, :, :, :, :, :) .* (-1).^(m(m<0));
%     if ~ip.is_log
%         int = exp(int);
%     end
end

function int = Integrand(src, m, k, cyl, rho_src, rho)
    u =  (src.CalProfile(rho_src) + (-1).^m .* src.CalProfile(-rho_src))/2; 
    switch cyl
        case 'J'
            int = log(u) + log(k) ...
                + BesselJLog(m, k*rho_src) + HankelHLog(m, k*rho);
        case 'H'
            int = log(u) + log(k) ...
                + BesselJLog(m, k*rho) + HankelHLog(m, k*rho_src);
        otherwise
            error('Wrong cylinder functions!')
    end
end
