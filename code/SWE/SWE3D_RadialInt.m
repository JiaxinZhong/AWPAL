% =========================================================================
% INTRO
%   - Calculate the integral involved in the radial component
%   - int_r1^r2 u(rs) * jn(k r<) * hn(k r>) * rs drs
%       - r< = min(r, rs)
%       - r> = max(r, rs)
% -------------------------------------------------------------------------
% INPUT
%   - n: order
%   - sph: the kernel. spherical function type 'j' or 'h'
%       - 'j':  hn(kr) * int_r1^r2 u(rs) j_n(k rs) k^2 rs drs
%       - 'h':  jn(kr) * int_r1^r2 u(rs) j_n(k rs) k^2 rs drs
%   - r1: the lower limit of the integral
%   - r2: the upeer limit of the integral
%   - r: the argument of the function outside the integral
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: n (order)
%   - 2 ~ N: r1 .* r2 .* r
% =========================================================================
function int = SWE3D_RadialInt(src, n, sph, r1, r2, r, varargin)

    validateattributes(n, {'numeric'}, {'column'});
    validatestring(sph, {'j', 'h'});
    validateattributes(r1, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    validateattributes(r2, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    validateattributes(r, {'numeric'}, {'size', [1, nan, nan, nan, nan]});
    
    ip = inputParser;
    ip.addParameter('int_num', []);
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % normalization 
    parse(ip, varargin{:});
    ip = ip.Results;

    % to ensure the convergence at 40 kHz
    if isempty(ip.int_num)
        if src.r < 0.2
            ip.int_num = 1e2;
        elseif src.r < 0.3
            ip.int_num = 1.5e2;
        else
            ip.int_num = 2e2;
        end
    end
    % ip.int_num = 2e2;

    [n_unique, ~, idx] = unique(n);

    int = GaussLegendreQuad(@(rs) ...
        Integrand(src, n_unique, sph, rs, r, ip.is_farfield), ...
        r1, r2, 'int_num', ip.int_num, 'is_log', true, 'dim', 10);
    int = exp(int(idx, :, :, :, :, :, :));
end

function int = Integrand(src, n, sph, rs, r, is_farfield)
%     u = src.CalProfile(rs);
    u = src.prf.val(rs);
    switch sph
        case 'j'
            rj = rs;
            rh = r;
        case 'h'
            rj = r;
            rh = rs;
        otherwise
            error('Wrong spherical functions!')
    end
    int = log(u) + log(rs) ...
        + SphBesselJ(n, src.wav.num*rj, 'is_log', true) ...
        + SphHankelH(n, src.wav.num*rh, 'is_log', true, 'arg_is_large', is_farfield);
    % int = log(u) + 2*log(src.wav.num) + log(rs) ...
    %     + SphBesselJ(n, src.wav.num*rj, 'is_log', true) ...
    %     + SphHankelH(n, src.wav.num*rh, 'is_log', true, 'arg_is_large', is_farfield);
%     int = u.* (src.wav.num)^2 .* (rs) ...
%         .* exp(SphBesselJ(n, src.wav.num*rj, 'is_log', true) ...
%         + SphHankelH(n, src.wav.num*rh, 'is_log', true, 'arg_is_large', is_farfield));
end
