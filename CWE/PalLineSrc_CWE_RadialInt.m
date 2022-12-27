% =========================================================================
% INTRO
%   - Calculate the integral involved in the radial component for audio 
%       sound when using the CWE
%   - int_rho1^rho2 Jn(krho<) * Hn(krho>) * R1 * R2 * ka^2 * rhov drhov
%       - rho< = min(rho, rhov)
%       - rho> = max(rho, rhov)
% -------------------------------------------------------------------------
% INPUT
%   - cyl: the kernel cylindrical function type 'J' or 'H'
%       - 'J':  Hn(krho) * int_rho1^rho2 ... j_n(k rhov) drhov
%       - 'H':  Jn(krho) * int_rho1^rho2 ... j_n(k rhov) drhov
%   - rho1: the lower limit of the integral
%   - rho2: the upeer limit of the integral
%   - rho: the argument of the function outside the integral
% DIMENSION
%   - 1: rho 
%   - 2: phi (holder)
%   - 3: ma
%   - 4: m1
%   - 5: dummy variable for integration
% =========================================================================
function R = PalLineSrc_CWE_RadialInt(...
    pal, ma_max, m1_max, rho1, rho2, rho, cyl, varargin)

    validatestring(cyl, {'J', 'H'});

    ip = inputParser;
    % number of points for the numerical integration
    ip.addParameter('int_num', 3e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % is_log = true: return the logarithm of the result
    % perform_sum = true: summation the result here 
    ip.addParameter('perform_sum', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    parse(ip, varargin{:});
    ip = ip.Results;

    R = GaussLegendreQuad(@(rhov) ...
        Integrand(pal, cyl, rhov, rho, ma_max, m1_max, ip.is_farfield), ...
        rho1, rho2, 'int_num', ip.int_num, 'is_log', true, 'dim', 5);
    R = exp(R);
    if ip.perform_sum
        R = sum(R, 4);
    end
end

function int = Integrand(pal, cyl, rhov, rho, ma_max, m1_max, is_farfield)
    ma = permute((-ma_max:ma_max).', [3,2,1]);
    m1 = permute((-m1_max:m1_max).', [4,2,3,1]);

    % dim: rho .* rhov -> 1 -> m1
    R1 = LineSrc_CWE_Radial(pal.src_low, rhov(:), m1_max);
    % dim: rho -> rhov -> m1
    R1 = reshape(R1, [size(permute(rhov, [1,5,3,4,2])), 2*m1_max+1]);
    % dim: rho -> 1 -> 1 -> m1 -> rhov
    R1 = permute(permute(R1, [1,5,3,4,2]), [1,2,4,3,5]);
    
    % dim: rho .* rhov -> 1 -> m2
    R2_buf = LineSrc_CWE_Radial(pal.src_high, rhov(:), m1_max + ma_max);
    % dim: rho -> rhov -> m2
    R2_buf = reshape(R2_buf, [size(permute(rhov, [1,5,3,4,2])), 2*(m1_max+ma_max)+1]);
    % dim: rho -> 1 -> m2 -> 1 -> rhov
    R2_buf = permute(R2_buf, [1,5,3,4,2]);
    % dim: rho -> 1 -> ma -> m1 -> rhov
    R2 = 0 * rhov .* m1 .* ma;
    for ma_now = -ma_max : ma_max
        for m1_now = -m1_max : m1_max
            R2(:, 1, ma_now+ma_max+1, m1_now+m1_max+1, :) ...
                = R2_buf(:, 1, (ma_now+m1_now)+1+ma_max+m1_max, 1, :);
        end
    end

    %% the integrand
    switch cyl
        case 'J'
            rhoJ = rhov;
            rhoH = rho;
        case 'H'
            rhoJ = rho;
            rhoH = rhov;
        otherwise
            error('Wrong cylindrical functions!');
    end

    int = BesselJ(ma(:), pal.audio.num*permute(rhoJ, [3,2,1,4,5]), 'is_log', true) ...
        + HankelH(ma(:), pal.audio.num*permute(rhoH, [3,2,1,4,5]), 'is_log', true, 'arg_is_large', is_farfield);
    int = permute(int, [3,2,1,4,5]) + log(conj(R1)) + log(R2) + 2*log(pal.audio.num) + log(rhov);
end
