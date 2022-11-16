% =========================================================================
% INTRO
% -------------------------------------------------------------------------
% INPUT
% DIMENSION
%   - 1: rho 
%   - 2: 1 (holder for phi)
%   - 3: ma
%   - 4: m1
%   - 5: rho_vsrc
% =========================================================================
function R = PalLineSrc_CWE_RadialJ(...
    pal, ma_max, m1_max, rho1, rho2, rho, varargin)

    validateattributes(rho1, {'numeric'}, {'column'});
    validateattributes(rho2, {'numeric'}, {'column'});
    validateattributes(rho, {'numeric'}, {'column'});

    ip = inputParser;
    % number of points for the numerical integration
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    parse(ip, varargin{:});
    ip = ip.Results;

    k1 = pal.ultra_low.num;
    k2 = pal.ultra_high.num;
    ka = pal.audio.num;
    a = pal.src_low.radius;
    
    m1 = permute((-m1_max:m1_max).', [4,2,3,1]);
    ma = permute((-ma_max:ma_max).', [3,2,1]);
%     m2 = m1+ma;

    [rho_vsrc, weight] = GaussLegendreQuadParam(ip.int_num, rho1.', rho2.', 'dim', 5);
    rho_vsrc = permute(rho_vsrc, [2, 1, 3, 4, 5]);
    weight = permute(weight, [2, 1, 3, 4, 5]);

    % dim: rho_vsrc => 1 => m1
    R_m1 = LineSrc_CWE_Radial(...
        pal.src_low, ...
        permute(rho_vsrc, [5,2,3,4,1]), ...
        m1_max,...
        'int_num', 5e2);
    % dim: 1 => 1 => 1 => m1 => rho_vsrc
    R_m1 = permute(permute(R_m1, [5,2,3,4,1]), [1, 2, 4, 3, 5]);
    
    R_m2 = 0 * rho_vsrc .* m1 .* ma;
    % dim: 1 => 1 => (0:m1_max+ma_max) => 1 => rho_vsrc
    R_m2_buf = permute(...
        LineSrc_CWE_Radial(...
        pal.src_high, ...
        permute(rho_vsrc, [5,2,3,4,1]), ...
        m1_max + ma_max, ...
        'int_num', 5e2), [5,2,3,4,1]);
    for i = -ma_max : ma_max
        for j = -m1_max : m1_max
            R_m2(1,1,i+ma_max+1,j+m1_max+1,:) ...
                = R_m2_buf(1,1,(i+j)+1+m1_max+ma_max,1,:);
        end
    end

    % dim: ma => 1 => 1 => 1 => rho_vsrc
    J = BesselJLog(ma(:), ka .* rho_vsrc);
    J = permute(J, [3, 2, 1, 4, 5]);
    if ip.is_farfield
        H = log(sqrt(2./pi./(ka.*rho)) .* exp(1i * (ka*rho - ma.*pi/2 - pi/4)));
    else
        % dim: ma => rho
        H = HankelHLog(ma(:), ka .* rho.');
        H = permute(permute(H, [3,2,1]), [2,1,3]);
    end

    JH = exp(J + H);
    JH(isinf(JH)) = 0;
    R = sum(conj(R_m1) .* R_m2 .* JH .* ka^2 .* rho_vsrc .* weight, [4, 5]);
end
