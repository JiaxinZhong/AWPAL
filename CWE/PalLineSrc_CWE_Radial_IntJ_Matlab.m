% =========================================================================
% INTRO
% -------------------------------------------------------------------------
% INPUT
% DIMENSION
%   - 1, 2: rho 
%   - 3: ma
%   - 4: m1
%   - 5: integration
% =========================================================================
function radial = PalLineSrc_CWE_Radial_IntJ_Matlab(...
    k1, k2, ka, a, lower_limit, upper_limit, ...
    varargin)

    validateattributes(k1, {'numeric'}, {'scalar'});
    validateattributes(k2, {'numeric'}, {'scalar'});
    validateattributes(ka, {'numeric'}, {'scalar'});
    validateattributes(lower_limit, {'numeric'}, {'scalar', '>=', 0});

    ip = inputParser;
    % number of points for the numerical integration
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('profile', 'uniform', @(x)any(validatestring(x, {'uniform', 'steerable'})));
    ip.addParameter('steer_angle', [pi/4;pi/2]);
    % ip.addParameter('steer_angle_lower', pi/2);
    % ip.addParameter('steer_angle_upper', pi/2);
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('int_method', 'deepest');
    ip.addParameter('m1', nan);
    ip.addParameter('ma', nan);
    ip.addParameter('ma_max', nan);
    ip.addParameter('m1_max', nan);
    parse(ip, varargin{:});
    ip = ip.Results;

    ma = ip.ma;
    m1 = ip.m1;
    if ~isnan(ip.ma_max) && ~isnan(ip.m1_max)
        ma = permute((-ip.ma_max:ip.ma_max).', [3,2,1]);
        m1 = permute((-ip.m1_max:ip.m1_max).', [4,2,3,1]);
    end

    m1_col = m1 + 0*ma;
    ma_col = ma + 0*m1;
    m2_col = m1_col + ma_col;
    m1_col = m1_col(:);
    m2_col = m2_col(:);
    ma_col = ma_col(:);

    radial = integral(@(rho_vsrc) Fun1(k1, k2, ka, ...
        m1_col, m2_col, ma_col, rho_vsrc), ...
        lower_limit, upper_limit, ...
        'ArrayValued', true, ...
        'RelTol', 1e-5);
    F1 = LineSrc_CWE_IntJ(k1, a, 0, a, 'order', m1_col, 'profile', 'uniform',...
        'is_log', false);
    F2 = LineSrc_CWE_IntJ(k2, a, 0, a, 'order', m2_col, 'profile', 'uniform', ...
        'is_log', false);

    radial = reshape(conj(F1) .* F2 .* radial, size(m1 .* ma));
end

function res = Fun1(k1, k2, ka, m1, m2, ma, rho_vsrc)
    res = conj(HankelH_Matlab(m1, k1*rho_vsrc)) ...
        .* HankelH_Matlab(m2, k2*rho_vsrc) ...
        .* BesselJ_Matlab(ma, ka*rho_vsrc) ...
        .* ka^2 .* rho_vsrc;
end


function res = Fun2(k1, k2, ka, m1, m2, ma, rho_vsrc)
    res = exp(conj(HankelHLog(m1, k1*rho_vsrc)) ...
        + HankelHLog(m2, k2*rho_vsrc) ...
        + BesselJLog(ma, ka*rho_vsrc)) ...
        .* ka^2 .* rho_vsrc;
end
