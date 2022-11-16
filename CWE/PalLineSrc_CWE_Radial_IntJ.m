% =========================================================================
% INTRO
% -------------------------------------------------------------------------
% INPUT
% DIMENSION
%   - 1, 2: rho 
%   - 3: ma
%   - 4: m1
%   - 5: rho_vsrc
% =========================================================================
function radial = PalLineSrc_CWE_Radial_IntJ(...
    pal, varargin)

    ip = inputParser;
    % number of points for the numerical integration
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('profile', 'uniform', @(x)any(validatestring(x, {'uniform', 'steerable'})));
    ip.addParameter('steer_angle', [pi/4;pi/2]);
    % ip.addParameter('steer_angle_lower', pi/2);
    % ip.addParameter('steer_angle_upper', pi/2);
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('int_method', 'deepest');
    % suggest to use log to avoid overflow
    ip.addParameter('is_log', true)
    ip.addParameter('m1', nan);
    ip.addParameter('ma', nan);
    ip.addParameter('ma_max', nan);
    ip.addParameter('m1_max', nan);
    ip.addParameter('lower_limit', nan);
    ip.addParameter('upper_limit', nan);
    parse(ip, varargin{:});
    ip = ip.Results;

    k1 = pal.ultra_low.num;
    k2 = pal.ultra_high.num;
    ka = pal.audio.num;
    a = pal.src_low.radius;
    
    ma = ip.ma;
    m1 = ip.m1;
    if ~isnan(ip.ma_max) && ~isnan(ip.m1_max)
        ma = permute((-ip.ma_max:ip.ma_max).', [3,2,1]);
        m1 = permute((-ip.m1_max:ip.m1_max).', [4,2,3,1]);
    end
    m2 = m1+ma;

    if (ip.lower_limit == 0) && (ip.upper_limit == a)
        % to do
    end
    if isinf(ip.upper_limit) && (ip.lower_limit == a)
        radial_H = PalLineSrc_CWE_RadialH(...
            k1, k2, ka, a, m1, ma, ...
            'is_log', true, ...
            'int_method', 'hybrid', ...
            'ma_max', ip.ma_max, ...
            'm1_max', ip.m1_max);
        radial_H2 = PalLineSrc_CWE_RadialH2(...
            k1, k2, ka, ...
            a, m1, ma, ...
            'is_log', true, ...
            'int_method', 'gauss', ...
            'ma_max', ip.ma_max, ...
            'm1_max', ip.m1_max);
        radial = LogAdd(radial_H, radial_H2) + log(1/2);
%         J1 = LineSrc_CWE_IntJ(k1, a, ...
%             0, a, ...
%             'max_order', ip.m1_max, ...
%             'is_log', true,...
%             'int_method', 'gauss', 'int_num', 150);
%         J1 = permute(J1, [4,2,3,1]);
%         radial = radial + conj(J1);
%         J2 = LineSrc_CWE_IntJ(k2, a, 0, a, ...
%             'max_order', ip.m1_max+ip.ma_max, ...
%             'is_log', true,...
%             'int_method', 'gauss', 'int_num', 150);
%         for i = -ip.ma_max:ip.ma_max
%             for j = -ip.m1_max:ip.m1_max
%                 radial(:,:,i+ip.ma_max+1,j+ip.m1_max+1) ...
%                     = radial(:,:,i+ip.ma_max+1,j+ip.m1_max+1) ...
%                     + J2(i+j + ip.ma_max + ip.m1_max + 1);
%             end
%         end
    end
end

function res = FunLog1(k1, k2, ka, m1, m2, ma, rho_vsrc)
    alpha = imag(k1+k2);
    H_m1 = HankelHScaled(m1(:), k1*conj(rho_vsrc), 'is_log', true);
    H_m2 = HankelHScaled(m2(:), k2*rho_vsrc, 'is_log', true);
    H_ma = HankelHScaled(ma(:), ka*rho_vsrc, 'is_log', true);
    IsInvalid(H_m1, 'var_name', 'H_m1');
    IsInvalid(H_m2, 'var_name', 'H_m2');
    IsInvalid(H_ma, 'var_name', 'H_ma');
    H_m1 = reshape(H_m1, size(m1 .* rho_vsrc));
    H_m2 = reshape(H_m2, size(m2 .* rho_vsrc));
    H_ma = reshape(H_ma, size(ma .* rho_vsrc));
    res = (conj(H_m1) + H_m2 + H_ma + 2*1i*ka*rho_vsrc-alpha.*rho_vsrc)...
        + log(ka^2 .* rho_vsrc);
end

function res = FunLog2(k1, k2, ka, m1, m2, ma, a, rho_vsrc)
    alpha = imag(k1+k2);
    K = (2*1i*ka + alpha) ./ (4*ka.^2 + alpha.^2);
    rho_vsrc = K.*rho_vsrc+a;
    res = FunLog1(k1, k2, ka, m1, m2, ma, rho_vsrc) + log(K);
end

function res = Fun1(k1, k2, ka, m1, m2, ma, rho_vsrc)
%     res = conj(HankelH_Matlab(m1, k1*conj(rho_vsrc), 'kind', 1, 'scaled', true)) ...
%         .* HankelH_Matlab(m2, k2*rho_vsrc, 'kind', 1, 'scaled', true) ...
%         .* HankelH_Matlab(ma, ka*rho_vsrc, 'kind', 1, 'scaled', true);
    alpha = imag(k1+k2);
    H_m1 = HankelHScaled(m1(:), k1*conj(rho_vsrc), 'is_log', true);
    H_m2 = HankelHScaled(m2(:), k2*rho_vsrc, 'is_log', true);
    H_ma = HankelHScaled(ma(:), ka*rho_vsrc, 'is_log', true);
    H_m1 = reshape(H_m1, size(m1 .* rho_vsrc));
    H_m2 = reshape(H_m2, size(m2 .* rho_vsrc));
    H_ma = reshape(H_ma, size(ma .* rho_vsrc));
    res = exp(conj(H_m1) + H_m2 + H_ma + 2*1i*ka*rho_vsrc-alpha.*rho_vsrc)...
        .* ka^2 .* rho_vsrc;

%     res = conj(HankelH_Matlab(m1, k1*conj(rho_vsrc), 'kind', 1, 'scaled', true)) ...
%         .* HankelH_Matlab(m2, k2*rho_vsrc, 'kind', 1, 'scaled', true) ...
%         .* HankelH_Matlab(ma, ka*rho_vsrc, 'kind', 1, 'scaled', true) ...
%         .* exp(2*1i*ka*rho_vsrc-alpha.*rho_vsrc) .* ka^2 .* rho_vsrc;
%     res = conj(HankelHScale(m1, k1*conj(rho_vsrc))) ...
%         .* HankelH_Matlab(m2, k2*rho_vsrc, 'kind', 1, 'scaled', true) ...
%         .* HankelH_Matlab(ma, ka*rho_vsrc, 'kind', 1, 'scaled', true);
end

function res = Fun2(k1, k2, ka, m1, m2, ma, rho_vsrc)
%     res = conj(HankelH_Matlab(m1, k1*conj(rho_vsrc), 'kind', 1, 'scaled', true)) ...
%         .* HankelH_Matlab(m2, k2*rho_vsrc, 'kind', 1, 'scaled', true) ...
%         .* HankelH_Matlab(ma, ka*rho_vsrc, 'kind', 1, 'scaled', true);
    alpha = imag(k1+k2);
%     H_m1 = HankelHScaled(m1(:), k1*conj(rho_vsrc), 'is_log', true);
%     H_m2 = HankelHScaled(m2(:), k2*rho_vsrc, 'is_log', true);
%     H_ma = HankelHScaled(ma(:), ka*rho_vsrc, 'is_log', true);
%     H_m1 = reshape(H_m1, size(m1 .* rho_vsrc));
%     H_m2 = reshape(H_m2, size(m2 .* rho_vsrc));
%     H_ma = reshape(H_ma, size(ma .* rho_vsrc));
%     res = exp(conj(H_m1) + H_m2 + H_ma)...
%         .* exp(2*1i*ka*rho_vsrc-alpha.*rho_vsrc) .* ka^2 .* rho_vsrc;

    m1_col = m1 + 0*ma;
    m1_col = m1_col(:);
    m2_col = m2(:);
    ma_col = ma + 0*m1;
    ma_col = ma_col(:);

    res = conj(HankelH_Matlab(m1_col, k1*conj(rho_vsrc), 'kind', 1, 'scaled', true)) ...
        .* HankelH_Matlab(m2_col, k2*rho_vsrc, 'kind', 1, 'scaled', true) ...
        .* HankelH_Matlab(ma_col, ka*rho_vsrc, 'kind', 1, 'scaled', true) ...
        .* exp(2*1i*ka*rho_vsrc-alpha.*rho_vsrc) .* ka^2 .* rho_vsrc;
    res = reshape(res, size(m2));
end
