% =========================================================================
% INTRO
%   - Calculate the radial part for a baffled line source in a 2D problem
%   - Evaluati int_0^a bar{u}(rho_s) J_n(k rho_<) H_n(k rho_>) k drho_src
% -------------------------------------------------------------------------
% INPUT
%   - M: maximum order
%   - k: wavenumber
%   - a: half-width (radius)
%   - rho: the polar coordinates
% -------------------------------------------------------------------------
% OUTPUT
%   - R: the radial component 
%   - R_prime: the derivative of the radial component
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: rho
%   - 2: 1 (holder for phi)
%   - 3: order
%   - 4: integration
% =========================================================================
function [R, R_prime] = LineSrc_CWE_Radial(src, rho, m_max, varargin)

    ip = inputParser;
    ip.addParameter('int_num', 2e2);
    ip.addParameter('profile', 'uniform', @(x)any(validatestring(x, {'uniform', 'steerable'})));
    % ip.addParameter('focus_dist', 0.2);
    ip.addParameter('steer_angle', [pi/4;pi/2]);
    % ip.addParameter('is_farfield', 0);
    % ip.addParameter('is_discrete', 0);
    % ip.addParameter('discrete_size', 0.01);
    % ip.addParameter('is_effective', 0);
    % ip.addParameter('effective_size', 0.01);
    ip.addParameter('is_cal_prime', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('array', []);
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % normalization 
    parse(ip, varargin{:});
    ip = ip.Results;

    R_prime = nan;
    
    k = src.wav.num;
    a = src.radius;

    m = permute((-m_max:m_max).', [3,2,1]);

    if ip.is_farfield
        % dim: m => rho
        int = LineSrc_CWE_RadialInt(src, m(:), 'J', 0, a, rho.', ...
            'int_num', ip.int_num, 'is_farfield', ip.is_farfield);
        R = permute(permute(int, [3,2,1]), [2,1,3]);
    else
        % origin points
        idx_origin = rho == 0;
        % interior points
        idx_int = (rho > 0) & (rho < a);
        % exterior points
        idx_ext = rho >= a;
        rho_origin = rho(idx_origin);
        rho_int = rho(idx_int);
        rho_ext = rho(idx_ext);   

        % radial component
        R = 0 * rho .* m;

        %% process origin points
        if ~isempty(rho_origin)
            % dim: m => rho_origin
            int = LineSrc_CWE_RadialInt(src, m(:), 'H', 0, a, rho_origin.', ...
                'int_num', ip.int_num);
            R(idx_origin, :, :) = permute(permute(int, [3,2,1]), [2,1,3]);
        end

        %% process interior points
        if ~isempty(rho_int)
            % dim: m => rho_origin
            int = LineSrc_CWE_RadialInt(src, m(:), 'J', 0, rho_int.', rho_int.', ...
                'int_num', ip.int_num) ...
                + LineSrc_CWE_RadialInt(src, m(:), 'H', rho_int.', a, rho_int.', ...
                'int_num', ip.int_num);
            R(idx_int, :, :) = permute(permute(int, [3,2,1]), [2,1,3]);
        end

        %% process exterior points
        if ~isempty(rho_ext)
            % dim: m => rho_origin
            int = LineSrc_CWE_RadialInt(src, m(:), 'J', 0, a, rho_ext.', ...
                'int_num', ip.int_num);
            R(idx_ext, :, :) = permute(permute(int, [3,2,1]), [2,1,3]);
        end
    end

    R_prime = nan;
    if ~ip.is_cal_prime
        return
    end
    R_prime = R * 0;

    %% the origin, interior, and exterior points
    idx_origin = rho == 0;
    idx_int = (rho > 0) & (rho < a);
    idx_ext = rho >= a;
    rho_origin = rho(idx_origin);
    rho_int = rho(idx_int);
    rho_ext = rho(idx_ext);

    %% process the origin points
    if sum(idx_origin(:)) > 0
        [JDH, JDH_exp] = BesselJPrimeHankelProduct( ...
            m_max, ...
            k*permute(rho_origin, [3,2,1]), k*permute(rho_src,[3,2,1,4]));
        JDH = JDH .* 10.^JDH_exp;
        JDH = permute(JDH, [3,2,1,4]);
        R_prime(repmat(idx_origin, [1,1,size(R_prime,3)])) ...
            = sum(u .* JDH .* weight .* k, 4);
    end

    %% process the interior points
    if sum(idx_int(:)) > 0
        % dim: int, rho2, 1, rho1
        [rho_src1, weight1] = GaussLegendreQuadParam(ip.int_num,  0, permute(rho_int,[4,2,3,1]));
        % dim: rho1, rho2, 1, int
        rho_src1 = permute(rho_src1, [4,2,3,1]);
        weight1 = permute(weight1, [4,2,3,1]);
        % dim: int, rho2, 1, rho1
        [rho_src2, weight2] = GaussLegendreQuadParam(ip.int_num,  permute(rho_int,[4,2,3,1]), a);
        % dim: rho1, rho2, 1, int
        rho_src2 = permute(rho_src2, [4,2,3,1]);
        weight2 = permute(weight2, [4,2,3,1]);

        [JHD, JHD_exp] = BesselJHankelPrimeProduct( ...
            m_max, ...
            k*permute(rho_src1,[3,2,1,4]), k*permute(rho_int, [3,2,1]));
        JHD = JHD .* 10.^JHD_exp;
        JHD = permute(JHD, [3,2,1,4]);
    
        [JDH, JDH_exp] = BesselJPrimeHankelProduct( ...
            m_max, ...
            k*permute(rho_int, [3,2,1]), k*permute(rho_src2,[3,2,1,4]));
        JDH = JDH .* 10.^JDH_exp;
        JDH = permute(JDH, [3,2,1,4]);
    
        R_prime(repmat(idx_int,[1,1,size(R_prime,3)])) ...
            = sum(u .* JHD .* weight1 .* k, 4) ...
            + sum(u .* JDH .* weight2 .* k, 4);
    end

    %% Process the exteior points
    if sum(idx_ext(:)) > 0
        [JHD, JHD_exp] = BesselJHankelPrimeProduct( ...
            m_max, ...
            k*rho_src, k*permute(rho_ext, [3,2,1]));
        JHD = JHD .* 10.^JHD_exp;
        JHD = permute(JHD, [3,2,1,4]);
    
        R_prime(repmat(idx_ext,[1,1,size(R_prime,3)])) = sum(u .* JHD .* weight .* k, 4);
    end
end

function u = SteerableProfile(m, k, steer_angle, rho)
    steer_angle = reshape(steer_angle, 1, 1, 1, 1, 1, 1, length(steer_angle));
%     delta_rho_pos = (a - rho) .* cos(steer_angle);
%     delta_rho_neg = (a + rho) .* cos(steer_angle);
    % u = sum(1i.^(-n) .* exp(-1i*k*delta_rho_pos) + 1i.^n .* exp(-1i*k*delta_rho_neg), 7);
%     u = sum(exp(-1i*k*delta_rho_pos) + (-1).^m .* exp(-1i*k*delta_rho_neg), 7);
    u = sum(exp(1i*real(k)*rho.*cos(steer_angle)) + (-1).^m .* exp(-1i*real(k).*rho.*cos(steer_angle)), 7) / 2;
%     u_phase = [reshape(delta_rho_neg, length(rho), length(steer_angle)); ...
%         reshape(delta_rho_pos, length(rho), length(steer_angle))];
end

%% piston
function u = UniformProfile(m)
    u = (1 + (-1).^m)/2;
end

%% Focusing PAL
function [u, u_phase] = FocusProfile(n, k, a, focus_dist, rho)
    delta_rho = sqrt(rho.^2+focus_dist.^2) - focus_dist.^2;
    % u = 1i.^(-n) .* exp(-1i*k*delta_rho) + 1i.^n .* exp(-1i*k*delta_rho);
    u = exp(-1i*k*delta_rho) + (-1).^n .* exp(-1i*k*delta_rho);
    u_phase = [delta_rho(:); delta_rho(:)];
end

function rho_new = DiscreteProfile(a0, a, rho)
    rho_new = (floor(rho/a0)+1/2)*a0;
end

function u_effective = EffectiveProfile(rho_src, discrete_size, effective_size, a)
    u_effective = 0 * rho_src;
    for i = 0:floor(a/discrete_size)
        rho_buf = rho_src - (i+1/2)*discrete_size;
        u_effective(rho_buf >= -effective_size/2 & rho_buf < effective_size/2) = 1;
    end

end

function int = Integrand(k, m, rho, x_src)
    rho_src = abs(x_src);
    phi = exp(1i*m.*phase(x_src));
    rho_min = min(rho, rho_src);
    rho_max = max(rho, rho_src);
    J = BesselJLog(m, k .* rho_min);
    H = HankelHLog(m, k .* rho_max);
    int = exp(J + H) .* phi;
end
