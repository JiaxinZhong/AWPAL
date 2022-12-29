% =========================================================================
% INTRO
%   - Calculate the integral involved in the radial component
%   - int_r1^r2 jn(kr<) * hn(kr>) * R1 * R2 * ka^3 * rv^2 drv
%       - r< = min(r, rv)
%       - r> = max(r, rv)
% -------------------------------------------------------------------------
% INPUT
%   - sph: the kernel function type 'j' or 'h'
%       - 'j':  hn(kr) * int_r1^r2 ... h_n(k rv) drv
%       - 'h':  jn(kr) * int_r1^r2 ... j_n(k rv) drv
%   - r1: the lower limit of the integral
%   - r2: the upeer limit of the integral
%   - r: the argument of the function outside the integral
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: r 
%   - 2: 1 (holder for theta)
%   - 3: 1 (holder for phi)
%   - 4: la
%   - 5: l1
%   - 6: l2
%   - 7: rv
% =========================================================================
function R = PalCircSrc_SWE_RadialInt(pal, la_max, l1_max, l2_max, r1, r2, r, sph, varargin)

%     validateattributes(r1, {'numeric'}, {'column'});
%     validateattributes(r2, {'numeric'}, {'column'});
%     validateattributes(r, {'numeric'}, {'column'});
    validatestring(sph, {'j', 'h'});

    ip = inputParser;
    ip.addParameter('int_method', 'direct', @(x)any(validatestring(x, {'direct', 'steepest_descent'})));
    % number of points for the numerical integration
    ip.addParameter('int_num', 3e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    % true: The spherical Hankel function is calculated using the limiting
    %   form at large arguments
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % calculate using the log
    ip.addParameter('is_log', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    parse(ip, varargin{:});
    ip = ip.Results;

    switch ip.int_method
        case 'direct'
            R = GaussLegendreQuad(@(rv) ...
                Integrand(pal, sph, rv, r, la_max, l1_max, l2_max, ip.is_farfield), ...
                r1, r2, 'int_num', ip.int_num, 'is_log', false, 'dim', 7);
            % R = exp(R);
        % todo
        case 'steepest_descent'
%             R = GaussLegendreQuad(@(rv) ...
%                 Integrand(pal, sph, r1+1i*rv, r, la_max, l1_max, l2_max, ip.is_farfield), ...
%                 0, inf, 'int_num', ip.int_num, 'is_log', true, 'dim', 7);
        otherwise
            error("Wrong integration methods!");
    end
end

function int = Integrand(pal, sph, rv, r, la_max, l1_max, l2_max, is_farfield)
    m1 = pal.src_low.prf.azimuth_order;
    m2 = pal.src_high.prf.azimuth_order;
    ma = m2 - m1;
    la = permute((0:la_max).', [4,2,3,1]);

    R12_int_num = 2e2;
    
    % dim: r .* rv -> 1 -> l1
    R1 = CircSrc_SWE_Radial(pal.src_low, rv(:), l1_max, 'int_num', R12_int_num);
    % dim: r -> rv -> l1
    R1 = reshape(R1, [size(permute(rv, [1,7,3,4,5,6,2])), l1_max+1]);
    % dim: r -> 1 -> 1 -> 1 -> l1 -> 1 -> rv
    R1 = permute(permute(R1, [1,7,3,4,5,6,2]), [1,2,5,4,3,6,7]);

    % dim: r .* rv -> 1 -> 1 -> l2
    R2 = CircSrc_SWE_Radial(pal.src_high, rv(:), l2_max, 'int_num', R12_int_num);
    % dim: r -> rv -> l2
    R2 = reshape(R2, [size(permute(rv, [1,7,3,4,5,6,2])), l2_max+1]);
    % dim: r -> 1 -> 1 -> 1 -> -> 1 -> l2 -> rv
    R2 = permute(permute(R2, [1,7,3,4,5,6,2]), [1,2,6,4,5,3,7]);
%     min(abs(imag(R1(:))))

    if 0
        max(rv(:));
        xx1 = squeeze(R1);
        figure;
        subplot(211);plot(real(xx1));
        subplot(212);plot(imag(xx1));
    end

    %% the integrand
    switch sph
        case 'j'
            rj = rv;
            rh = r;
        case 'h'
            rj = r;
            rh = rv;
        otherwise
            error('Wrong spherical functions!')
    end
    int = SphBesselJ(2*la+abs(ma), ...
        pal.audio.num*permute(rj, [6,2,3,4,5,1,7]), ...
        'is_log', true) ...
        + SphHankelH(2*la+abs(ma), ...
        pal.audio.num*permute(rh, [6,2,3,4,5,1,7]), ...
        'is_log', true, 'arg_is_large', is_farfield);
%     if max(rv(:)) > 10
%         max(rv(:))
%         xx1 = squeeze(exp(int));
%         figure;
%         subplot(211);plot(real(xx1));
%         subplot(212);plot(imag(xx1));
%     end
    % int = SphBesselJ(2*la+abs(ma), ...
        % pal.audio.num*permute(rj, [6,2,3,4,5,1,7])) ...
        % .* SphHankelH(2*la+abs(ma), ...
        % pal.audio.num*permute(rh, [6,2,3,4,5,1,7]), ...
        % 'arg_is_large', is_farfield);
    int = permute(exp(int), [6,2,3,4,5,1,7]) ...
        .* conj(R1) .* R2 .* (pal.audio.num).^3 .* (rv.^2);
end
