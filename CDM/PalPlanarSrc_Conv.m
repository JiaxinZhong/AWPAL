% ==========================================================================
% - Calculate PAL sound pressure using the Berktay's convolution model
% INPUT
%   - pal: PAL wave info
%   - theta, phi: zenithal and azimuthal coordinates of field points
% ==========================================================================
function dir = PalPlanarSrc_Conv(pal, fp, varargin)

    ip = inputParser();
    ip.addParameter('int_num', 500);
    ip.addParameter('is_norm_dB', true);
    ip.addParameter('type', 'conventional')
    ip.parse(varargin{:});
    ip = ip.Results;

    k1 = pal.ultra_low.num;
    k2 = pal.ultra_high.num;
    ka = pal.audio.num;
    alpha = (imag(k1) + imag(k2))/2;
    
    [theta_vsrc, theta_vsrc_weight] = GaussLegendreQuadParam(ip.int_num, 0, pi, 'dim', 8);
    [phi_vsrc, phi_vsrc_weight] = GaussLegendreQuadParam(ip.int_num, 0, 2*pi, 'dim', 9);

    dir1 = pal.src_low.CalDirectivity(theta_vsrc, phi_vsrc);
    dir2 = pal.src_high.CalDirectivity(theta_vsrc, phi_vsrc);

%     di_westervelt = sqrt(1-1i*ka/alpha.*(sin((phi-phi_vsrc)/2)).^2)...
%         .* sqrt(1-1i*ka/alpha.*(sin((phi+phi_vsrc)/2)).^2);
    gam = cos(fp.theta) .* cos(theta_vsrc) + sin(fp.theta) .* sin(theta_vsrc) ...
        .* cos(fp.phi - phi_vsrc);
    di_westervelt = alpha ./ (k2 - conj(k1) - ka .* gam);
    
    switch ip.type
        case 'direct'
            dir = sum(di_westervelt .*conj(dir1) .* (dir2).* phi_vsrc_weight ...
                .* theta_vsrc_weight .* sin(theta_vsrc), [8,9]);
        case 'improved'
            dira = pal.src_audio.CalDirectivity(fp.theta, fp.phi);
            dir = dira .* sum(di_westervelt .*conj(dir1) .* dir2.* phi_vsrc_weight ...
                .* theta_vsrc_weight.* sin(theta_vsrc), [8,9]);
        case 'Westervelt'
            dir = pal.absorp_coef ...
                ./ (pal.ultra_high.num - conj(pal.ultra_low.num) ...
                - pal.audio.num .* cos(fp.theta)) + 0*fp.phi; 
        otherwise
            error('Wrong type!');
    end

    if ip.is_norm_dB
        dir = 20*log10(abs(dir));
        dir = dir - max(dir(:));
    end
end
