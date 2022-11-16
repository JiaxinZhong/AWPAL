% ==========================================================================
% - Calculate PAL sound pressure using the Berktay's convolution model
% ==========================================================================
function [dir, dira] = PalLineSrc_Conv(...
        pal, phi, varargin)

    ip = inputParser();
    ip.addParameter('int_num', 2e3);
    ip.addParameter('type', 'conventional')
    % Normalize in decibels
    ip.addParameter('is_norm_dB', false)
    ip.parse(varargin{:});
    ip = ip.Results;

    [phi_vsrc, weight] = GaussLegendreQuadParam(ip.int_num, 0, 2*pi, 'dim', 3);

    dir1 = pal.src_low.CalDirectivity(phi_vsrc);
    dir2 = pal.src_high.CalDirectivity(phi_vsrc);
    di_westervelt = pal.absorp_coef ./ (pal.ultra_high.num - conj(pal.ultra_low.num) ...
        - pal.audio.num .* cos(phi-phi_vsrc));
    
    switch ip.type
        case 'direct'
            dira = nan;
            dir = sum(di_westervelt .*conj(dir1) .* (dir2).* weight, 3);
        case 'improved'
            dira = pal.src_audio.CalDirectivity(phi);
%             dira2 = pal.src_audio.CalDirectivity(phi - phi_vsrc);
            dir = sum(di_westervelt .* conj(dir1) .* (dir2) .* weight, 3).*dira;
        case 'Westervelt'
            dir = pal.absorp_coef ...
                ./ (pal.ultra_high.num - conj(pal.ultra_low.num) ...
                - pal.audio.num .* cos(phi-pi/2));
        otherwise
            error('Wrong convolution model!')
    end

    if ip.is_norm_dB
        dir = 20*log10(abs(dir));
        dir = dir - max(dir(:));
    end
end
