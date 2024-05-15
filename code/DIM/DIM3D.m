% =========================================================================
% INTRO
%   - Calculate the ultrasound pressure radiated by rectangular source
%       using the direct integration method (DIM)
%   - Modified based on LineSrc_DIM.m
% -------------------------------------------------------------------------
% INPUT
%   - fp.x, the x-coordinate of the field point
%   - fp.y, the y-coordinate of the field point
% DIMENSION
%   - [fp.x <=> fp.y <=> fp.z] => rho_src => phi_src
% =========================================================================
function [prs, vel] = DIM3D(src, fp, varargin)

    ip = inputParser();
    ip.addParameter('is_cal_vel', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
    % The coordinates used for the integration
    ip.addParameter('int_coord', 'Cartesian', @(x)any(validatestring(x, {'Cartesian', 'polar'})));
    ip.parse(varargin{:});
    ip = ip.Results;
    
    switch ip.int_coord
        case 'Cartesian'
            % todo
        case 'polar'
            [rho_src, weight_rho_src] = GaussLegendreQuadParam(ip.int_num, 0, src.r, 'dim', 4);
            [phi_src, weight_phi_src] = GaussLegendreQuadParam(ip.int_num, 0, 2*pi, 'dim', 5);
%             rho_src = permute(rho_src, [4,2,3,1]);
%             weight_rho_src = permute(weight_rho_src, [4,2,3,1]);
%             phi_src = permute(phi_src, [5,2,3,4,1]);
%             weight_phi_src = permute(weight_phi_src, [5,2,3,4,1]);
            x_src = rho_src .* cos(phi_src);
            y_src = rho_src .* sin(phi_src);
            dist = sqrt((fp.x-x_src).^2 + (fp.y-y_src).^2 + fp.z.^2);
            prs = (src.wav.num) / (1i*2*pi) .* sum(sum(src.prf.val(rho_src) .* exp(1i*src.wav.num.*dist) ./ dist .* rho_src .* weight_phi_src, 5) .* weight_rho_src, 4);
    end
    
    %% Calculate velocity
    if ~ip.is_cal_vel
        vel = nan;
        return
    end
    switch ip.int_coord
        case 'Cartesian'
            % tod
        case 'polar'
            vel.x = 1/ (2*pi) .* sum(sum(src.prf.val(rho_src) ...
                .* (fp.x - x_src) .* (1 - 1i*src.wav.num.*dist) ./ dist.^3 ...
                .* exp(1i*src.wav.num.*dist) .* rho_src .* weight_phi_src, 5) ...
                .* weight_rho_src, 4);
            vel.y = 1/ (2*pi) .* sum(sum(src.prf.val(rho_src) ...
                .* (fp.y - y_src) .* (1 - 1i*src.wav.num.*dist) ./ dist.^3 ...
                .* exp(1i*src.wav.num.*dist) .* rho_src .* weight_phi_src, 5) ...
                .* weight_rho_src, 4);
            vel.z = 1/ (2*pi) .* sum(sum(src.prf.val(rho_src) ...
                .* (fp.z - 0) .* (1 - 1i*src.wav.num.*dist) ./ dist.^3 ...
                .* exp(1i*src.wav.num.*dist) .* rho_src .* weight_phi_src, 5) ...
                .* weight_rho_src, 4);
    end
end
