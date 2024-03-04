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
function [prs, vel] = DIM3D_CircSrc(src, fp, varargin)

    ip = inputParser();
    % normalization scheme
    ip.addParameter('is_cal_vel', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;
    
    [rho_src, weight_rho_src] = GaussLegendreQuadParam(ip.int_num, 0, src.r);
    [phi_src, weight_phi_src] = GaussLegendreQuadParam(ip.int_num, 0, 2*pi);
    rho_src = permute(rho_src, [4,2,3,1]);
    weight_rho_src = permute(weight_rho_src, [4,2,3,1]);
    phi_src = permute(phi_src, [5,2,3,4,1]);
    weight_phi_src = permute(weight_phi_src, [5,2,3,4,1]);

    x_src = rho_src .* cos(phi_src);
    y_src = rho_src .* sin(phi_src);

    dist = sqrt((fp.x-x_src).^2 + (fp.y-y_src).^2 + fp.z.^2);
    prs = (src.wav.num) / (1i*2*pi) .* sum(sum(src.prf.val(rho_src) .* exp(1i*src.wav.num.*dist) ./ dist .* rho_src .* weight_phi_src, 5) .* weight_rho_src, 4);

    %% Calculate velocity
    if ~ip.is_cal_vel
        vel = nan;
        return
    end
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


%         case 'Gauss'
%             vel.x = GaussLegendreQuad(@(ys) ...
%                 IntegrandX(src, fp, ys), ...
%                 -src.ay, src.ay, ...
%                 'int_num', ip.int_num, ...
%                 'is_log', false, 'dim', 10);
%             vel.y = GaussLegendreQuad(@(ys) ...
%                 IntegrandY(src, fp, ys), ...
%                 -src.ay, src.ay, ...
%                 'int_num', ip.int_num, ...
%                 'is_log', false, 'dim', 10);
%         otherwise
%             error("Wrong integration methods!")
%     end

%     vel.x = ampl * vel.x * 1i / 2;
%     vel.y = ampl * vel.y * 1i / 2;

    %% convert the velocity under the source coordinate to the global coordinate
%     vel = vel.Rotate(-src.dir.phi);

    %% Calculate the Lagrangian density
    % if ~ip.is_cal_lag
        % return
    % end
    % lag = rho0/2 * (abs(vel.x).^2 + abs(vel.y).^2) - abs(prs).^2 / (2*rho0*c0^2);
end

% function int = Integrand(src, fp, ys, is_trunc_back)
%     fp_new = fp.Translate(src.pos, 'is_inverse', false);
%     fp_new = fp_new.Rotate(src.dir.phi);
%     int = src.CalProfile(ys) ...
%         .* besselh(0, src.wav.num .* sqrt(fp_new.x.^2 + (fp_new.y - ys).^2));
%     if is_trunc_back
%         int = int .* (fp_new.x >= 0);
%     end
% end
% 
% function int_x = IntegrandX(src, fp, ys)
%     fp_new = fp.Translate(src.pos, 'is_inverse', false);
%     fp_new = fp_new.Rotate(src.dir.phi);
%     dist = sqrt(fp_new.x.^2 + (fp_new.y - ys).^2);
%     int_x = src.CalProfile(ys) .* besselh(1, src.wav.num .* dist) ./ dist .* src.wav.num .* fp_new.x;
% end
% 
% function int_y = IntegrandY(src, fp, ys)
%     fp_new = fp.Translate(src.pos, 'is_inverse', false);
%     fp_new = fp_new.Rotate(src.dir.phi);
%     dist = sqrt(fp_new.x.^2 + (fp_new.y - ys).^2);
%     int_y = src.CalProfile(ys) .* besselh(1, src.wav.num .* dist) ./ dist .* src.wav.num .* (fp_new.y - ys);
% end
