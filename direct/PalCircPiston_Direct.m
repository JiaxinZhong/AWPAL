% =========================================================================
% INTRO
%   - Calculate the sound pressure radiated by a baffled circular PAL
%   - Direct method: direct numerical integrations
%   - The amplitude of the on-surface pressure is assumed to be 1
% -------------------------------------------------------------------------
% INPUT
%   - k, the wavenumber
%   - a, the piston radius
%   - x, the x-coorodinate of the field point
%   - y, the y-coorodinate of the field point
%   - z, the z-coordinate of the field point
% =========================================================================
function prs = PalCircPiston_Direct(ka, k1, k2, a, x, y, z, varargin)
    
    ip = inputParser();
    ip.addParameter('gauss_num', 5e1, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.parse(varargin{:});

    validateattributes(ka, {'numeric'}, {'scalar'});
    validateattributes(k1, {'numeric'}, {'scalar'});
    validateattributes(k2, {'numeric'}, {'scalar'});
    validateattributes(a, {'numeric'}, {'scalar', '>', 0});
    validateattributes(x, {'numeric'}, {'size', [nan,nan,nan,1,1,1]});
    validateattributes(y, {'numeric'}, {'size', [nan,nan,nan,1,1,1]});
    validateattributes(z, {'numeric'}, {'size', [nan,nan,nan,1,1,1]});

    ip = ip.Results;
    rho0 = 1.21;
    c0 = 343;
    beta = 1.2;

    % dimension: 1, 1, 1, r_v, theta_v, phi_v
    [r_v, weight_r_v] = GaussLegendreQuadParam(ip.gauss_num, 0, 5);
    [theta_v, weight_theta_v] = GaussLegendreQuadParam(ip.gauss_num, 0, pi);
    [phi_v, weight_phi_v] = GaussLegendreQuadParam(ip.gauss_num, 0, 2*pi);
    r_v = permute(r_v, [4,2,3,1]);
    weight_r_v = permute(weight_r_v, [4,2,3,1]);
    theta_v = permute(theta_v, [5,2,3,4,1]);
    weight_theta_v = permute(weight_theta_v, [5,2,3,4,1]);
    phi_v = permute(phi_v, [6,2,3,4,5,1]);
    weight_phi_v = permute(weight_phi_v, [6,2,3,4,5,1]);

    x_v = r_v .* sin(theta_v) .* cos(phi_v) + x;
    y_v = r_v .* sin(theta_v) .* sin(phi_v) + y;
    z_v = r_v .* cos(theta_v) + z;

    prs1 = CircPiston_Direct(k1, a, squeeze(x_v), squeeze(y_v), squeeze(z_v));
    prs1 = reshape(prs1, size(x_v));
    prs2 = CircPiston_Direct(k2, a, squeeze(x_v), squeeze(y_v), squeeze(z_v));
    prs2 = reshape(prs2, size(x_v));
    int = conj(prs1) .* prs2 ...
        .* exp(1i* ka .* r_v) .* ka.^2 .* r_v .* sin(theta_v);
    prs = -beta / (4 * pi * rho0 * c0^2) .* sum(sum(sum(int .* weight_phi_v, 6) .* weight_theta_v, 5) .* weight_r_v, 4);

end


