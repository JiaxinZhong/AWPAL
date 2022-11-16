% =========================================================================
% INTRO
%   - Calculate the sound pressure radiated by a baffled circular piston
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
function prs = CircPiston_Direct(k, a, x, y, z, varargin)
    
    ip = inputParser();
    ip.addParameter('gauss_num', 5e1, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('profile_type', 'uniform');
    ip.addParameter('steering_angle', 0);
    ip.parse(varargin{:});
    ip = ip.Results;

    validateattributes(k, {'numeric'}, {'scalar'});
    validateattributes(a, {'numeric'}, {'scalar', '>', 0});
    validateattributes(x, {'numeric'}, {'size', [nan,nan,nan,1,1]});
    validateattributes(y, {'numeric'}, {'size', [nan,nan,nan,1,1]});
    validateattributes(z, {'numeric'}, {'size', [nan,nan,nan,1,1]});

    % dimension: 1, 1, 1, rho_src, phi_src
    [rho_src, weight_rho_src] = GaussLegendreQuadParam(ip.gauss_num, 0, a);
    [phi_src, weight_phi_src] = GaussLegendreQuadParam(ip.gauss_num, 0, 2*pi);
    rho_src = permute(rho_src, [4,2,3,1]);
    weight_rho_src = permute(weight_rho_src, [4,2,3,1]);
    phi_src = permute(phi_src, [5,2,3,4,1]);
    weight_phi_src = permute(weight_phi_src, [5,2,3,4,1]);

    x_src = rho_src .* cos(phi_src);
    y_src = rho_src .* sin(phi_src);

    % the velocity profile 
    switch ip.profile_type
        case 'uniform'
            profile = 1;
        case 'steer'
            profile = exp(1i*real(k).*x_src.*sin(ip.steering_angle));
    end

    dist = sqrt((x-x_src).^2 + (y-y_src).^2 + z.^2);
    prs = 1/(1i*2*pi) .* sum(sum(profile .* exp(1i*k.*dist) ./ dist .* k.* rho_src .* weight_phi_src, 5) .* weight_rho_src, 4);
end

