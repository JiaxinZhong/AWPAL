% =========================================================================
% INTRO
%   - Calculate the sound pressure radiated by a baffled circular piston
%   - The amplitude of the on-surface pressure is assumed to be 1
% -------------------------------------------------------------------------
% INPUT
%   - k, the wavenumber
%   - a, the piston radius
%   - rho, the polar radius of the field point
%   - z, the z coordinate of the field point
% =========================================================================
function prs = CircSrc_GBE(k, a, rho, z, varargin)

    validateattributes(rho, {'numeric'}, {'2d'});
    validateattributes(z, {'numeric'}, {'2d'});
    validateattributes(k, {'numeric'}, {'scalar'});
    validateattributes(a, {'numeric'}, {'scalar'});

    ip = inputParser();
    ip.addParameter('coef', 'Huang1999JASA');
    ip.addParameter('A', nan, @(x)validateattributes(x, {'numeric'}, {'vector'}));
    ip.addParameter('B', nan, @(x)validateattributes(x, {'numeric'}, {'vector'}));
    ip.parse(varargin{:});
    ip = ip.Results;
    if isnan(ip.A)
        % import GBE coefficients
        [ip.A, ip.B] = GBE_Coef('set', ip.coef);
    end
    % dimension: 1, 1, AB
    ip.A = reshape(ip.A, 1, 1, length(ip.A));
    ip.B = reshape(ip.B, 1, 1, length(ip.B));

    % Rayleigh distance 
    R = k*a.^2/2;

    % dimensionless coordinates
    rho_norm = rho ./ a;
    z_norm = z ./ R;

    % calculate the sound pressure
    prs = sum(ip.A ./ (1+1i*ip.B.*z_norm) .* ...
        exp(1i*k.*z - ip.B.*rho_norm.^2 ./ (1+1i*ip.B.*z_norm)), 3);
end
