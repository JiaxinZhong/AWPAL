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
    
    rho0 = 1.21;
    c0 = 343;
    beta = 1.2;

    prs = integral3(@(r_vsrc, theta_vsrc, phi_vsrc) ...
        Integrand(ka, k1, k2, a, x, y, z, r_vsrc, theta_vsrc, phi_vsrc), ...
        0, Inf, 0, pi, 0, 2*pi, ...
        'AbsTol', 0, 'RelTol', 1e-2);
    prs = -beta / (4 * pi * rho0 * c0^2) .* prs;

end

function int = Integrand(ka, k1, k2, a, x, y, z, r_vsrc, theta_vsrc, phi_vsrc)

    x_vsrc = r_vsrc .* sin(theta_vsrc) .* cos(phi_vsrc) + x;
    y_vsrc = r_vsrc .* sin(theta_vsrc) .* sin(phi_vsrc) + y;
    z_vsrc = r_vsrc .* cos(theta_vsrc) + z;

    int = conj(CircPiston_Direct(k1, a, x_vsrc, y_vsrc, z_vsrc)) ...
        .* CircPiston_Direct(k2, a, x_vsrc, y_vsrc, z_vsrc) ...
        .* exp(1i* ka .* r_vsrc) .* ka.^2 .* r_vsrc .* sin(theta_vsrc);
end
