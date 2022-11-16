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
    
    prs = k / (2*pi*1i) .* integral(@(phis) ...
        integral(@(rhos) Integrand(k, x, y, z, rhos, phis), ...
        0, a, 'ArrayValued', true, 'AbsTol', 0, 'RelTol', 1e-3), ...
        0, 2*pi, 'ArrayValued', true, 'AbsTol', 0, 'RelTol', 1e-3);

end

function int = Integrand(k, x, y, z, rhos, phis)
    xs = rhos .* cos(phis);
    ys = rhos .* sin(phis);
    dist = sqrt((x-xs).^2 + (y-ys).^2 + z.^2);
    int = exp(1i*k.*dist) ./ dist .* rhos;
end
