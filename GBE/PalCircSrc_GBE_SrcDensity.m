% =========================================================================
% INTRO
%   - Calculate the sound pressure radiated by a baffled circular piston
% -------------------------------------------------------------------------
% INPUT
%   - k, the wavenumber
%   - a, the piston radius
%   - rho, the polar radius of the field point
%   - z, the z coordinate of the field point
% =========================================================================
function q = PalCircSrc_GBE_SrcDensity(pal, a, rhov, zv, varargin)

    ip = inputParser();
    ip.addParameter('coef', 'Liu2008JASA_Table2');
    ip.addParameter('int_num', 200);
    ip.addParameter('A', nan, @(x)validateattributes(x, {'numeric'}, {'vector'}));
    ip.addParameter('B', nan, @(x)validateattributes(x, {'numeric'}, {'vector'}));
    ip.parse(varargin{:});
    ip = ip.Results;
    if isnan(ip.A)
        % import GBE coefficients
        [ip.A, ip.B] = GBE_Coef('set', ip.coef);
    end
    % dimension: 1, 1, AB
    A1 = reshape(ip.A, 1, 1, length(ip.A));
    B1 = reshape(ip.B, 1, 1, length(ip.B));
    A2 = reshape(ip.A, 1, 1, 1, length(ip.A));
    B2 = reshape(ip.B, 1, 1, 1, length(ip.B));

    % Rayleigh distance 
    k1 = pal.ultra_low.num;
    k2 = pal.ultra_high.num;
    R1 = k1*a.^2/2;
    R2 = k2*a.^2/2;

    % calculate the sound pressure
    zv1 = zv ./ R1;
    zv2 = zv ./ R2;
    A = conj(A1) .* A2 ./ conj(1+1i*B1.*zv1) ./ (1+1i*B2.*zv2);
    B = conj(B1 ./ (1 + 1i*B1 .* zv1)) + B2 ./ (1+1i*B2.*zv2);
    rhov_tilde = rhov ./ a;

    gauss = sum(A .* exp(-B .* rhov_tilde), [3,4]);

    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    % prs = prs .* beta .* pal.audio.angfreq ./ (2*1i*rho0*c0^3);

    q = gauss .* exp(1i.*(k2-conj(k1)).*zv);
end


% function int = Integrand(zv, z, rho, k1, k2, ka, A1, A2, B1, B2, R1, R2, Ra)
% end
