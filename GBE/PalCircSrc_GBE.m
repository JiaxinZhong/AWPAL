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
function prs = PalCircSrc_GBE(pal, a, rho, z, varargin)

    validateattributes(rho, {'numeric'}, {'2d'});
    validateattributes(z, {'numeric'}, {'2d'});
    validateattributes(a, {'numeric'}, {'scalar'});

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
    R1 = pal.ultra_low.num*a.^2/2;
    R2 = pal.ultra_high.num*a.^2/2;
    Ra = pal.audio.num * a.^2 / 2;

    % calculate the sound pressure
    prs = 0;
    interval = [0; 1; 2; 3; 4];
    for i = 1:length(interval)-1
        prs = prs + GaussLegendreQuad(...
            @(zv) Integrand(zv, z, rho/a, ...
            pal.ultra_low.num, pal.ultra_high.num, pal.audio.num, ...
            A1, A2, B1, B2, R1, R2, Ra), ...
            interval(i), interval(i+1),...
            'int_num', ip.int_num, 'dim', 5);
    end

    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    prs = prs .* beta .* pal.audio.angfreq ./ (2*1i*rho0*c0^3);
end


function int = Integrand(zv, z, rho, k1, k2, ka, A1, A2, B1, B2, R1, R2, Ra)
    zv1 = zv ./ R1;
    zv2 = zv ./ R2;
    za = abs(z-zv)./Ra;
    A = conj(A1) .* A2 ./ conj(1 + 1i * B1 .* zv1) ./ (1+1i*B2 .* zv2);
    B = conj(B1 ./ (1 + 1i*B1 .* zv1)) + B2 ./ (1+1i*B2.*zv2);

    gauss = sum(A./(1+1i*B.*za) .* exp(-B./(1+1i*B.*za).*rho), [3,4]);
    int = gauss.*exp(1i.*(k2-conj(k1)).*zv+1i*ka.*abs(z-zv));
end
