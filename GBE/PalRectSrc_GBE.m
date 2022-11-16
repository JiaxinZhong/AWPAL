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
function prs = PalRectSrc_GBE(pal, ax, ay, x, y, z, varargin)

    % validateattributes(a, {'numeric'}, {'scalar'});

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
    A1 = reshape(ip.A, 1, 1, 1, length(ip.A));
    B1 = reshape(ip.B, 1, 1, 1, length(ip.B));
    A2 = reshape(ip.A, 1, 1, 1, 1, length(ip.A));
    B2 = reshape(ip.B, 1, 1, 1, 1, length(ip.B));

    % Rayleigh distance 
    R1x = pal.ultra_low.num * ax.^2 / 2;
    R1y = pal.ultra_low.num * ay.^2 / 2;
    R2x = pal.ultra_high.num *ax.^2 / 2;
    R2y = pal.ultra_high.num * ay.^2 / 2;
    Rax = pal.audio.num * ax.^2/2;
    Ray = pal.audio.num * ay.^2/2;

    % calculate the sound pressure
    prs = 0;
    interval = [0; 1; 2; 4];
    for i = 1:length(interval)-1
        prs = prs + GaussLegendreQuad(...
            @(zv) Integrand(zv, x/ax, y/ay, z, ...
            pal.ultra_low.num, pal.ultra_high.num, pal.audio.num, ...
            A1, A2, B1, B2, R1x, R1y, R2x, R2y, Rax, Ray), ...
            interval(i), interval(i+1),...
            'int_num', ip.int_num, 'dim', 6);
    end

    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    prs = prs .* beta .* pal.audio.angfreq ./ (2*1i*rho0*c0^3);
end


function int = Integrand(zv, x, y, z, k1, k2, ka, A1, A2, B1, B2, R1x, R1y, R2x, R2y, Rax, Ray)
    zv1x = zv ./ R1x;
    zv1y = zv ./ R1y;
    zv2x = zv ./ R2x;
    zv2y = zv ./ R2y;
    zax = abs(z-zv)./Rax;
    zay = abs(z-zv)./Ray;
    Ax = conj(A1) .* A2 ./ conj(sqrt(1 + 1i * B1 .* zv1x)) ./ sqrt(1+1i*B2.*zv2x);
    Bx = conj(B1 ./ (1 + 1i*B1 .* zv1x)) + B2 ./ (1+1i*B2.*zv2x);
    Ay = conj(A1) .* A2 ./ conj(sqrt(1 + 1i * B1 .* zv1y)) ./ sqrt(1+1i*B2.*zv2y);
    By = conj(B1 ./ (1 + 1i*B1 .* zv1y)) + B2 ./ (1+1i*B2.*zv2y);

    gauss_x = sum(Ax ./ sqrt(1+1i*Bx.*zax) .* exp(-Bx.*x.^2./(1+1i*Bx.*zax)), [4,5]);
    gauss_y = sum(Ay ./ sqrt(1+1i*By.*zay) .* exp(-By.*y.^2./(1+1i*By.*zay)), [4,5]);

    int = gauss_x .* gauss_y .* exp(1i.*(k2-conj(k1)).*zv+1i*ka.*abs(z-zv));
end
