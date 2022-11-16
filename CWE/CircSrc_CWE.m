% ===========================================================================
% INTRO
%   - Calculate the radiation from a baffled circular source
%   - Method: CWE (cylindrical wave expansion)
%   - Default on-surface pressure: p0 = rho0c0v0 = 1
% -------------------------------------------------------------------------
% INPUT
%   - k, the wavenumber
%   - a, the radius of the source
%   - rho, the polar coordinates of the field points
%   - z, the z coordinates of the field points
% DIMENSION
%   - 1 and 2: rho .* z
%   - 3: int
% OUTPUT
%   - prs, sound pressure
% ===========================================================================
function prs = CircSrc_CWE(k, a, rho, z, varargin)

    validateattributes(k, {'numeric'}, {'scalar'});
    validateattributes(a, {'numeric'}, {'scalar'});
    validateattributes(rho, {'numeric'}, {'size', [nan, nan, 1]});
    validateattributes(z, {'numeric'}, {'size', [nan, nan, 1]});

    ip = inputParser();
    % the number of points for the numerical integration
    ip.addParameter('int_num', 1e2);
    ip.parse(varargin{:});
    ip = ip.Results;

    [node1, weight1] = GaussLegendreQuadParam(ip.int_num, 0, pi/2);
    [node2, weight2] = GaussLegendreQuadParam(ip.int_num, 0, inf);
    node1 = permute(node1, [3,2,1]);
    weight1 = permute(weight1, [3,2,1]);
    node2 = permute(node2, [3,2,1]);
    weight2 = permute(weight2, [3,2,1]);

    krho1 = real(k) .* sin(node1);
    kz1 = sqrt(k^2 - krho1.^2);
    krho2 = real(k) .* cosh(node2);
    kz2 = sqrt(k^2 - krho2.^2);

    int1 = besselj(0, krho1 .* rho) .* besselj(1, krho1 .* a) ...
        .* exp(1i*kz1.*z) ./ kz1 .* real(k) .* cos(node1);
    int2 = besselj(0, krho2 .* rho) .* besselj(1, krho2 .* a) ...
        .* exp(1i*kz2.*z) ./ kz2 .* real(k) .* sinh(node2);
    int1(isnan(int1)) = 0;
    int2(isnan(int2)) = 0;
%     int1 = 0;
    prs = k*a .* (sum(int1 .* weight1, 3) + sum(int2 .* weight2, 3));

end
