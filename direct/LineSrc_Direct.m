% =========================================================================
% INTRO
%   - Calculate the ultrasound pressure radiated by a baffled line source
% -------------------------------------------------------------------------
% INPUT
%   - x, the x-coordinate of the field point
%   - y, the y-coordinate of the field point
% =========================================================================
function [prs, vel_x, vel_y, lag] = LineSrc_Direct(src, x, y, varargin)

    ip = inputParser();
    % normalization scheme
    ip.addParameter('src_norm', 'vel', @(x)any(validatestring(x, {'prs', 'vel', 'surf_vel'})));
    ip.addParameter('is_cal_vel', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_cal_lag', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    % the wavenumber
    k = src.wav.num;
    % the half width of the line source
    a = src.radius;

    prs = integral(@(xs) Integrand(k, x, y, xs), -a, a, ...
        'ArrayValued', true, 'AbsTol', 0, 'RelTol', 1e-3);

    rho0 = 1.21;
    c0 = 343;
    % normalization scheme
    switch ip.src_norm
        case 'prs'
            ampl = 1;
        case 'vel'
            ampl = rho0*c0;
        case 'surf_vel'
            ampl = rho0 *c0 / (2*a);
    end
    prs = prs .* ampl / 2;

    vel_x = nan;
    vel_y = nan;
    if ~ip.is_cal_vel
        return
    end

    vel_x = integral(@(xs) IntegrandX(k, x, y, xs), -a, a, ...
        'ArrayValued', true, 'AbsTol', 0, 'RelTol', 1e-3);
    vel_y = integral(@(xs) IntegrandY(k, x, y, xs), -a, a, ...
        'ArrayValued', true, 'AbsTol', 0, 'RelTol', 1e-3);

    vel_x = 1i * ampl * vel_x / rho0 / c0 / 2;
    vel_y = 1i * ampl * vel_y / rho0 / c0 / 2;

    lag = nan;
    if ~ip.is_cal_lag
        return
    end
    lag = rho0/2 * (abs(vel_x).^2 + abs(vel_y).^2) - abs(prs).^2 / (2*rho0*c0^2);
end

function int = Integrand(k, x, y, xs)
    int = k .* besselh(0, k .* sqrt((x-xs).^2 + y.^2));
end

function int_x = IntegrandX(k, x, y, xs)
    dist = sqrt((x-xs).^2 + y.^2);
    int_x = besselh(1, k .* dist) ./ dist .* k .* (x-xs);
end

function int_y = IntegrandY(k, x, y, xs)
    dist = sqrt((x-xs).^2 + y.^2);
    int_y = besselh(1, k .* dist) ./ dist .* k .* y;
end
