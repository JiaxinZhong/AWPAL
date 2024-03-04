% =========================================================================
% INTRO
%   - Log of the spherical Bessel function of first kind
% -------------------------------------------------------------------------
% INPUT
%   - n: order
% -------------------------------------------------------------------------
% DIMENSION
%   - n before z
% =========================================================================

function [j, j_prime] = SphBesselJ(n, z, varargin)

    validateattributes(n, {'numeric'}, {'>=', 0});

    CheckDim('preceding', n, z);

    ip = inputParser;
%     ip.addParameter('z_large', 1e3);
    ip.addParameter('is_log', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    [j, j_prime] = SphBesselJLog(n, z, 'is_cal_derivative', ip.is_cal_derivative);
    if ~ip.is_log
        j = exp(j);
        j_prime = exp(j_prime);
    end
end

