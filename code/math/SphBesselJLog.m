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

function [j, j_prime] = SphBesselJLog(n, z, varargin)

    validateattributes(n, {'numeric'}, {'>=', 0});

    ip = inputParser;
%     ip.addParameter('z_large', 1e3);
    % ip.addParameter('is_log', false, 
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    N = max(abs(n));

    n_full = (0:N).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    [j, j_prime] = BesselJLog(n_full, z_row, 'nu0', .5, ...
        'is_cal_derivative', ip.is_cal_derivative);

    if ip.is_cal_derivative
        j_prime = log(-.5*sqrt(pi/2./z_row.^3).*exp(j) + sqrt(pi/2./z_row).*exp(j_prime));
        j_prime = j_prime(n+1,:);
        j_prime = j_prime(:, idx_z_row);
        j_prime = reshape(j_prime, size(n .* z));    
    end

    j = j + 1/2 .* log(pi./2./z_row);

    j(1,z_row==0) = log(1);
    j(2:end, z_row == 0) = log(0);

    j = j(n+1,:);

    j = j(:, idx_z_row);
    j = reshape(j, size(n .* z));
end

