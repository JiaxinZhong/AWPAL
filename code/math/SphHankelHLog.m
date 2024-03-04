% =========================================================================
% INTRO
%   - The Log of Bessel function of first kind
% INPUT
%   - n: order
% DIMENSION
%   - n before z
% =========================================================================

function [h, h_prime] = SphHankelHLog(n, z, varargin)

    validateattributes(n, {'numeric'}, {'>=', 0});

    ip = inputParser;
%     ip.addParameter('z_large', 1e3);
    % Using the limiting form when the argument is very large
    ip.addParameter('arg_is_large', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    N = max(abs(n));

    n_full = (0:N).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    [h, h_prime] = HankelH(n_full, z_row, 'nu0', .5, 'is_log', true, ...
        'arg_is_larg', ip.arg_is_large, 'is_cal_derivative', ip.is_cal_derivative);

    if ip.is_cal_derivative
        h_prime = log(-.5*sqrt(pi/2./z_row.^3).*exp(h) + sqrt(pi/2./z_row).*exp(h_prime));
        h_prime = h_prime(n+1,:);
        h_prime = h_prime(:, idx_z_row);
        h_prime = reshape(h_prime, size(n .* z));
    end

    h = h + 1/2 .* log(pi./2./z_row);

%     h(1,z_row==0) = log(1);
%     h(2:end, z_row == 0) = log(0);

    h = h(n+1,:);
    h = h(:, idx_z_row);
    h = reshape(h, size(n .* z));
end

