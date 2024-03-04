% legacy
% =========================================================================
% INTRO
%   - The Log of Bessel function of first kind
% INPUT
%   - n: order
% =========================================================================

function [J, J_prime] = BesselJLog(n, z, varargin)

    CheckDim('preceding', n, z);

    ip = inputParser;
    ip.addParameter('z_large', 1e3);
    ip.addParameter('nu0', 0, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<', 1}));
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    J_prime = [];

    N = max(abs(n));

    n_full = (0:N).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    J = 0 * n_full .* z_row;

    idx_z_large = abs(z_row) > ip.z_large;
    idx_z_small = abs(z_row) <= ip.z_large;
    z_small = z_row(idx_z_small);
    z_large = z_row(idx_z_large);

    if sum(idx_z_small(:)) > 0
        J(:,idx_z_small) = BesselJLog_Rothwell(N, z_small, 'nu0', ip.nu0, ...
            'is_cal_derivative', ip.is_cal_derivative);
    end
    if sum(idx_z_large(:)) > 0
        J(:,idx_z_large) = BesselJ_Asym(n_full+ip.nu0, z_large, ...
            'approx_order', N+5, 'is_log', true);
    end
    if ip.is_cal_derivative
        % M+1 to avoid M = 0
        [~, J_prime] = BesselJLog_Rothwell(N+1, z_row, 'nu0', ip.nu0, ...
            'is_cal_derivative', ip.is_cal_derivative);
    end

    if ip.is_cal_derivative
        J_prime = J_prime(abs(n)+1, :);
        J_prime = reshape(J_prime, size(n.*z));
    end
    
    J(1, z_row==0) = log(1);
    J(2:end, z_row==0) = log(0);
    J = J(abs(n)+1,:);
    if ~isempty(n(n<0))
        J(n<0,:) = J(n<0, :) + log(-1) .* n(n<0);
    end

    J = J(:, idx_z_row);
    J = reshape(J, size(n .* z));

end

