% =========================================================================
% INTRO
%   - The Log of Hankel function of first kind
% INPUT
%   - m: order
% =========================================================================
function [H, H_prime] = HankelH(m, z, varargin)

    CheckDim('preceding', m, z);

    validateattributes(m, {'numeric'}, {'size', [nan, 1]});
    
    ip = inputParser;
    ip.addParameter('z_large', 1e3);
    ip.addParameter('nu0', 0, @(x)validateattributes(x, {'numeric'}, ...
        {'scalar', '>=', 0, '<', 1}));   
    ip.addParameter('kind', 1, ...
        @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 1, '<=', 2}));
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_log', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % Using the limiting form when the argument is very large
    ip.addParameter('arg_is_large', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    M = max(abs(m));

    n_full = (0:M).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    H = 0 * n_full .* z_row;
    H_prime = 0 * n_full .* z_row;

    if ip.arg_is_large
        H = HankelH_Asym(n_full+ip.nu0, z_row, ...
            'approx_order', 0, 'is_log', true, 'kind', ip.kind);
    else
        idx_z_large = abs(z_row) > ip.z_large;
        idx_z_small = abs(z_row) <= ip.z_large;
        z_small = z_row(idx_z_small);
        z_large = z_row(idx_z_large);

        if sum(idx_z_small(:)) > 0
            H(:,idx_z_small) = HankelHLog_Rothwell(M, z_small, 'nu0', ip.nu0, ...
                'kind', ip.kind);
        end
        if sum(idx_z_large(:)) > 0
            H(:,idx_z_large) = HankelH_Asym(n_full+ip.nu0, z_large, ...
                'approx_order', M+5, 'is_log', true, 'kind', ip.kind);
        end
        if ip.is_cal_derivative
            % M+1 to avoid M = 0
            [~, H_prime] = HankelHLog_Rothwell(M+1, z_row, 'nu0', ip.nu0, ...
                'kind', ip.kind, 'is_cal_derivative', ip.is_cal_derivative);
        end
    end
    
    H = H(:, idx_z_row);
    H = H(abs(m)+1, :);
    if ip.is_cal_derivative
        H_prime = H_prime(abs(m)+1, :);
    end
    % negative order
    if ~isempty(m(m<0))
        H(m<0,:) = H(m<0, :) + log(-1) .* m(m<0);
        if ip.is_cal_derivative
            H_prime(m<0,:) = H_prime(m<0, :) + log(-1) .* m(m<0);
        end
    end
    H = reshape(H, size(m .* z));
    if ip.is_cal_derivative
        H_prime = reshape(H_prime, size(m .* z));
    end

    if ~ip.is_log
        H = exp(H);
        H_prime = exp(H_prime);
    end
end

