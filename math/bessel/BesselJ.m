% =========================================================================
% INTRO
%   - Bessel function of the first kind
% -------------------------------------------------------------------------
% INPUT
%   - m: the order
%   - z: argument
% -------------------------------------------------------------------------
% DIMENSION
%   - m before z
% =========================================================================

function J = BesselJ(m, z, varargin)

    CheckDim('preceding', m, z);

    ip = inputParser;
    ip.addParameter('z_large', 1e3);
    ip.addParameter('is_log', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('nu0', 0, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<', 1}));
    ip.parse(varargin{:});
    ip = ip.Results;

    M = max(abs(m));
    m_col = m(:);

    n_full = (0:M).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    J = 0 * n_full .* z_row;

    idx_z_large = abs(z_row) > ip.z_large;
    idx_z_small = abs(z_row) <= ip.z_large;
    z_small = z_row(idx_z_small);
    z_large = z_row(idx_z_large);

    if sum(idx_z_small(:)) > 0
        J(:,idx_z_small) = BesselJLog_Rothwell(M, z_small, 'nu0', ip.nu0);
    end
    if sum(idx_z_large(:)) > 0
        J(:,idx_z_large) = BesselJ_Asym(n_full+ip.nu0, z_large, ...
            'approx_order', M+5, 'is_log', true);
    end
    
    J(1, z_row==0) = log(1);
    J(2:end, z_row==0) = log(0);
    J = J(abs(m_col)+1,:);
    if ~isempty(m_col(m_col<0))
        J(m_col<0,:) = J(m_col<0, :) + log(-1) .* m_col(m_col<0);
    end

    J = J(:, idx_z_row);
    J = reshape(J, size(m .* z));

    if ~ip.is_log
        J = exp(J);
    end
end

