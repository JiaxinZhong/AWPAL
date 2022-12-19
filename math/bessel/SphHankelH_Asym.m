% =========================================================================
% INTRO
%   - The asymptotic expansions of the spherical Hankel functions
% =========================================================================

function h = SphHankelH_Asym(m, z, varargin)

    CheckDim('preceding', m, z);

    ip = inputParser;
    ip.addParameter('is_log', false);
    ip.addParameter('approx_order', 0);
    ip.parse(varargin{:})
    ip = ip.Results;

    M = max(abs(m(:)));
    m_full_col = (0:M).';
    m_col = m(:);
    z_row = z(:).';

    h = HankelH_Asym(m_full_col + 0.5, z, 'approx_order', ip.approx_order, ...
        'is_log', true) + 1/2 .* log(pi/2./z_row);
    h = h(m_col+1, :);
    h = reshape(h, size(m .* z));
    if ~ip.is_log
        h = exp(h);
    end
end
