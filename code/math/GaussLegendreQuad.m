% =========================================================================
% INTRO
%   - Calculate the integral using the Gauss-Legendre quadrature
% -------------------------------------------------------------------------
% INPUT
%   - int_num: the number of points for the integration
%   - a, b: the lower and upper limit 
%   Dimension:  Gauss => a .* b
% -------------------------------------------------------------------------
% Output
%	node, weight	--- node and weights
% =========================================================================
function int = GaussLegendreQuad(integrand, a, b, varargin)


    ip = inputParser;
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    % true: integrand and int are log
    ip.addParameter('is_log', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('dim', nan, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

    % the dimension where stores the nodes and weights
    if isnan(ip.dim)
        ip.dim = max(ndims(a), ndims(b)) + 1;
    end

    [node, weight] = GaussLegendreQuadParam(ip.int_num, a, b, 'dim', ip.dim);

    if ip.is_log
        buf = integrand(node) + log(weight);
%         IsInvalid(buf);
        max_exp = max(real(buf), [], ip.dim);
        int = max_exp + log(sum(exp(buf-max_exp), ip.dim));
        int(isnan(int)) = -inf;
    else
        int = sum(integrand(node) .* weight, ip.dim);
    end
end
