% =========================================================================
% INTRO
%   - Calculate the radial component for a baffled circular source using
%       the spherical wave expansion method
% -------------------------------------------------------------------------
% INPUT
%   - M: maximum order
%   - k: wavenumber
%   - a: half-width (radius)
%   - rho: the polar coordinates
% -------------------------------------------------------------------------
% OUTPUT
%   - R: the radial component 
%   - R_prime: the derivative of the radial component
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: r
%   - 2: 1 (holder for theta)
%   - 3: 1 (holder for phi)
%   - 4: order l
%   - 5: integration
% =========================================================================
function [R, R_prime] = CircSrc_SWE_Radial(src, r, ell_max, varargin)

    ip = inputParser;
    ip.addParameter('int_num', []);
    ip.addParameter('is_cal_prime', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % normalization 
    parse(ip, varargin{:});
    ip = ip.Results;

    a = src.radius;

    ell = permute((0:ell_max).', [4, 2, 3, 1]);

    if ip.is_farfield
        % dim: ell -> r
        R = CircSrc_SWE_RadialInt(src, 2*ell(:)+abs(src.prf.azimuth_order), ...
            'j', 0, a, r.', 'int_num', ip.int_num, 'is_farfield', ip.is_farfield);
        % dim: r -> 1 -> 1 -> ell
        R = permute(permute(R, [4,2,3,1]), [2,1,3,4]);
    else
        % origin points
        idx_origin = (r == 0);
        r_origin = r(idx_origin);
        % interior points
        idx_int = (r > 0) & (r < a);
        r_int = r(idx_int);
        % exterior points
        idx_ext = r >= a;
        r_ext = r(idx_ext);   

        % store radial component
        R = 0 * r .* ell;

        %% process origin points
        if ~isempty(r_origin)
            % dim: ell -> r_origin
            int = CircSrc_SWE_RadialInt(src, 2*ell(:)+abs(src.prf.azimuth_order), ...
                'h', 0, a, r_origin.', 'int_num', ip.int_num);
            % dim: r_origin -> 1 -> 1 -> ell
            R(idx_origin, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
        end

        %% process interior points
        if ~isempty(r_int)
            % dim: ell -> r_int
            int = CircSrc_SWE_RadialInt(src, 2*ell(:)+abs(src.prf.azimuth_order), ...
                'j', 0, r_int.', r_int.', 'int_num', ip.int_num) ...
                + CircSrc_SWE_RadialInt(src, 2*ell(:)+abs(src.prf.azimuth_order), ...
                'h', r_int.', a, r_int.', 'int_num', ip.int_num);
            % dim: r_int -> 1 -> 1 -> ell
            R(idx_int, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
        end

        %% process exterior points
        if ~isempty(r_ext)
            % dim: ell -> r_ext
            int = CircSrc_SWE_RadialInt(src, 2*ell(:)+abs(src.prf.azimuth_order), ...
                'j', 0, a, r_ext.', 'int_num', ip.int_num);
            % dim: r_ext -> 1 -> 1 -> ell
            R(idx_ext, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
        end
    end

    R_prime = nan;
end
