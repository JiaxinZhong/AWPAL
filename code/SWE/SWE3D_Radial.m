% =========================================================================
% INTRO
%   - Calculate the radial component for a baffled acoustic source using
%       the spherical wave expansion (SWE) method
% -------------------------------------------------------------------------
% INPUT
%   - src: acoustic source info
%       - src.radius: half-width (radius)
%   - r: radial coordinate
%   - ell_max: maximum order
% -------------------------------------------------------------------------
% OUTPUT
%   - R: the radial component 
%   - R_prime: the derivative of the radial component
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: r
%   - 2: 1 (holder for theta)
%   - 3: 1 (holder for phi)
%   - 4: order ell
%   - 5: integration
% =========================================================================
function [R, R_prime] = SWE3D_Radial(src, r, ell_max, varargin)

    ip = inputParser;
    ip.addParameter('int_num', []);
    ip.addParameter('is_cal_prime', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    % normalization 
    parse(ip, varargin{:});
    ip = ip.Results;

    ell = permute((0:ell_max).', [4, 2, 3, 1]);
    if isfield(src.prf, 'azimuth_order')
        ell_m_col =  2*ell(:)+abs(src.prf.azimuth_order);
    else
        ell_m_col =  2*ell(:);
    end
    % store radial component
    R = 0 * r .* ell;

%     if src.r_is_discrete % the source is discrete in the radial direction
%         for i = 1:length(src.pos.r)
%             src_r = src.pos.r(i);
%             r_row = r.';
%             src_r_big = max(src_r, r_row);
%             src_r_small = min(src_r, r_row);
%             R_tmp = exp(2*log(src.wav.num) + log(src_r) ...
%                 + SphBesselJ(ell_m_col, src.wav.num*src_r_small, ...
%                     'is_log', true) ...
%                 + SphHankelH(ell_m_col, src.wav.num*src_r_big, ...
%                 'is_log', true, 'arg_is_large', ip.is_farfield));
%             R = R + permute(permute(R_tmp, [4,2,3,1]), [2,1,3,4]);
%         end
%     else % the source is continuous in the radial direction
        if ip.is_farfield
            % dim: ell -> r
            R = SWE3D_RadialInt(src, ell_m_col, ...
                'j', 0, src.r, r.', 'int_num', ip.int_num, 'is_farfield', ip.is_farfield);
            % dim: r -> 1 -> 1 -> ell
            R = permute(permute(R, [4,2,3,1]), [2,1,3,4]);
        else
            % origin points
            idx_origin = (r == 0);
            r_origin = r(idx_origin);
            % interior points
            idx_int = (r > 0) & (r < src.r);
            r_int = r(idx_int);
            % exterior points
            idx_ext = r >= src.r;
            r_ext = r(idx_ext);   

            % process origin points
            if ~isempty(r_origin)
                % dim: ell -> r_origin
                int = SWE3D_RadialInt(src, ell_m_col, ...
                    'h', 0, src.r, r_origin.', 'int_num', ip.int_num);
                % dim: r_origin -> 1 -> 1 -> ell
                R(idx_origin, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
            end

            % process interior points
            if ~isempty(r_int)
                % dim: ell -> r_int
                int = SWE3D_RadialInt(src, ell_m_col, ...
                    'j', 0, r_int.', r_int.', 'int_num', ip.int_num) ...
                    + SWE3D_RadialInt(src, ell_m_col, ...
                    'h', r_int.', src.r, r_int.', 'int_num', ip.int_num);
                % dim: r_int -> 1 -> 1 -> ell
                R(idx_int, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
            end

            % process exterior points
            if ~isempty(r_ext)
                % dim: ell -> r_ext
                int = SWE3D_RadialInt(src, ell_m_col, ...
                    'j', 0, src.r, r_ext.', 'int_num', ip.int_num);
                % dim: r_ext -> 1 -> 1 -> ell
                R(idx_ext, 1, 1, :) = permute(permute(int, [4,2,3,1]), [2,1,3,4]);
            end
        end
%     end

    % dim: ell, 1, 1, 1
    Y0 = SphHarmonic(2*ell(:)+abs(src.prf.phi_m), src.prf.phi_m, pi/2, 0);
    % dim: 1, 1, 1, ell
    Y0 = permute(Y0, [4,2,3,1]);
    R = 4*pi * R .* Y0;

    %% Calculate the prime of the radial component
    R_prime = [];
    if ip.is_cal_prime
        
    end
end
