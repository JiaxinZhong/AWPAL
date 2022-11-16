% =========================================================================
% INTRO
%   - Calculate the profile for a baffled line source in a 2D problem
%       when using the CWE method
% -------------------------------------------------------------------------
% INPUT
%   - m: order
%   - k: wavenumber
%   - a: half-width (radius)
%   - rho_src: the polar coordinates of source point
%   - profile: the profile type
% -------------------------------------------------------------------------
% OUTPUT
%   - u: the profile
% -------------------------------------------------------------------------
% DIMENSION
%   - 1: m
%   - 2: rho_src
% =========================================================================
function u = LineSrc_CWE_Profile(m, k, a, rho_src, prf, varargin)

%     ip = inputParser;
%     ip.addParameter('steer_angle', nan);
%     parse(ip, varargin{:});
%     ip = ip.Results;

    % validateattributes(m, {'numeric'}, {'size', [nan,1]});

    u = (LineSrcProfile(a, k, rho_src, prf) ...
        + (-1).^m .* LineSrcProfile(a, k, -rho_src, prf))/2; 
end
