% =========================================================================
% INTRO
%   - Calculate the sound pressure generated by a 3D source using
%       the spherical wave expansion (SWE)
% -------------------------------------------------------------------------
% INPUT
%   - src: data of the circular source 
%   - fp.r, fp.theta, fp.phi: the radial, zenithal, and azimuthal
%       coordinates of the field points
% OUTPUT
%   - prs: normalized to rho0*c0*v0
% -------------------------------------------------------------------------
% DIMENSION
%   fp.r -> fp.theta -> fp.phi -> ell -> int
% =========================================================================
function [prs, vel] = SWE3D(src, fp, varargin)

    ip = inputParser();
    ip.addParameter('int_num', 2e2, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 2}));
    ip.addParameter('is_cal_vel', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('is_farfield', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.addParameter('max_order', nan, ...
        @(x)validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
    ip.addParameter('is_norm', false)
    % on-surface pressure amplitude
    % ip.addParameter('p0', 1, @(x)validateattributes(x, {'numeric'}, {'scalar'}));
    ip.parse(varargin{:});
    ip = ip.Results;

%     if isnan(ip.max_order)
%         switch src.shape
%             case 'sphere'
                ip.max_order = ceil(src.r*real(src.wav.num)*1.2);
%             case 'point'
%                 ip.max_order = ceil(max(src.pos.r)*real(src.wav.num)*1.2);
%         end
%     end

    %% Calculate angular components: spherical harmonics
    % dim: 1, 1, 1, ell
    ell = permute((0:ip.max_order).', [4,2,3,1]); % order
    % dim: ell, 1, 1, 1
%     Y0 = SphHarmonic(2*ell(:)+abs(prf_azimuth_order), ...
%         prf_azimuth_order, pi/2, 0);
    % dim: 1, 1, 1, ell
%     Y0 = permute(Y0, [4,2,3,1]);
    % dim: ell, 1, theta, phi
    Y = SphHarmonic(2*ell(:)+abs(src.prf.phi_m), ...
        src.prf.phi_m, permute(fp.theta,[1,3,2]), ...
        permute(fp.phi, [1,2,4,3]));
    % dim: 1, theta, phi, ell
    Y = permute(permute(permute(Y, [1,3,2,4]), [1,2,4,3]), [4,2,3,1]);

    %% Calculate radial components
    R = SWE3D_Radial(src, fp.r, ip.max_order, ...
        'int_num', ip.int_num, ...
        'is_farfield', ip.is_farfield);

    %% Calculate the sound pressure
%     prs = 4*pi * sum(Y0 .* Y .* R, 4);
    prs = sum(Y .* R, 4);

    if ip.is_norm
        prs = prs ./ max(abs(prs(:)));
    end

    %% calculate velocity
    % todo
    vel.r = [];
    vel.theta = [];
    vel.phi = [];
    if ip.is_cal_vel
        ;
    end


end