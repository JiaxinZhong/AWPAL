% =========================================================================
% INTRO
%   - The Log of Hankel function of first kind using the Rothwell's method
% INPUT
%   - M: maximum order
% =========================================================================

function [H, H_prime] = HankelHLog_Rothwell(M, z, varargin)

    ip = inputParser;
    ip.addParameter('kind', 1, @(x)validateattributes(x, {'numeric'}, {'scalar', '>=', 1, '<=', 2}));
    ip.addParameter('nu0', 0, @(x)validateattributes(x, {'numeric'}, ...
        {'scalar', '>=', 0, '<', 1}));
    ip.addParameter('is_cal_derivative', false, ...
        @(x)validateattributes(x, {'logical'}, {'scalar'}));
    ip.parse(varargin{:})
    ip = ip.Results;

    H_prime = [];

    m = (0:M).';
    z_row0 = z(:).';
    [z_row, ~, idx_z_row] = unique(z_row0);

    J = 0 * m .* z_row;
    H = 0 * m .* z_row;
    R = 0 * m .* z_row;
            if ip.is_cal_derivative
    H_prime = 0 * m .* z_row;
            end

%     RH = 0 * m .* z_row;

    R(end, :) = ContFracLentz(@(m, z_row) RatioA(m, z_row), ...
        @(m, z_row) RatioB(m, z_row, M, ip.nu0), z_row,...
        'itr_max', max(10,ceil(max(abs(z_row)*3))));

    for mm = M:-1:1
        R(mm, :) = 2 * (mm + ip.nu0 - 1) ./ z_row - 1./R(mm+1, :);
    end

    J(1, :) = log(besselj(ip.nu0, z_row, 1)) + abs(imag(z_row));
    for mm = 1:M
        J(mm+1,:) = J(mm,:) - log(R(mm+1,:));
    end

    H(1,:) = log(besselh(ip.nu0, ip.kind, z_row, 1)) + (-1).^(ip.kind+1) .*1i.*z_row;
    % H(1,:) = log(besselh(ip.nu0, 1, z_row, 1)) + 1i.*z_row;
    for mm = 1:M
        RH = 1./(1./R(mm+1,:) +(-1)^ip.kind *2*1i./(pi .* z_row) .* exp(-J(mm,:) - H(mm,:)));
        H(mm+1,:) = H(mm,:) - log(RH);
        if ip.is_cal_derivative
            H_prime(mm+1, :) = H(mm+1,:) + log(RH - (mm+ip.nu0) ./ z_row);
        end
    end
    
    if ip.is_cal_derivative
        H_prime(1, :) = log(-exp(H(2,:)));
    end
    

%     J(1,z_row==0) = log(1);
%     J(2:end,z_row==0) = log(0);

%     J = J(:, idx_z_row);
    H = H(:, idx_z_row);
%     J = reshape(J, size(m .* z_row0));
    H = reshape(H, size(m .* z_row0));
    if ip.is_cal_derivative
        H_prime = H_prime(:, idx_z_row);
        H_prime = reshape(H_prime, size(m .* z_row0));
    end
end

function res = RatioA(m, z)
    res = 1;
end

function res = RatioB(m, z, M, nu0)
    res = 2 * (-1).^m .* (M + m + nu0) ./ z;
end
