% INTRO
%   - Closed form solution of a baffled rigid circular piston
% OUTPUT
%   - prs: pressure, normalized to p0 = rho0*c0*v0
%   - vel: particle velocity, normalized to v0
function [prs, vel] = ClosedForm_CircPiston(src, fp)

    prs = -2*1i * sin(src.wav.num .* src.r / 2 .* (sqrt(1+fp.z.^2./src.r.^2) - fp.z./src.r)) ...
        .* exp(1i* src.wav.num .* src.r / 2 .* (sqrt(1+fp.z.^2./src.r.^2) + fp.z./src.r));
    vel.x = 0;
    vel.y = 0;
    vel.z = exp(1i * src.wav.num .* fp.z) ...
        - fp.z ./ src.r .* (1 + fp.z.^2 ./ src.r.^2).^(-1/2) ...
        .* exp(1i .* src.wav.num .* src.r .* sqrt(1 + fp.z.^2 ./ src.r.^2));
    
end