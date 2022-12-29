%% Test convergence
clear all

prf = SrcProfile('name', 'uniform');
% prf = SrcProfile('name', 'quadratic', 'order', 1);
src = CircSrc('radius', .2, 'prf', prf);
pal = PalSrc('audio_freq', 4e3, 'ultra_freq', 40e3, 'src', src);

r = 1;
r1 = 2;
r2 = 5;

la_max = 70;
l1_max = ceil(real(pal.ultra_low.num)*pal.src_ultra.radius*1.2);
l2_max = ceil(real(pal.ultra_high.num)*pal.src_ultra.radius*1.2);

tic
sph = 'j';
R0 = PalCircSrc_SWE_RadialInt(pal, la_max, l1_max, l2_max, r1, r2, r, sph,...
    'is_farfield', true);
toc
if 1
   %% Coefficient
   tic
    la = permute((0:la_max).', [4,2,3,1]);
    l1 = permute((0:l1_max).', [5,2,3,4,1]);
    l2 = permute((0:l2_max).', [6,2,3,4,5,1]);
    m1 = pal.src_low.prf.azimuth_order;
    m2 = pal.src_high.prf.azimuth_order;
    ma = m2 - m1;
	wigner = 0 * l1 .* l2 .* la;
	for i1 = 1:length(l1)
		for i2 = 1:length(l2)
			for ia = 1:length(la)
% 				wigner(1, 1, 1, ia, i1, i2) = Wigner3j000(...
%                     2*l1(i1) + abs(m1), 2*l2(i2) + abs(m2), 2*la(ia) + abs(ma))...
%                     .* Wigner3j(...
%                     2*l1(i1) + abs(m1), 2*l2(i2) + abs(m2), 2*la(ia) + abs(ma),...
%                     -m1, m2, -ma);
                wigner(1, 1, 1, ia, i1, i2) = Wigner3j000(...
                    2*l1(i1) + abs(m1), 2*l2(i2) + abs(m2), 2*la(ia) + abs(ma))...
                    .* Wigner3j(...
                    2*l1(i1) + abs(m1), 2*l2(i2) + abs(m2), 2*la(ia) + abs(ma),...
                    -m1, m2, -ma);
			end
		end
    end
    toc

    % dim: l1
    Y1 = SphHarmonic(2*l1(:) + abs(m1), m1, pi/2, 0);
    % dim: 1 -> 1 -> 1 -> 1 -> l1
    Y1 = permute(Y1, [5,2,3,4,1]);
    % dim: l2
    Y2 = SphHarmonic(2*l2(:) + abs(m2), m2, pi/2, 0);
    % dim: 1 -> 1 -> 1 -> 1 -> 1 -> l2
    Y2 = permute(Y2, [6,2,3,4,5,1]);
    A = (-1).^m2 .* Y1 .* Y2 .* wigner ...
        .* sqrt((4*l1+2*abs(m1)+1).*(4*l2+2*abs(m2)+1).*(4*la+2*abs(ma)+1)./4/pi);
    
    R = R0.*A;
end
R = squeeze(R);

la_idx = 44;
R_la = squeeze(R(la_idx, :, :));
A_la = squeeze(A(1,1,1,la_idx, :, :));
wigner_la = squeeze(wigner(1,1,1,la_idx,:,:));

fig = Figure;
pcolor((0:l1_max), (0:l2_max), (abs(real(R_la.'))));
fig.Init;
% fig2 = Figure;
% pcolor((0:l1_max), (0:l2_max), imag(R_la.'));
% fig2.Init;