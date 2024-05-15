% DIMENSION
%   - [fp.x <=> fp.y <=> fp.z] => vsrc.x => vsrc.y => vsrc.z

function [prs_tot, prs_cum, prs_local] = PalDIM3D(pal, fp, varargin)
    ip = inputParser;

    ip.addParameter('x_vsrc_int_num', 1e2);
    ip.addParameter('y_vsrc_int_num', 1e2);
    ip.addParameter('z_vsrc_int_num', 1e2);
    ip.addParameter('rho_vsrc_int_num', 1e2);
    ip.addParameter('phi_vsrc_int_num', 1e2);

    ip.addParameter('ultra_int_num', 1e2);

    ip.addParameter('is_incl_local', false);

    % integration range 
    ip.addParameter('x_vsrc', [-4,4]);
    ip.addParameter('y_vsrc', [-4,4]);
    ip.addParameter('z_vsrc', [-4,4]);
    ip.addParameter('rho_vsrc', [0,4]);
    ip.addParameter('phi_vsrc', [0,2*pi]);

    % The coordinates used for the integration
    ip.addParameter('int_coord', 'Cartesian', @(x)any(validatestring(x, {'Cartesian', 'spherical', 'cylindrical'})));
    ip.addParameter('ultra_int_coord', 'Cartesian', @(x)any(validatestring(x, {'Cartesian', 'polar'})));

    ip.parse(varargin{:});
    ip = ip.Results;

    % initial parameters
    beta = 1.2;
    rho0 = 1.21;
    c0 = 343;
    pal.audio.wav.angfreq = 2*pi*pal.audio.wav.freq;
    pal.audio.wav.num = pal.audio.wav.angfreq / c0;

    fprintf('Run PalDIM3D.m\n');
    prs_cum = 0;
    switch ip.int_coord
        case 'Cartesian'
            for ix = 1:length(ip.x_vsrc) - 1
                for iy = 1:length(ip.y_vsrc) - 1
                    for iz = 1:length(ip.z_vsrc) - 1
                        t0 = tic;
                        fprintf('Integrating vsrc: xyz = [%gm, %gm] x [%gm, %gm] x [%gm, %gm]...', ...
                            ip.x_vsrc(ix), ip.x_vsrc(ix+1), ip.y_vsrc(iy), ip.y_vsrc(iy+1), ip.z_vsrc(iz), ip.z_vsrc(iz+1));
                        
                        [vsrc.x, weight.x] = GaussLegendreQuadParam(ip.x_vsrc_int_num, ip.x_vsrc(ix), ip.x_vsrc(ix+1), 'dim', 4);
                        [vsrc.y, weight.y] = GaussLegendreQuadParam(ip.y_vsrc_int_num, ip.y_vsrc(iy), ip.y_vsrc(iy+1), 'dim', 5);                
                        [vsrc.z, weight.z] = GaussLegendreQuadParam(ip.z_vsrc_int_num, ip.z_vsrc(iz), ip.z_vsrc(iz+1), 'dim', 6);
        
                        prs_cum = prs_cum - beta .* pal.audio.wav.angfreq^2 / (4 * pi * rho0 * c0^4) ...
                                    .* sum(sum(sum(Integrand(pal, fp, vsrc, ip) .* weight.z, 6) .* weight.y, 5) .* weight.x, 4);
                
                        fprintf(' Finished. Elpased time: %g s.\n', toc(t0));
                    end
                end
            end
        case 'cylindrical' % todo
            for irho = 1:length(ip.rho_vsrc) - 1
                for iphi = 1:length(ip.phi_vsrc) - 1
                    for iz = 1:length(ip.z_vsrc) - 1
                        t0 = tic;
                        fprintf('Integrating vsrc: rho, phi, z = [%gm, %gm] x [%g, %g] x [%gm, %gm]...', ...
                            ip.rho_vsrc(irho), ip.rho_vsrc(irho+1), ip.phi_vsrc(iphi), ip.phi_vsrc(iphi+1), ip.z_vsrc(iz), ip.z_vsrc(iz+1));
                        
                        [vsrc.rho, weight.rho] = GaussLegendreQuadParam(ip.rho_vsrc_int_num, ip.rho_vsrc(irho), ip.rho_vsrc(irho+1), 'dim', 4);
                        [vsrc.phi, weight.phi] = GaussLegendreQuadParam(ip.phi_vsrc_int_num, ip.phi_vsrc(iphi), ip.phi_vsrc(iphi+1), 'dim', 5);                
                        [vsrc.z, weight.z] = GaussLegendreQuadParam(ip.z_vsrc_int_num, ip.z_vsrc(iz), ip.z_vsrc(iz+1), 'dim', 6);
        
                        [vsrc.x, vsrc.y] = Polar2Cart(vsrc.rho, vsrc.phi);
                        prs_cum = prs_cum - beta .* pal.audio.wav.angfreq^2 / (4 * pi * rho0 * c0^4) ...
                                    .* sum(sum(sum(Integrand(pal, fp, vsrc, ip) .* weight.z, 6) .* weight.phi, 5) .* vsrc.rho .* weight.rho, 4);
                
                        fprintf(' Finished. Elpased time: %g s.\n', toc(t0));
                    end
                end
            end
        case 'spherical'
    end

    prs_local = 0;
    if ip.is_incl_local
        % todo
    end
    prs_tot = prs_cum + prs_local;
end

function int = Integrand(pal, fp, vsrc, ip)    
    % calculate ultrasound pressure field
    vsrc0 = vsrc;
    vsrc0.x = permute(vsrc0.x, [4,5,6,1,2,3]);
    vsrc0.y = permute(vsrc0.y, [4,5,6,1,2,3]);
    vsrc0.z = permute(vsrc0.z, [4,5,6,1,2,3]);
    prs1 = DIM3D(pal.ultra1, vsrc0, 'int_num', ip.ultra_int_num, 'int_coor', ip.ultra_int_coord);
    prs2 = DIM3D(pal.ultra2, vsrc0, 'int_num', ip.ultra_int_num, 'int_coor', ip.ultra_int_coord);
    prs1 = permute(prs1, [4,5,6,1,2,3]);
    prs2 = permute(prs2, [4,5,6,1,2,3]);

    dist = sqrt((fp.x - vsrc.x).^2 + (fp.y - vsrc.y).^2 + (fp.z - vsrc.z).^2);
    int = conj(prs1) .* prs2 .* exp(1i * pal.audio.wav.num .* dist) ./ dist;
end
