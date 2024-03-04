function [x, y, z] = Sph2Cart(r, theta, phi)
    
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);

end