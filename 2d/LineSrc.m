% A line source in the 2D problem
classdef LineSrc < SoundSrc
    properties
        radius % radius of the circular source
    end

    properties (Dependent)
        elem_num % total number of array elements
    end

    properties (Constant)
        shape = 'line';
    end

    methods
        function obj = LineSrc(varargin)
            var_list = {'radius'};
            var_here = {};
            cnt_here = 0;
            var_parent = {};
            cnt_parent = 0;
            
            for i = 1:length(varargin)/2
                if cell2mat(strfind(var_list, varargin{2*i-1})) == 1
                    cnt_here = cnt_here + 1;
                    var_here{2*cnt_here-1} = varargin{2*i-1};
                    var_here{2*cnt_here} = varargin{2*i};
                else
                    cnt_parent = cnt_parent + 1;
                    var_parent{2*cnt_parent-1} = varargin{2*i-1};
                    var_parent{2*cnt_parent} = varargin{2*i};
                end
            end

            ip = inputParser();
            ip.addParameter('radius', []);
            ip.parse(var_here{:});
            ip = ip.Results;

            obj = obj@SoundSrc(var_parent{:});
            obj.radius = ip.radius;
        end

        function u = CalProfile(obj, ys)
            a = obj.radius;
            k = real(obj.wav.num);
            switch obj.prf.name
                case 'uniform'
                    u = 1;
                case 'steerable'
                    % u = exp(1i*k*xs.*cos(obj.prf.phi));
                    u = exp(1i*k*ys.*sin(obj.prf.phi));
                case 'cosine'
                    u = (cos(pi*ys/2./a)).^obj.prf.order / 2 * pi;
                case 'cosine_steerable'
                    u = (cos(pi*xs/2./a)).^obj.prf.order ...
                        .* exp(1i*k*xs.*cos(obj.prf.phi));
                case 'hanning'
                    u = (cos(pi*xs/2./a)).^2;
                case 'hamming'
                    u = (0.54 + 0.46 .* cos(pi*xs./a));
                case 'blackman'
                    u = 0.42 + 0.5 .* cos(pi.*xs./a) + 0.08 .* cos(2*pi*xs./a);
                case 'triangular'
                    u = 1 - abs(xs) ./ a;
                otherwise
                    error('Wrong profile!')
            end
        end

        function dir = CalDirectivity(obj, phi)
            ky = real(obj.wav.num) .* sin(phi);
            switch obj.prf.name
                case 'uniform'
                    dir = sinc(ky .* obj.radius/pi);
                case 'steerable'
                    dir = sinc((ky - real(obj.wav.num) .* sin(obj.prf.phi)) .* obj.radius/pi);
                case 'cosine'
                    switch obj.prf.order
                        case 1
%                             dir = 2*obj.radius .* cos(ky.*obj.radius) ./ (1 - (2*ky.*obj.radius/pi).^2);
                            dir = obj.radius .* (sinc(ky.*obj.radius/pi+1/2) + sinc(ky.*obj.radius/pi-1/2));
                        case 2
                            dir = 2*sinc(ky .* obj.radius/pi) ...
                                + (sinc(ky*obj.radius/pi + 1) + sinc(ky.*obj.radius/pi-1));
                        otherwise
                            error('Wrong cosine type!')
                    end
                case 'cosine_steerable'
                    ky = ky - real(obj.wav.num) .* cos(obj.prf.phi);
                    switch obj.prf.order
                        case 1
                            dir = 2*obj.radius .* cos(ky.*obj.radius) ./ (1 - (2*ky.*obj.radius/pi).^2);
                        case 2
                            dir = 2*sinc(ky .* obj.radius/pi) ...
                                + (sinc(ky*obj.radius/pi + 1) + sinc(ky.*obj.radius/pi-1));
                        otherwise
                            error('Wrong cosine type!')
                    end
                case 'hanning'
                    dir = sinc(ky.*obj.radius/pi) ./ (1-(ky.*obj.radius/pi).^2);
                case 'hamming'
                    dir = sinc(ky.*obj.radius/pi) .* (1.08 - 0.16*(ky.*obj.radius/pi).^2) ...
                        ./ (1-(ky.*obj.radius/pi).^2) / 1.08;
                case 'blackman'
                    dir = sinc(ky.*obj.radius/pi) .* (0.84 + (ky.*obj.radius/pi).^2 ...
                        ./ (1 - (ky.*obj.radius/pi).^2) - .16.*(ky.*obj.radius/2/pi).^2 ...
                        ./ (1 - (ky.*obj.radius/2/pi).^2)) / 0.84;
                case 'triangular'
                    dir = (sinc(ky.*obj.radius/2/pi)).^2;
                otherwise
                    error('Wrong profile type!')
            end
        end
    end
end
