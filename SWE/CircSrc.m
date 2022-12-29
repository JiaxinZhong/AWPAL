% a circular source in a 3D problem
classdef CircSrc < SoundSrc
    properties
        radius  % radius of the circular source
    end

    properties (Dependent)
        elem_num % total number of array elements
    end

    properties (Constant)
        shape = 'circle';
    end

    methods
        function obj = CircSrc(varargin)
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

        function u = CalProfile(obj, rs)
            switch obj.prf.name
                case 'uniform' 
                    u = 1;
                case 'quadratic'
                    u = (obj.prf.order + 1) .* (1 - (rs/obj.radius).^2).^obj.prf.order;
                case 'Zernike'
                    u = ZernikeRadial(obj.prf.degree, obj.prf.azimuth_order, rs/obj.radius);
                otherwise
                    error('Wrong profile type!')
            end
        end

        function dir = CalDirectivity(obj, theta, phi)
            theta = theta + 0*phi;
            phi = phi + 0*theta;
            switch obj.prf.name
                case 'uniform'
                    dir = Jinc(real(obj.wav.num)*obj.radius.*sin(theta));
                case 'steerable'
                    dir = Jinc(real(obj.wav.num) * obj.radius ...
                        .* sqrt((cos(obj.prf.phi) .* sin(obj.prf.theta) ...
                        - cos(phi) .* sin(theta)).^2 ...
                        + (sin(obj.prf.phi) .* sin(obj.prf.theta) ...
                        - sin(phi) .* sin(theta)).^2));
                case 'quadratic'
                    dir = exp((obj.prf.order + 1) .* log(2) ...
                        + gammaln(obj.prf.order+2) ...
                        - (obj.prf.order+1) .* log(obj.wav.num.*obj.radius.*sin(theta)) ...
                        + BesselJLog(obj.prf.order+1, obj.wav.num .* obj.radius.*sin(theta)));
                    dir(theta==0) = 1;
                otherwise
                    error('Wrong profile type!')
            end
        end

    end
end
