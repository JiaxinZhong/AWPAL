% a rectangular source
classdef RectSrc < SoundSrc
    properties
        ax  % half-width in x directions
        ay  % half-width in y directions
    end

    properties (Constant)
        shape = 'rectangle';
    end

    methods
        function obj = RectSrc(varargin)
            var_list = {'ax', 'ay'};
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
            ip.addParameter('ax', []);
            ip.addParameter('ay', []);
            ip.parse(var_here{:});
            ip = ip.Results;

            obj = obj@SoundSrc(var_parent{:});
            obj.ax = ip.ax;
            obj.ay = ip.ay;
        end

        function dir = CalDirectivity(obj, theta, phi)
            theta = theta + 0*phi;
            phi = phi + 0*theta;
            k = real(obj.wav.num);
            kx = k*sin(theta).*cos(phi);
            ky = k*sin(theta).*sin(phi);
            switch obj.prf.name
                case 'uniform'
                    dir = sinc(kx.*obj.ax/pi) .* sinc(ky.*obj.ay/pi);
                case 'steerable'
                    kx0 = k*sin(obj.prf.theta).*cos(obj.prf.phi);
                    ky0 = k*sin(obj.prf.theta).*sin(obj.prf.phi);
                    dir = sinc((kx - kx0).*obj.ax/pi) .* sinc((ky - ky0).*obj.ay/pi);
                otherwise
                    error('Wrong profile type!')
            end
        end
    end
end

