classdef PointSrc < SoundSrc

%     properties (Dependent)
%         n_elem % number of elements
%     end

    properties (Constant)
        shape = 'point'
    end

    methods
        function obj = PointSrc(varargin)
            
%             ip = inputParser();
%             ip.addParameter('wav', []);
%             ip.addParameter('freq', []);
%             ip.addParameter('pos', []);
%             ip.addParameter('Q', 1);
%             ip.addParameter('dir', Point3D('phi', 0, 'theta', 0));
%             ip.parse(varargin{:});
%             ip = ip.Results;

            obj = obj@SoundSrc(varargin{:});
        end

%         function val = get.n_elem(obj)
%             val = numel(obj.pos);
%         end
    end
end