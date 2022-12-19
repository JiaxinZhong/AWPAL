classdef SrcProfile < handle
    properties
        name
        degree
        azimuth_order
        order
        theta
        phi
    end

    methods
        function obj = SrcProfile(varargin)
            ip = inputParser();
            ip.addParameter('name', 'uniform');
            ip.addParameter('azimuth_order', 0);
            ip.addParameter('degree', 0);
            ip.addParameter('order', 0);
            ip.addParameter('theta', 0);
            ip.addParameter('phi', 0);
            ip.parse(varargin{:});
            ip = ip.Results;
            
            obj.name = ip.name;
            obj.azimuth_order = ip.azimuth_order;
            obj.order = ip.order;
            obj.theta = ip.theta;
            obj.phi = ip.phi;
            obj.phi = ip.degree;
        end
    end
end
