% sound source
classdef SoundSrc < handle
    properties
        wav     % the wave info
        pos     % position of the source
        dir     % direction of the source
        Q       % source strength
        pos_centroid % position of the centroid
    end

    properties (Dependent)
    end

    methods
        function obj = SoundSrc(varargin)
            ip = inputParser();
            ip.addParameter('wav', []);
            ip.addParameter('freq', []);
            ip.addParameter('pos', []);
            ip.addParameter('pos_centroid', []);
            ip.addParameter('Q', 1);
            ip.addParameter('dir', Point3D('phi', 0, 'theta', 0));
            ip.parse(varargin{:});
            ip = ip.Results;

            if ~isempty(ip.freq)
                obj.wav = SoundWave('freq', ip.freq);
            else
                obj.wav = ip.wav;
            end
            obj.pos = ip.pos;
            obj.Q = ip.Q;
            obj.dir = ip.dir;
            obj.pos_centroid = ip.pos_centroid;
        end

    end
end
