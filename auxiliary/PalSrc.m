classdef PalSrc < handle
    properties
        % precal
        ultra_high
        ultra_low
        ultra
        audio
        src_low
        src_high
        src_audio
        src_ultra
    end
    

    properties (Dependent)
        absorp_coef
    end
    
    methods 
        function obj = PalSrc(varargin)

            ip = inputParser();
            % for parent
            ip.addParameter('audio_freq', []);
            ip.addParameter('ultra_freq', 40e3);
            ip.addParameter('src', []);
            ip.parse(varargin{:});
            ip = ip.Results;

            obj.audio = SoundWave('freq', ip.audio_freq);
            obj.ultra_high = SoundWave('freq', ip.ultra_freq + ip.audio_freq/2);
            obj.ultra_low = SoundWave('freq', ip.ultra_freq - ip.audio_freq/2);
            obj.ultra = SoundWave('freq', ip.ultra_freq);
            
            obj.src_low = copy(ip.src);
            obj.src_low.wav = obj.ultra_low;
            
            obj.src_high = copy(ip.src);
            obj.src_high.wav = obj.ultra_high;
            
            obj.src_ultra = copy(ip.src);
            obj.src_ultra.wav = obj.ultra;
            
            obj.src_audio = copy(ip.src);
            obj.src_audio.wav = obj.audio;
            obj.src_audio.prf.order = 2*obj.src_audio.prf.order;
            
%             if ~isempty(ip.src)
%                 switch ip.src.shape
%                     case 'circle'
%                         prf_audio = SrcProfile('name', ip.src.prf.name, ...
%                             'order', 2*ip.src.prf.order, ...
%                             'azimuth_order', ip.src.prf.azimuth_order, ...
%                             'theta', ip.src.prf.theta, ...
%                             'phi', ip.src.prf.phi);
%                         obj.src_audio = CircSrc(...
%                             'radius', ip.src.radius, ...
%                             'wav', obj.audio, ...
%                             'prf', prf_audio);
%                     case 'line'
%                         prf_audio = SrcProfile('name', ip.src.prf.name, ...
%                             'order', 2*ip.src.prf.order, ...
%                             'azimuth_order', ip.src.prf.azimuth_order, ...
%                             'phi', ip.src.prf.phi);
%                         obj.src_audio = LineSrc(...
%                             'radius', ip.src.radius, ...
%                             'wav', obj.audio, ...
%                             'prf', prf_audio);
%                 end
%             end
            
        end

        function absorp_coef = get.absorp_coef(obj)
            absorp_coef = imag(obj.ultra_low.num + obj.ultra_high.num);
        end
    end
end
