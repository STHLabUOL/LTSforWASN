function [y] = resampler_sinc_on(x, WindowSize, SRO_ppm, Delay)
% Resample single channel audio data x with sub-herz precision using
% precise SINC interpolation of data and sampling theorem.
%
% 05.11.2016 Johannes Deuse-Kleinsteuber
%       Added param init, if they are not passed by a call
%       Eleminated a for-loop for Window sizes > 50
%
% :param x: signal (vector)
% :param WindowSize: Size of window (scalar e.g. 100)
% :param SRO_ppm: SRO_ppm of sampling-rate (scalar in ppm)
% :param Delay: Constant time offset of signal (scalar in signal taps)
% :return y: resampled data

    y = 0*x;
    x_len_min1 = length(x)-1;
    ratio =  (1+SRO_ppm/10^6);
    % Delay = Delay/fsbase; AC
    
    hh=hann(2*WindowSize+1);
    % better: parfor (if supported)
    for k=0:x_len_min1
        k_ratio_rnd = round(ratio*k);
        IndexFirst = max([0, k_ratio_rnd - WindowSize]);
        IndexLast = min([x_len_min1, k_ratio_rnd + WindowSize]);

        m = IndexFirst:IndexLast;
        y(k+1) = sum(x(m+1).*sinc(k*ratio-m-Delay).*hh(end-length(m)+1:end)');
    end
end