function [sampleSeries] = seriesFramesToSamples(frameSeries, frameSize, frameShift, lenSmp)

    % Translates series of values where each value corresponds to one frame
    % of some origin signal to a series where each value corresponds to 
    % one sample of that signal (values are repeated)
    %
    % Input:
    %   frameSeries: Series of values where each corresponds to one frame
    %   frameSize: used frame size
    %   frameShift: used frame shift
    %   lenSmp: length of origin signal in samples
    % Output:
    %   sampleSeries: Equivalent series of values, each corresponding to
    %                   one sample of the origin signal.

    sigLen_smp = lenSmp;
    sampleSeries = zeros(1, sigLen_smp);

    for ii = 1:length(frameSeries)
        from = (ii-1)*frameShift +1;
        to = from + frameShift;
        if ii == length(frameSeries)
            to = from + frameSize;
        end
        sampleSeries(from:to) = frameSeries(ii);
    end

end




