function [sampleSeries] = seriesFramesToSamples(frameSeries, frameSize, frameShift, lenSmp)

    %
    % Translates series of values where each sample represents one frame
    % to a corresponding sample-wise series by repeating values.
    %
    %
    % Note: This was created specifically for evaluation 11.10.22
    % and probably needs a rework.

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




