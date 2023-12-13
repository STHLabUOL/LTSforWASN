import numpy as np

def seriesFramesToSamples(frameSeries, frameSize, frameShift, lenSmp):

    sigLen_smp = lenSmp
    sampleSeries = np.zeros((sigLen_smp,))

    for ff, val in enumerate(frameSeries):
        start = ff*frameShift
        end = start + frameShift + 1
        if ff+1 == np.size(frameSeries):
            end = start + frameSize
        
        sampleSeries[start:end] = val

    return sampleSeries