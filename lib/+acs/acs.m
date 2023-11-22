function [ACS, preACS, preACS_unrestricted, dACS] = acs(GCSD2_avg, mapping, output_type)

    % Compute acoustic coherence state (ACS) quantities
    % based on averged secondary generalized cross spectral density (GCSD2)
    %
    % Input:
    %   GCSD2_avg: Averaged GCSD2 up until Nyq. Freq. Bin
    %   mapping: which mapping preACS->ACS to use. ("EXP")
    %   output_type: "cont" (continuous) or "bin" (binary)
    %               binary output is transformed via 
    %               thresholding (rules may differ depending on mapping)
    %
    % Output:
    %   ACS: Acoustic coherence state
    %   preACS: preliminary ACS, before mapping is applied
    %   preACS_unrestricted: preACS before temporal restrictions  are
    %                       applied
    %   dACS: differential ACS

    preACS_unrestricted = mean(abs(GCSD2_avg), 1); % no SFM
    preACS = preACS_unrestricted; %redundant
    dACS = 0; %redundant
    
    % Obtain ACS
    ACS = acs.acs__mapping(preACS_unrestricted, mapping, output_type);   

end