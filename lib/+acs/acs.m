function [ACS, preACS, preACS_unrestricted, dACS] = acs(GCSD2_avg, mapping, output_type)

    % Compute ACS quantity
    %
    % Input:
    % GCSD2_avg: Averaged GCSD2 up until Nyq. Freq. Bin
    % mapping: which mapping preACS->ACS to use. ("EXP")
    % output_type: "cont" (continuous) or "bin" (binary)
    %   binary output is transformed via thresholding (rules may differ depending on mapping)

    preACS_unrestricted = mean(abs(GCSD2_avg), 1); % no SFM
    preACS = preACS_unrestricted; %redundant
    dACS = 0; %redundant
    
    % Obtain ACS
    ACS = acs.acs__mapping(preACS_unrestricted, mapping, output_type);   

end