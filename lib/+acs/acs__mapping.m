function ACS = acs__mapping(preACS, mapping, output_type)

    % Transform preliminary ACS to final ACS quantity
    % using specified mapping
    %
    % Input: 
    %   mapping: "EXP"
    %   output_type: "cont" or "bin".
    % Output:
    %   ACS: final acoustic coherence state


    %% MAPPING FUNCTIONS
    function out = map_EXP(in, output_type)
        % Based on exp. curve via percentiles
        % (note: exp. curve is defined in log-domain)
        a = -0.636; b = -3.362; c = 23.4; %95th percentiles, chi0=0.0712
        %a = -0.6523; b= -5.0616; c = 24.3381; %based on libri-speech db
        basic_map = @(x) a - b.*exp(-c.*x);
        chi_max = 0.25; % is mapped to 1
        maxLogSROEstErr = 0; %1ppm
        out = basic_map(in) - maxLogSROEstErr;
        out = out./abs(basic_map(chi_max));
        out= -out; % ACS inversely prop. to SRO-Est-Err
        out(out < 0) = 0; 
        out(out > 1) = 1;
        if output_type == "bin"
            out(out > 0) = 1;
        end
    end
  

    %% OBTAIN ACS 
    if mapping == "EXP"
        ACS = map_EXP(preACS, output_type);
    else
        ACS = 0;
        error('Invalid Mapping Type.');
    end

end