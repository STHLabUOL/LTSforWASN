function parDXCPPhaT = DXCP_PhaT_DefaultParams()

    % Returns default parameters for DXCP-PhaT.
    % Reduces clutter when using DXCP-PhaT module.

    fs_Hz = 16e3;

    parDXCPPhaT.B_sec = 5;                  % accumulation time / frame shift in sec
    parDXCPPhaT.FFTsize = 2^13;             % FFT size (also frame size)
    parDXCPPhaT.FFTshift = 2^11;            % frame shift (FFTsize/2)
    
    parDXCPPhaT.B_NumFr = floor(parDXCPPhaT.B_sec*fs_Hz/parDXCPPhaT.FFTshift);  % number of frames acc. to desired B
    parDXCPPhaT.Lambda = 50;    % maximum lag of the second xcorr
    parDXCPPhaT.window = blackman(parDXCPPhaT.FFTsize, 'periodic');
    
    parDXCPPhaT.alpha1 = .5;               % smoothing constant of the first CSD (.53 for FFTsize = 2^13 & FFTshift = 2^11)
    % parDXCPPhaT.alpha2 = 0;                 % smoothing constant of the second CSD (.998)
    parDXCPPhaT.SetCSD2avg_NumFr = parDXCPPhaT.B_NumFr/2;               % settling time of CSD2 averaging (B_NumFr/2)
    parDXCPPhaT.Z_12_abs_min = 1e-12;        % minimum value of |X1*conj(X2)| to avoid devision by 0 in GCC-PhaT
    % parDXCPPhaT.al_d12 = .995;               % smoothing constant of SRO-compensated CCF-1 (used to estimate d12)
    parDXCPPhaT.p_upsmpFac = 4; % upsampling factor
    
    parDXCPPhaT.flag_plot = 0;
    parDXCPPhaT.alpha2 = .99;
    parDXCPPhaT.flag_avgXcSec = 1;   % 1 - recursive averaging, 0 - simple averaging (w.o. alpha2)
    parDXCPPhaT.alpha3 = .99;
    parDXCPPhaT.flag_avgCSDcomp = 1; % 1 - recursive averaging, 0 - simple averaging (w.o. alpha2)
    
    parDXCPPhaT.flag_STOest = 1;        % 1 - perform STO estimation, or 0 - not
    parDXCPPhaT.flag_calcFigOut = 1;    % 1 - calculate GCCFs for plots or 0 - mot
    parDXCPPhaT.flag_genWeight = 1;     % 1 - PhaT weighting or 0 - no weighting (conventional cross-sorrelation)

end
