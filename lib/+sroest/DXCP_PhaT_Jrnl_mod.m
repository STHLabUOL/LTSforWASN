function [SROest_ppm, SSOest_smp, GCCF1, GCCF2, GCCF2_avg, GCCF2_avg_upsmp, ...
    GCCF1_sroComp, GCCF1_sroComp_avg, GCSD2_avg_hist, GCSD2_hist] = DXCP_PhaT_Jrnl_mod(sig_z, parDXCPPhaT)

% Implementation of Double-Cross-Correlation Processor (DXCP) developed for
% journal paper 'Blind Synchronization in Acoustic Sensor Networks Using
% Double-Cross-Correlation Prozessor with Phase Transform (DXCP-PhaT)'.
% Generalized cross-correlation functions (CCFs) are calculated in STFT
% domain via generalized cross-spectral densities (GCSDs) in the spirit of
% the Generalized Cross-Correlation with Phase Transform (GCC-PhaT).
%   Sampling rate offset (SRO) in ppm (1e-6) and sampling time offset (SSO)
% in samples are estimated.
%   In this function DXCP_PhaT_DJ() some profitable implementation aspects
% are incorporated, which were obtained during application of DXCP-PhaT
% approach on Raspberry Pis used for development of a demonstrator deployed
% in real-world experiments. Abbreviation DJ means Demonstrator-motivated
% implemented for Journal paper.
%   Extensions:
%     - compensation of frequency dependent phase delay caused by recursive
%   averaging (mean calc.) of the primary CSD with alpha1 (AC, 16.01.2020)
% 
% :input sig_z:                     Two asynchronous signals [sig_lenx2]
% :param parDXCPPhaT:               Parameters of DXCPrec approach (see below)
% 
% :return SROest_ppm:               SRO estimates in ppm (vector)
% :return SSOest_smp:               SSO estimates in samples (vector)
% :return GCCF1:                    First GCCF (matrix)
% :return GCCF2:                    Second GCCF (matrix)
% :return GCCF2_avg:                Final averaged GCCF-2 (vector)
% :return GCCF2_avg_upsmp:          Final upsampled GCCF2_avg (vector)
% :return GCCF1_sroComp:            SRO-compensated GCCF-2 (matrix)
% :return GCCF1_sroComp_avg:        SRO-comp. averaged GCCF-2 (vector)
%
% written by Aleksej Chinaev (August 30, 2019)

    % external parameters
    Lambda = parDXCPPhaT.Lambda;     % maximum lag of the second xcorr
    B_NumFr = parDXCPPhaT.B_NumFr;   % number of frames acc. to desired B
    window = parDXCPPhaT.window;     % analysis window
    
    FFTsize = parDXCPPhaT.FFTsize;   % FFT size (frame length)
    FFTshift = parDXCPPhaT.FFTshift; % frame shift
    
    alpha1 = parDXCPPhaT.alpha1;     % smoothing constant of the first CSD
    alpha2 = parDXCPPhaT.alpha2;     % smoothing constant of the second CSD
    alpha3 = parDXCPPhaT.alpha3;     % smoothing constant of the SRO-compensated CSD
    
    SetCSD2avg_NumFr = parDXCPPhaT.SetCSD2avg_NumFr;    % settling time of CSD2 averaging
    Z_12_abs_min = parDXCPPhaT.Z_12_abs_min;            % to avoid devision by 0 in GCC-PhaT
    p_upsmpFac = parDXCPPhaT.p_upsmpFac; % upsampling factor
    
    flag_avgXcSec = parDXCPPhaT.flag_avgXcSec;
    flag_avgCSDcomp = parDXCPPhaT.flag_avgCSDcomp;
    
    flag_STOest = parDXCPPhaT.flag_STOest;              % 1 - perform STO estimation, or 0 - not
    flag_calcFigOut = parDXCPPhaT.flag_calcFigOut;      % 1 - calculate GCCFs for plots or 0 - mot
    flag_genWeight = parDXCPPhaT.flag_genWeight;        % 1 - PhaT weighting or 0 - no weighting (conventional cross-sorrelation)
    
    % internal parameters
    B = B_NumFr*FFTshift;
    Cont_NumFr = B_NumFr+1;                                 % number of frames in DXCP container
    DXCPPhaT_Cont_GCSD_avg = zeros(FFTsize, Cont_NumFr);    % DXCP container
    
    % further useful parameters
    [len_sig, num_chan] = size(sig_z);
    Upsilon = FFTsize/2-1;          % maximum lag of the first xcorr
    FrNum_L = floor((len_sig-FFTsize)/FFTshift)+1;          % number of signal frames for frame shift
    FFT_Nyq = FFTsize/2+1;
    lambda_vec_upsmp = -Lambda:1/p_upsmpFac:(Lambda+1-1/p_upsmpFac);
    
    % output variables
    SROest_ppm = zeros(1,FrNum_L);
    SSOest_smp = zeros(1,FrNum_L);
    GCCF1 = zeros(FrNum_L, 2*Upsilon+1);
    GCCF2 = zeros(FrNum_L, 2*Lambda+1);
    GCCF1_sroComp = zeros(FrNum_L, 2*Upsilon+1);
    GCCF1_sroComp_avg = zeros(1, 2*Upsilon+1);    
    
    % some variables for averaging
    SROppm_act = 0;
    GCSD_PhaT_avg = zeros(FFTsize, 1);
    GCSD2_avg = zeros(FFTsize, 1);
    GCSD2_avg_hist = zeros(FFTsize, FrNum_L); %maybe transpose
    GCSD2_hist = zeros(FFTsize, FrNum_L); % without rec. avg
    GCSD1_sroComp_avg = zeros(FFTsize, 1);
    
    % frame-based processing of the input signal
    index_vec = 1:FFTsize;
    ell_GCCF2avg = 1;
    ell_avg_sh = 1;
    for ell=1:FrNum_L
        
        % 1) Cutting and windowing the current frames acc. to eq. (5) in [1]
        z_12_win = sig_z(index_vec,:).*(window*ones(1,num_chan));
        
        % 2) Calculate generalized (normalized) GCSD with Phase Transform (GCSD-PhaT)
        % and apply recursive averaging
        Z_12 = fft(z_12_win, FFTsize, 1);
        Z_12_act = Z_12(:,1).*conj(Z_12(:,2));
        if flag_genWeight
            % PhaT weighting
            Z_12_act_abs = abs(Z_12_act);
            Z_12_act_abs(Z_12_act_abs<Z_12_abs_min) = Z_12_abs_min; % avoid division by 0
            GenWeights = 1./Z_12_act_abs;
        else % conventional cross-correlation
            GenWeights = 1;
        end
        GCSD_PhaT_act = GenWeights.*Z_12_act;
        GCSD_PhaT_avg = alpha1*GCSD_PhaT_avg + (1-alpha1)*GCSD_PhaT_act;
        
        % compensation of frequency dependent phase delay caused by recursive averaging with alpha1
        theta_recAvg1 = 2*pi/FFTsize*(0:(FFTsize-1))'*SROppm_act*1e-6*FFTshift;
        termCompRecAvg1 = exp(-1i*atan(alpha1*sin(theta_recAvg1)./(1-alpha1*cos(theta_recAvg1))));
        GCSD_PhaT_avgCmpRecAvg1 = GCSD_PhaT_avg.*termCompRecAvg1;
        
        % Calculation of the first generalized CCF (GCCF-1) in time domain
        if flag_calcFigOut
            GCCF1_avg_big = fftshift(real(ifft(GCSD_PhaT_avgCmpRecAvg1, FFTsize))');
            GCCF1(ell,:) = GCCF1_avg_big((FFT_Nyq-Upsilon):(FFT_Nyq+Upsilon));
        end
        
        % 3) Fill of DXCP-PhaT-container with Cont_NumFr number of past GCCF1_avg
        if ell<=Cont_NumFr
            DXCPPhaT_Cont_GCSD_avg(:,ell) = GCSD_PhaT_avgCmpRecAvg1;
        else
            DXCPPhaT_Cont_GCSD_avg(:,1:(Cont_NumFr-1)) = DXCPPhaT_Cont_GCSD_avg(:,2:Cont_NumFr); %move to left (overwrite oldest)
            DXCPPhaT_Cont_GCSD_avg(:,Cont_NumFr) = GCSD_PhaT_avgCmpRecAvg1; %insert new
        end
        
        % 4) As soon as DXCP-container is filled, compute a second CSD based on last and first
        % vectors GCCF1_avg from DXCP-container similar to (10) in [1] and do averaging
        if ell>=Cont_NumFr
            % Calculate the second GCSD
            GCSD2_act = DXCPPhaT_Cont_GCSD_avg(:,end).*conj(DXCPPhaT_Cont_GCSD_avg(:,1));
            GCSD2_hist(:,ell) = GCSD2_act;
            % Get GCCF-2 in time domain (as an additional output)
            if flag_calcFigOut
                GCCF2_act_big = fftshift(real(ifft(GCSD2_act, FFTsize))');
                GCCF2_act = GCCF2_act_big((FFT_Nyq-Lambda):(FFT_Nyq+Lambda));
                GCCF2(ell,:) = GCCF2_act;
            end
            %
            if flag_avgXcSec
                % recursive averaging over the second CSD or xcorr
                GCSD2_avg = alpha2*GCSD2_avg + (1-alpha2)*GCSD2_act;
            else
                % simple averaging over the whole signal
                GCSD2_avg = (ell_GCCF2avg-1)*GCSD2_avg/ell_GCCF2avg + GCSD2_act/ell_GCCF2avg;
                ell_GCCF2avg = ell_GCCF2avg + 1;
            end
            % Calculate averaged CCF-2 in time domain
            GCCF2_avg_big = fftshift(real(ifft(GCSD2_avg, FFTsize))');
            GCCF2_avg = GCCF2_avg_big((FFT_Nyq-Lambda):(FFT_Nyq+Lambda));
            %NK: Save averaged GCCF2 instead
            GCCF2(ell,:) = GCCF2_avg;
            % NK: Add: GCSD2_avg history
            GCSD2_avg_hist(:,ell) = GCSD2_avg;
        end
        
        % 5) Parabolic interpolation (13) with (14) with maximum search as in [1]
        % and calculation of the current SRO estimate acc. to (15) in [1]
        if ell>=Cont_NumFr+SetCSD2avg_NumFr     % As soon as useful SRO estimates are available
            % parabolic interpolation with upsampling
            GCCF2_avg_upsmp = resample(GCCF2_avg,p_upsmpFac,1);
            % figure; plot(-Lambda:Lambda, GCCF2_avg,'o'); hold on; plot(lambda_vec_upsmp, GCCF2_avg_upsmp,'sig_z'); grid on;
            [~, idx_max] = max(GCCF2_avg_upsmp);
            if (idx_max == 1) || (idx_max == (2*Lambda+1)*p_upsmpFac)
                DelATSest_ell_frac = 0;
            else
                sup_pnts = GCCF2_avg_upsmp((idx_max-1):(idx_max+1));  % supporting points y(x) for x={-1,0,1}
                % fractional maxima x_max=-b/2/a of y(x) = a*x^2 + b*x + c
                DelATSest_ell_frac = (sup_pnts(3)-sup_pnts(1))/2/(2*sup_pnts(2)-sup_pnts(3)-sup_pnts(1));
            end
            % resulting real-valued x_max of y(x)
            DelATSest_ell = lambda_vec_upsmp(idx_max) + DelATSest_ell_frac/p_upsmpFac;
            SROppm_act = DelATSest_ell/B*1e6;
            SROest_ppm(ell) = SROppm_act;   % in ppm
        end
        
        % 6) SSO estimation based on the shifted CCF-1
        if flag_STOest==1 && ell>=Cont_NumFr+SetCSD2avg_NumFr
            % a) phase shifting of GCSD-1 to remove SRO-induced time offset and phase of the recursive averaging
            timeOffset_SROinduced = SROppm_act*1e-6*FFTshift*(ell-1);
            termCompSROinducedShift = exp(1i*2*pi/FFTsize*timeOffset_SROinduced*(0:(FFTsize-1))');
            GCSD1_sroComp = GCSD_PhaT_avgCmpRecAvg1.*termCompSROinducedShift;
%             theta_recAvg1 = 2*pi/FFTsize*(0:(FFTsize-1))'*SROppm_act*1e-6*FFTshift;
%             termCompRecAvg1 = exp(-1i*atan(alpha1*sin(theta_recAvg1)./(1-alpha1*cos(theta_recAvg1))));
%             GCSD1_sroComp = GCSD_PhaT_avgCmpRecAvg1.*termCompSROinducedShift.*termCompRecAvg1;
            % b) averaging over time and zero-phase filtering within the frame (if necessary)
            if flag_avgCSDcomp
                GCSD1_sroComp_avg = alpha3*GCSD1_sroComp_avg + (1-alpha3)*GCSD1_sroComp;
            else
                GCSD1_sroComp_avg = (ell_avg_sh-1)*GCSD1_sroComp_avg/ell_avg_sh + GCSD1_sroComp/ell_avg_sh;
                ell_avg_sh = ell_avg_sh+1;
            end
            % c) go into the time domain via calculation of shifted GCC-1
            GCCF1_sroComp_big = fftshift(real(ifft(GCSD1_sroComp_avg, FFTsize))');
            GCCF1_sroComp_avg = GCCF1_sroComp_big((FFT_Nyq-Upsilon):(FFT_Nyq+Upsilon));
            % d) Maximum search over averaged filtered shifted GCC-1 (with real-valued SSO estimates)
            GCCF1_smShftAvg_forMax = GCCF1_sroComp_avg;
            % GCCF1_smShftAvg_forMax = abs(GCCF1_sroComp_avg); % removed because unusual for the following GCC-PhaT
            [~, idx_max] = max(GCCF1_smShftAvg_forMax);
            % D12_estGeoRev_act = Upsilon_vec(idx_max); % TO-DO reellwertig statt ganzzahling
            if (idx_max == 1) || (idx_max == 2*Upsilon+1)
                SSOsmp_est_act_frac = 0;
            else
                sup_pnts = GCCF1_smShftAvg_forMax((idx_max-1):(idx_max+1)); % supporting points y(x) for x={-1,0,1}
                % fractional maxima x_max=-b/2/a of y(x) = a*x^2 + b*x + c
                SSOsmp_est_act_frac = (sup_pnts(3)-sup_pnts(1))/2/(2*sup_pnts(2)-sup_pnts(3)-sup_pnts(1));
            end
            SSOsmp_est_act = idx_max - Upsilon - 1 + SSOsmp_est_act_frac;   % resulting real-valued x_max
            %
            GCCF1_sroComp(ell,:) = GCCF1_sroComp_avg;
            SSOest_smp(ell) = SSOsmp_est_act;
        end
        
        index_vec = index_vec + FFTshift;  % Update indices for next signal frame 
    end
end
