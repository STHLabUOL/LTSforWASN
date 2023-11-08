function [SROest_ppm, SROrem_ppm_vec, SROrem_ppm_raw_vec, ATDres_smp_vec, SSOest_smp, GCCF1, GCCF2, GCCF2_avg,...
    GCCF2_avg_upsmp, GCCF1_sroComp, GCCF1_sroComp_avg, sig_z2_res, GCCF2_avg_time, GCCF2_avg_upsmp_time, ...
    RTFs, IACS] = DXCP_PhaT_CL1_with_SAD(sig_z, parDXCPPhaT_cl, OnlineSAD)
% Uses OnlineSAD object to obtain pseudo-ACS (1 or 0) per frame.
% OnlineSAD object uses oracle information
%
% Online implementation of Double-Cross-Correlation Processor (DXCP)
% presented in the publication 'A Double-Cross-Correlation Processor (DXCP)
% for Blind Sampling Rate Offset Estimation in Acoustic Sensor Networks'
% submitted to the ICASSP-2019 [1].
% 
% :input sig_z:                     Two asynchronous signals [sig_lenx2]
% :param parDXCPPhaT_cl:            Parameters of DXCPrec approach (see below)
% :input OnlineIACS                 OnlineIACS Object for iacs calc
% 
% :return SROest_ppm:               SRO estimates in ppm (vector)
% :return SSOest_smp:               SSO estimates in samples (vector)
% :return GCCF1:                    First GCCF (matrix)
% :return GCCF2:                    Second GCCF (matrix)
% :return GCCF2_avg:                Final averaged GCCF-2 (vector)
% :return GCCF2_avg_upsmp:          Final upsampled GCCF2_avg (vector)
% :return GCCF1_sroComp:            SRO-compensated GCCF-2 (matrix)
% :return GCCF1_sroComp_avg:        SRO-comp. averaged GCCF-2 (vector
% :return SROrem_ppm_vec:                Estimates of ATD difference (vector)
% :return sig_z2_res:                    Resampled signal chanel 2 (vector)
    
    % external parameters
    WinSizeSINC = parDXCPPhaT_cl.WinSizeSINC;
    
    FFTsize = parDXCPPhaT_cl.FFTsize;   % FFT size (frame length)
    FFTshift = parDXCPPhaT_cl.FFTshift; % frame shift
    B_NumFr = parDXCPPhaT_cl.B_NumFr;   % number of frames acc. to desired B_smp
    Lambda = parDXCPPhaT_cl.Lambda;     % maximum lag of the second xcorr
    window = parDXCPPhaT_cl.window;     % analysis window
    
    alpha1 = parDXCPPhaT_cl.alpha1;     % smoothing constant of the first CCF
    alpha2 = parDXCPPhaT_cl.alpha2;     % smoothing constant of the second CCF
    alpha3 = parDXCPPhaT_cl.alpha3;     % smoothing constant of the SRO-compensated CSD
        
    SetCSD2avg_NumFr = parDXCPPhaT_cl.SetCSD2avg_NumFr;    % settling time of CSD2 averaging
    Z_12_abs_min = parDXCPPhaT_cl.Z_12_abs_min;            % to avoid devision by 0 in GCC-PhaT
    p_upsmpFac = parDXCPPhaT_cl.p_upsmpFac;                % upsampling factor
    
    flag_STOest = parDXCPPhaT_cl.flag_STOest;              % 1 - perform STO estimation, or 0 - not
    flag_calcFigOut = parDXCPPhaT_cl.flag_calcFigOut;      % 1 - calculate GCCFs for plots or 0 - mot
    flag_genWeight = parDXCPPhaT_cl.flag_genWeight;        % 1 - PhaT weighting or 0 - no weighting (conventional cross-sorrelation)
    flag_avgCSDcomp = parDXCPPhaT_cl.flag_avgCSDcomp;
    flag_avgXcSec = parDXCPPhaT_cl.flag_avgXcSec;
    
    fs_Hz = parDXCPPhaT_cl.fs_Hz;   % for calculation of RTFs
    
    ForwardControl = parDXCPPhaT_cl.ForwardControl;
    L_hold = parDXCPPhaT_cl.L_hold;
    K_Nomin = parDXCPPhaT_cl.K_Nomin;
    K_Denom = parDXCPPhaT_cl.K_Denom;
    
    % internal parameters
    B_smp = B_NumFr*FFTshift;
    Upsilon = FFTsize/2-1;                                  % maximum lag of the first xcorr
    Cont_NumFr = B_NumFr+1;                                 % number of frames in DXCP container
    DXCPPhaT_Cont_GCSD_avg = zeros(FFTsize, Cont_NumFr);    % DXCP container
    
     % further useful parameters
    [len_sig, num_chan] = size(sig_z);
    FrNum_L = floor((len_sig-FFTsize)/FFTshift)+1;          % number of signal frames for frame shift
    FFT_Nyq = FFTsize/2+1;
    lambda_vec_upsmp = -Lambda:1/p_upsmpFac:(Lambda+1-1/p_upsmpFac);
    
    % output variables
    SROest_ppm = zeros(1,FrNum_L);
    SROrem_ppm_vec = zeros(1,FrNum_L);
    SROrem_ppm_raw_vec = zeros(1,FrNum_L); 
    IACS = zeros(1,FrNum_L);
    ATDres_smp_vec = zeros(1,FrNum_L);
    SSOest_smp = zeros(1,FrNum_L);
    GCCF1 = zeros(FrNum_L, 2*Upsilon+1);
    GCCF2 = zeros(FrNum_L, 2*Lambda+1);
    GCCF1_sroComp = zeros(FrNum_L, 2*Upsilon+1);
    GCCF1_sroComp_avg = zeros(1, 2*Upsilon+1);
    sig_z2_res = zeros(size(sig_z,1),1);
    GCCF2_avg_time=zeros(FrNum_L, 2*Lambda+1);
    GCCF2_avg_upsmp_time=zeros(FrNum_L,p_upsmpFac*(2*Lambda+1));
    RTF_vec = zeros(1,2*FrNum_L);
    RTFs = zeros(1,3);

    % some variables for averaging
    SROppm_act = 0;
    GCSD1_avg = zeros(FFTsize, 1);
    GCSD2_avg = zeros(FFTsize, 1);
    GCSD1_sroComp_avg = zeros(FFTsize, 1);
    
    % frame-based processing of the input signal
    index_vec = 1:FFTsize;
    index_vec_sinc = 1:(FFTsize+WinSizeSINC);
%     idx_vec_sinc = (-WinSizeSINC+1):(FFTsize+WinSizeSINC);

    ell = 1;
    ell_GCCF2avg = 1;
    ell_GCSD1avg = 1;
    idx_ring_start = 1;
    
    del_sinc_int = 0;
%     del_sinc_frac = 0;
%     SSOsmp_est_act =0;
    SROppm_act_control_vec = zeros(1,50);
    SROppm_rem_control_vec = zeros(1,50);
    flag_forwContr_done = 0;
    SROppm_init=0;
    ATDres_smp=0;
%     d12_corr=0;
%     t1=0;t2=0;t3=0;t4=0;t5=0;t0=0;t31=0;FFTs=0;
    %     for ell=1:FrNum_L
    tic;
    while max(index_vec)<len_sig && max(index_vec_sinc-del_sinc_int)<len_sig
        %tic
        
        %% Synchronization (va moving buffer and resampling)
        % Sync-1: calculation of the current SRO and an accumulating time drift (ATD)
        SROppm_forCorr = SROppm_act;
        SROfact_corr = -1/(1+SROppm_forCorr*1e-6);
        SROppm_corr = SROfact_corr*SROppm_forCorr;
        % efficient arbitrary resampling of the 2. channel
        if ell==1
            RESsize = FFTsize;  % of the whole signal frame
        else
            RESsize = FFTshift; % only of the last FFTshift frame samples
        end
        idx_vec_res = (index_vec(FFTsize-RESsize+1)-WinSizeSINC):(index_vec(FFTsize)+WinSizeSINC);
        % calculate control values for resampling (w.o. STO or TDOA compensation)
        if ell==1 % whole frame in the beginning
            ATDres_smp = SROppm_corr*1e-6*(FFTsize-1)/2;
        elseif ell==2 % efficient resampling
            ATDres_smp = ATDres_smp + SROppm_corr*1e-6*(FFTsize+FFTshift)/2;
        else % further frames
            ATDres_smp = ATDres_smp + SROppm_corr*1e-6*FFTshift;
        end
        ATDres_smp_vec(ell) = ATDres_smp;
        %ATDres_smp = ATDres_smp + SROppm_corr*1e-6*FFTshift-d12_corr+d12_corr_old;
        ATDres_int = round(ATDres_smp);
        ATDres_frac = ATDres_smp - ATDres_int;
        % indexes for a buffer movable on async-signal 
        idx_vec_res_z2 = idx_vec_res+ATDres_int;
        % avoid overflow at the end of asynchronous signal
        if max(idx_vec_res_z2)>len_sig
            idx_vec_res_z2 = (len_sig-length(idx_vec_res)+1):len_sig;
        end
        
        % Sync-2: Moving buffer (compensation for integer part of accumulating time drift)
        if ell==1
            % extend with zeros in the beginning of asynchronous signal
            z_2_ell_resIn = [zeros(WinSizeSINC,1); sig_z(idx_vec_res_z2((WinSizeSINC+1):end), 2)];
        else
            z_2_ell_resIn = sig_z(idx_vec_res_z2,2);
        end
        
        % Sync-3: resampling (compensation for SRO and a fractional part of ATD)
        z_2_ell_resOut = multirate.resampler_sinc_on(z_2_ell_resIn', WinSizeSINC, SROppm_corr, -ATDres_frac);
        %z_2_ell_resOut = asn.freqresample(z_2_ell_resIn', 2^10, 2^9, SROppm_corr)';
        
        % store a resampled signal and the current sigal segments for DXCP-PhaT
        sig_z2_res(index_vec((FFTsize-RESsize+1):FFTsize)) = z_2_ell_resOut((1:RESsize)+WinSizeSINC);
        z_1 = sig_z(index_vec,1); z_2_res = sig_z2_res(index_vec);
        
        %t0=mean([t0,toc]);
        RTF_vec(2*(ell-1)+1) = toc/(FFTshift/fs_Hz);
        %tic
        
        % 1) Cutting current frames and window them acc. to eq. (5) in [1]
        z_12_win = [z_1 z_2_res].*(window*ones(1,num_chan));
        
        % 2) Calculate generalized (normaized) GCSD with Phase Transform (GCSD-PhaT)
        % and apply recursive averaging
        Z_12 = fft(z_12_win, FFTsize, 1);
        %FFTs=FFTs+2;
        Z_12_act = Z_12(:,1).*conj(Z_12(:,2));
        %t1 = mean([t1,toc]);
        %tic
        
        if flag_genWeight
            % PhaT weighting
            Z_12_act_abs = abs(Z_12_act);
            Z_12_act_abs(Z_12_act_abs<Z_12_abs_min) = Z_12_abs_min; % avoid division by 0
            GenWeights = 1./Z_12_act_abs;
        else % conventional cross-correlation
            GenWeights = 1;
        end
        GCSD1_avg = alpha1*GCSD1_avg + (1-alpha1)*GenWeights.*Z_12_act;
        
        if flag_STOest==1
            % compensation of frequency dependent phase delay caused by recursive averaging with alpha1
            theta_recAvg1 = 2*pi/FFTsize*(0:(FFTsize-1))'*SROppm_act*1e-6*FFTshift;
            termCompRecAvg1 = exp(-1i*atan(alpha1*sin(theta_recAvg1)./(1-alpha1*cos(theta_recAvg1))));
            GCSD_PhaT_avgCmpRecAvg1 = GCSD1_avg.*termCompRecAvg1;
        end
        
        % Calculation of the first generalized CCF (GCCF-1) in time domain
        if flag_calcFigOut
            GCCF1_avg_big = fftshift(real(ifft(GCSD_PhaT_avgCmpRecAvg1, FFTsize))');
            GCCF1(ell,:) = GCCF1_avg_big((FFT_Nyq-Upsilon):(FFT_Nyq+Upsilon));
        end
        %t2 = mean([t2,toc]);
        
        % 3) Filling of DXCP-container with Cont_NumFr number of past XcorrFirst_sm
        DXCPPhaT_Cont_GCSD_avg(:,idx_ring_start) = GCSD1_avg;
        if(idx_ring_start == Cont_NumFr)
            idx_ring_end = 1;
        else
            idx_ring_end = idx_ring_start+1;
        end
        
        % 4) As soon as DXCP-container is filled, compute second CSD based
        % on last and first vectors of DXCP-container and perform time averaging
        if ell>=Cont_NumFr
            %tic
            % Calculate the second GCSD
%             GCSD2_act = DXCPPhaT_Cont_GCSD_avg(:,end).*conj(DXCPPhaT_Cont_GCSD_avg(:,1));
            GCSD2_act = DXCPPhaT_Cont_GCSD_avg(:,idx_ring_start).*conj(DXCPPhaT_Cont_GCSD_avg(:,idx_ring_end));
            % Get GCCF-2 in time domain (as an additional output)
            if flag_calcFigOut
                GCCF2_act_big = fftshift(real(ifft(GCSD2_act, FFTsize))');
                GCCF2_act = GCCF2_act_big((FFT_Nyq-Lambda):(FFT_Nyq+Lambda));
                GCCF2(ell,:) = GCCF2_act;
            end
            % Averaging of the second GCSD over time
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
            GCCF2_avg_time(ell,:) = GCCF2_avg;
           % t4 = mean([t4,toc]);
        end
        
        % Compute current IACS
        OnlineSAD.process(ell);
        IACS(ell) = OnlineSAD.SAD;

        %Ringbuffer variable handling
        idx_ring_start=idx_ring_start+1;
        if mod(ell,Cont_NumFr)==0
            idx_ring_start=1;
        end
        
        % 5) Polynomial interpolation (13) with (14) with maximum search as in
        % [1] and calculation of the remaining current SRO estimate sim. to (15) in [1]
        if ell>=Cont_NumFr+SetCSD2avg_NumFr % As soon as useful SRO estimates available
            %tic
            % parabolic interpolation with upsampling
            GCCF2_avg_upsmp = resample(GCCF2_avg,p_upsmpFac,1);
            GCCF2_avg_upsmp_time(ell,:)=GCCF2_avg_upsmp;
            [~, idx_max] = max(GCCF2_avg_upsmp);
            if (idx_max == 1) || (idx_max == (2*Lambda+1)*p_upsmpFac)
                DelATSest_ell_frac = 0;
            else
                sup_pnts = GCCF2_avg_upsmp((idx_max-1):(idx_max+1));  % supporting points y(x) for x={-1,0,1}
                % fractional maxima x_max=-b/2/a of y(x) = a*x^2 + b*x + c
                DelATSest_ell_frac = (sup_pnts(3)-sup_pnts(1))/2/(2*sup_pnts(2)-sup_pnts(3)-sup_pnts(1));
            end
            % resulting real-valued x_max of y(x)
            DelATSest_act = lambda_vec_upsmp(idx_max) + DelATSest_ell_frac/p_upsmpFac;
            SROppm_rem_act_raw = DelATSest_act/B_smp*1e6;    
            SROrem_ppm_raw_vec(ell) = SROppm_rem_act_raw;
            % --- apply IACS ---                            
            if ell > 1
                %SROppm_rem_act = IACS(ell)*SROppm_rem_act_raw + (1-IACS(ell))*SROrem_ppm_vec(ell-1);
                SROppm_rem_act = IACS(ell)*SROppm_rem_act_raw; %scale, set to 0 when IACS=0
            else
                SROppm_rem_act = SROppm_rem_act_raw;
            end
            
            SROrem_ppm_vec(ell) = SROppm_rem_act;
            %t5 = mean([t5,toc]);
        end
        
        % IMC-controller
        if ell>=Cont_NumFr+SetCSD2avg_NumFr % As soon as useful SRO estimates available
            % save SRO estimates for step detection and for use in IMC control
            SROppm_rem_vec = [SROppm_rem_act, SROppm_rem_control_vec(1:end-1)];
            SROppm_rem_control_vec = [SROppm_rem_act, SROppm_rem_control_vec(1:end-1)];
            %feedforward-control with step detection  
            if (ForwardControl && abs(SROppm_rem_vec(1)-SROppm_rem_vec(2)) > 10^6/(B_NumFr*FFTshift) && flag_forwContr_done==0)
                %SROppm_init = SROppm_rem_act + SROppm_init;   % update rough SRO estimation
                SROppm_init = SROppm_rem_act + SROppm_act;
                ell_init = ell;                         % new Step frame deonoted 
                flag_forwContr_done = 1;                 % SRO Step detected use forward control
            end
            %wait for stabilization of second CCF than feedback-control
            if ( flag_forwContr_done==1 && ell<ell_init+L_hold )
                SROppm_act_control_vec = zeros(1,50);   % Reset IMC controller
                SROppm_rem_control_vec = zeros(1,50);   % Reset IMC controller
                SROppm_act = SROppm_init;               % hold rough SRO for 12.8 seconds
            else
                flag_forwContr_done=0;                   % Reset feed forward control and step detection
                %feedback-control
                SROppm_act_control = K_Nomin*SROppm_rem_control_vec(1:length(K_Nomin))' - K_Denom(2:length(K_Denom))*SROppm_act_control_vec(1:length(K_Denom)-1)';
                SROppm_act_control_vec = [SROppm_act_control, SROppm_act_control_vec(1:end-1)]; % save IMC output for use in IMC control
                SROppm_act = SROppm_act_control+SROppm_init;    % update controlled SRO
            end
            SROest_ppm(ell) = SROppm_act;               % Save SRO in array for output
        end
        
        % 6) SSO estimation based on the shifted CCF-1
        if flag_STOest==1 && ell>=Cont_NumFr+SetCSD2avg_NumFr
            % a) phase shifting of GCSD-1 to remove SRO-induced time offset and phase of the recursive averaging
            timeOffset_SROinduced = SROppm_act*1e-6*FFTshift*(ell-1);
            termCompSROinducedShift = exp(1i*2*pi/FFTsize*timeOffset_SROinduced*(0:(FFTsize-1))');
            GCSD1_sroComp = GCSD_PhaT_avgCmpRecAvg1.*termCompSROinducedShift;
            
            % b) averaging over time and zero-phase filtering within the frame (if necessary)
            if flag_avgCSDcomp
                GCSD1_sroComp_avg = alpha3*GCSD1_sroComp_avg + (1-alpha3)*GCSD1_sroComp;
            else
                GCSD1_sroComp_avg = (ell_GCSD1avg-1)*GCSD1_sroComp_avg/ell_GCSD1avg + GCSD1_sroComp/ell_GCSD1avg;
                ell_GCSD1avg = ell_GCSD1avg + 1;
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
        RTF_vec(2*(ell-1)+2) = toc/(FFTshift/fs_Hz);

        index_vec_sinc = index_vec_sinc + FFTshift;
        index_vec = index_vec + FFTshift;  % Update indices for next signal frame
        ell = ell + 1;
        
    end
    RTFdiff_vec = diff(RTF_vec);
    RTFdiff_vec_ell = diff(RTF_vec(1:2:end));
    RTFs(1) = mean(RTFdiff_vec(2*(Cont_NumFr+floor(SetCSD2avg_NumFr))+1:2:end));   % RTF of resampling method (1)
    RTFs(2) = mean(RTFdiff_vec(2*(Cont_NumFr+floor(SetCSD2avg_NumFr))+2:2:end));   % RTF of sro estimator (2)
    RTFs(3) = mean(RTFdiff_vec_ell(Cont_NumFr+floor(SetCSD2avg_NumFr)+1:end));     % RTF of one frame of (1) and (2)
end