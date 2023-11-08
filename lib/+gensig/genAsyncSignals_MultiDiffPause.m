function [signals, setup] = genAsyncSignals_MultiDiffPause(parSource, parRoom, parNoise, parDiffuse, parAsync)
% Generate asynchronous noisy signals z_p for a single/multi-source room-
% acoustic setup with P receivers (microphones) and Ns acoustic sources.
% Mic-signals are distorted by both diffuse noise and self-noise:
% 
% z_p = DDC(y_p) & y_p = sum_{i=1}^Ns(s_i*h_{i,p}) + w_p + v_p
% 
% where DDC - digital-to-digital converter, y_p - noisy synchronous signal
% of the p-th mic, s_i - an acoustic source signal of i-th source,
% h_{i,p} - a room-impuls response between i-th source and p-th mic,
% Ns - a number of acoustic sources, w_p - an addivite diffuse noise on
% p-th mic, v_p - a self-noise of the p-th mic.
%
% INPUT:
%
%   parSource               (struct) parameters for source signal
% 
%       .fs_Hz              (double) reference sampling frequency in Hz
% 
%       .len_s              (double) signal length in seconds, if the 
%                           signal is too long it will be cut, if it is too
%                           short it will be repeated to match the desired 
%                           length
%
%       .file               (1 x Ns cell) specified (single-channel) .wav file, 
%                           not required if .preset is set
%
%       .preset             (1 x Ns cell) preset source signals:
%                               'speech' male and female e.g. from TIMIT
%                               'white' Gaussian white noise
%                           setting .preset will overwrite .file
%
%   parRoom                 (struct) parameters for room size and acoustics
%
%       .roomSize_m         (3 x 1 double) z/y/z dim [m] of the room, not
%                           required if .preset is set
%
%       .srcPos_m           (3 x 1 double) z/y/z pos [m] of the source, not
%                           required if .preset is set
%
%       .micPos_m           (3 x P double) z/y/z pos [m] of P mics, not
%                           required if .preset is set
%   
%       .preset             (int) choose one of four (1..4) presets (P = 2
%                           overwrites .roomSize_m, .srcPos_m, .micPos_m)
%
%       .T60_s              (double) reverberation time (T60_s) in seconds
%   
%       .RIRlen_smp         (int) length of generated RIR (samples)
%
%       .plot               (int) toggle to plot room configuration:
%                               0: do not plot (default)
%                               1: plot room setup
%
%   parNoise                (struct) parameters additive independent noise
%
%       .file               (string) specified (P-channel) .wav file, 
%                           not required if .preset is set
%
%       .preset             (string) preset source signals:
%                               'babble' babble from Noisex
%                               'factory1' from Noisex
%                               'factory2' from Noisex
%                               'white' Gaussian white noise
%                           setting .preset will overwrite .file
%
%       .SNR_dB             (double) global SNR level [dB], global noise 
%                           power is kept fixed and SNR is linearly 
%                           averaged across all mics
% 
%   parDiffuse             (struct) parameters additive diffuse noise
%
%       .file               (string) specified (P-channel) .wav file, 
%                           not required if .preset is set
% 
%       .preset             (string) preset source signals:
%                               e.g. 'babble' babble from Noisex
% 
%       .SNR_dB             (double) global SNR level of diffuse noise [dB]
% 
%       .spatialType        (string) type of diffuse noise
%                               either 'spherical' or 'cylindrical'
% 
%       .spatialFFT         (int) FFT-length to generate diffuse matrix
%   
%   parAsync                (struct) parameters for simulated SRO
%       
%       .SRO_ppm            (P x 1 double) SROs measured in parts per
%                           million (ppm) for each channel relative to
%                           reference sampling rate fs_Hz (channels with
%                           value zero are not resampled)
%
%       .winLenSINC         (int) length of sinc window length
%
%       .SSO_interval_smp   (1 x 2) interval for drawing real-valued values
%                           of simulated sampling start offset (SSO). Note,
%                           (channels with SRO_ppm = 0 get SSO_smp = 0).
%
% OUTPUT:
% 
%   signals                   (struct) generated signals
%
%       .s                  (len_smp x Ns) source signal
%
%       .s_h                (len_smp x P) noise-free reverberant mic signals:
%                           s_h_p = sum_{i=1}^{num_src}(s_i * h_i_p)
%   
%       .w                  (len_smp x P) diffuse-noise signals
%   
%       .v                  (len_smp x P) self-noise signals
%
%       .y                  (len_smp x P) noisy reverberant synchronous
%                           mic signals: x_p_sync = s_h_p + v
% 
%       .z                  (len_smp x P) noisy reverberant asynchronous
%                           mic signals: z_p = DDC(y_p)
%
%   setup                   (struct) information about simulation setup
%   
%       .fs_Hz              (double) original source fs in Hz
%
%       .roomSize_m         (3 x 1 double) z/y/z dim [m] of the room
%
%       .srcPos_m           (3 x 1 double) z/y/z pos [m] of the source
%
%       .micPos_m           (3 x P double) z/y/z pos [m] of P mics
% 
%       .T60_s              (double) T60 time in seconds
%   
%       .rir_h              (Ns x P x RIRlen_smp) room impulse responses
%
%       .SRO_ppm            (1 x P) simulated sampling rate offset in ppm
%                           (fixed and defined by user)
%
%       .SSO_smp            (1 x P) simulated sampling start offset in
%                           samples (random and real-valued)
% 
%       .TDOAmatrix_smp     (P x P x Ns) matrix with time difference of
%                           arrival (TDOA) values [smp]
% 
%       .SNR_dB_eff         (double) mean effective SNR [dB]
%
% written by Philipp Thuene (June 12, 2018)
% - modified by Aleksej Chinaev (February 5, 2020)
% (adding noise signals, resampling for input signal length, diffuse noise)

    %% PARAMTERS: BASIC
    
    fs_Hz = parSource.fs_Hz;
    len_smp = parSource.len_s * fs_Hz;
    smp_add = parAsync.winLenSINC + ceil(max(parAsync.SRO_ppm)*1e-6*(len_smp-1));
    len_smp_add = len_smp + smp_add;


    %% COMPOSE CLEAN SIGNAL AND MASK
    
    Ns = length(parSource.file); % number of acoustic sources
    s = zeros(len_smp_add, Ns);
    s_mask = zeros(len_smp_add,Ns);
    %len_coh_smp = 60*fs_Hz;
    %src_loudness = parSource.loudness;
    for idx_src=1:Ns   
        % load clean sig
        s_i = audioread(parSource.file{idx_src});
        % set mask based on activity range and loudness
        activity = parSource.activity{idx_src};
        for ii = 1:size(activity, 1)
            from = floor(activity(ii, 1)*fs_Hz)+1;
            to = floor(activity(ii, 2)*fs_Hz);
            s_mask(from:to, idx_src) = parSource.loudness(idx_src);
        end
        % guarantee sufficient sig length 
        len_read = length(s_i);
        if len_read < len_smp_add
            s_i = repmat(s_i,ceil(len_smp_add/len_read),1); % take care of it
        end
        s(:,idx_src) = s_i(1:len_smp_add).*s_mask(:,idx_src);
    end
    

    %% GENREATE RIRs FOR P CHANNELS

    % shorthands
    roomSize_m = parRoom.roomSize_m;
    micPos_m = parRoom.micPos_m;
    srcPos_m = parRoom.srcPos_m;
    
    % calculate TDOA_smp(s): implemented here only for 2 mics !!!
    soundSpeed_c_mps = 340; % in m/s (340 for 15°C)
    % TDOA_smp = zeros(Ns,1);
    P = length(micPos_m(1,:));
    TDOAmatrix_smp = zeros(P,P,Ns);
    for idx_src=1:Ns
%         distMic1_m = sqrt(sum((micPos_m(:,1)-srcPos_m(:,idx_src)).^2));
%         distMic2_m = sqrt(sum((micPos_m(:,2)-srcPos_m(:,idx_src)).^2));
%         TDOA_smp(idx_src,1) = (distMic1_m-distMic2_m)*fs_Hz/soundSpeed_c_mps;
        for p=2:P
            for q=1:(p-1)
                distMic_p_m = sqrt(sum((micPos_m(:,p)-srcPos_m(:,idx_src)).^2));
                distMic_q_m = sqrt(sum((micPos_m(:,q)-srcPos_m(:,idx_src)).^2));
                TDOAmatrix_smp(p,q,idx_src) = (distMic_q_m-distMic_p_m)*fs_Hz/soundSpeed_c_mps;
            end
        end
    end

    % generate RIR h
    len_h = parRoom.RIRlen_smp;
    T60_s = parRoom.T60_s;
    h = zeros(Ns, P, len_h);
    H_cell = gensig.createRIR_funcAC(roomSize_m, srcPos_m, micPos_m, T60_s, soundSpeed_c_mps, len_h, fs_Hz);
    for idx_src=1:Ns
        for p=1:P
            h(idx_src,p,:) = H_cell{idx_src,p}(1:len_h);
        end
    end
%     h = h/max(max(max(abs(h)))); % normalization to max amplitude of h

    
    %% GENERATE REVERBERANT SIGNALS VIA RIRs
    
    s_h = zeros(len_smp_add, P);
    for p=1:P
        for idx_src=1:Ns
            s_h_p_src = filter(squeeze(h(idx_src, p, :)), 1, s(:, idx_src));
            s_h(:,p) = s_h(:,p) + s_h_p_src(1:len_smp_add);
        end
    end

    
    %% COMPOSE SELF NOISE SIGNALS

    n_in = audioread(parNoise.file);
    n_in = n_in(:,1);
    v = zeros(len_smp_add,P);

    % construct noise signals
    startPoints_n  = floor(length(n_in)/P)*(0:(P-1))+1;
    for idx_chan=1:P
        v(:,idx_chan) = n_in(startPoints_n(idx_chan)+(0:(len_smp_add-1)));
        v(:,idx_chan) = v(:,idx_chan) - mean(v(:,idx_chan));
    end
    

    %% GENERATE DIFF. NOISE SIGNALS
        
    spatialType = parDiffuse.spatialType;      % spherically diffuse 
    spatialFFT = parDiffuse.spatialFFT;
    
    % Generate diffuse noise excitation
    w_in = audioread(parDiffuse.file);  % read excitation diffuse signal
    len_read_w_in = length(w_in);
    
    w_in_channels = zeros(len_smp_add,P);
    if len_read_w_in < P*len_smp_add
        w_in = repmat(w_in,ceil(P*len_smp_add/len_read_w_in),1); % take care of it
    end
    w_in = w_in(1:P*len_smp_add);
    w_in = w_in - mean(w_in);
    w_in_channels(:) = w_in;

    DistMat = zeros(P);
    for p=1:P
        for q=1:P
            DistMat(p,q) = norm(micPos_m(:,p) - micPos_m(:,q));
        end
    end
    % define coherence (arbitrary geometry)
    Coh = gensig.defineCoherenceArbitr_func(spatialType, P, DistMat, soundSpeed_c_mps, fs_Hz, spatialFFT);
    % mix signals (toolbox call, modified)
    w = gensig.mix_signals(w_in_channels, Coh, 'eigen');     % never really tried 'cholesky'
    
    
    %% SCALE SELF-NOISE AND DIFF. NOISE (ACCORDING TO SNR) AND ADD TO REVERBERANT SIGNALS
    
    SNR_lin_selfNoise = 10^(parNoise.SNR_dB/10);
    SNR_lin_diffuseNoise = 10^(parDiffuse.SNR_dB/10);
    pow_s_h = zeros(1,P);
    for p=1:P
        pow_s_h(p) = var(s_h(:,p));
        v(:,p) = v(:,p)/sqrt(var(v(:,p)));  % self noise with power 1
        w(:,p) = w(:,p)/sqrt(var(w(:,p)));  % diffuse noise with power 1
%         % scale self noise
%         fact_selfNoise_p = sqrt(pow_s_h(p)/SNR_lin_selfNoise);
%         v(:,p) = fact_selfNoise_p * v(:,p);
%         % scale diffuse noise
%         fact_diffuseNoise_p = sqrt(pow_s_h(p)/SNR_lin_diffuseNoise);
%         w(:,p) = fact_diffuseNoise_p * w(:,p);
    end
    % scale self noise
    fact_selfNoise = sqrt(mean(pow_s_h)/SNR_lin_selfNoise);
    v = fact_selfNoise * v;
    % scale diffuse noise
    fact_diffuseNoise = sqrt(mean(pow_s_h)/SNR_lin_diffuseNoise);
    w = fact_diffuseNoise * w;
    
    % calculate effective SNR w.r.t. diffuse and mic-self noises
    v_plus_w = v + w;
    pow_v_plus_w = zeros(1,P);
    for p=1:P
        pow_v_plus_w(p) = var(v_plus_w(:,p));
    end
    SNR_dB_eff = 10*log10(mean(pow_s_h)/mean(pow_v_plus_w));
    
    % mic signal with mic-self noise (still synchronous for later evaluation)
    % synchronous reverberant signal with the additional noises (diffuse, self-mic and transient)
    y = s_h + w + v;
    
    
    %% GENERATE ASYNC. SIGNALS VIA DIGITAL-TO-DITIGAL CONVERSION (DDC)
    
    z = y;  % only channels with SRO_ppm \not= 0 are resampled
    SSO_smp = zeros(1,P);
    for p=2:P   % first channel is the reference channel with SRO_ppm = SSO_smp = 0
        % draw SSO value from uniform distribution
        SSO_smp(p) = parAsync.SSO_interval_smp(1) + diff(parAsync.SSO_interval_smp)*rand(1);
        SSO_smp_int = round(SSO_smp(p));
        SSO_smp_frac = SSO_smp(p) - SSO_smp_int;
        if SSO_smp_int>=0
            y_p_sso = [y((SSO_smp_int+1):len_smp_add,p); zeros(SSO_smp_int,1)];
        else
            y_p_sso = [zeros(abs(SSO_smp_int),1); y(1:(len_smp_add-abs(SSO_smp_int)),p)];
        end
        z(:,p) = (multirate.resampler_sinc(y_p_sso.', parAsync.winLenSINC, parAsync.SRO_ppm(p), -SSO_smp_frac)).';
    end
    
    
    %% CUT TO DESIRED SIGNAL LENGTH

    signals.s = s(1:len_smp,:);
    signals.h = h;
    signals.s_h = s_h(1:len_smp,:);
    signals.s_mask = s_mask(1:len_smp,:);
    signals.w = w(1:len_smp,:);
    signals.v = v(1:len_smp,:);
    signals.y = y(1:len_smp,:);
    signals.z = z(1:len_smp,:);
    

    %% SET OUTPUTS

    setup.fs_Hz = fs_Hz;
    setup.roomSize_m = roomSize_m;
    setup.srcPos_m = srcPos_m;
    setup.micPos_m = micPos_m;
    setup.T60_s = T60_s;
    setup.rir_h = h;
    setup.SRO_ppm = parAsync.SRO_ppm;
    setup.SSO_smp = SSO_smp;
    setup.SNR_dB_eff = SNR_dB_eff;
%     setup.TDOA_smp = TDOA_smp;
    setup.TDOAmatrix_smp = TDOAmatrix_smp;
    
end