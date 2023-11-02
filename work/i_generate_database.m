% Generate async microphone signals in simulated environment
%
%


%% INIT

clear variables; %close all; %clc;
cl_srt = clock;
conf = Config();


%% PARAMETERS: MAIN

fs_Hz = 16000;      % reference sampling rate [Hz]
sigLen_s = 300;      % signal length [seconds] (previously 90)
T60_ms = 200;       % reverberation time T_{60} in ms
SRO_ppm_range = [0 100];     % desired sampling rate offset [ppm]
STO_smp = 0;     % desired start time offset [smp]
SNR_dB_selfNoise = 20;  % signal-to-noise ratio [dB] for mic. noise
SNR_dB_diffuseNoise_range = [10 20]; 

% randomiser params
n_setups = 30; % number of random setups

globSrcPause_start_range = [90, 110];
globSrcPause_end_range = [190, 210];
locSrcA_rel_start_range = [0.1, 0.2]; % range relative start point
locSrcA_rel_end_range = [0.6, 0.7]; 
locSrcB_rel_start_range = [0.2, 0.3]; % range relative start point
locSrcB_rel_end_range = [0.8, 0.9]; 

minDistFromWall = .5;
minDistFromGlobSrc = 3;
distBetweenLocSrc= linspace(1, 8, n_setups); % uniform
micDistFromLocSrc_range = [.1, .3];
FFTsize = 2^13;
FFTshift = 2^11;

sigPath = conf.pathTIMIT;

%% PARAMETERS: SOURCE(S)

% parameters of the source signal
parSource.fs_Hz = fs_Hz;            % reference sampling frequency
parSource.len_s = sigLen_s;         % signal duration in seconds
parSource.numSrcs = 3;              % 3. fixed. number of acoustic sources
preset_source = 'speech';          % 'speech'/'white'
parSource.preset(1) = {preset_source};
parSource.preset(2) = {preset_source};
parSource.preset(3) = {preset_source}; % files drawn randomly 
parSource.loudness_glob = 5;
parSource.loudness = [
    parSource.loudness_glob 
    .5 
    .5
];


%% PARAMETERS: ROOM
% (including positions of sources and mics)
% to be randomised: srcPos_m, micPos_m

parRoom.plot = 0;
dist = @(vec1, vec2) sqrt(sum((vec1-vec2).^2, 2));
% parRoom.soundSpeed_c_mps = 340; % in meter per second (340 m/s for 15°C)
% general geometrical parameters
parRoom.roomSize_m = [6, 8, 3].'; % 50 m^2 (SINC apartment)
parRoom.T60_ms = T60_ms;           % T60 time [ms] % 500
parRoom.T60_s = T60_ms*1e-3;
parRoom.RIRlen_smp = ceil(T60_ms*1e-3*fs_Hz/2);
parRoom.dist_wall_m = minDistFromWall;
parRoom.Pos_Reg = [minDistFromWall, parRoom.roomSize_m(1)-minDistFromWall;
                    minDistFromWall, parRoom.roomSize_m(2)-minDistFromWall;
                    minDistFromWall, parRoom.roomSize_m(3)-minDistFromWall;
                   ]; %positions region


%% PARAMETSR: NOISE

% self-noise
parNoise.SNR_dB = SNR_dB_selfNoise;
parNoise.preset = 'white'; % 'white'/'pink'
parNoise.file = [sigPath '../simulated_noise_white_16kHz_16min.wav'];

% diffuse noise
parDiffuse.preset = 'white';
parDiffuse.file = [sigPath '../simulated_noise_white_16kHz_16min_2.wav'];
parDiffuse.spatialType = 'spherical';      % 'white'/'spherical'/'cylindrical'
parDiffuse.spatialFFT = 1024;


%% PARAMETERS: ASYNC

parAsync.SSO_absMax = STO_smp;  % maximum absolute value of sampling start offset (SSO)
parAsync.SSO_interval_smp = parAsync.SSO_absMax*[1 1]; % interval SSO
parAsync.winLenSINC = 256;      % one-sided window length of sinc resampler


%% LOOP: GENERATE RANDOM SETUPS AND CORRESP. ASYNC SIGNALS

% Prepare target folder
[year, month, day] = ymd(datetime);
[hour, minute, second] = hms(datetime);
target_dirname = [conf.pathDB_root, num2str(year) '_' num2str(month) '_' num2str(day) '_' num2str(hour) '_' num2str(minute) '_' num2str(round(second))];
if isfolder(target_dirname)
    error('Target directory for new results already exists!')
else
    mkdir(target_dirname);
    mkdir([target_dirname '/mat_files']);
    mkdir([target_dirname '/wav_files']);
end

% Export workspace (save all variables)
save([target_dirname '/parameters.mat']);

% Generate signals...
for nn = 1:n_setups
            
    % draw diff. noise SNR
    parDiffuse.SNR_dB = unifrnd(SNR_dB_diffuseNoise_range(1), SNR_dB_diffuseNoise_range(2));

    % draw SRO
    SRO_ppm = unifrnd(SRO_ppm_range(1), SRO_ppm_range(2));
    parAsync.SRO_ppm = [0 SRO_ppm]; % first mic is reference

    % draw source activity
    globSrcPause_start = unifrnd(globSrcPause_start_range(1), globSrcPause_start_range(2));
    globSrcPause_end = unifrnd(globSrcPause_end_range(1), globSrcPause_end_range(2));
    parSource.activity{1} = [0, globSrcPause_start; globSrcPause_end, sigLen_s];
    globSrcPause_len = globSrcPause_end - globSrcPause_start;
    parSource.activity{2} = [globSrcPause_start+globSrcPause_len*unifrnd(locSrcA_rel_start_range(1), locSrcA_rel_start_range(2)), ...
        globSrcPause_start+globSrcPause_len*unifrnd(locSrcA_rel_end_range(1), locSrcA_rel_end_range(2))];
    parSource.activity{3} = [globSrcPause_start+globSrcPause_len*unifrnd(locSrcB_rel_start_range(1), locSrcB_rel_start_range(2)), ...
        globSrcPause_start+globSrcPause_len*unifrnd(locSrcB_rel_end_range(1), locSrcB_rel_end_range(2))];

    % draw source positions
    disp('Drawing source positions...'); tic;
    randomRoomPos = @() [unifrnd(minDistFromWall, parRoom.roomSize_m(1)-minDistFromWall), unifrnd(minDistFromWall, parRoom.roomSize_m(2)-minDistFromWall), unifrnd(minDistFromWall, parRoom.roomSize_m(3)-minDistFromWall)];
    pos_globSrc = randomRoomPos();
    pos_locSrcA = pos_globSrc;
    pos_locSrcB = pos_globSrc;
    while dist(pos_globSrc, pos_locSrcA) < minDistFromGlobSrc || dist(pos_globSrc, pos_locSrcB) < minDistFromGlobSrc
        pos_locSrcA = gensig.randomPositionForFixedDistance(parRoom.roomSize_m, distBetweenLocSrc(nn), minDistFromWall);
        pos_locSrcB = gensig.randomPositionFixedDistance(parRoom.roomSize_m, pos_locSrcA, distBetweenLocSrc(nn), minDistFromWall);
    end
    parRoom.srcPos_m = [pos_globSrc; pos_locSrcA; pos_locSrcB]';
    disp('Done!'); toc;

    %draw mic positions
    parRoom.micPos_m = [gensig.randomMicPosNearLocSrc(pos_locSrcA, micDistFromLocSrc_range, parRoom.roomSize_m', minDistFromWall); ... 
                         gensig.randomMicPosNearLocSrc(pos_locSrcB, micDistFromLocSrc_range, parRoom.roomSize_m', minDistFromWall)]';

    % draw source files (make sure no file appears more than once)
    src_file_numbers = 0;
    while isequal(src_file_numbers, 0) || length(unique(src_file_numbers)) ~= length(src_file_numbers)
        src_file_numbers = round(unifrnd(1, 50, 1, 3)); % timit files 1.wav ... 50.wav
    end
    parSource.file(1) = {[sigPath num2str(src_file_numbers(1)) '.wav']}; % global src 
    parSource.file(2) = {[sigPath num2str(src_file_numbers(2)) '.wav']}; % local src A
    parSource.file(3) = {[sigPath num2str(src_file_numbers(3)) '.wav']}; % local src B

    % Generate signals...
    strSetup = ['setup_' num2str(nn)];
    target_filename = [target_dirname '/mat_files/' strSetup '.mat'];
    disp('Generate asynchronous mic signals and save them into file:'); disp([target_filename ' ...']);
    [signals, setup] = gensig.genAsyncSignals_MultiDiffPause(parSource, parRoom, parNoise, parDiffuse, parAsync);
    
    % save matfile
    sig_z = signals.z; sig_y = signals.y; sig_s_mask = signals.s_mask; %sig_s = signals.s; %sig_s_h = signals.s_h;
    TDOA_smp = squeeze(setup.TDOAmatrix_smp(2,1,:)); setup.TDOA_smp = TDOA_smp;
    save(target_filename, 'sig_z','sig_y','sig_s_mask','setup','parSource', 'parRoom', 'parNoise', 'parDiffuse', 'parAsync');
    
    % save audio
    audiowrite([target_dirname '/wav_files/z_' strSetup '.wav'], sig_z/(max(max(abs(sig_z)))+1e-3), fs_Hz);

end



