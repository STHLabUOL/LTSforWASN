% Estimate SRO with DXCP PhaT 
% Estimate prelim. ACS
%
%


%% INIT

clear variables; %close all; %clc;
conf = Config();


%% PARAMETERS: GENERAL

dbPath = conf.pathDB_test;%conf.pathDB_train;
load([dbPath, 'parameters.mat'], 'fs_Hz', 'sigLen_s', 'n_setups');
sigLen_smp = sigLen_s*fs_Hz;


%% PARAMETERS: DXCP-PhaT

parDXCPPhaT = sroest.DXCP_PhaT_DefaultParams();
sigLen_frames = floor((sigLen_smp - parDXCPPhaT.FFTsize)/parDXCPPhaT.FFTshift) + 1;


%% ESTIMATE SRO AND PRE-ACS

% Prepare Output (save results matrix-style)
SRO_Est = zeros(n_setups, sigLen_frames); 
chi = zeros(n_setups, sigLen_frames); 

% Loop Iteration parameters and load corresponding signals
parfor nn = 1:n_setups
    %Load signals
    filename = [dbPath, 'mat_files/setup_', num2str(nn), '.mat'];
    data = load(filename, 'sig_z');
    sig_z = data.sig_z;
    tic;
    %Estimate SRO (open-loop)
    disp(['Processing for setup ' num2str(nn)]);
    [SROest_ppm, ~, ~, ~, ~, ~, ~, ~, GCSD2_avg_hist, ~] = sroest.DXCP_PhaT_Jrnl_mod(sig_z, parDXCPPhaT);       
    SRO_Est(nn, :) = SROest_ppm;
    % Compute chi a.k.a. prelim. ACS
    [~, ~, preACS, ~] = acs.acs(GCSD2_avg_hist(1:4097,:), "EXP", "bin");
    chi(nn, :) = preACS;
    toc;     
end    


%% SAVE RESULTS

if ~isfolder([dbPath, 'results'])
    mkdir([dbPath, 'results']);
end
if ~isfile([dbPath '/results/est_results.mat'])
    save([dbPath '/results/est_results.mat'], 'SRO_Est', 'chi');
else
    error('Warning: Did not save results. A file with the same name already exists!');
end
