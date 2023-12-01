% Evaluate (CL-)DXCPP Perfomance: Native vs. ACS-supported vs.
% VAD-supported
%
% Only consider the (mostly incoherent) middle segment
%
%

%% INIT

clear variables; clc; %close all;
conf = Config();


%% SET PARAMETERS: GENERAL

dbPath = conf.pathDB_test;
load([dbPath, 'parameters.mat'], 'fs_Hz', 'sigLen_s', 'n_setups');

sigLen_smp = sigLen_s*fs_Hz;
sigLen_frames = floor((sigLen_smp-2^13)/2^11) + 1;
include_CL_DXPP = 1;

frames_add = 100; % how many frames to add to range at the end. 100frames ~ DXCP-PhaT Time constant
ACSmapping = "EXP";


%% LOAD SRO AND PIACS DATA

load([dbPath 'results/est_results.mat'], 'SRO_Est', 'chi');


%% SET PARAMETERS: DXCP-PhaT

parDXCPPhaT = sroest.DXCP_PhaT_DefaultParams();
sigLen_frames = floor((sigLen_smp - parDXCPPhaT.FFTsize)/parDXCPPhaT.FFTshift) + 1;


%% PARAMETERS: CL-DXCP-PHAT

Tf = '8'; % 1,3 or 8 as text string
ControllerStructure = 'PIT1'; % exact_IMC, approx_IMC, PIT1, PI
parDXCPPhaTcl = parDXCPPhaT;
parDXCPPhaTcl.fs_Hz = fs_Hz;   % for calculation of RTFs2
parDXCPPhaTcl.WinSizeSINC = 16;    % one-sided window size for sinc (16/64/256)
% controller parameters
parDXCPPhaTcl.ForwardControl = false; % Enable forward control
parDXCPPhaTcl.L_hold = 100; % Forward control duration in frames
[parDXCPPhaTcl.K_Nomin, parDXCPPhaTcl.K_Denom]=...
    sroest.Controller_Choice(ControllerStructure, Tf, parDXCPPhaT.B_NumFr);


%% COMPUTE EVAL RESULTS: MAX-EST-ERR, SSNR, AMSC

% Prepare results
% -- OpenLoop
SRO_maxErr_default = zeros(1, n_setups);
SRO_maxErr_ACS = zeros(1, n_setups);
SRO_maxErr_VAD = zeros(1, n_setups);
AMSC_default = zeros(1, n_setups);
AMSC_ACS = zeros(1, n_setups);
AMSC_VAD = zeros(1, n_setups);
% -- ClosedLoop
SRO_maxErr_cl_default = zeros(1, n_setups);
SRO_maxErr_cl_ACS = zeros(1, n_setups);
SRO_maxErr_cl_VAD = zeros(1, n_setups);
AMSC_cl_default = zeros(1, n_setups);
AMSC_cl_ACS = zeros(1, n_setups);
AMSC_cl_VAD = zeros(1, n_setups);

% Start eval...
t = tic;
parfor nn = 1:n_setups

    disp(['Now computing for setup ' num2str(nn)]);

    % Load reference data
    data = load([dbPath 'mat_files/setup_' num2str(nn) '.mat'], 'parAsync', 'sig_s_mask', 'sig_z', 'sig_y');
    parAsync = data.parAsync; sig_s_mask = data.sig_s_mask; sig_z = data.sig_z; sig_y = data.sig_y;
    % Get range of interest: Middle Segment 
    % -- sample range
    select_smp_raw = sig_s_mask(:,1) == 0;
    idx_start_smp = find(select_smp_raw == 1, 1);
    idx_end_smp = floor(idx_start_smp + sum(select_smp_raw) + frames_add*parDXCPPhaT.FFTshift);
    select_smp = zeros(1, sigLen_smp);
    select_smp(idx_start_smp:idx_end_smp) = 1;
    select_smp = boolean(select_smp);
    % -- frame range (translated)
    relative_length = sum(select_smp_raw)/length(sig_s_mask(:,1));
    relative_start = idx_start_smp/length(sig_s_mask(:,1));
    idx_start = floor(sigLen_frames*relative_start);
    idx_end = idx_start + floor(sigLen_frames*relative_length) + frames_add;
    select = zeros(1, sigLen_frames);
    select(idx_start:idx_end) = 1;
    select = boolean(select);

    % /// OL: DEFAULT
    SRO_Est_ = squeeze(SRO_Est(nn, :));
    % MaxErr
    SRO_maxErr_default(nn) = max(abs(SRO_Est_(select) - parAsync.SRO_ppm(2)));
    % SSNR
    SRO_Est_fine = multirate.seriesFramesToSamples(SRO_Est_, parDXCPPhaT.FFTsize, parDXCPPhaT.FFTshift, sigLen_smp);
    y_synced = multirate.resampler_sinc_dyn(sig_z(:,2)', parDXCPPhaTcl.WinSizeSINC, -SRO_Est_fine, 0, [idx_start_smp, idx_end_smp]);
    % AMSC
    AMSC_default(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(select_smp), hann(2^13), 2^11));

    % /// OL: ACS-CONTROLLED
    % Obtain controlled SRO estimate
    chi_ = squeeze(chi(nn, :));   
    ACS = acs.acs__mapping(chi_, ACSmapping, "bin");
    SRO_Est_ACS = smooth.recursiveAvgDyn(SRO_Est_, 1-ACS, 2);
    % MaxErr
    SRO_maxErr_ACS(nn) = max(abs(SRO_Est_ACS(select) - parAsync.SRO_ppm(2)));
    % SSNR
    SRO_Est_fine = multirate.seriesFramesToSamples(SRO_Est_ACS, parDXCPPhaT.FFTsize, parDXCPPhaT.FFTshift, sigLen_smp);
    y_synced = multirate.resampler_sinc_dyn(sig_z(:,2)', parDXCPPhaTcl.WinSizeSINC, -SRO_Est_fine, 0, [idx_start_smp, idx_end_smp]);
    % AMSC
    AMSC_ACS(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(select_smp), hann(2^13), 2^11));

    % /// OL: VAD-CONTROLLED
    % Get VAD mask
    VADoracle_mic1 = (sig_s_mask(:,1)+sig_s_mask(:,2)) > 0; % either global and or local source active
    VADoracle_mic2 = (sig_s_mask(:,1)+sig_s_mask(:,3)) > 0; % either global and or local source active
    VADoracel_smp = VADoracle_mic1 & VADoracle_mic2;% VAD1 and VAD2
    VAD_oracle = zeros(1, sigLen_frames);
    % -- convert sample- to frame-mask
    index_vec = 1:parDXCPPhaT.FFTsize;
    for ell=1:sigLen_frames
        VAD_oracle(ell) = round(mean(VADoracel_smp(index_vec)));
        index_vec = index_vec + parDXCPPhaT.FFTshift; 
    end
    % -- shift by DXCPP Time Constant to allow fair comparison
    VAD_oracle = multirate.shift_signal(VAD_oracle, 100);
    % Obtain controlled estimate
    SRO_Est_VAD = smooth.recursiveAvgDyn(SRO_Est_, 1-VAD_oracle, 2);
    % MaxErr
    SRO_maxErr_VAD(nn) = max(abs(SRO_Est_VAD(select) - parAsync.SRO_ppm(2)));
    % SSNR
    SRO_Est_fine = multirate.seriesFramesToSamples(SRO_Est_VAD, parDXCPPhaT.FFTsize, parDXCPPhaT.FFTshift, sigLen_smp);
    y_synced = multirate.resampler_sinc_dyn(sig_z(:,2)', parDXCPPhaTcl.WinSizeSINC, -SRO_Est_fine, 0, [idx_start_smp, idx_end_smp]);
    % AMSC
    AMSC_VAD(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(select_smp), hann(2^13), 2^11));

    if include_CL_DXPP
        % Note: Input signals are truncated at end of selected range to
        % avoid unneccessary processing. this leads to slightly
        % different indexing compared to OL

        % /// CL: DEFAULT
        % Estimate...
        [SROest_ppm, ~, ~, ~, ~, ~, ~, ~, ~, ~, y_synced, ~, ~, ~] = sroest.DXCP_PhaT_CL1(sig_z(1:idx_end_smp,:), parDXCPPhaTcl);
        % -- compensate RTO
        RTO_smp = sum(SROest_ppm - parAsync.SRO_ppm(2))*1e-6*parDXCPPhaT.FFTshift;
        y_synced = multirate.shift_signal(y_synced, round(RTO_smp));
        % MaxErr
        SRO_maxErr_cl_default(nn) = max(abs(SROest_ppm(idx_start:end) - parAsync.SRO_ppm(2)));
        % AMSC
        AMSC_cl_default(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(idx_start_smp:end), hann(2^13), 2^11));

        % /// CL: VAD-CONTROLLED
        % Estimate...
        OnlineSAD_Inst = acs.OnlineSAD(VAD_oracle);
        [SROest_ppm, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, y_synced, ~, ~, ~, ~] = sroest.DXCP_PhaT_CL1_with_SAD(sig_z(1:idx_end_smp,:), parDXCPPhaTcl, OnlineSAD_Inst); 
        % -- compensate RTO
        RTO_smp = sum(SROest_ppm - parAsync.SRO_ppm(2))*1e-6*parDXCPPhaT.FFTshift;
        y_synced = multirate.shift_signal(y_synced, round(RTO_smp));
        % MaxErr
        SRO_maxErr_cl_VAD(nn) = max(abs(SROest_ppm(idx_start:end) - parAsync.SRO_ppm(2)));
        % AMSC
        AMSC_cl_VAD(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(idx_start_smp:end), hann(2^13), 2^11));

        % /// CL: ACS-CONTROLLED
        % Estimate...
        n_frames_acs_overwrite = floor((60*fs_Hz)/parDXCPPhaTcl.FFTshift);
        OnlineACS_Inst = acs.OnlineACS(ACSmapping, "bin", n_frames_acs_overwrite);
        [SROest_ppm, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, y_synced, ~, ~, ~, ~, ~, ~] = sroest.DXCP_PhaT_CL1_with_IACS(sig_z(1:idx_end_smp,:), parDXCPPhaTcl, OnlineACS_Inst); 
        % -- compensate RTO
        RTO_smp = sum(SROest_ppm - parAsync.SRO_ppm(2))*1e-6*parDXCPPhaT.FFTshift;
        y_synced = multirate.shift_signal(y_synced, round(RTO_smp));
        % MaxErr
        SRO_maxErr_cl_ACS(nn) = max(abs(SROest_ppm(idx_start:end) - parAsync.SRO_ppm(2)));
        % AMSC
        AMSC_cl_ACS(nn) = mean(mscohere(sig_y(select_smp, 2), y_synced(idx_start_smp:end), hann(2^13), 2^11));


    end
end
toc(t);


%% SAVE RESULTS

[year, month, day] = ymd(datetime);
[hour, minute, second] = hms(datetime);
target_filename = [dbPath 'results/eval_results_' num2str(year) '_' num2str(month) '_' num2str(day) '_' num2str(hour) '_' num2str(minute) '.mat'];
save(target_filename); 


%% DISPLAY AVERAGE VALUES

disp("---------- Average AMSC Values:")
disp("Avg. AMSC for OL-DXCP-Phat (plain, +SAD, +ACS): ")
disp([mean(AMSC_default), mean(AMSC_VAD), mean(AMSC_ACS)])

disp("Avg. AMSC for CL-DXCP-Phat (plain, +SAD, +ACS): ")
disp([mean(AMSC_cl_default), mean(AMSC_cl_VAD), mean(AMSC_cl_ACS)]);

disp("---------- Average AE_max Values:")
disp("Avg. AE_max for OL-DXCP-Phat [ppm] (plain, +SAD, +ACS): ")
disp([mean(SRO_maxErr_default), mean(SRO_maxErr_VAD), mean(SRO_maxErr_ACS)])

disp("Avg. AE_max for CL-DXCP-Phat [ppm] (plain, +SAD, +ACS): ")
disp([mean(SRO_maxErr_cl_default), mean(SRO_maxErr_cl_VAD), mean(SRO_maxErr_cl_ACS)])


%% PLOTS: MAX-ERR

figpos = [100, 100, 530, 170];
figure('Renderer', 'painters', 'Position', figpos); 

% -- OL
subplot(2, 1, 1);
hold on; grid on; box on;
%set(gca, 'YGrid', 'on', 'XGrid', 'none');
set(gca, 'MinorGridLineStyle', 'none');
colors = ["#0072BD", '#D95319', "#000"];
markerSize = 6;
legend_str = ["regular DXCP-PhaT", "+SAD", "+ACS"];
xlim([0, n_setups+1]);
ylim([1.5e-2, 5e2]);
[~, sorted_idcs] = sort(SRO_maxErr_default);
for ii = 1:n_setups
    nn = sorted_idcs(ii);
    plot(ii, SRO_maxErr_default(nn), 'o', 'MarkerFaceColor', colors(1), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(ii, SRO_maxErr_VAD(nn), 's', 'MarkerFaceColor', colors(2), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(ii, SRO_maxErr_ACS(nn), 'v', 'MarkerFaceColor', colors(3), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    if ii > 1
        legend_str = [legend_str, '', '', '', ''];
    end
end
yline(1, '--k', 'LineWidth', 2);
plot(1:n_setups, SRO_maxErr_default(sorted_idcs), 'color', colors(1));
plot(1:n_setups, SRO_maxErr_VAD(sorted_idcs), 'color', colors(2));
plot(1:n_setups, SRO_maxErr_ACS(sorted_idcs), 'color', colors(3));
l = legend(legend_str, 'location', 'northwest', 'NumColumns', 3, 'FontSize', 14, 'EdgeColor', 'none');
l.ItemTokenSize(1) = 20;
set(gca, 'YScale', 'log');
set(gca,'xtick', 1:n_setups);
%set(gca,'xtick', []);
set(gca,'xticklabels', sorted_idcs);
set(gca, 'ytick', [0.1 1 10 ])
set(gca,'yticklabels', [0.1 1 10]);
ax = gca;
ax.XAxis.FontSize = 12;
xtickangle(90);
ylabel('AE$_{\varepsilon,max}$ [ppm]', 'Interpreter', 'Latex', 'FontSize', 14);

% -- CL
%figpos2 = figpos + [0, 0, 0, 25]; %extra height for xlabel
%figure('Renderer', 'painters', 'Position', figpos2); 
subplot(2, 1, 2)
hold on; grid on; box on;
%set(gca, 'YGrid', 'on', 'XGrid', 'off');
set(gca, 'MinorGridLineStyle', 'none');
colors = ["#0072BD", '#D95319', "#000"];
markerSize = 6;
legend_str = ["regular CL-DXCP-PhaT", "+SAD", "+ACS"];
xlim([0, n_setups+1]); 
ylim([1.5e-2, 5e2]);
[~, sorted_idcs] = sort(SRO_maxErr_cl_default);
for ii = 1:n_setups
    nn = sorted_idcs(ii);
    plot(ii, SRO_maxErr_cl_default(nn), 'o', 'MarkerFaceColor', colors(1), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(ii, SRO_maxErr_cl_VAD(nn), 's', 'MarkerFaceColor', colors(2), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(ii, SRO_maxErr_cl_ACS(nn), 'v', 'MarkerFaceColor', colors(3), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    if ii > 1
        legend_str = [legend_str, '', '', '', ''];
    end
end
yline(1, '--k', 'LineWidth', 2);
plot(1:n_setups, SRO_maxErr_cl_default(sorted_idcs), 'color', colors(1));
plot(1:n_setups, SRO_maxErr_cl_VAD(sorted_idcs), 'color', colors(2));
plot(1:n_setups, SRO_maxErr_cl_ACS(sorted_idcs), 'color', colors(3));
l = legend(legend_str, 'location', 'northwest', 'NumColumns', 3, 'FontSize', 14, 'EdgeColor', 'none');
l.ItemTokenSize(1) = 20;
set(gca, 'YScale', 'log');
set(gca,'xtick', 1:n_setups);
set(gca,'xticklabels', sorted_idcs);
set(gca,'xticklabels', sorted_idcs);
set(gca, 'ytick', [0.1 1 10 ])
set(gca,'yticklabels', [0.1 1 10]);
ax = gca;
ax.XAxis.FontSize = 12;
xtickangle(90);
ylabel('AE$_{\varepsilon,max}$ [ppm]', 'Interpreter', 'Latex', 'FontSize', 14);
xlabel('setup id');