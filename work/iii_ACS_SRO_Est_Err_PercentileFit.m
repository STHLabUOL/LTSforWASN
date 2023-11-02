% Investigate correlation: pACS vs. SRO Est. Error
% 
%


%% INIT

clear variables; %close all; clc;
conf = Config();


%% PARAMETERS: GENERAL

dbPath = conf.pathDB_train;
load([dbPath, 'parameters.mat'], 'fs_Hz', 'sigLen_s', 'n_setups');
sigLen_smp = sigLen_s*fs_Hz;
sigLen_frames = floor((sigLen_smp-2^13)/2^11) + 1;


%% LOAD SRO AND PACS DATA

load([dbPath 'results/est_results.mat']);


%% COLLECT SRO-EST-ERR AND PRE-ACS SAMPLES FROM ALL SETUPS

chi_all = [];
sro_absErr_all = [];

for nn = 1:n_setups

    % Load ground-truth values
    load([dbPath 'mat_files/setup_' num2str(nn) '.mat'], 'parAsync', 'parDiffuse');

    % Grab SRO and ACS Estimates
    SRO_Est_ = squeeze(SRO_Est(nn, 100:end));
    chi_ = squeeze(chi(nn, 100:end));

    % Compute SRO Est. Error
    SRO_absErr = abs(SRO_Est_ - parAsync.SRO_ppm(2));

    % Samplewise ACS vs. SRO-AbsErr
    chi_all = [chi_all, chi_];
    sro_absErr_all = [sro_absErr_all, SRO_absErr];

end


%% GET SRO-EST-ERR PERCENTILES WITHIN DISCRETE CHI-BINS AND FIT EXP-CURVE

% params
binwidth = 0.01;
nth_prctile = 95;
omit_start = 0; % omit samples from start for better fitting
omit_end = 0; % omit samples from end for better fitting
lift_factor = 1; % scale fitting-samples -> bias log-curve

% set bins and prepare results
chi_bin_edges = min(chi_all):binwidth:max(chi_all);
chi_bin_edges = chi_bin_edges((1+omit_start):end-omit_end);
chi_bin_centers = chi_bin_edges(2:end) - binwidth/2;
sro_absErr_values = zeros(1, length(chi_bin_centers));
debug_valuesPerBin = zeros(1, length(chi_bin_centers));
% calc percentiles
for ii = 2:length(chi_bin_edges)
    select = chi_all > chi_bin_edges(ii-1) & chi_all <= chi_bin_edges(ii);
    sro_absErr_values(ii-1) = prctile(sro_absErr_all(select), nth_prctile)*lift_factor;
    debug_valuesPerBin(ii-1) = length(sro_absErr_all(select));
end

% curve fitting (clear up details about StartPoint argument)
g = fittype('a-b*exp(-c*x)');
f0 = fit(chi_bin_centers', log10(sro_absErr_values'), g,'StartPoint',[[ones(size(chi_bin_centers')), -exp(-chi_bin_centers')]\log10(sro_absErr_values'); 1]);
fit_coefs = coeffvalues(f0);
fit_x = linspace(chi_bin_centers(1), chi_bin_centers(end), 100000);
fit_y = 10.^f0(fit_x);
% estimate zero-crossing
chi0 = fit_x(find(abs(log10(fit_y)) == min(abs(log10(fit_y)))));
disp(['Curve zero crossing (SRO-Est-Err = 1ppm) approx. at chi0=' num2str(chi0)]);


%% PLOT BI. HISTOGRAM

yNumBins = 150;
yEdges = linspace(log10(min(sro_absErr_all)), log10(max(sro_absErr_all)), yNumBins);
[hcounts, xEdges, yEdges] =  histcounts2(chi_all, log10(sro_absErr_all), chi_bin_edges, yEdges, 'Normalization', 'probability');
yCenters = yEdges(2:end) - (yEdges(2)-yEdges(1))/2;
xCenters = xEdges(2:end) - (xEdges(2)-xEdges(1))/2;
figure; grid on;
for ii = 1:size(hcounts, 1)
   yvals = movmean(hcounts(ii,:), 5);
   fill3(xCenters(ii)*ones(1, length(yCenters)), 10.^yCenters, yvals, 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);  
   grid on;
   hold on;
end
plot3(chi_bin_centers, sro_absErr_values, zeros(size(sro_absErr_values)), 'x', 'Color', '#A2142F', 'MarkerSize', 8);
plot3(fit_x, fit_y, zeros(size(fit_y)), 'LineWidth', 3, 'Color', '#A2142F'); 
xlabel('$\chi$');
ylabel('$|\varepsilon - \hat{\varepsilon}|$ [ppm]');
zlabel('Rel. Freq.');
xlim([chi_bin_edges(1), chi_bin_edges(end)]);
view(80, 60);
title('Lin. Interpolated and Smoothed Histograms');
set(gca, 'YScale', 'log');
ylim([1e-3, 20]);

