% Data assimilation master file
% Author: Jean-Sebastien Spratt -- jspratt@caltech.edu
% Edits: 5/28/26 Eleanor Anderson-Zych -- eanderzy@umich.edu
%       - to read in acoustic cavitation data
% This file is a wrapper which runs a variety of data assimilation methods
% on laser-induced cavitation radius vs time data. It uses the code from
% Spratt et al. (2020). Depending on which data assimilation method and
% physical model is used, parameters from different sections below must be
% specified

clear all
clc

%% Data import
% see import_data_exp.m file for details on data loading. The following is
% meant for data in the same format as RofT data from Selda's experiments.
% If the format is different, the import_data_exp.m file will need to be
% modified accordingly

data_type     = 'exp';
%data_set = 'SoftPA_nobeads';
%data_filepath = (['example_data/']);
%data_filename = 'A_E2_002.mat';
data_filepath = 'data/opt_param_data/';
knn_filename  = 'IMR_FullDataset_Results_03ag_NH_2025-07-31-10-57-31.mat';

% Load KNN results once: exp_data{i} has .time (µs) and .radius (µm) per bubble;
% results.G(i), results.mu(i), results.R0(i) are per-bubble optimal parameters.
% Verify field names match your actual results struct before running.
knn_data  = load([data_filepath, knn_filename]);
N_bubbles = length(knn_data.exp_data);

num_peaks = 1; % single growth+collapse cycle per acoustic experiment

%% Data assimilation parameters

method = 'En4D'; % data assimilation method ('En4D','EnKS',('EnKF'))

% Fixed parameters across all experiments
%G_guess = 500;
%mu_guess = 0.05;
%R0_guess = 1e-6;
G1_guess        = 0;
alpha_guess     = 0.05;
lambda_nu_guess = 0.01;
% Per-experiment G_guess, mu_guess, R0_guess are set inside the loop below
q = 48; % Number of ensemble members
std = 0.01; % expected standard deviation of measurements;
init_scheme = 2; % leave as 2, initializes ensemble with truth + noise

% The following are ending criteria for iterative optimization:
epsilon = 1e-5; % threshold difference in norm of state vector between steps
max_iter = 10; % max # of iterations until end optimization (5 to 10 is usually good)

% IEnKS only:
if strcmp(method, 'EnKS')
    l = 3; %lag of smoother
end

% Note: Beta coefficients are fixed as equal here. They can be modified in
% the main_En4D_peaks and main_mda_peaks files if needed, if more weight
% needs to be attributed to earlier or later points in assimilation

%% Modeling parameters
model = 'neoHook'; % 'neoHook','nhzen','sls','linkv','fung','fung2','fungexp','fungnlvis'
NT = 30; % Amount of nodes inside the bubble
NTM = 30; % Amount of nodes outside the bubble
Pext_type = 'ga'; % Type of external forcing
ST = 0.072; % (N/m) Liquid Surface Tension

Tgrad = 0; % Thermal effects inside bubble
Tmgrad = 0; % Thermal effects outside bubble
Cgrad = 0; % Vapor diffusion effects
comp = 1; % Activates the effect of compressibility (0=Rayleigh-Plesset, 1=Keller-Miksis)
disp_timesteps = 1; % displays timesteps in En4D run (for debugging)

% following should not be changed (untested):
disptime = 0; % 1 = display simulation time
Dim = 0; % 1 = output variables in dimensional form

%% Covariance Inflation parameters
% The following is only used for the IEnKS. Default values provided below
% See Spratt et al. (2020) section 3.3 for details on theta, lambda

CI_scheme = 2; % 1 is scalar alpha CI, 2 is RTPS (BEST)
CI_theta = 0.7; % multiplicative covariance parameter (0.5 < theta < 0.95)
CI_add = 0; % Set to 1 for additive covariance (else 0)
beta = 1.02; % additive covariance parameter (lambda in paper) (1.005 < beta < 1.05)
alpha = 0.005; % random noise in forecast step (only for EnKF)

%% Spread of parameters in the ensemble
% Material parameter fractional spreads derived from KNN population statistics (CV = std/mean)
Rspread         = 0.01;
Uspread         = 0.1;
Pspread         = 0.1;
Sspread         = 0.1;
tauspread       = 0.1;
Cspread         = 0.001;
Tmspread        = 0.0005;
Brspread        = 0.01;
fohspread       = 0.01;
Caspread        = knn_data.weighted_stats.G.std  / knn_data.weighted_stats.G.mean;
Respread        = knn_data.weighted_stats.mu.std / knn_data.weighted_stats.mu.mean;
Despread        = 0;
alphaspread     = 0;
lambda_nuspread = 0;
R0spread        = knn_data.weighted_stats.R0.std / knn_data.weighted_stats.R0.mean;
%{
Rspread = 0;
Uspread = 0;
Pspread = 0;
Sspread = 0;
tauspread = 0;
Cspread = 0;
Tmspread = 0;
Brspread = 0;
fohspread = 0;
Caspread = 0;
Respread = 0;
Despread = 0; % set to 0 if not used in model
alphaspread = 0; % set to 0 if not used in model
lambda_nuspread = 0; % set to 0 if not used in model
%}

Caspread = min(Caspread, 0.25);
Respread = min(Respread, 0.25);
R0spread = min(R0spread, 0.25);

% Loop over all experiments
est_params = [];
R = std^2; % Measurement error covariance
%{
for bubble_idx = 1%1:N_bubbles
    fprintf('=== Bubble %d / %d ===\n', bubble_idx, N_bubbles);

    % Per-experiment priors from KNN-NRMSE search
    G_guess  = knn_data.results_all.G(bubble_idx);   % Pa
    mu_guess = knn_data.results_all.mu(bubble_idx);   % Pa·s
    R0_guess = knn_data.results_all.R0(bubble_idx);   % m

    visco_params = struct('G',G_guess,'G1',G1_guess,'mu',mu_guess, ...
        'alpha',alpha_guess,'lambda_nu',lambda_nu_guess);

    if strcmp(method,'En4D')
        main_En4D_peaks
    elseif strcmp(method,'EnKS')
        main_mda_peaks
    end

    save(sprintf('./test_bubble%02d_%s_%s_q%d_G%.0f_mu%.4f_%s.mat', ...
        bubble_idx, method, model, q, G_guess, mu_guess, ...
        datestr(now,'yyyy-mm-dd_HH-MM')), '-v7.3')
end
%}

%% Results collection 

% Pre-allocate results storage
DA_results = struct();
DA_results.G      = zeros(N_bubbles, 1);
DA_results.mu     = zeros(N_bubbles, 1);
DA_results.R0     = zeros(N_bubbles, 1);
DA_results.NRMSE  = zeros(N_bubbles, 1);
DA_results.status = zeros(N_bubbles, 1); % 1=converged, 0=failed

for bubble_idx = 1:N_bubbles
    fprintf('=== Bubble %d / %d ===\n', bubble_idx, N_bubbles);

    G_guess  = knn_data.results_all.G(bubble_idx);
    mu_guess = knn_data.results_all.mu(bubble_idx);
    R0_guess = knn_data.results_all.R0(bubble_idx);

    visco_params = struct('G',G_guess,'G1',G1_guess,'mu',mu_guess, ...
        'alpha',alpha_guess,'lambda_nu',lambda_nu_guess);

    try
        if strcmp(method,'En4D')
            main_En4D_peaks
        elseif strcmp(method,'EnKS')
            main_mda_peaks
        end

        % Extract final estimates using same logic as plot_En4D_exp.m
        Uc = sqrt(P_inf/rho);
        Ca_final = exp(x_est(2*NT+NTM+7, jj-1));
        Re_final = exp(x_est(2*NT+NTM+8, jj-1));
        R0_final = exp(x_est(2*NT+NTM+9, jj-1));

        G_final  = P_inf / Ca_final;
        mu_final = (P_inf * R0_final) / (Re_final * Uc);

        % Compute NRMSE with final parameters
        xf_final = [x1, squeeze(mean(E2, 2))];
        R_sim    = xf_final(1, 1:length(yth));
        NRMSE_final = sqrt(mean((R_sim - yth).^2)) / mean(yth);

        % Store
        DA_results.G(bubble_idx)      = G_final;
        DA_results.mu(bubble_idx)     = mu_final;
        DA_results.R0(bubble_idx)     = R0_final;
        DA_results.NRMSE(bubble_idx)  = NRMSE_final;
        DA_results.status(bubble_idx) = 1;

        fprintf('  G=%.1f Pa, mu=%.4f Pa.s, R0=%.3e m, NRMSE=%.4f\n', ...
            G_final, mu_final, R0_final, NRMSE_final)

    catch ME
        fprintf('  FAILED: %s\n', ME.message)
        DA_results.status(bubble_idx) = 0;
    end

    save(sprintf('./bubble%02d_%s_%s_q%d_%s.mat', ...
        bubble_idx, method, model, q, datestr(now,'yyyy-mm-dd_HH-MM')), '-v7.3')
end

%% Weighted statistics - NRMSE-weighted mean and std for each parameter

ok  = DA_results.status == 1;
nrmse_ok = DA_results.NRMSE(ok);

% Weights: inverse NRMSE, normalized to sum to 1
% Lower NRMSE = better fit = higher weight
raw_weights = 1 ./ (nrmse_ok + 1e-6);
w = raw_weights / sum(raw_weights);

G_ok  = DA_results.G(ok);
mu_ok = DA_results.mu(ok);
R0_ok = DA_results.R0(ok);

% Weighted mean
DA_results.weighted_stats.G.mean  = sum(w .* G_ok);
DA_results.weighted_stats.mu.mean = sum(w .* mu_ok);
DA_results.weighted_stats.R0.mean = sum(w .* R0_ok);

% Weighted std: sqrt(sum(w*(x-mean)^2))
DA_results.weighted_stats.G.std  = sqrt(sum(w .* (G_ok  - DA_results.weighted_stats.G.mean).^2));
DA_results.weighted_stats.mu.std = sqrt(sum(w .* (mu_ok - DA_results.weighted_stats.mu.mean).^2));
DA_results.weighted_stats.R0.std = sqrt(sum(w .* (R0_ok - DA_results.weighted_stats.R0.mean).^2));

% Median and IQR as robust alternatives
DA_results.weighted_stats.G.median  = median(G_ok);
DA_results.weighted_stats.mu.median = median(mu_ok);
DA_results.weighted_stats.R0.median = median(R0_ok);

DA_results.weighted_stats.G.iqr  = iqr(G_ok);
DA_results.weighted_stats.mu.iqr = iqr(mu_ok);
DA_results.weighted_stats.R0.iqr = iqr(R0_ok);

%% Print summary
fprintf('\n========== DA Results Summary ==========\n')
fprintf('Bubbles converged: %d / %d\n', sum(ok), N_bubbles)
fprintf('\nWeighted mean +/- weighted std:\n')
fprintf('  G  = %.1f +/- %.1f Pa\n',   DA_results.weighted_stats.G.mean,  DA_results.weighted_stats.G.std)
fprintf('  mu = %.4f +/- %.4f Pa.s\n', DA_results.weighted_stats.mu.mean, DA_results.weighted_stats.mu.std)
fprintf('  R0 = %.3e +/- %.3e m\n',    DA_results.weighted_stats.R0.mean, DA_results.weighted_stats.R0.std)
fprintf('\nMedian [IQR]:\n')
fprintf('  G  = %.1f [%.1f] Pa\n',   DA_results.weighted_stats.G.median,  DA_results.weighted_stats.G.iqr)
fprintf('  mu = %.4f [%.4f] Pa.s\n', DA_results.weighted_stats.mu.median, DA_results.weighted_stats.mu.iqr)
fprintf('  R0 = %.3e [%.3e] m\n',    DA_results.weighted_stats.R0.median, DA_results.weighted_stats.R0.iqr)
fprintf('========================================\n')

%% Compare DA vs KNN estimates
fprintf('\n=== DA vs KNN comparison (converged bubbles) ===\n')
knn_G  = knn_data.results_all.G(ok);
knn_mu = knn_data.results_all.mu(ok);
knn_R0 = knn_data.results_all.R0(ok);

pct_diff_G  = 100*(G_ok  - knn_G)  ./ knn_G;
pct_diff_mu = 100*(mu_ok - knn_mu) ./ knn_mu;
pct_diff_R0 = 100*(R0_ok - knn_R0) ./ knn_R0;

fprintf('Mean pct diff G:  %.1f%%\n', mean(pct_diff_G))
fprintf('Mean pct diff mu: %.1f%%\n', mean(pct_diff_mu))
fprintf('Mean pct diff R0: %.1f%%\n', mean(pct_diff_R0))

%% Save final results
save('DA_results_all.mat', 'DA_results', '-v7.3')
fprintf('\nResults saved to DA_results_all.mat\n')

%% Optional: plot parameter distributions
figure(200)
clf
subplot(3,1,1)
histogram(G_ok/1e3, 20)
xlabel('G [kPa]'); ylabel('count'); grid on
xline(DA_results.weighted_stats.G.mean/1e3, 'r-', 'linewidth', 2)
title('Shear modulus distribution')

subplot(3,1,2)
histogram(mu_ok, 20)
xlabel('\mu [Pa·s]'); ylabel('count'); grid on
xline(DA_results.weighted_stats.mu.mean, 'r-', 'linewidth', 2)
title('Viscosity distribution')

subplot(3,1,3)
histogram(R0_ok*1e6, 20)
xlabel('R_0 [\mum]'); ylabel('count'); grid on
xline(DA_results.weighted_stats.R0.mean*1e6, 'r-', 'linewidth', 2)
title('Equilibrium radius distribution')
%% Plotting
%{
if strcmp(method,'En4D')
    plot_En4D_exp
elseif strcmp(method,'EnKS')
    plot_Kalman_exp
end
%}
