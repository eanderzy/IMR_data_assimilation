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
N_bubbles = 1;%length(knn_data.exp_data);

num_peaks = 1; % single growth+collapse cycle per acoustic experiment

%% Data assimilation parameters

method = 'En4D'; % data assimilation method ('En4D','EnKS',('EnKF'))

% Fixed parameters across all experiments
%G_guess = 500;
%mu_guess = 0.05;
%R0_guess = 1e-6;
G1_guess        = 1e9;
alpha_guess     = 0.5;
lambda_nu_guess = 0.1;
% Per-experiment G_guess, mu_guess, R0_guess are set inside the loop below
q = 48; % Number of ensemble members
std = 0.01; % expected standard deviation of measurements;
init_scheme = 2; % leave as 2, initializes ensemble with truth + noise

% The following are ending criteria for iterative optimization:
epsilon = 1e-5; % threshold difference in norm of state vector between steps
max_iter = 5; % max # of iterations until end optimization (5 to 10 is usually good)

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

Tgrad = 1; % Thermal effects inside bubble
Tmgrad = 0; % Thermal effects outside bubble
Cgrad = 1; % Vapor diffusion effects
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

%% Loop over all experiments
est_params = [];
R = std^2; % Measurement error covariance

for bubble_idx = 1:N_bubbles
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

    save(sprintf('./bubble%02d_%s_%s_q%d_G%.0f_mu%.4f_%speaks%s.mat', ...
        bubble_idx, method, model, q, G_guess, mu_guess, num_peaks, ...
        datestr(now,'yyyy-mm-dd_HH-MM')), '-v7.3')
end

%% Plotting
%{
if strcmp(method,'En4D')
    plot_En4D_exp
elseif strcmp(method,'EnKS')
    plot_Kalman_exp
end
%}
