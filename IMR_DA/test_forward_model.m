% test_forward_model.m
% Runs the bubble forward model once with KNN-optimal parameters for bubble 1
% and plots dimensionless radius and forcing pressure vs time.
% Run this from the IMR_DA directory.

clear; clc;
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'IMR-vanilla', 'functions'));
addpath(fullfile(script_dir, 'IMR-vanilla', 'src'));

%% Parameters — must match DA_master.m
data_filepath = 'data/opt_param_data/';
knn_filename  = 'IMR_FullDataset_Results_10ag_NH_2025-07-31-11-30-09.mat';
knn_data      = load([data_filepath, knn_filename]);

bubble_idx = 1;
G_guess    = knn_data.results_all.G(bubble_idx);
mu_guess   = knn_data.results_all.mu(bubble_idx);
R0_guess   = knn_data.results_all.R0(bubble_idx);

model    = 'neoHook';
NT       = 240;
NTM      = 240;
Pext_type     = 'ga';
Pext_Amp_Freq = [-24e6 1e6 5e-6 3.7];
Tgrad    = 0;
Tmgrad   = 0;
Cgrad    = 0;
comp     = 1;
disptime = 0;
num_peaks = 1;
data_type = 'exp';

G1_guess        = 1e9;
alpha_guess     = 0.5;
lambda_nu_guess = 0.1;
est_params      = [];

visco_params = struct('G', G_guess, 'G1', G1_guess, 'mu', mu_guess, ...
    'alpha', alpha_guess, 'lambda_nu', lambda_nu_guess);

%% Load experimental data (sets tspan, n, yth, peak_time)
import_data_exp

%% Initialize state vector (sets x0_true, N, and all dimensionless groups)
R0 = R0_guess;
initialize

%% Build vars cell (same order as main_En4D_peaks.m)
tspan_star = tspan / t0;
vars = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    Cgrad comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};

%% Run forward model once over the full time span
xi = x0_true';
xi(3) = log(xi(3));   % P is stored as log(P) in state vector

fprintf('Running forward model: G=%.0f Pa, mu=%.4f Pa·s, R0=%.2e m\n', ...
    G_guess, mu_guess, R0_guess);
fprintf('t0 = %.3e s, tspan_star = %.1f, dt_star = %.1f, w_star = %.4f\n', ...
    t0, tspan_star, dt_star, w_star);

[xf, ~] = f(0, tspan_star, xi, vars, []);

fprintf('Forward model complete.\n');
