% plot_bubble_dynamics.m
% Runs funIMRsolver (same ODE physics as f.m / DA forward model) once with
% user-specified material parameters and plots key bubble dynamics outputs.
% Run from the project root: IMR_data_assimilation/

clear; close all; clc;

addpath(genpath('functions'));
addpath(genpath('src'));

%% ---- Parameters -------------------------------------------------------

model  = 'neoHook';     % 'neoHook' | 'linkv' | 'sls' | 'nhzen' | 'fung'
G      = 1e3;          % shear modulus (Pa)
mu     = 0.01;          % viscosity (Pa·s)
G1     = 0;           % dashpot modulus (Inf → KV limit)
alpha  = 0;             % Fung nonlinearity (0 → neoHook)
lambda_nu = 0;          % nonlinear viscosity exponent

matprop = struct('G',G,'mu',mu,'G1',G1,'alpha',alpha,'lambda_nu',lambda_nu);

R0     = 1e-6;         % initial bubble radius (m)
tspan  = 50e-6;          % simulation duration (s)
NT     = 30;           % nodes in bubble (>=500 for accuracy)
NTM    = 30;            % nodes in medium

Pext_type     = 'ga';   % Rayleigh-type collapse from out-of-equilibrium IC
                        % ~226 for 11 kPa PA gel; tune to match eqR
Pext_Amp_Freq = [-24e6 5e-6 2*pi*1e6 3.7];

Tgrad    = 0;   % 1 = temperature gradient in bubble
Tmgrad   = 0;   % 0 = cold liquid assumption
Cgrad    = 0;   % 1 = vapor–gas diffusion
Dim      = 0;   % 1 = dimensional outputs
comp     = 1;   % 1 = Keller-Miksis; 0 = Rayleigh-Plesset
disptime = 0;

% REq: equilibrium radius ratio R_eq/R0 (0 = auto-compute inside solver)
REq = 0;

IMRsolver_RelTolX = 1e-7;

%% ---- Run solver -------------------------------------------------------

tic;
[t, R, U, P, S, T, C, Tm, tdel, Tdel, Cdel, Tau, pf] = funIMRsolver( ...
    model, matprop, tspan, R0, NT, NTM, ...
    Pext_type, Pext_Amp_Freq, disptime, Tgrad, Tmgrad, Cgrad, Dim, comp, REq, IMRsolver_RelTolX,0);


%% ---- Plot -------------------------------------------------------------

figure('Name','Bubble Dynamics','Position',[100 100 900 700]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- Radius ---
nexttile;
plot(t, R/R0, 'b', 'LineWidth',1.5);
xlabel('Time (µs)'); ylabel('R/R0');
title('Bubble Radius');
grid on;

% --- Internal pressure ---
nexttile;
Pmt = IMRcall_parameters(R0,G,G1,mu);
P_inf = Pmt(19);
plot(t, pf./1e6, 'k', 'LineWidth',1.5);
xlabel('Time (µs)'); ylabel('P (MPa)');
title('Bubble Pressure');
grid on;

