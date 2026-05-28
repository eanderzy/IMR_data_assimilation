% This file exports data in vector yth, it must be re-written for each data
% set to ensure the data is formatted correctly

%dataset = 8; %between 1-20 for water case
%addpath ./IMR-vanilla/functions
%addpath ./IMR-vanilla/src
% load('/scratch/jspratt/EnKS/exp_data/11kPa_PA/RofTdata.mat')
%% old data format
%{
load([data_filepath,data_filename]);
Rexp = Rnew(dataset,:)*1e-6;
[R0,max_index] = max(Rexp);

yth = Rexp(max_index:end)./R0;
t = t(max_index:end);

tspan = t(end)-t(1);
n = length(t)-1;

if exist('l') == 0
    l = n;
end
timesteps = n+1;

% Find peak_time
peak_indices = find(islocalmin(yth));
peak_time = peak_indices(num_peaks);
%}

%% new data format
exp_data = load([data_filepath,data_filename]);
Rnew = exp_data.R;
t = exp_data.t;

% fixing shape if needed
if length(Rnew(1,:)) == 1
    Rnew = Rnew';
end
if length(t(1,:)) == 1
    t = t';
end

Rexp = Rnew*1e-6;
% CHANGE: The two lines below set R0 = Rmax and crop the data to start at
% the maximum radius. For acoustic forcing with full growth+collapse:
%   1. Replace R0 with the prescribed equilibrium/stress-free radius R0_guess
%      (defined in DA_master.m from your KNN results).
%   2. Remove the max_index crop so yth spans the full time series from nucleation.
% Replace with:
%   R0 = R0_guess;                 % prescribed stress-free radius (m)
%   yth = Rexp ./ R0;              % normalize by stress-free radius, not Rmax
%   (keep t unchanged - do not crop to max_index)
% Note: yth will now start at values << 1 (growth phase) and peak at Rmax/R0 >> 1.

% for IMR
%[R0,max_index] = max(Rexp);       
%yth = Rexp(max_index:end)./R0;   
%t = t(max_index:end);             

R0=R0_guess;
yth = Rexp ./ R0;


% delet NaNs
kk = 1;
for jj = 1:length(yth)
    if isnan(yth(jj))
    else
        yth(kk) = yth(jj);
        t(kk) = t(jj);
        kk = kk + 1;
    end
end
yth = yth(1:kk-1);
t = t(1:kk-1);
%

tspan = t(end)-t(1);
n = length(t)-1;

if exist('l') == 0
    l = n;
end
timesteps = n+1;

% Find peak_time
% CHANGE: islocalmin finds collapse minima in yth. With the full growth+collapse
% window yth is NOT normalized to start at 1, so the first local minimum is the
% first collapse. num_peaks=1 (set in DA_master.m) selects the end of the
% single growth+collapse cycle as the DA window boundary.
% If yth has spurious early minima before growth peaks, add a depth threshold:
%   peak_indices = find(islocalmin(yth) & yth < 0.2*max(yth));
peak_indices = find(islocalmin(yth));
peak_time = peak_indices(num_peaks);

% Exceptionally, for test run without collapse points
%{
if method == 'EnKS'
    
    if num_peaks > 2 % checking for indexing error present in some runs
        collapse_1_idx = peak_indices(2);
    else
        collapse_1_idx = peak_indices(1);
    end
    
    yth = [yth(1:collapse_1_idx-2),yth(collapse_1_idx+2:peak_time-2), ...
        yth(peak_time+2:end)];
    t = [t(1:collapse_1_idx-2),t(collapse_1_idx+2:peak_time-2), ...
        t(peak_time+2:end)];
    
    tspan = t(end) - t(1);
    n = length(t) - 1;
    timesteps = n+1;
end
%}