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

%% new data format — acoustic cavitation from KNN results file
% knn_data (loaded in DA_master.m) and bubble_idx must be set before this script runs.
% Old format (laser-cav, .R/.t struct):
%   exp_data = load([data_filepath,data_filename]);
%   Rnew = exp_data.R;  t = exp_data.t;
bubble_table = knn_data.exp_data{bubble_idx};
t    = bubble_table.time(:);
Rexp = bubble_table.radius(:); 

t = t-t(1); % time shift
R0  = R0_guess;    % m, from knn_data.results_all.R0(bubble_idx) set in DA_master.m
yth = (Rexp ./ R0)';     % normalised radius (peaks above 1 at maximum expansion)
sprintf('max yth: %0.5f',max(yth))

% delete NaNs
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

% Use the full experimental record as the assimilation window.
% With sparse acoustic data (no clear collapse/rebound), using all points is correct.
% peak_time = length(yth) → l = n in main_En4D_peaks, one window over all data.
peak_time = length(yth);

% Old collapse-finding logic (laser-cav / rebounding bubble case):
%peak_indices = find(islocalmin(yth) & yth < 0.3 * max(yth));
%if isempty(peak_indices)
%    [~, peak_time] = min(yth);
%else
%    peak_time = peak_indices(num_peaks);
%end

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