clear all % clear workspace
rng(23); % set random number seed for reproducability
dbstop if error % stop execution if there is an error

% construct the appropriate path depending on the system this is run on
if ispc 
    root = 'L:/';
elseif ismac
    root = '/Volumes/labs/';
elseif isunix 
    root = '/media/labs/';
end

% read in list of subjects
group_list = readtable("./group_list.csv");
subjects = cellstr(group_list.record_id)';
% first get behavioral file from task run without resistance, then form
% task run with resistance
loaded_unloaded = {'unloaded', 'loaded'};

% set results directory
resdir = './modeling_results';
% get timestamp
timestamp = datestr(datetime('now'), 'mm_dd_yy_THH-MM-SS');

% This will read all (relevant) files from the directory pointed to
[big_table, subj_mapping] = merge_horizon_data(root, subjects, group_list, loaded_unloaded);
outpath_beh = sprintf([resdir '/beh_%s_HC_iMUD.csv'], timestamp);
writetable(big_table, outpath_beh);

% Reads in the above 'outpath_beh' file and fit this file
% outpath_beh = './modeling_results/beh_12_24_24_T11-08-48_HC_iMUD.csv';
fits = fit_model(outpath_beh);
fits.id = subj_mapping(:, 1);
outpath_fits = sprintf('./modeling_results/HC_MUD_fits_%s.csv', timestamp);
writetable(fits, outpath_fits);

