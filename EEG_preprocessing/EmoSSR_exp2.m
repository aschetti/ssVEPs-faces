
%% EmoSSR_exp2

close all; clear all; clc;

% experiment settings
randseed = rng(9001); % seed for pseudo-random number generator (it's over 9000!)
addpath(genpath('E:\Experiments\EmoSSR_faces\toolboxes_functions\')); % add path to EEGLAB + plugins & various functions
expname = 'exp2'; % experiment name
pathexp = ['E:\Experiments\EmoSSR_faces\' expname '\']; % main directory
pathdata = [pathexp 'EEG\data\']; % where to get the .bdf files
pathanalysis = [pathexp 'EEG\analysis\']; % where to store the logfile of artifact rejection
Avg.channs = 'E:\Experiments\EmoSSR_faces\toolboxes_functions\leipzig68_Christopher.epf'; % channel locations
Avg.channs64 = 'E:\Experiments\EmoSSR_faces\toolboxes_functions\leipzig64_Christopher.locs'; % channel locations for plot_topo.m
begin_epoch = 0; % begin epoch (in seconds)
end_epoch = 7.5; % end epoch (in seconds)
% conditions:
% 1: upright angry, female identity
% 2: upright neutral, female identity
% 3: upright irregular, female identity
% 4: inverted angry, female identity
% 5: inverted neutral, female identity
% 6: inverted irregular, female identity
% 11: upright angry, male identity
% 12: upright neutral, male identity
% 13: upright irregular, male identity
% 14: inverted angry, male identity
% 15: inverted neutral, male identity
% 16: inverted irregular, male identity
Avg.num = [1 2 3 4 5 6 11 12 13 14 15 16]; % indices conditions (all identities)
events = {'1' '2' '3' '4' '5' '6' '11' '12' '13' '14' '15' '16'}; % triggers events (all identities)
Avg.numfem = [1 2 3 4 5 6]; % indices female identity
Avg.nummale = [11 12 13 14 15 16]; % indices male identity

Avg.pathin_low_emotion_variability = [pathexp 'EEG\analysis\total\low_emotion_variability\']; % where to store the analysis of face identity with low emotion variability (i.e., female)
Avg.pathin_high_emotion_variability = [pathexp 'EEG\analysis\total\high_emotion_variability\']; % where to store the analysis of face identity with high emotion variability (i.e., male)
Avg.pathout_low_emotion_variability = [pathexp 'EEG\analysis\mean\low_emotion_variability\']; % where to save the grand average of face identity with low emotion variability (i.e., female)
Avg.pathout_high_emotion_variability = [pathexp 'EEG\analysis\mean\high_emotion_variability\']; % where to save the grand average of face identity with high emotion variability (i.e., male)

logfile = []; elapsed = []; % preallocate logfile and matrix with elapsed time (in seconds) of single-subject preprocessing

filenames = dir([pathdata '*.bdf']); % read file names in folder path (puts it in a structure called filenames)

%%

% start preprocessing
for isub = 1:numel(filenames) % loop through files
    
    % if the file name contains "part 2", skip it and continue from the next file
    if ~isempty(strfind(filenames(isub).name, 'part2'))
        disp('***********************')
        disp(['Skipping ' filenames(isub).name(1:end-4) ' (it has already been merged)'])
        disp('***********************')
        continue
    end
    
    tic % start timer
    
    % file name
    if ~isempty(strfind(filenames(isub).name, 'vpn08')) || ~isempty(strfind(filenames(isub).name, 'vpn23'))
        dispname = filenames(isub).name(1:end-10);
    else
        dispname = filenames(isub).name(1:end-4);
    end
    
    % display participant number
    disp('***********************')
    disp(['Processing ' dispname '...'])
    disp('***********************')

    % import .bdf files
    if ~isempty(strfind(dispname, 'vpn08')) || ~isempty(strfind(dispname, 'vpn23')) % EXCEPTION: data of participants 8 and 23 are split in 2 parts
        % load all 68 channels (64 EEG + 4 EOG), average reference
        EEG1 = pop_biosig([pathdata dispname '_part1.bdf'], 'channels', [1:68], 'ref', [1:64], 'refoptions', {'keepref' 'on'}); % (part 1)
        EEG2 = pop_biosig([pathdata dispname '_part2.bdf'], 'channels', [1:68], 'ref', [1:64], 'refoptions', {'keepref' 'on'}); % (part 2)
        EEG = pop_mergeset(EEG1, EEG2); % merge EEG data
        clear EEG1 EEG2  % clear temporary EEG files
    else
        EEG = pop_biosig([pathdata dispname '.bdf'], 'channels', [1:68], 'ref', [1:64], 'refoptions', {'keepref' 'on'});
    end
    
    % resampling
    if EEG.srate ~= 256 % if sampling rate is not 256 Hz...
        EEG = pop_resample(EEG, 256); % ... resample
    end
    
    %     pop_eegplot(EEG, 1, 1, 1); % check data
    
    EEG = pop_chanedit(EEG, 'load', {Avg.channs, 'filetype', 'besa (elp)'}); % assign channel locations
    %     topoplot([],EEG.chanlocs, 'style', 'blank', 'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo); % see channel locations
    EEG = pop_rmbase(EEG, 'Warning', 'off'); % remove DC offset
    
    %     pop_eegplot(EEG, 1, 1, 1); % check data
    
    % extract epochs
    EEG = pop_epoch(EEG, events, [begin_epoch end_epoch]);
    EEG = pop_selectevent(EEG, 'type', cellfun(@str2num, events), 'deleteevents', 'on'); % eliminate non-condition triggers (e.g., 253) that can cause errors later on
    
    %     pop_eegplot(EEG, 1, 1, 1); % check data
    
    % for log file
    for irej = 1:numel(EEG.epoch)
        trigtotepochs{irej, :} = EEG.epoch(1, irej).eventtype; % create vector with the sum of total epochs for each condition
    end
    if ~isempty(strfind(dispname, 'vpn08')) || ~isempty(strfind(dispname, 'vpn23')) % probably due to merging EEG datasets, the cell2mat function does ont work on these participants... use cellfun
        trigtotepochs = cellfun(@str2num, trigtotepochs)'; % convert this vector from cell to number (otherwise the "sum" function won't work) and transpose
    else
        trigtotepochs = cell2mat(trigtotepochs)'; % convert this vector from cell to number (otherwise the "sum" function won't work) and transpose
    end
    for k = 1:numel(Avg.num)
        totepochs(:, k) = sum(trigtotepochs == Avg.num(k)); % find all triggers belonging to each condition and sum them
    end
    
    %%%%%%%%%%%%%%%%%%
    %%% ARTIFACT REJECTION %%%
    %%%%%%%%%%%%%%%%%%
    
    % channel interpolation & artifact rejection using FASTER
    cfg = []; % create Fieldtrip-like structure
    cfg.datachan = 1:64; % select scalp electrodes
    %     cfg.timerange = [begin_epoch end_epoch-0.5];
    % The field above allows us to select the time range (within the epoch) that should be scanned for artifacts.
    % Because Christian would use the last 0.5s only as buffer for his spatial filter, that segment could be kept noisy.
    % However, cutting out the last 0.5s creates a huge spike that distorts the signal so much that it even screws up the grand averages.
    % Therefore, I prefer to discard a few more trials but avoid such a problem.
    cfg.thresh = [3 3 0 3 3 12]; % see help eegF_FASTER for a description of each number
    [gen_bad_chans, EEG, trials2remove_FASTER] = eegF_FASTER(cfg, EEG); % run eegF_FASTER function
    EEG = pop_select(EEG, 'notrial', trials2remove_FASTER); % remove bad epochs
    
    % (following https://sccn.ucsd.edu/wiki/Chapter_01:_Rejecting_Artifacts)
    % detect abrupt spikes or flat activity based on kurtosis
    [EEG, ~, ~, nrej_rejkurt] = pop_rejkurt(EEG, ...
        1, ... % data to reject on (1: electrode data)
        [1:64], ... % channels (exclude external channels)
        5, ... % single-channel kurtosis limit (in standard deviations)
        5, ... % all-channel kurtosis limit (in standard deviations)
        0, ... % do not superpose rejection marks on previously marks stored in the dataset
        1, ... % rejection of marked trials (0: do not reject but store the marks; 1: reject)
        1, ... % visualization type (calls eegplot)
        [], ... % topcommand (deprecated)
        0); % plots (1: on; 0: off)
    
    % reject epochs using spectral estimates
    % blinks (high power at low frequencies)
    [EEG, trials2remove_rejspecBlinks] = pop_rejspec(EEG, ...
        1, ... % data to reject on (1: electrode data)
        'elecrange', [1:64], ... % channels (exclude external channels)
        'method', 'fft', ... % method to compute spectrum
        'threshold', [-50 50], ... % threshold limits (in dB)
        'freqlimits', [0 2], ... % frequency limits (in Hz)
        'eegplotplotallrej', 1, ... % do not superpose rejection marks on previously marks stored in the dataset
        'eegplotreject', 1); % rejection of marked trials (0: do not reject but store the marks; 1: reject)
    
    % muscle artifacts (high power at high frequencies)
    [EEG, trials2remove_rejspecMuscle] = pop_rejspec(EEG, ...
        1, ... % data to reject on (1: electrode data)
        'elecrange', [1:64], ... % channels (exclude external channels)
        'method', 'fft', ... % ['fft'|'multitaper'] method to compute spectrum.
        'threshold', [-100 25], ... % threshold limits (in dB)
        'freqlimits', [20 40], ... % frequency limits (in Hz)
        'eegplotplotallrej', 1, ... % do not superpose rejection marks on previously marks stored in the dataset
        'eegplotreject', 1); % rejection of marked trials (0: do not reject but store the marks; 1: reject)
        
    % for log file
    for irej = 1:numel(EEG.epoch)
        trigepochsleft{irej, :} = EEG.epoch(1, irej).eventtype; % create vector with the sum of epochs left for each condition (after artifact rejection)
    end
    if ~isempty(strfind(dispname, 'vpn08')) || ~isempty(strfind(dispname, 'vpn23')) % probably due to merging EEG datasets, the cell2mat function does not work on these participants... use cellfun
        trigepochsleft = cellfun(@str2num, trigepochsleft)'; % convert this vector from cell to number (otherwise the "sum" function won't work) and transpose
    else
        trigepochsleft = cell2mat(trigepochsleft)'; % convert this vector from cell to number (otherwise the "sum" function won't work) and transpose
    end
    for j = 1:numel(Avg.num)
        epochsleft(:, j) = sum(trigepochsleft == Avg.num(j)); % find all triggers belonging to each condition and sum them
        rejepxcond(:, j) = 100 - (epochsleft(:, j) * 100) / totepochs(:, j); % calculate percentage of rejected epochs separately for each condition
    end
    
    % save log file in a structure variable
    temp_logfile = struct('subject', dispname, ... % participant number
        'interp_chans', gen_bad_chans', ... % interpolated channels
        'perc_rejepochs', rejepxcond); % percentage of rejected epochs per condition
    logfile = [logfile; temp_logfile];
    
    EEG = pop_select(EEG, 'channel', 1:64); % eliminate external channels
    EEG = pop_reref(EEG, [], 'refstate', 0); % re-reference to average (the eegF_FASTER function re-references to Cz, hence the need for re-referencing here)
    
    %     pop_eegplot(EEG, 1, 1, 1); % check data
    
    % save preprocessed data
    % keep only trials with female identity
    EEG2 = EEG;
    EEG2 = pop_selectevent(EEG2, 'type', Avg.numfem);
    pop_saveset(EEG2, 'filepath', Avg.pathin_low_emotion_variability, 'filename', [dispname '_elist_be_artrej_lowVar.set']);
    
    % keep only trials with male identity
    EEG3 = EEG;
    EEG3 = pop_selectevent(EEG3, 'type', Avg.nummale);
    pop_saveset(EEG3, 'filepath', Avg.pathin_high_emotion_variability, 'filename', [dispname '_elist_be_artrej_highVar.set']);
    
    % stop timer
    temp_elapsed = struct('subject', dispname, ... % participant number
        'elapsed_time', toc); % interpolated channels
    elapsed = [elapsed, temp_elapsed];
    
    % clear variables before starting the preprocessing of a new participant
    clear EEG EEG2 EEG3 trigtotepochs totepochs trigepochsleft epochsleft rejepxcond temp_logfile temp_elapsed
    
end

% save log file
cd(pathanalysis) % set current directory
save exp2_log_preproc.mat logfile
save exp2_elapsed_time.mat elapsed

% create a .csv file with summary of number of interpolated channels and percentage of rejected epochs per condition
% (if you want to know which electrodes have been interpolated, check log_preproc.mat)
ssj = {}; temp = {}; temp_elec = [];
for itemp = 1:numel(logfile) % loop through participants
    ssj{itemp, 1} = logfile(itemp).subject; % participant number
    temp(itemp, :) = strsplit(num2str(logfile(itemp).perc_rejepochs)); % percentage of rejected epochs per condition
    temp_elec{itemp, 1} = num2str(numel(logfile(itemp).interp_chans)); % number of interpolated electrodes
end
summary = [ssj temp_elec temp]; % merge data
head_summary = {'ssj' 'interp_chans' 'upright_angry_lowVar' 'upright_neutral_lowVar' 'upright_irregular_lowVar' 'inverted_angry_lowVar' 'inverted_neutral_lowVar' 'inverted_irregular_lowVar' 'upright_angry_highVar' 'upright_neutral_highVar' 'upright_irregular_highVar' 'inverted_angry_highVar' 'inverted_neutral_highVar' 'inverted_irregular_highVar'}; % file header
summary = [head_summary; summary]; % merge header and data
cell2csv([pathanalysis expname '_logfile_interp_artrej.csv'], summary); % save as .csv

%%
