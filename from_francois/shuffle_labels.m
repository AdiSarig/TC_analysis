%==========================================================================
%
% Shuffle labels from behavior in order to run permutation analyses.
%
% François Stockart <francois.stockart@univ.grenobles_alpes.fr>
% 16 jun 2022
%
%========================================================================

%% Start up

% clear
addpath functions plots

%band = 'low'; sh = 40; b = 0;
band = 'smnhg'; sh = 1; b = -1;
region = 'wholebrain';

conf = getconfig();
addpath(conf.dir.fieldtrip);
ft_defaults;

nsub = length(conf.subjects);

% Specify the number of permutations that you want here!!!
n_perm = 1000;

%%

for s = 1:nsub
    
    % Get subject id
    subj = conf.subjects{s};
    
    % Skip subjects that are not in conf.subjectspreproc
    if all(strcmp(conf.subjectsanalysis,subj)==0)
        fprintf(' [ step 3.1 ] skipping %s\n',subj);
        continue
    else
        fprintf(' [ step 3.1 ] Processing %s\n',subj);
    end
    
    if subj == "sub-guia"
        ses = 1:2;  % Add here subjects for which there's two sessions
    else
        ses = 1;
    end
    
    for ises = ses
        
        %% Load behav
        stage1_dir = [conf.dir.bidsproc '/ieegproc/' subj '/ses-' num2str(ises) '/ieeg'];
        behav_filename = [stage1_dir(1:end-4) '/beh/' subj '_events.tsv'];
        events_ = readtable(behav_filename,'FileType','text','Delimiter','\t');
        sel = find(strcmp(events_.type,'start'))+3 ;
        events = events_(sel,:);
        
        %% Exclude bad trials

        cfg = [];
        cfg.trials = true(size(events,1),1);
        nev = size(events,1); % Number of events before trial exclusion

        % Identify trials with artifacts detected during visualization.
        if subj ~= "sub-guia"
            cfg.trials(conf.sub(s).badtrials) = 0; 
        else
            cfg.trials(conf.sub(s).badtrials.(strcat('ses', num2str(ises)))) = 0;
        end
        % Identify trials with too long or short RTs.
        cfg.trials(events.rt1 > conf.RT_immediate(2) & events.task == "staircase") = 0;
        cfg.trials(events.rt1 < conf.RT_immediate(1) & events.task == "staircase") = 0;
        cfg.trials(events.rt1 > conf.RT_delayed(2) & events.task == "psycho") = 0;
        cfg.trials(events.rt1 < conf.RT_delayed(1) & events.task == "psycho") = 0;

        % Exclude trials
        events = events(cfg.trials == 1,:);

        % Print number of bad trials
        sprintf("For subject %s, %i bad trials are discarded, accounting for %fpc of trials.\n", ...
            subj, nnz(~cfg.trials), nnz(~cfg.trials)/nev*100)
        
        %% Select events by task 
        ev.sta = events(events.task == "staircase",:);
        ev.sta_tp = ev.sta(ev.sta.trialtype ~= "catch",:);
        ev.psy = events(events.task == "psycho",:);
        ev.psy_tp = ev.psy(ev.psy.trialtype ~= "catch" & ev.psy.trialtype ~= "faceface" & ev.psy.resp ~= 2,:);
        ev.psy_ta = ev.psy(ev.psy.trialtype == "catch" & ev.psy.resp ~= 2,:);
        ev.str = events(events.task == "stream",:);
        
        %% Shuffle labels
        lab.sta_tp = nan(size(ev.sta_tp,1), n_perm);
        lab.psy_tp = nan(size(ev.psy_tp,1), n_perm);
        lab.psy_ta = nan(size(ev.psy_ta,1), n_perm);
        lab.str = nan(size(ev.str,1), n_perm);
        str_high = ev.str.trialtype == "high";
        rng('default');
        for p = 1:n_perm
            lab.sta_tp(:,p) = ev.sta_tp.resp(randperm(size(ev.sta_tp,1)));
            lab.psy_tp(:,p) = ev.psy_tp.resp(randperm(size(ev.psy_tp,1)));
            lab.psy_ta(:,p) = ev.psy_ta.resp(randperm(size(ev.psy_ta,1)));
            lab.str(:,p) = str_high(randperm(size(str_high,1)));
        end
        
        %% Save labels
        
        if ~exist([conf.dir.bidsproc 'permutations'], 'dir')
            mkdir([conf.dir.bidsproc 'permutations']);
        end
        if ~exist([conf.dir.bidsproc 'permutations/shuffled_labels'], 'dir')
            mkdir([conf.dir.bidsproc 'permutations/shuffled_labels']);
        end
        
        fprintf(['Saving shuffled labels for ' subj ' in %s.\n'], [conf.dir.bidsproc 'permutations/shuffled_labels']);
        save([conf.dir.bidsproc 'permutations/shuffled_labels/' subj num2str(ises)], 'lab');
        
    end
end