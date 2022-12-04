clear
% close all

addpath ../
addpath ../functions/
addpath plots/
conf = getconfig();

band = 'smnhg';
regions = {'fusiform','ifg','dlpfc'};

nsub = length(conf.subjects);

% Load the data
fprintf(' [ analysis ] Loading %s\n',[conf.dir.bidsproc 'decode/decode_max.mat']);
load([conf.dir.bidsproc 'decode/decode_max.mat']);

% create tables for analysis
for reg = 1:length(regions)
    
    s_names = fieldnames(max_hit);
    
    % Prepare empty tables
    table_max_tp = [];
    table_max_conf = [];
    table_max_conf_hit = [];
    table_max_ta = [];
    table_max_str = [];
    
    for s = 1:length(s_names)
        
        subj = char(s_names(s));
        
        % Skip subjects without electrodes in current ROI.
        if ~isfield(max_hit.(subj), regions{reg})
            continue
        end
        
        % Put data in table
        len.hit = length(max_hit.(subj).(regions{reg}));
        len.miss = length(max_miss.(subj).(regions{reg}));
        table_max_tp = [table_max_tp; table(repmat(string(subj),len.hit+len.miss,1),...
            [ones(len.hit,1); zeros(len.miss,1)],...
            [max_hit.(subj).(regions{reg}); max_miss.(subj).(regions{reg})],...
                'VariableNames', {'sub' 'resp' 'max'})];
            
        len.hit_hc = length(max_hit_hc.(subj).(regions{reg}));
        len.hit_lc = length(max_hit_lc.(subj).(regions{reg}));
        len.miss_hc = length(max_miss_hc.(subj).(regions{reg}));
        len.miss_lc = length(max_miss_lc.(subj).(regions{reg}));
        table_max_conf = [table_max_conf; table(repmat(string(subj),len.hit+len.miss,1),...
            [ones(len.hit,1); zeros(len.miss,1)], [ones(len.hit_hc,1);...
            zeros(len.hit_lc,1); ones(len.miss_hc,1); zeros(len.miss_lc,1)],...
            [max_hit_hc.(subj).(regions{reg}); max_hit_lc.(subj).(regions{reg});...
            max_miss_hc.(subj).(regions{reg}); max_miss_lc.(subj).(regions{reg})],...
                'VariableNames', {'sub' 'resp' 'conf' 'max'})];
        table_max_conf_hit = [table_max_conf_hit; table(repmat(string(subj),len.hit,1),...
            [ones(len.hit_hc,1); zeros(len.hit_lc,1)],...
            [max_hit_hc.(subj).(regions{reg}); max_hit_lc.(subj).(regions{reg})],...
                'VariableNames', {'sub' 'conf' 'max'})];
            
        if isfield(max_fa, subj)
            len.fa = length(max_fa.(subj).(regions{reg}));
            len.crej = length(max_crej.(subj).(regions{reg}));
            table_max_ta = [table_max_ta; table(repmat(string(subj),len.fa+len.crej,1),...
                [ones(len.fa,1); zeros(len.crej,1)],...
                [max_fa.(subj).(regions{reg}); max_crej.(subj).(regions{reg})],...
                    'VariableNames', {'sub' 'resp' 'max'})];
        end
        
        if isfield(max_strlow, subj)
            len.strlow = length(max_strlow.(subj).(regions{reg}));
            len.strhigh = length(max_strhigh.(subj).(regions{reg}));
            table_max_str = [table_max_str; table(repmat(string(subj),len.strlow+len.strhigh,1),...
                [repmat("low",len.strlow,1); repmat("high",len.strhigh,1)],...
                [max_strlow.(subj).(regions{reg}); max_strhigh.(subj).(regions{reg})],...
                    'VariableNames', {'sub' 'intensity' 'max'})];
        end
            
    end
    
    % Run models
    lme_tp_rim.(regions{reg}) = fitlme(table_max_tp, 'max ~ resp + (1 | sub)');
    lme_tp_rem.(regions{reg}) = fitlme(table_max_tp, 'max ~ resp + (resp | sub)');
    lme_conf_rim.(regions{reg}) = fitlme(table_max_conf, 'max ~ resp*conf + (1 | sub)');
    lme_conf_rem.(regions{reg}) = fitlme(table_max_conf, 'max ~ resp*conf + (conf | sub)');
    lme_conf_hit_rim.(regions{reg}) = fitlme(table_max_conf_hit, 'max ~ conf + (1 | sub)');
    lme_conf_hit_rem.(regions{reg}) = fitlme(table_max_conf_hit, 'max ~ conf + (conf | sub)');
    lme_ta_rim.(regions{reg}) = fitlme(table_max_ta, 'max ~ resp + (1 | sub)');
    lme_ta_rem.(regions{reg}) = fitlme(table_max_ta, 'max ~ resp + (resp | sub)');
    lme_str_rim.(regions{reg}) = fitlme(table_max_str, 'max ~ intensity + (1 | sub)');
    lme_str_rem.(regions{reg}) = fitlme(table_max_str, 'max ~ intensity + (intensity | sub)');
end
