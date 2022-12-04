function [results, params, masking_subtracted_data] = get_lmm_data(params, data, data_avg)

%% Masking
cond_names = params.cond_names.masking_no_blank;

% select subjects
exclude_subjects = params.exclude_subjects.masking;
existing_subs = [data_avg.Masking.V_f.sub_num];
all_subs = ones(1, length(existing_subs));
for sub = exclude_subjects
    all_subs = all_subs & existing_subs ~= sub;
end
subs_masking = find(all_subs);
params.subs_masking = subs_masking;

[masking_subtracted_data] = subtract_masks(cond_names,subs_masking, data.Masking, data_avg.Masking, false, true);

% select time points
p3_latency = params.p3_latency;
[~, p3_t_min] = min(abs(masking_subtracted_data(1).V_f.time{1, 1} - p3_latency(1)));
[~, p3_t_max] = min(abs(masking_subtracted_data(1).V_f.time{1, 1} - p3_latency(2)));

van_latency = params.van_latency.masking;
[~, van_t_min] = min(abs(masking_subtracted_data(1).V_f.time{1, 1} - van_latency(1)));
[~, van_t_max] = min(abs(masking_subtracted_data(1).V_f.time{1, 1} - van_latency(2)));

% select channels
p3_channel_vec = params.p3_channel_vec;
van_channel_vec = params.van_channel_vec;

for ind = 1:length(params.p3_channels)
    p3_channel_vec(ind) = find(strcmp(data_avg.Masking.V_f(1).EEG.label, params.p3_channels{ind}));
end
for ind = 1:length(params.van_channels)
    van_channel_vec(ind) = find(strcmp(data_avg.Masking.V_f(1).EEG.label, params.van_channels{ind}));
end

[results_masking] = add_data_one_paradigm(data.Masking, data_avg.Masking, cond_names,...
    'masking', subs_masking, p3_t_min, p3_t_max, van_t_min, van_t_max,...
    p3_channel_vec, van_channel_vec, [], [], masking_subtracted_data);

%% DCF
% cond names, time points and channels are the same as in masking

% select subjects
exclude_subjects = params.exclude_subjects.dcf;
existing_subs = [data_avg.color_fusion.V_f.sub_num];
all_subs = ones(1, length(existing_subs));
for sub = exclude_subjects
    all_subs = all_subs & existing_subs ~= sub;
end
subs_dcf = find(all_subs);
params.subs_dcf = subs_dcf;

% select time points
p3_latency = params.p3_latency;
[~, p3_t_min] = min(abs(data_avg.color_fusion.V_f(1).EEG.time - p3_latency(1)));
[~, p3_t_max] = min(abs(data_avg.color_fusion.V_f(1).EEG.time - p3_latency(2)));

van_latency = params.van_latency.dcf;
[~, van_t_min] = min(abs(data_avg.color_fusion.V_f(1).EEG.time - van_latency(1)));
[~, van_t_max] = min(abs(data_avg.color_fusion.V_f(1).EEG.time - van_latency(2)));

results_dcf = add_data_one_paradigm(data.color_fusion, data_avg.color_fusion, cond_names,...
    'DCF', subs_dcf, p3_t_min, p3_t_max, van_t_min, van_t_max,...
    p3_channel_vec, van_channel_vec, [], [], []);

%% IB
% time points and channels are the same as in masking
cond_names = params.cond_names.ib_no_noise;
noise_cond = params.cond_names.ib_noise;

% select subjects (only those who didn't notice the stimuli in phase 1)
based_on_exclusion_criteria = params.exclude_subjects.ib;
noticers = params.exclude_subjects.ib_noticers;
exclude_subjects = [based_on_exclusion_criteria];
exclude_subjects = unique(exclude_subjects);
existing_subs = [data_avg.Inattentional_Blindness.P1_f.sub_num];
all_subs = ones(1, length(existing_subs));
for sub = exclude_subjects
    all_subs = all_subs & (existing_subs ~= sub);
end
subs_ib = find(all_subs);
params.subs_ib = subs_ib;

% select time points
p3_latency = params.p3_latency;
[~, p3_t_min] = min(abs(data_avg.Inattentional_Blindness.P1_f(1).EEG.time - p3_latency(1)));
[~, p3_t_max] = min(abs(data_avg.Inattentional_Blindness.P1_f(1).EEG.time - p3_latency(2)));

van_latency = params.van_latency.ib;
[~, van_t_min] = min(abs(data_avg.Inattentional_Blindness.P1_f(1).EEG.time - van_latency(1)));
[~, van_t_max] = min(abs(data_avg.Inattentional_Blindness.P1_f(1).EEG.time - van_latency(2)));

results_ib = add_data_one_paradigm(data.Inattentional_Blindness, data_avg.Inattentional_Blindness, cond_names,...
    'IB', subs_ib, p3_t_min, p3_t_max, van_t_min, van_t_max,...
    p3_channel_vec, van_channel_vec, noise_cond, noticers, []);


%% save results
results.p300 = [results_masking.p300; results_dcf.p300; results_ib.p300];
results.van = [results_masking.van; results_dcf.van; results_ib.van];

if ~params.debug
    if length(van_channel_vec) < 13 % masking post hoc analysis
        van_table_name = 'results_van_masking_posthoc.csv';
    else
        van_table_name = 'results_van.csv';
    end
    results_p300 = struct2table(results.p300);
    results_van  = struct2table(results.van);
    writetable(results_p300, 'results_p300.csv');
    writetable(results_van, van_table_name);
    
    IB_between_results_p300 = struct2table(results_ib.IB_between_p300);
    IB_between_results_van  = struct2table(results_ib.IB_between_VAN);
    writetable(IB_between_results_p300, 'IB_between_results_p300.csv');
    writetable(IB_between_results_van, 'IB_between_results_van.csv');
    
    noise_results_p300 = struct2table(results_ib.noise_p300);
    noise_results_van  = struct2table(results_ib.noise_van);
    writetable(noise_results_p300, 'noise_results_p300.csv');
    writetable(noise_results_van, 'noise_results_van.csv');

    noise_results_p300_for_lmm = struct2table(results_ib.noise_p300_for_lmm);
    noise_results_van_for_lmm  = struct2table(results_ib.noise_van_for_lmm);
    writetable(noise_results_p300_for_lmm, 'noise_results_p300_for_lmm.csv');
    writetable(noise_results_van_for_lmm, 'noise_results_van_for_lmm.csv');
end


end

