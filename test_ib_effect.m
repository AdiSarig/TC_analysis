function test_ib_effect(params, data_avg)

% ib_test_latency = [0.2, 0.4]; % P7 visual inspection
ib_test_latency = [0.4, 0.6]; % p300
% ib_test_latency = [0.16, 0.19]; % van masking
% ib_test_latency = [0.24, 0.32]; % van dcf
% ib_test_channels = {'P7'};
ib_test_channels = {'C3', 'C1', 'Cz', 'C2', 'C4',...
               'CP3', 'CP1', 'CPz', 'CP2', 'CP4',...
               'P3', 'P1', 'Pz', 'P2', 'P4'};
% ib_test_channels = {'P5', 'P7', 'PO3', 'PO7',...
%                 'O1', 'POz', 'Oz', 'Iz', 'O2'...
%                 'P6', 'P8', 'PO4', 'PO8'};

%% Masking
cond_names = params.cond_names.masking_no_blank;

% select IB subjects
exclude_subjects = unique([params.exclude_subjects.masking, params.exclude_subjects.ib]);
existing_subs = [data_avg.Masking.V_f.sub_num];
ib_subs = [];
noticer_subs = [];
for sub = existing_subs
    if any(sub==exclude_subjects)
        continue
    elseif any(sub==params.exclude_subjects.ib_noticers)
        noticer_subs = [noticer_subs, find(sub == existing_subs)];
    else
        ib_subs = [ib_subs, find(sub == existing_subs)];
    end
end

% subtract masks
for ind = 1:length(cond_names)
    sub_count = 1;
    for sub = [ib_subs, noticer_subs]
        masking_subtracted_data(sub_count).sub = existing_subs(sub);
        if any(sub == ib_subs)
            masking_subtracted_data(sub_count).is_ib = true;
        else
            masking_subtracted_data(sub_count).is_ib = false;
        end
        % select the compatible mask for each condition
        switch cond_names{ind}(1)
            case 'V'
                current_blank = data_avg.Masking.V_b(sub).EEG;
            case 'I'
                current_blank = data_avg.Masking.IV_b(sub).EEG;
        end
        masking_subtracted_data(sub_count).(cond_names{ind}) = data_avg.Masking.(cond_names{ind})(sub).EEG.avg - current_blank.avg;
        sub_count = sub_count+1;
    end
end

% select time points
[~, ib_test_t_min] = min(abs(data_avg.Masking.V_f(1).EEG.time - ib_test_latency(1)));
[~, ib_test_t_max] = min(abs(data_avg.Masking.V_f(1).EEG.time - ib_test_latency(2)));

% select channels
for ind = 1:length(ib_test_channels)
    ib_test_channels_vec(ind) = find(strcmp(data_avg.Masking.V_f(1).EEG.label, ib_test_channels{ind}));
end

% save into table
for ind = 1:length(cond_names)
    for sub = 1:length(masking_subtracted_data)
        selected_avg = mean(masking_subtracted_data(sub).(cond_names{ind})(ib_test_channels_vec,ib_test_t_min:ib_test_t_max), 'all');
        masking_subtracted_data(sub).([cond_names{ind} '_avg']) = selected_avg;
    end
end

f_diff = num2cell([masking_subtracted_data.V_f_avg] - [masking_subtracted_data.IV_f_avg]);
[masking_subtracted_data.f_diff] = f_diff{:};
h_diff = num2cell([masking_subtracted_data.V_h_avg] - [masking_subtracted_data.IV_h_avg]);
[masking_subtracted_data.h_diff] = h_diff{:};

f_diff_ibs = [masking_subtracted_data([masking_subtracted_data.is_ib]).f_diff];
f_diff_noticers = [masking_subtracted_data(~[masking_subtracted_data.is_ib]).f_diff];
[H, p] = ttest2(f_diff_ibs,f_diff_noticers);
disp('Masking:');
if H
    disp('A difference was found between face visibility of IBs and noticers')
    disp(p)
else
    disp('No difference was found between face visibility of IBs and noticers')
    disp(p)
end

h_diff_ibs = [masking_subtracted_data([masking_subtracted_data.is_ib]).h_diff];
h_diff_noticers = [masking_subtracted_data(~[masking_subtracted_data.is_ib]).h_diff];
[H, p] = ttest2(h_diff_ibs,h_diff_noticers);
if H
    disp('A difference was found between house visibility of IBs and noticers');
    disp(p)
else
    disp('No difference was found between house visibility of IBs and noticers');
    disp(p)
end


%% DCF
% cond names, time points and channels are the same as in masking

% select IB subjects
exclude_subjects = unique([params.exclude_subjects.dcf, params.exclude_subjects.ib]);
existing_subs = [data_avg.color_fusion.V_f.sub_num];
ib_subs = [];
noticer_subs = [];
for sub = existing_subs
    if any(sub==exclude_subjects)
        continue
    elseif any(sub==params.exclude_subjects.ib_noticers)
        noticer_subs = [noticer_subs, find(sub == existing_subs)];
    else
        ib_subs = [ib_subs, find(sub == existing_subs)];
    end
end

% save into summary struct
for ind = 1:length(cond_names)
    for sub = 1:length([ib_subs noticer_subs])
        selected_avg = mean(data_avg.color_fusion.(cond_names{ind})(sub).EEG.avg(ib_test_channels_vec,ib_test_t_min:ib_test_t_max), 'all');
        summary_struct(sub).sub = data_avg.color_fusion.(cond_names{ind})(sub).sub_num;
        if any(sub == ib_subs)
            summary_struct(sub).is_ib = true;
        else
            summary_struct(sub).is_ib = false;
        end
        summary_struct(sub).([cond_names{ind} '_avg']) = selected_avg;
    end
end

f_diff = num2cell([summary_struct.V_f_avg] - [summary_struct.IV_f_avg]);
[summary_struct.f_diff] = f_diff{:};
h_diff = num2cell([summary_struct.V_h_avg] - [summary_struct.IV_h_avg]);
[summary_struct.h_diff] = h_diff{:};

f_diff_ibs = [summary_struct([summary_struct.is_ib]).f_diff];
f_diff_noticers = [summary_struct(~[summary_struct.is_ib]).f_diff];
[H, p] = ttest2(f_diff_ibs,f_diff_noticers);
disp('DCF:');
if H
    disp('A difference was found between face visibility of IBs and noticers')
    disp(p)
else
    disp('No difference was found between face visibility of IBs and noticers')
    disp(p)
end

h_diff_ibs = [summary_struct([summary_struct.is_ib]).h_diff];
h_diff_noticers = [summary_struct(~[summary_struct.is_ib]).h_diff];
[H, p] = ttest2(h_diff_ibs,h_diff_noticers);
if H
    disp('A difference was found between house visibility of IBs and noticers');
    disp(p)
else
    disp('No difference was found between house visibility of IBs and noticers');
    disp(p)
end

end