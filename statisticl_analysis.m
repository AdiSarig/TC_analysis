% Statistical analysis - EEG data

% MUST have before running:
% field-trip toolbox
% mvpa-light toolbox

% ------Must modify before running--------
% The original field-trip toolbox was modified for this analysis. If the
% toolbox was updated or this analysis is running on another computer the
% toolbox must be modified again.
% Before line 364 in ft_statistics_montecarlo add:
% if strcmp(cfg.statistic, 'mvpa')
%     cfg.dim = length(statobs);
%     cfg.connectivity = true(1);
% end
%
% then should be the rest of the default function:
% [stat, cfg] = clusterstat(cfg, statrnd, statobs);
% ...
% ----------------------------------------

%% Parameters
load_masking = true;
load_IB = true;
load_DCF = true;
debug = false;
save_to_csv = true;
skip_permutation = true;
data_loaded = false;
run_masking_posthoc = false;

[params] = init_params(load_masking, load_IB, load_DCF, debug, save_to_csv);

%% Load data
if ~data_loaded
    tic
    [data,data_avg] = load_data(params);
    running_time.loading_data = toc;
end

%% Linear mixed model
% Prepare data for linear mixed model (later performed in JASP)
% 3(paradigm: IB/masking/DCF) X 2(visibility: visible/invisible) 
% X 2(stimulus type: face/house)
tic
[~, params] = get_lmm_data(params, data, data_avg);
running_time.lmm_data = toc;

% for post hoc masking analysis
if run_masking_posthoc
    params.van_channels = {'P5', 'P7', 'P9', 'PO7',...
                           'P6', 'P8', 'P10', 'PO8'};
    params.van_channel_vec = zeros(length(params.van_channels),1);
    get_lmm_data(params, data, data_avg);
end

% test the difference between IBs and noticers in DCF and masking
test_ib_effect(params, data_avg);

%% Permutation
% Masking/DCF/IB: 2(visiblity: visible/invisible) X 2(stimulus type: face/house)

% exclude also noticers
exclude_subjects = [params.exclude_subjects.ib, params.exclude_subjects.ib_noticers];
exclude_subjects = unique(exclude_subjects);
existing_subs = [data_avg.Inattentional_Blindness.P1_f.sub_num];
all_subs = ones(1, length(existing_subs));
for sub = exclude_subjects
    all_subs = all_subs & (existing_subs ~= sub);
end
params.subs_ib_no_noticers = find(all_subs);

if ~skip_permutation
    tic
    % Masking
    perm_res.masking = run_permutation_one_paradigm('masking', params.subs_masking, params.cond_names.masking_no_blank, data.Masking, data_avg.Masking);
    
    % DCF
    perm_res.dcf = run_permutation_one_paradigm('dcf', params.subs_dcf, params.cond_names.dcf_no_blank, data.color_fusion, data_avg.color_fusion);
    
    % IB
    perm_res.ib = run_permutation_one_paradigm('ib', params.subs_ib_no_noticers, params.cond_names.ib_no_noise, data.Inattentional_Blindness, data_avg.Inattentional_Blindness);
    
    running_time.permutation_analysis = toc;

    %% permutation results
    perm_cond = {'masking', 'dcf', 'ib'};
    effects = {'visibility_stat', 'stim_type_stat', 'visibility_within_stim_type_stat'};
    for cond = perm_cond
        for effect = effects
            if any(perm_res.(cond{:}).(effect{:}).mask)
                fprintf('An effect was found for paradigm %s, condition %s\n', cond{:}, effect{:});
            else
                fprintf('No effect was found for paradigm %s, condition %s\n', cond{:}, effect{:});
            end
        end
    end
end
%% Decoding

% Set generalization parameters
TC_subs = params.TC_subs;
existing_subs_masking = [data_avg.Masking.V_f.sub_num];
existing_subs_dcf = [data_avg.color_fusion.V_f.sub_num];
existing_subs_ib = [data_avg.Inattentional_Blindness.P1_f.sub_num];

[masking_subtracted_data, problematic_subs] = subtract_masks(params.cond_names.masking_no_blank,...
    params.subs_masking, data.Masking, data_avg.Masking, false, false);

%
tic
% Masking
[faces_vis_stat.masking, houses_vis_stat.masking] = run_decoding_one_paradigm('masking', masking_subtracted_data, params.subs_masking); % params.subs_masking

% DCF
[faces_vis_stat.dcf, houses_vis_stat.dcf] = run_decoding_one_paradigm('dcf', data.color_fusion, params.subs_dcf);

% IB
[faces_vis_stat.ib, houses_vis_stat.ib] = run_decoding_one_paradigm('ib', data.Inattentional_Blindness, params.subs_ib_no_noticers);

running_time.decoding_paradigms_separately = toc;

% train on two, test on one
tic
for ind = 1:length(TC_subs)
    sub_ind_masking = find(existing_subs_masking == TC_subs(ind));
    sub_ind_dcf = find(existing_subs_dcf == TC_subs(ind));
    sub_ind_ib = find(existing_subs_ib == TC_subs(ind));
    if params.debug
        sub_ind_masking = 1;
        sub_ind_dcf = 1;
        sub_ind_ib = 1;
    end

    cfg = [];
    cfg.toilim    = [-0.5 1];
    color_fusion.V_f = ft_redefinetrial(cfg, data.color_fusion.V_f(sub_ind_dcf).EEG);
    color_fusion.V_h = ft_redefinetrial(cfg, data.color_fusion.V_h(sub_ind_dcf).EEG);
    color_fusion.IV_f = ft_redefinetrial(cfg, data.color_fusion.IV_f(sub_ind_dcf).EEG);
    color_fusion.IV_h = ft_redefinetrial(cfg, data.color_fusion.IV_h(sub_ind_dcf).EEG);

    IB.V_f = ft_redefinetrial(cfg, data.Inattentional_Blindness.P2_f(sub_ind_ib).EEG);
    IB.V_h = ft_redefinetrial(cfg, data.Inattentional_Blindness.P2_h(sub_ind_ib).EEG);
    IB.IV_f = ft_redefinetrial(cfg, data.Inattentional_Blindness.P1_f(sub_ind_ib).EEG);
    IB.IV_h = ft_redefinetrial(cfg, data.Inattentional_Blindness.P1_h(sub_ind_ib).EEG);

    masking.V_f = ft_redefinetrial(cfg, masking_subtracted_data(sub_ind_masking).V_f);
    masking.V_h = ft_redefinetrial(cfg, masking_subtracted_data(sub_ind_masking).V_h);
    masking.IV_f = ft_redefinetrial(cfg, masking_subtracted_data(sub_ind_masking).IV_f);
    masking.IV_h = ft_redefinetrial(cfg, masking_subtracted_data(sub_ind_masking).IV_h);

    % train: masking+DCF, test: IB
    cfg = [];
    faces_vis.test_ib{1, ind} =  ft_appenddata(cfg, masking.V_f, color_fusion.V_f);
    houses_vis.test_ib{1, ind} =  ft_appenddata(cfg, masking.V_h, color_fusion.V_h);
    faces_vis.test_ib{2, ind} =  ft_appenddata(cfg, masking.IV_f, color_fusion.IV_f);
    houses_vis.test_ib{2, ind} =  ft_appenddata(cfg, masking.IV_h, color_fusion.IV_h);
    faces_vis.test_ib{3, ind} =  IB.V_f;
    houses_vis.test_ib{3, ind} =  IB.V_h;
    faces_vis.test_ib{4, ind} =  IB.IV_f;
    houses_vis.test_ib{4, ind} =  IB.IV_h;
       
    % train: masking+IB, test: DCF
    cfg = [];
    faces_vis.test_dcf{1, ind} =  ft_appenddata(cfg, masking.V_f, IB.V_f);
    houses_vis.test_dcf{1, ind} =  ft_appenddata(cfg, masking.V_h, IB.V_h);
    faces_vis.test_dcf{2, ind} =  ft_appenddata(cfg, masking.IV_f, IB.IV_f);
    houses_vis.test_dcf{2, ind} =  ft_appenddata(cfg, masking.IV_h, IB.IV_h);
    faces_vis.test_dcf{3, ind} =  color_fusion.V_f;
    houses_vis.test_dcf{3, ind} =  color_fusion.V_h;
    faces_vis.test_dcf{4, ind} =  color_fusion.IV_f;
    houses_vis.test_dcf{4, ind} =  color_fusion.IV_h;
    
    % train: DCF+IB, test: masking
    cfg = [];
    faces_vis.test_masking{1, ind} =  ft_appenddata(cfg, color_fusion.V_f, IB.V_f);
    houses_vis.test_masking{1, ind} =  ft_appenddata(cfg, color_fusion.V_h, IB.V_h);
    faces_vis.test_masking{2, ind} =  ft_appenddata(cfg, color_fusion.IV_f, IB.IV_f);
    houses_vis.test_masking{2, ind} =  ft_appenddata(cfg, color_fusion.IV_h, IB.IV_h);
    faces_vis.test_masking{3, ind} =  masking.V_f;
    houses_vis.test_masking{3, ind} =  masking.V_h;
    faces_vis.test_masking{4, ind} =  masking.IV_f;
    houses_vis.test_masking{4, ind} =  masking.IV_h;    
end

faces_vis_stat.test_ib = cross_decoding(faces_vis.test_ib);
houses_vis_stat.test_ib = cross_decoding(houses_vis.test_ib);
faces_vis_stat.test_dcf = cross_decoding(faces_vis.test_dcf);
houses_vis_stat.test_dcf = cross_decoding(houses_vis.test_dcf);
faces_vis_stat.test_masking = cross_decoding(faces_vis.test_masking);
houses_vis_stat.test_masking = cross_decoding(houses_vis.test_masking);

running_time.decoding_train_2_test_1 = toc;

% plot_decoding_results;

save('decoding_svm', 'faces_vis_stat', 'houses_vis_stat');

% print running time
running_sections = fieldnames(running_time);
for ind = 1:length(running_sections)
    current_duration = duration(0,0,running_time.(running_sections{ind}));
    fprintf('Running %s took %s (hh:mm:ss)\n', running_sections{ind}, current_duration);
end

%% temporal generalization