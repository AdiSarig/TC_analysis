function [perm_res] = run_permutation_one_paradigm(paradigm, sub_inds, cond_names, data, data_avg)

% Prepare data for main effects and interaction
visibility = cell(2, length(sub_inds));
stim_type  = cell(2, length(sub_inds));
visibility_within_stim_type = cell(2, length(sub_inds));
f_visibility = cell(2, length(sub_inds));
h_visibility = cell(2, length(sub_inds));

if strcmp(paradigm, 'masking')
    [masking_subtracted_data] = subtract_masks(cond_names,sub_inds, data, data_avg, true, true);
end

for ind = 1:length(sub_inds)
    switch paradigm
        case 'masking'
            data_v_f = masking_subtracted_data(sub_inds(ind)).V_f;
            data_v_h = masking_subtracted_data(sub_inds(ind)).V_h;
            data_iv_f = masking_subtracted_data(sub_inds(ind)).IV_f;
            data_iv_h = masking_subtracted_data(sub_inds(ind)).IV_h;
        case 'dcf'
            cfg = [];
            cfg.resamplefs = 200;
            data_v_f = ft_resampledata(cfg, data.V_f(sub_inds(ind)).EEG);
            data_v_h = ft_resampledata(cfg, data.V_h(sub_inds(ind)).EEG);
            data_iv_f = ft_resampledata(cfg, data.IV_f(sub_inds(ind)).EEG);
            data_iv_h = ft_resampledata(cfg, data.IV_h(sub_inds(ind)).EEG);
        case 'ib'
            cfg = [];
            cfg.resamplefs = 200;
            data_v_f = ft_resampledata(cfg, data.P2_f(sub_inds(ind)).EEG);
            data_v_h = ft_resampledata(cfg, data.P2_h(sub_inds(ind)).EEG);
            data_iv_f = ft_resampledata(cfg, data.P1_f(sub_inds(ind)).EEG);
            data_iv_h = ft_resampledata(cfg, data.P1_h(sub_inds(ind)).EEG);
    end
    
    % main effect - visibility
    cfg = [];
    merged_data = ft_appenddata(cfg, data_v_f, data_v_h);
    visibility{1,ind} = ft_timelockanalysis(cfg, merged_data);
    merged_data = ft_appenddata(cfg, data_iv_f, data_iv_h);
    visibility{2,ind} = ft_timelockanalysis(cfg, merged_data);
    
    % main effect - stimulus type
    merged_data = ft_appenddata(cfg, data_v_f, data_iv_f);
    stim_type{1,ind} = ft_timelockanalysis(cfg, merged_data);
    merged_data = ft_appenddata(cfg, data_v_h, data_iv_h);
    stim_type{2,ind} = ft_timelockanalysis(cfg, merged_data);
    
    % interaction
    cfg = []; % average trials
    data_v_f = ft_timelockanalysis(cfg, data_v_f);
    data_iv_f = ft_timelockanalysis(cfg, data_iv_f);
    data_v_h = ft_timelockanalysis(cfg, data_v_h);
    data_iv_h = ft_timelockanalysis(cfg, data_iv_h);
    
    cfg.parameter = 'avg';
    cfg.operation = 'subtract';
    visibility_within_stim_type{1, ind} = ft_math(cfg, data_v_f, data_iv_f);
    visibility_within_stim_type{2, ind} = ft_math(cfg, data_v_h, data_iv_h);
    
    % simple effects
    f_visibility{1, ind} = data_v_f;
    f_visibility{2, ind} = data_iv_f;
    h_visibility{1, ind} = data_v_h;
    h_visibility{2, ind} = data_iv_h;
end

% Run permutation analysis
perm_res.visibility_stat = permutation_erp(visibility);
perm_res.stim_type_stat  = permutation_erp(stim_type);
perm_res.visibility_within_stim_type_stat  = permutation_erp(visibility_within_stim_type);
perm_res.f_visibility = permutation_erp(f_visibility);
perm_res.h_visibility = permutation_erp(h_visibility);

end

