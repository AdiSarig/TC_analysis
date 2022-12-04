
cond_names = {'V_b', 'IV_b'};

for idx = 1:length(data.Masking.V_f)
    for ind = 1:length(cond_names)
        data_seg = data.Masking.(cond_names{ind})(idx).EEG;
        % Average trials
        cfg = [];
        data_avg.Masking.(cond_names{ind})(idx).EEG = ft_timelockanalysis(cfg, data_seg);
%         data_avg.Masking.(cond_names{ind})(idx).sub_num = sub_num;
    end
end

