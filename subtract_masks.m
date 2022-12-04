function [masking_subtracted_data, problematic_subs] = subtract_masks(cond_names,sub_inds, data, data_avg, is_resample, use_iv_b)
% subtract mask related activity from each trial separately
% figure();
problematic_subs = [];
for ind = 1:length(cond_names)
    for sub = sub_inds
        num_trials = length(data.(cond_names{ind})(sub).EEG.trial);
        
        % select the compatible mask for each condition
        switch cond_names{ind}(1)
            case 'V'
                current_blank = data_avg.V_b(sub).EEG;
            case 'I'
                if use_iv_b
                    current_blank = data_avg.IV_b(sub).EEG;
                else
                    cfg = [];
                    cfg.offset    = -(0.075*data.(cond_names{ind})(sub).EEG.fsample); % shift in 75 ms times the sampling rate
                    current_blank = ft_redefinetrial(cfg, data_avg.V_b(sub).EEG);
%                     if strcmp(cond_names{ind}(4),'f')
%                         [~, t_min] = min(abs(data_avg.IV_b(sub).EEG.time - (-0.5)));
%                         [~, t_max] = min(abs(data_avg.IV_b(sub).EEG.time - 1));
%                         [~, t_min_shifted] = min(abs(current_blank.time - (-0.5)));
%                         [~, t_max_shifted] = min(abs(current_blank.time - 1));
% 
%                         deviation_from_true_b = current_blank.avg(:,t_min_shifted:t_max_shifted) - data_avg.IV_b(sub).EEG.avg(:,t_min:t_max);
%                         
%                         plot(current_blank.time(t_min_shifted:t_max_shifted),mean(deviation_from_true_b));
%                         hold on
%                         if any(mean(deviation_from_true_b) > 0.05)
%                             problematic_subs = [problematic_subs, data_avg.IV_f(sub).sub_num];
%                             fprintf('*********\nsub %d has large deviations between the shifted mask to the actual one\n***********\n', data_avg.IV_f(sub).sub_num);
%                         end
%                     end
                end
        end
%         cfg = [];
%         cfg.toilim    = [-0.5 1];
%         current_blank = ft_redefinetrial(cfg, current_blank);
%         current_data = ft_redefinetrial(cfg, data.(cond_names{ind})(sub).EEG);
        current_data = data.(cond_names{ind})(sub).EEG;

        if is_resample
            cfg = [];
            cfg.resamplefs = 200;
            current_blank = ft_resampledata(cfg, current_blank);
            current_data = ft_resampledata(cfg, current_data);
        end
        
        % subtract mask from each trial        
        trial_cell = cell(1, num_trials);
        for trial = 1:num_trials
            trial_cell{trial} = current_data.trial{trial} - current_blank.avg; 
        end
        
        current_data.trial = trial_cell;
        masking_subtracted_data(sub).(cond_names{ind}) = current_data;
    end
end


end

