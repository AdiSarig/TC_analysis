function [data,data_avg] = load_data(params)
% load data for all paradigms

for name = params.paradigm_names
    switch name{:}
        case 'Inattentional_Blindness'
            cond_names = params.cond_names.ib;
            cond_code = params.cond_code.ib;
        case 'Masking'
            cond_names = params.cond_names.masking;
            cond_code = params.cond_code.masking;
        case 'color_fusion'
            cond_names = params.cond_names.dcf;
            cond_code = params.cond_code.dcf;
    end
    for ind = 1:length(cond_names)
        current_condition = sprintf('*_Baseline Correction_%s.vhdr', cond_names{ind});
        files = dir(fullfile('..',name{:},'Analysis','Analyzer','Export',current_condition));
        for idx = 1:length(files)
            % extract subject number
            sub_num = regexp(files(idx).name,'\d*','Match');
            sub_num = str2double(sub_num{1});
            if ~strcmp(files(idx).name(1:4),'Main') % change overlapping sub num
                sub_num = sub_num + 1000;
            end
            
            % Load segmented data as continueous data
            cfg = [];
            cfg.dataset = fullfile(files(idx).folder, files(idx).name);
            cfg.channel = {params.neighbours.label}';
            data_cont = ft_preprocessing(cfg);
            
            % Parse data into original segments
            cfg = [];
            cfg.dataset = fullfile(files(idx).folder, files(idx).name);
            cfg.trialdef.eventtype  = 'Stimulus';
            cfg.trialdef.eventvalue = cond_code.(cond_names{ind});
            cfg.trialdef.prestim    = params.prestim;
            cfg.trialdef.poststim   = params.poststim;
            cfg                     = ft_definetrial(cfg);
            
            % save segmented data into general structure
            data_seg = ft_redefinetrial(cfg, data_cont);
            data.(name{:}).(cond_names{ind})(idx).EEG = data_seg;
            
            % Average trials
            cfg = [];
            data_avg.(name{:}).(cond_names{ind})(idx).EEG = ft_timelockanalysis(cfg, data_seg);
            data_avg.(name{:}).(cond_names{ind})(idx).sub_num = sub_num;

            if params.save_to_csv
                % save data to csv
                cfg = [];
                cfg.toilim = [-0.2, 0.8];
                data_seg = ft_redefinetrial(cfg, data_seg); % first, use only timewindow of interest

                cfg = [];
                cfg.resamplefs = 200;
                data_seg = ft_resampledata(cfg, data_seg); % second, resample the data

                all_trials = cell2mat(data_seg.trial);
                num_trials = 1:length(data_seg.trial);
                num_trials_vec = repelem(num_trials, length(data_seg.time{1, 1}));
                all_trials_with_nums = [num_trials_vec; all_trials];
                file_name = sprintf('data_in_csv%c%s_%s_sub_%d.csv', filesep,name{:},cond_names{ind},sub_num);
                writematrix(all_trials_with_nums,file_name);
            end
            
            if params.debug && idx == 2
                break
            end
        end
    end   
end

% save('saved_data.mat','data','-v7.3');

end

