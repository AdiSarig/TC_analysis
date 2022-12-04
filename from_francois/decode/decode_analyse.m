%========================================================
%
%   Decoding analyses
%   
%========================================================

% Need to run decode_prep.m before.

%% Set up

clear
close all

addpath ../
addpath ../functions/
addpath plots/
conf = getconfig();

band = 'smnhg';
regions = {'fusiform', 'spc', 'ifg', 'dlpfc'};

nsub = length(conf.subjects);

thresh = 5;

compute_gen = 0; % Switch to compute generalization or not

% if ~exist('figs','dir')
%     mkdir('figs');
% end
% short = 0;
% if short == 1
%     lw = 1;
%     figsz = [100,400,150,120];
% else
%     lw = 2;
%     figsz = [100,400,1000,850];
% end
ci = @(x) 1.96*std(x)./sqrt(size(x,1));

%% Compute generalization

for s = 1:nsub

    subj = conf.subjects{s};
    
    % Skip subjects that are not in conf.subjectsanalysis
    if all(strcmp(conf.subjectsanalysis,subj)==0)
        fprintf(' [ step 3.2 ] skipping %s\n',subj);
        continue
    else
        fprintf(' [ step 3.2 ] Processing %s\n',subj);
    end 
        
    if subj == "sub-guia"
        ses = 1:2;  % Add here subjects for which there's two sessions
    else
        ses = 1;
    end
    
    for ises = ses
        
        % Load data
        dir = [conf.dir.bidsproc '/analysis/' subj '/ses-' num2str(ises)];
        filename = [dir '/data_by_condition_resampled'];
        fprintf('Loading resampled data by condition for %s from %s\n',subj,dir);
        load(filename);
        
        for reg = 1:length(regions) 
            
            % Load electrodes
            elecfile = [conf.dir.bids '/' subj '/ses-' num2str(ises) '/ieeg/' subj '_space-mni_electrodes.tsv'];
            fprintf('|- Selecting electrodes from %s\n',elecfile);
            elecpostbl = readtable(elecfile,'FileType','text','Delimiter','\t');
            [elpos, sel] = getelec(elecpostbl,data.StaH.label,regions{reg},conf.sub(s).anat);
            fprintf('   |- Found %d electrodes for region *%s*\n',length(sel),regions{reg});
            
            if isempty(elpos)
                fprintf('Skipping region *%s* for %s because no electrodes were found.\n',regions{reg}, subj);
                continue
            elseif length(elpos) < thresh
                fprintf('Skipping region *%s* for %s because only %i electrodes were found.\n',regions{reg}, subj, length(elpos));
                continue
            end
                        
            fprintf([' [ analysis ] loading decoding results from ' regions{reg} ' for ' subj '.\n']);
            dec = load([conf.dir.bidsproc 'decode/decoder_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim','clas','clas_max','clas_max_all');
            dec_unl = load([conf.dir.bidsproc 'decode/decoder_unlocked_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim','clas');

            % get behavior and data excluding faceface and 1+ trials
            tr = string(data.PsyH_and_M.events.trialtype) ~= "faceface" & data.PsyH_and_M.events.resp ~= 2;
            ev = data.PsyH_and_M.events(tr,:);
            hg.tp_tl = data.PsyH_and_M.trial(sel,tr,:); % target-present - target-locked
            hg_max = squeeze(max(hg.tp_tl, [], 3));
            
            % Binarize confidence
            conf_binarized = nan(size(ev,1),1);
            for t = 1:size(ev,1)
                if (~contains(subj, conf.subjectslowconfdistrib) && ev.conf(t) == 3) ||... % for subjects that have more "1" responses than "3".
                        (contains(subj, conf.subjectslowconfdistrib) && ev.conf(t) ~= 1)
                    conf_binarized(t) = 1; % High conf
                else
                    conf_binarized(t) = 0; % Low conf
                end
            end

            % Get data size and labels
            [nel,ntr_tp,nt_tl] = size(hg.tp_tl);
            label.resp_tp = ev.resp;
            label.conf = conf_binarized;
            label.resp_ta = data.PsyFA_and_CR.events.resp;
            
            time_axis_tl = conf.win_disp(1):1/conf.resample_freq:conf.win_disp(2);
            time_axis_tl = time_axis_tl(1:end-1);
            
            % Get data for fixation cross-locked window, excluding the visual 
            % response to the mask.
            ind_fl = data.PsyH_and_M_unl.time{1} > conf.win_fix_locked(1);
            hg.tp_fl = data.PsyH_and_M_unl.trial(sel,tr,ind_fl); % target-present - fixation-locked
            hg.ta = data.PsyFA_and_CR.trial(sel,:,ind_fl); % target-absent - fixation-locked
            
            nt_fl = size(hg.tp_fl,3);
            ntr_ta = size(hg.ta,2);
            
            time_axis_fl = conf.win_fix_locked(1):1/conf.resample_freq:conf.win_fix_locked(2);
            time_axis_fl = time_axis_fl(1:end-1);
            
            % get data in the no report task
            if isfield(data, 'Str')
                hg.str = data.Str.trial(sel,:,:);
                ntr_str = size(hg.str,2);
            end
            
%             % Plot the high gamma response in the fixation-locked window.
%             % The goal is to see if their is an initial visual response to
%             % the onset of the masks (for some electrodes it lasts till 600ms).
%             
%             figure(349+reg)
%             nexttile()
%             plot(time_axis_fl,squeeze(mean(hg.tp_fl,2)))
%             ylabel('HGA (a.u.)');
%             xlabel("time (sec)");
%             title([subj ' (' num2str(nel) ' elecs, ' num2str(ntr_tp) ' trials)'])
%             sgtitle(['Fixation-locked activity in ' regions{reg}])

            % Set events identity
            id.high = ev.trialtype == "high";
            id.adapt = ev.trialtype == "adapt";
            id.hit = ev.resp == 1;
            id.miss = ev.resp == 0;
            id.hithigh = ev.trialtype == "high" & ev.resp == 1;
            id.hitadapt = ev.trialtype == "adapt" & ev.resp == 1;
            id.misshigh = ev.trialtype == "high" & ev.resp == 0;
            id.missadapt = ev.trialtype == "adapt" & ev.resp == 0;
            id.hc = conf_binarized == 1;
            id.lc = conf_binarized == 0;
            id.hithc = conf_binarized == 1 & ev.resp == 1;
            id.hitlc = conf_binarized == 0 & ev.resp == 1;
            id.misshc = conf_binarized == 1 & ev.resp == 0;
            id.misslc = conf_binarized == 0 & ev.resp == 0;
            id.fa = data.PsyFA_and_CR.events.resp == 1;
            id.crej = data.PsyFA_and_CR.events.resp == 0;
            if isfield(data, 'Str')
                id.strlow = data.Str.events.trialtype == "low";
                id.strhigh = data.Str.events.trialtype == "high";
            end
            
            % Load shuffled labels for permutations in target-absent trials
            fprintf(['Loading shuffled labels for ' subj ' in %s.\n'], [conf.dir.bidsproc 'permutations/shuffled_labels']);
            load([conf.dir.bidsproc 'permutations/shuffled_labels/' subj num2str(ises)], 'lab');
            
            shuf_id.fa = lab.psy_ta == 1;
            shuf_id.crej = lab.psy_ta == 0;
            
            %% Compute generalisation (decoder trained on one time point on another time point) for stim-locked epochs
            if compute_gen == 1
                generalization_tl = nan(ntr_tp,nt_tl,nt_tl);
%                 generalization_max = nan(ntr_tp,nt_fl);
                generalization_fl = nan(ntr_tp,nt_fl,nt_fl);
                gen_tl_from_max = nan(ntr_tp,nt_tl);
                gen_fl_from_max = nan(ntr_tp,nt_fl);
                gen_ta_from_max = nan(ntr_ta,nt_fl);
                if isfield(data, 'Str')
                    gen_str_from_max = nan(ntr_str,nt_tl);
                else
                    gen_str_from_max = [];
                end
                fprintf(' [ analysis ] computing decoder generalization\n');
                fprintf('    ');
                for xv = 1:dec.clas.cv.NumTestSets
                    idtest = dec.clas.cv.test(xv);
                    % For decoder trained on stim-locked epochs
                    fprintf('\b\b\b\b[%2d]',xv);
                    for k = 1:nt_tl
                        w = dec.clas.cl_xval{k,xv}.Coeffs(2,1).Linear;
                        b = dec.clas.cl_xval{k,xv}.Coeffs(2,1).Const;
                        f = squeeze(hg.tp_tl(:,idtest,:));
                        % so dimension 2 is "validation" and dimension 3 is "training"
                        generalization_tl(idtest,:,k) = reshape(w.'*f(:,:)+b,[nt_tl,sum(idtest)]).';
                    end
                    % For decoder trained on fixation-locked epochs
                    tmp = 1:size(dec_unl.tstim,2);
                    for k = tmp(ind_fl)
                        w = dec_unl.clas.cl_xval{k,xv}.Coeffs(2,1).Linear;
                        b = dec_unl.clas.cl_xval{k,xv}.Coeffs(2,1).Const;
                        f = squeeze(hg.tp_fl(:,idtest,:));
                        generalization_fl(idtest,:,k) = reshape(w.'*f(:,:)+b,[nt_fl,sum(idtest)]).';
                    end
                    
                    % For decoder trained on the max of each trial
                    idtest = dec.clas_max.cv.test(xv);
                    w_max = dec.clas_max.cl_xval{xv}.Coeffs(2,1).Linear;
                    b_max = dec.clas_max.cl_xval{xv}.Coeffs(2,1).Const;
%                     f = squeeze(hg_max(:,idtest));
%                     generalization_max(idtest,:) = w_max.'*f+b_max;
                    f_all = squeeze(hg.tp_tl(:,idtest,:));
                    gen_tl_from_max(idtest,:) = reshape(w_max.'*f_all+b_max,[nt_tl,sum(idtest)]).';
                    % test on fixation-locked window for different
                    % conditions
                    f_all = squeeze(hg.tp_fl(:,idtest,:));
                    gen_fl_from_max(idtest,:) = reshape(w_max.'*f_all+b_max,[nt_fl,sum(idtest)]).';
                end
                
                % Test decoder trained on max on target-absent trials and 
                % no report trials - this decoder was trained on all trials
                if iscell(dec.clas_max_all.cl_xval)
                    w_max = dec.clas_max_all.cl_xval{1}.Coeffs(2,1).Linear;
                    b_max = dec.clas_max_all.cl_xval{1}.Coeffs(2,1).Const;
                else
                    w_max = dec.clas_max_all.cl_xval.Coeffs(2,1).Linear;
                    b_max = dec.clas_max_all.cl_xval.Coeffs(2,1).Const;
                end
                if sum(label.resp_ta) > 0
                    for idtest = 1:ntr_ta
                        f = squeeze(hg.ta(:,idtest,:));
                        gen_ta_from_max(idtest,:) = reshape(w_max.'*f+b_max,[nt_fl,1]).';
                    end
                end
                if isfield(data, 'Str')
                    for idtest = 1:ntr_str
                        f = squeeze(hg.str(:,idtest,:));
                        gen_str_from_max(idtest,:) = reshape(w_max.'*f+b_max,[nt_tl,1]).';
                    end
                end

                % compute AUROC for all pairs of timepoints for high and
                % low intensity trials separately (stim-locked)
                fprintf(' [ analysis ] computing decoder generalization AUROC (stim-locked).\n');
                fprintf('     ');
                genaroc_tl = nan(nt_tl,nt_tl);
                genaroc_high_tl = nan(nt_tl,nt_tl);
                genaroc_adapt_tl= nan(nt_tl,nt_tl);
                for k1 = 1:nt_tl     
                    fprintf('\b\b\b\b\b[%3d]',k1);
                    for k2 = 1:nt_tl
                        [~,~,~,genaroc_tl(k1,k2)] = perfcurve(label.resp_tp,generalization_tl(:,k1,k2),1);
                        [~,~,~,genaroc_high_tl(k1,k2)] = perfcurve(label.resp_tp(id.high),generalization_tl(id.high,k1,k2),1);
                        [~,~,~,genaroc_adapt_tl(k1,k2)] = perfcurve(label.resp_tp(id.adapt),generalization_tl(id.adapt,k1,k2),1);
                    end
                end
                
                % compute AUROC based on max
                fprintf(' [ analysis ] Computing decoder generalization AUROC (stim-locked) based on trial max.\n');
                genaroc_max_del = nan(nt_tl,1);
                genaroc_max_str = nan(nt_tl,1);
                for k1 = 1:nt_tl
                    [~,~,~,genaroc_max_del(k1)] = perfcurve(label.resp_tp,gen_tl_from_max(:,k1),1);
                    if isfield(data, 'Str')
                        [~,~,~,genaroc_max_str(k1)] = perfcurve(id.strhigh,gen_str_from_max(:,k1),1);
                    end
                end
                
                % Compute unlocked AUROC (diagonal only)
                fprintf(' [ analysis ] computing decoder generalization AUROC (fixation-locked).\n');
                fprintf('     ');
                genaroc_fl = nan(nt_fl,1);
                for k1 = 1:nt_fl     
                    [~,~,~,genaroc_fl(k1)] = perfcurve(label.resp_tp,generalization_fl(:,k1,k1),1);
                end
                
                % Save 
                fprintf(' [ analysis ] saving in %s\n',[conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
                if isfield(data, 'Str')
                    save([conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
                        'genaroc_tl', 'genaroc_fl', 'genaroc_adapt_tl', 'genaroc_high_tl', 'generalization_tl',...
                        'generalization_fl', 'gen_fl_from_max',...
                        'gen_ta_from_max', 'gen_str_from_max', 'gen_tl_from_max'); 
                else
                    save([conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
                        'genaroc_tl', 'genaroc_fl', 'genaroc_adapt_tl', 'genaroc_high_tl', 'generalization_tl',...
                        'generalization_fl', 'gen_fl_from_max',...
                        'gen_ta_from_max', 'gen_tl_from_max'); 
                end

            else
                fprintf(' [ analysis ] Loading generalization from %s\n',[conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
                load([conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
                    'genaroc_tl', 'genaroc_fl', 'genaroc_adapt_tl', 'genaroc_high_tl', 'generalization_tl',...
                    'generalization_fl', 'gen_fl_from_max',...
                    'gen_ta_from_max', 'gen_str_from_max', 'gen_tl_from_max');
            end
                
            % Plot results
            figure(52+reg); 
            set(gcf, 'Position', get(0, 'Screensize'));
            gen_adapt.(regions{reg}){s} = gen_matrix(conf, genaroc_adapt_tl, "adapt");
            title([subj num2str(ises) ' (' num2str(nel) ' elecs, ' num2str(sum(id.adapt)) ' trials).'])
            sgtitle('Time generalization for adapt trials')
            print([conf.dir.bidsproc 'figs/decoding/time_gen_adapt_' regions{reg}], '-djpeg');
            
            figure(73+reg); 
            set(gcf, 'Position', get(0, 'Screensize'));
            gen_high.(regions{reg}){s} = gen_matrix(conf, genaroc_high_tl, "high");
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs, ' num2str(sum(id.high)) ' trials).'])
            sgtitle('Time generalization for high trials')
            print([conf.dir.bidsproc 'figs/decoding/time_gen_high_' regions{reg}], '-djpeg');
            
            % Plot results alongside decoding diagonal
            figure(22+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile()
            p1 = plot(time_axis_tl, diag(genaroc_high_tl), 'Color', conf.color.hit2); hold on
            p2 = plot(time_axis_tl, diag(genaroc_adapt_tl), 'Color', conf.color.hit1);
            ylabel("AUROC")
            xlabel("time (sec)");
            plot(xlim(),[0.5 0.5], 'k');
            legend([p1, p2], {'High intensity', 'Low intensity'})
            title([subj  num2str(ises)])
            sgtitle(['AUROC over time (diagonal) in ' regions{reg}]);
            print([conf.dir.bidsproc 'figs/decoding/time_gen_diag_' regions{reg}], '-djpeg');
            
%             figure(167+reg); 
%             nexttile(); hold on
%             p1 = plot(time_axis_tl, diag(genaroc_tl), 'Color', conf.color.hit);
%             ylabel('AUROC');
%             xlabel("time (sec)");
%             title([subj ' (' num2str(nel) ' elecs, ' num2str(ntr_tp) ' trials)'])
%             plot(xlim(),[0.5 0.5], 'k');
%             p2 = plot(xlim(),[genaroc_max genaroc_max], 'k--');
%             legend([p1, p2], {'Classifier over time', 'Classifier based on max'}, 'Location', 'southeast')
%             sgtitle(['Target-locked decoding in ' regions{reg}])

%             %% Look at whether the classifier outputs correlate with relative intensity
%             
%             tic
%             mn.hithigh = nan(nt_tl,1);
%             mn.hitadapt = nan(nt_tl,1);
%             mn.misshigh = nan(nt_tl,1);
%             mn.missadapt = nan(nt_tl,1);
%             ci.hithigh = nan(nt_tl,2);
%             ci.hitadapt = nan(nt_tl,2);
%             ci.misshigh = nan(nt_tl,2);
%             ci.missadapt = nan(nt_tl,2);
%             for it = 1:nt_tl-1
%                 mn.hithigh(it) = mean(generalization_tl(id.hithigh,it,it));
%                 mn.hitadapt(it) = mean(generalization_tl(id.hitadapt,it,it));
%                 mn.misshigh(it) = mean(generalization_tl(id.misshigh,it,it));
%                 mn.missadapt(it) = mean(generalization_tl(id.missadapt,it,it));
%                 ci.hithigh(it,:) = bootci(1000, @mean, generalization_tl(id.hithigh,it,it));
%                 ci.hitadapt(it,:) = bootci(1000, @mean, generalization_tl(id.hitadapt,it,it));
%                 ci.misshigh(it,:) = bootci(1000, @mean, generalization_tl(id.misshigh,it,it));
%                 ci.missadapt(it,:) = bootci(1000, @mean, generalization_tl(id.missadapt,it,it));
%             end
%             toc
%             
%             % Relative intensity
% %             intens_diff = nan(nt_tl,1);
% %             intens_diff_hit = nan(nt_tl,1);
% %             intens_diff_miss = nan(nt_tl,1);
% %             for it = 1:nt_tl
% %                 intens_diff(it) = mean(generalization_tl(id.high,it,it)) - mean(generalization_tl(id.adapt,it,it));
% %                 
% %                 intens_diff_hit(it) = mean(generalization_tl(id.hithigh,it,it)) - mean(generalization_tl(id.hitadapt,it,it));
% %                 intens_diff_miss(it) = mean(generalization_tl(id.misshigh,it,it)) - mean(generalization_tl(id.missadapt,it,it));
% %             end
%             
%             % Plot results alongside decoding diagonal
%             figure()
%             nexttile()
%             plot(time_axis_tl, diag(genaroc_tl)); hold on
%             title('AUROC over time (diagonal)');
%             xlabel("time (sec)");
%             plot(xlim(),[0.5 0.5], 'k');
% 
%             nexttile(); hold on
%             p1 = plot(time_axis_tl, mn.hitadapt, 'Color', conf.color.hit1);
%             p2 = plot(time_axis_tl, mn.hithigh, 'Color', conf.color.hit2);
%             plot(time_axis_tl, ci.hitadapt, 'LineStyle', '--', 'Color', conf.color.hit1)
%             plot(time_axis_tl, ci.hithigh, 'LineStyle', '--', 'Color', conf.color.hit2)
% %             a = area(time_axis_tl, [ci.hithigh(:,1) ci.hithigh(:,2)], 'EdgeAlpha', 0);
% %             a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
% %             a(2).FaceColor = conf.color.hit2; a(2).FaceAlpha = 0.1;
% %             a = area(time_axis_tl, [ci.hitadapt(:,1) ci.hitadapt(:,2)], 'EdgeAlpha', 0);
% %             a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
% %             a(2).FaceColor = conf.color.hit1; a(2).FaceAlpha = 0.1;
%             plot(xlim(),[0 0], 'k');
%             plot([0 0],ylim(), 'k');
%             title('Decoding predictions by intensity for hit trials');
%             ylabel("Predictor value");
%             xlabel("time (sec)");
%             legend([p1, p2], {'Adapt', 'High'})
%             
%             nexttile(); hold on
%             p1 = plot(time_axis_tl, mn.missadapt, 'Color', conf.color.miss2);
%             p2 = plot(time_axis_tl, mn.misshigh, 'Color', conf.color.miss1);
%             plot(time_axis_tl, ci.missadapt, 'LineStyle', '--', 'Color', conf.color.miss2)
%             plot(time_axis_tl, ci.misshigh, 'LineStyle', '--', 'Color', conf.color.miss1)
% %             a = area(time_axis_tl, [ci.misshigh(:,1) ci.misshigh(:,2)], 'EdgeAlpha', 0);
% %             a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
% %             a(2).FaceColor = conf.color.miss1; a(2).FaceAlpha = 0.1;
% %             a = area(time_axis_tl, [ci.missadapt(:,1) ci.missadapt(:,2)], 'EdgeAlpha', 0);
% %             a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
% %             a(2).FaceColor = conf.color.miss2; a(2).FaceAlpha = 0.1;
%             plot(xlim(),[0 0], 'k');
%             plot([0 0],ylim(), 'k');
%             title('Decoding predictions by intensity for miss trials');
%             ylabel("Predictor value");
%             xlabel("time (sec)");
%             legend([p1, p2], {'Adapt', 'High'})
%             
%             sgtitle([regions{reg} ' in ' subj ' (' num2str(nel) ' elecs, ' num2str(ntr_tp) ' trials).'])
            
%             %% Look at whether the classifier outputs correlate with relative intensity
%             
%             tic
%             mn.hithc = nan(nt_tl,1);
%             mn.hitlc = nan(nt_tl,1);
%             mn.misshc = nan(nt_tl,1);
%             mn.misslc = nan(nt_tl,1);
%             ci.hithc = nan(nt_tl,2);
%             ci.hitlc = nan(nt_tl,2);
%             ci.misshc = nan(nt_tl,2);
%             ci.misslc = nan(nt_tl,2);
%             for it = 1:nt_tl-1
%                 mn.hithc(it) = mean(generalization_tl(id.hithc,it,it));
%                 mn.hitlc(it) = mean(generalization_tl(id.hitlc,it,it));
%                 mn.misshc(it) = mean(generalization_tl(id.misshc,it,it));
%                 mn.misslc(it) = mean(generalization_tl(id.misslc,it,it));
%                 ci.hithc(it,:) = bootci(1000, @mean, generalization_tl(id.hithc,it,it));
%                 ci.hitlc(it,:) = bootci(1000, @mean, generalization_tl(id.hitlc,it,it));
%                 ci.misshc(it,:) = bootci(1000, @mean, generalization_tl(id.misshc,it,it));
%                 ci.misslc(it,:) = bootci(1000, @mean, generalization_tl(id.misslc,it,it));
%             end
%             toc
%             
%             % Confidence computations with Spearman correlation over three
%             % levels
%             rho = nan(1,nt_tl);
%             p_corr = nan(1,nt_tl);
%             
%             for it = 1:nt_tl
%                 [rho(it), p_corr(it)] = corr(ev.conf, ...
%                     generalization_tl(:,it,it), 'type', 'Spearman');
%             end
%             
%             % Plot results alongside decoding diagonal
%             figure()
%             nexttile()
%             plot(time_axis_tl, diag(genaroc_tl)); hold on
%             title('AUROC over time (diagonal)');
%             plot(xlim(),[0.5 0.5]);
%             
%             nexttile(); hold on
%             p1 = plot(time_axis_tl, mn.hitlc, 'Color', conf.color.hit1);
%             p2 = plot(time_axis_tl, mn.hithc, 'Color', conf.color.hit2);
%             plot(time_axis_tl, ci.hitlc, 'LineStyle', '--', 'Color', conf.color.hit1)
%             plot(time_axis_tl, ci.hithc, 'LineStyle', '--', 'Color', conf.color.hit2)
%             plot(xlim(),[0 0], 'k');
%             plot([0 0],ylim(), 'k');
%             title('Decoding predictions by confidence for hit trials');
%             ylabel("Predictor value");
%             xlabel("time (sec)");
%             legend([p1, p2], {'Low conf', 'High conf'})
%             
%             nexttile(); hold on
%             p1 = plot(time_axis_tl, mn.misslc, 'Color', conf.color.miss2);
%             p2 = plot(time_axis_tl, mn.misshc, 'Color', conf.color.miss1);
%             plot(time_axis_tl, ci.misslc, 'LineStyle', '--', 'Color', conf.color.miss2)
%             plot(time_axis_tl, ci.misshc, 'LineStyle', '--', 'Color', conf.color.miss1)
%             plot(xlim(),[0 0], 'k');
%             plot([0 0],ylim(), 'k');
%             title('Decoding predictions by confidence for miss trials');
%             ylabel("Predictor value");
%             xlabel("time (sec)");
%             legend([p1, p2], {'Low conf', 'High conf'})
%             
% %             nexttile()
% %             plot(time_axis_tl, rho); hold on
% %             plot(xlim(),[0 0]);
% %             if sum(p_corr<0.05) > 0
% %                 plot(time_axis_tl(p_corr<0.05), 0, 'k.');
% %             end
% %             title('Correlation with confidence (Spearman) over three levels');
%             
%             sgtitle([regions{reg} ' in ' subj ' (' num2str(nel) ' elecs, ' num2str(ntr_tp) ' trials).'])

            %% Decoding over time with the decoding trained on the max
            
            % Decoding for hit vs miss (target-locked)
            figure(161+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            p1 = plot(time_axis_tl, mean(gen_tl_from_max(id.hit,:)), 'Color', conf.color.hit);
            p2 = plot(time_axis_tl, mean(gen_tl_from_max(id.miss,:)), 'Color', conf.color.miss);
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.hit,:))'-ci(gen_tl_from_max(id.hit,:))' 2*ci(gen_tl_from_max(id.hit,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.hit; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.miss,:))'-ci(gen_tl_from_max(id.miss,:))' 2*ci(gen_tl_from_max(id.miss,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.miss; a(2).FaceAlpha = 0.1;
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.hit,:))'-ci(gen_tl_from_max(id.hit,:))'...
%                 mean(gen_tl_from_max(id.hit,:))'+ci(gen_tl_from_max(id.hit,:))'], 'Linestyle', '--', 'Color', conf.color.hit)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.miss,:))'-ci(gen_tl_from_max(id.miss,:))'...
%                 mean(gen_tl_from_max(id.miss,:))'+ci(gen_tl_from_max(id.miss,:))'], 'Linestyle', '--', 'Color', conf.color.miss)
            ylabel('w*hg')
            ylim([min(mean(gen_tl_from_max(id.miss,:)))-1 max(mean(gen_tl_from_max(id.hit,:)))+1])
            xlabel("Time from target onset (sec)")
            legend([p1 p2], {'Hit', 'Miss'}, 'Location', 'northwest')
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
            sgtitle(['Decoding (trained on max) hit vs miss for ' regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/train_max_test_HvsM_' regions{reg}], '-djpeg');

            figure(703+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            p1 = plot(time_axis_tl, mean(gen_tl_from_max(id.hithigh,:)), 'Color', conf.color.hit2);
            p2 = plot(time_axis_tl, mean(gen_tl_from_max(id.hitadapt,:)), 'Color', conf.color.hit1);
            p3 = plot(time_axis_tl, mean(gen_tl_from_max(id.misshigh,:)), 'Color', conf.color.miss1);
            p4 = plot(time_axis_tl, mean(gen_tl_from_max(id.missadapt,:)), 'Color', conf.color.miss2);
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.hithigh,:))'-ci(gen_tl_from_max(id.hithigh,:))' 2*ci(gen_tl_from_max(id.hithigh,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.hit2; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.hitadapt,:))'-ci(gen_tl_from_max(id.hitadapt,:))' 2*ci(gen_tl_from_max(id.hitadapt,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.hit1; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.misshigh,:))'-ci(gen_tl_from_max(id.misshigh,:))' 2*ci(gen_tl_from_max(id.misshigh,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.miss1; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.missadapt,:))'-ci(gen_tl_from_max(id.missadapt,:))' 2*ci(gen_tl_from_max(id.missadapt,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.miss2; a(2).FaceAlpha = 0.1;
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.hithigh,:))'-ci(gen_tl_from_max(id.hithigh,:))'...
%                 mean(gen_tl_from_max(id.hithigh,:))'+ci(gen_tl_from_max(id.hithigh,:))'], 'Linestyle', '--', 'Color', conf.color.hit2)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.hitadapt,:))'-ci(gen_tl_from_max(id.hitadapt,:))'...
%                 mean(gen_tl_from_max(id.hitadapt,:))'+ci(gen_tl_from_max(id.hitadapt,:))'], 'Linestyle', '--', 'Color', conf.color.hit1)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.misshigh,:))'-ci(gen_tl_from_max(id.misshigh,:))'...
%                 mean(gen_tl_from_max(id.misshigh,:))'+ci(gen_tl_from_max(id.misshigh,:))'], 'Linestyle', '--', 'Color', conf.color.miss1)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.missadapt,:))'-ci(gen_tl_from_max(id.missadapt,:))'...
%                 mean(gen_tl_from_max(id.missadapt,:))'+ci(gen_tl_from_max(id.missadapt,:))'], 'Linestyle', '--', 'Color', conf.color.miss2)
            ylabel('w*hg')
            xlabel("Time from target onset (sec)")
            ylim([min(mean(gen_tl_from_max(id.miss,:)))-1 max(mean(gen_tl_from_max(id.hit,:)))+1])
            legend([p1 p2 p3 p4], {'Hit high', 'Hit adapt', 'Miss high', 'Miss adapt'}, 'Location', 'northwest')
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
            sgtitle(['Decoding (trained on max) hit vs miss*intensity for ' regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/train_max_test_HvsMxintensity_' regions{reg}], '-djpeg');
            
            figure(222+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            p1 = plot(time_axis_tl, mean(gen_tl_from_max(id.hithc,:)), 'Color', conf.color.hit2);
            p2 = plot(time_axis_tl, mean(gen_tl_from_max(id.hitlc,:)), 'Color', conf.color.hit1);
            p3 = plot(time_axis_tl, mean(gen_tl_from_max(id.misshc,:)), 'Color', conf.color.miss1);
            p4 = plot(time_axis_tl, mean(gen_tl_from_max(id.misslc,:)), 'Color', conf.color.miss2);
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.hithc,:))'-ci(gen_tl_from_max(id.hithc,:))' 2*ci(gen_tl_from_max(id.hithc,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.hit2; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.hitlc,:))'-ci(gen_tl_from_max(id.hitlc,:))' 2*ci(gen_tl_from_max(id.hitlc,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.hit1; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.misshc,:))'-ci(gen_tl_from_max(id.misshc,:))' 2*ci(gen_tl_from_max(id.misshc,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.miss1; a(2).FaceAlpha = 0.1;
            a = area(time_axis_tl, [mean(gen_tl_from_max(id.misslc,:))'-ci(gen_tl_from_max(id.misslc,:))' 2*ci(gen_tl_from_max(id.misslc,:))'], 'EdgeAlpha', 0);
            a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
            a(2).FaceColor = conf.color.miss2; a(2).FaceAlpha = 0.1;
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.hithc,:))'-ci(gen_tl_from_max(id.hithc,:))'...
%                 mean(gen_tl_from_max(id.hithc,:))'+ci(gen_tl_from_max(id.hithc,:))'], 'Linestyle', '--', 'Color', conf.color.hit2)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.hitlc,:))'-ci(gen_tl_from_max(id.hitlc,:))'...
%                 mean(gen_tl_from_max(id.hitlc,:))'+ci(gen_tl_from_max(id.hitlc,:))'], 'Linestyle', '--', 'Color', conf.color.hit1)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.misshc,:))'-ci(gen_tl_from_max(id.misshc,:))'...
%                 mean(gen_tl_from_max(id.misshc,:))'+ci(gen_tl_from_max(id.misshc,:))'], 'Linestyle', '--', 'Color', conf.color.miss1)
%             plot(time_axis_tl, [mean(gen_tl_from_max(id.misslc,:))'-ci(gen_tl_from_max(id.misslc,:))'...
%                 mean(gen_tl_from_max(id.misslc,:))'+ci(gen_tl_from_max(id.misslc,:))'], 'Linestyle', '--', 'Color', conf.color.miss2)
            ylabel('w*hg')
            xlabel("Time from target onset (sec)")
            ylim([min(mean(gen_tl_from_max(id.miss,:)))-1 max(mean(gen_tl_from_max(id.hit,:)))+1])
            legend([p1 p2 p3 p4], {'Hit high conf', 'Hit low conf', 'Miss high conf', 'Miss low conf'}, 'Location', 'northwest')
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
            sgtitle(['Decoding (trained on max)  hit vs miss*conf for ' regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/train_max_test_HvsMxconf_' regions{reg}], '-djpeg');
            
            % Decoding in no report
            if isfield(data, 'Str')
                figure(832+reg)
                set(gcf, 'Position', get(0, 'Screensize'));
                nexttile(); hold on
                p1 = plot(time_axis_tl, mean(gen_str_from_max(id.strhigh,:)), 'Color', conf.color.crej2);
                p2 = plot(time_axis_tl, mean(gen_str_from_max(id.strlow,:)), 'Color', conf.color.crej1);
                a = area(time_axis_tl, [mean(gen_str_from_max(id.strhigh,:))'-ci(gen_str_from_max(id.strhigh,:))' 2*ci(gen_str_from_max(id.strhigh,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.crej2; a(2).FaceAlpha = 0.1;
                a = area(time_axis_tl, [mean(gen_str_from_max(id.strlow,:))'-ci(gen_str_from_max(id.strlow,:))' 2*ci(gen_str_from_max(id.strlow,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.crej1; a(2).FaceAlpha = 0.1;
%                 plot(time_axis_tl, [mean(gen_str_from_max(id.strhigh,:))'-ci(gen_str_from_max(id.strhigh,:))'...
%                     mean(gen_str_from_max(id.strhigh,:))'+ci(gen_str_from_max(id.strhigh,:))'], 'Linestyle', '--', 'Color', conf.color.crej2)
%                 plot(time_axis_tl, [mean(gen_str_from_max(id.strlow,:))'-ci(gen_str_from_max(id.strlow,:))'...
%                     mean(gen_str_from_max(id.strlow,:))'+ci(gen_str_from_max(id.strlow,:))'], 'Linestyle', '--', 'Color', conf.color.crej1)
                ylabel('w*hg')
                xlabel("Time from target onset (sec)")
                ylim([min(mean(gen_str_from_max(id.strlow,:)))-1 max(mean(gen_str_from_max(id.strhigh,:)))+1])
                legend([p1 p2], {'High', 'Low'}, 'Location', 'northwest')
                title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
                sgtitle(['Decoding (trained on max) in no report for ' regions{reg}])
                print([conf.dir.bidsproc 'figs/decoding/train_max_test_no_report_' regions{reg}], '-djpeg');
            end
            
            %% Figures for ASSC presentation
                    
            if subj == "sub-guia" && ises == 1
                
                y_scale = [0.9 4];
                
                % Test on delayed
                figure();  hold on
                set(gcf,'Position',[100 100 650 500])
                p1 = plot(time_axis_tl, mean(gen_tl_from_max(id.hit,:)), 'Color', conf.color.hit);
                p2 = plot(time_axis_tl, mean(gen_tl_from_max(id.miss,:)), 'Color', conf.color.miss);
                a = area(time_axis_tl, [mean(gen_tl_from_max(id.hit,:))'-ci(gen_tl_from_max(id.hit,:))' 2*ci(gen_tl_from_max(id.hit,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.hit; a(2).FaceAlpha = 0.1;
                a = area(time_axis_tl, [mean(gen_tl_from_max(id.miss,:))'-ci(gen_tl_from_max(id.miss,:))' 2*ci(gen_tl_from_max(id.miss,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.miss; a(2).FaceAlpha = 0.1;
                line([0 0], y_scale, 'Color', 'k', 'Linestyle', '--');
                ylim(y_scale)
                xlim([-0.25 1.2])
                ylabel('Weights * HGA', 'Fontsize', 24)
                xlabel("Time from face onset (sec)", 'Fontsize', 24)
                legend([p1 p2], {'Hit', 'Miss'}, 'Location', 'northwest', 'Fontsize', 18)
                ax = gca;
                ax.FontSize = 18;
%                 ax.BoxStyle='back';
%                 box off
                print([conf.dir.bidsproc 'figs/decoding/ASSC_' subj '_train_max_test_HvsM_' regions{reg}], '-djpeg');
                
                % Test on no report
                figure(); hold on
                set(gcf,'Position',[100 100 650 500])
                p1 = plot(time_axis_tl, mean(gen_str_from_max(id.strhigh,:)), 'Color', conf.color.crej2);
                p2 = plot(time_axis_tl, mean(gen_str_from_max(id.strlow,:)), 'Color', conf.color.crej1);
                a = area(time_axis_tl, [mean(gen_str_from_max(id.strhigh,:))'-ci(gen_str_from_max(id.strhigh,:))' 2*ci(gen_str_from_max(id.strhigh,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.crej2; a(2).FaceAlpha = 0.1;
                a = area(time_axis_tl, [mean(gen_str_from_max(id.strlow,:))'-ci(gen_str_from_max(id.strlow,:))' 2*ci(gen_str_from_max(id.strlow,:))'], 'EdgeAlpha', 0);
                a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
                a(2).FaceColor = conf.color.crej1; a(2).FaceAlpha = 0.1;
                line([0 0], y_scale, 'Color', 'k', 'Linestyle', '--');
                ylim(y_scale)
                xlim([-0.25 1.2])
                ylabel('Weights * HGA', 'Fontsize', 24)
                xlabel("Time from face onset (sec)", 'Fontsize', 24)
                legend([p1 p2], {'High', 'Low'}, 'Location', 'northwest', 'Fontsize', 18)
                ax = gca;
                ax.FontSize = 18;
                print([conf.dir.bidsproc 'figs/decoding/ASSC_' subj '_train_max_test_no_report_' regions{reg}], '-djpeg');
            end
            
            %% Find the max value for each trial when mutliplied by the weights of the max decoder
            
            % make a unique identifier accounting for subject's id + ses
            subjses = [subj(5:end) num2str(ises)];
            
            [max_hit.(subjses).(regions{reg}), ind_hit.(subjses).(regions{reg})] = max(gen_fl_from_max(id.hit,:),[],2);
            [max_miss.(subjses).(regions{reg}), ind_miss.(subjses).(regions{reg})] = max(gen_fl_from_max(id.miss,:),[],2);
            [max_hit_lc.(subjses).(regions{reg}), ind_hit_lc.(subjses).(regions{reg})] = max(gen_fl_from_max(id.hitlc,:),[],2);
            [max_hit_hc.(subjses).(regions{reg}), ind_hit_hc.(subjses).(regions{reg})] = max(gen_fl_from_max(id.hithc,:),[],2);
            [max_miss_lc.(subjses).(regions{reg}), ind_miss_lc.(subjses).(regions{reg})] = max(gen_fl_from_max(id.misslc,:),[],2);
            [max_miss_hc.(subjses).(regions{reg}), ind_miss_hc.(subjses).(regions{reg})] = max(gen_fl_from_max(id.misshc,:),[],2);
            if sum(label.resp_ta) > 0
                [max_fa.(subjses).(regions{reg}), ind_fa.(subjses).(regions{reg})] = max(gen_ta_from_max(id.fa,:),[],2);
                [max_crej.(subjses).(regions{reg}) ,ind_crej.(subjses).(regions{reg})] = max(gen_ta_from_max(id.crej,:),[],2);
            end
            if isfield(data, 'Str')
                [max_strlow.(subjses).(regions{reg}), ind_strlow.(subjses).(regions{reg})] = max(gen_str_from_max(id.strlow,:),[],2);
                [max_strhigh.(subjses).(regions{reg}), ind_strhigh.(subjses).(regions{reg})] = max(gen_str_from_max(id.strhigh,:),[],2);
            end
            
            % Some statistics. REPLACE TTEST WITH LM ACCOUNTING FOR STIM
            % INTENSITY
            [~,p_tp] = ttest2(max_hit.(subjses).(regions{reg}), max_miss.(subjses).(regions{reg}));
            [~,p_hit_conf] = ttest2(max_hit_lc.(subjses).(regions{reg}), max_hit_hc.(subjses).(regions{reg}));
            [~,p_miss_conf] = ttest2(max_miss_lc.(subjses).(regions{reg}), max_miss_hc.(subjses).(regions{reg}));
            if sum(label.resp_ta) > 0
                [~,p_ta] = ttest2(max_fa.(subjses).(regions{reg}), max_crej.(subjses).(regions{reg}));
            end
            if isfield(data, 'Str')
                [~,p_str] = ttest2(max_strlow.(subjses).(regions{reg}), max_strhigh.(subjses).(regions{reg}));
            end
            [rho_hit, p_corr_hit] = corr(ev.face_onset(id.hit), ind_hit.(subjses).(regions{reg}), 'type', 'Spearman');
            [rho_miss, p_corr_miss] = corr(ev.face_onset(id.miss), ind_miss.(subjses).(regions{reg}), 'type', 'Spearman');
            
            % Plot the results
            % Max decoding for hit vs miss (fixation-locked)
            figure(683+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            x_hit = ones(1,sum(id.hit)) + linspace(-0.2,0.2,sum(id.hit));
            x_miss = 2*ones(1,sum(id.miss)) + linspace(-0.2,0.2,sum(id.miss));
            scatter(x_hit, max_hit.(subjses).(regions{reg}), [], conf.color.hit, 'filled');
            scatter(x_miss, max_miss.(subjses).(regions{reg}), [], conf.color.miss, 'filled');
            errorbar([mean(max_hit.(subjses).(regions{reg})) mean(max_miss.(subjses).(regions{reg}))],...
                [ci(max_hit.(subjses).(regions{reg})) ci(max_miss.(subjses).(regions{reg}))], 'Color', 'k', 'LineWidth',2);
            ylabel('max(w*hg)')
            xlabel("Stimulus type")
            set(gca,'xtick',[1, 2])
            set(gca,'xticklabel',{'Hit', 'Miss'})
            offset = (max(max_hit.(subjses).(regions{reg}))-min(max_miss.(subjses).(regions{reg})))/10;
            axis([0.7 2.3 min(max_miss.(subjses).(regions{reg}))-offset max(max_hit.(subjses).(regions{reg}))+offset])
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
            if p_tp < 0.05
                line([1 2], [max(max_hit.(subjses).(regions{reg}))-offset max(max_hit.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([1 1], [max(max_hit.(subjses).(regions{reg}))-offset*2 max(max_hit.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([2 2], [max(max_hit.(subjses).(regions{reg}))-offset*2 max(max_hit.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
            end
            if p_tp < 0.001
                text(1.4, max(max_hit.(subjses).(regions{reg})), '***', 'Fontsize', 16);
            elseif p_tp < 0.01
                text(1.44, max(max_hit.(subjses).(regions{reg})), '**', 'Fontsize', 16);
            elseif p_tp < 0.05
                text(1.48, max(max_hit.(subjses).(regions{reg})), '*', 'Fontsize', 16);
            end
            sgtitle(['Max decoding per trial in fixation-locked window for ' regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/max_tp_' regions{reg}], '-djpeg');
            
            % Max decoding by conf
            figure(612+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            x_hit_lc = ones(1,sum(id.hitlc)) + linspace(-0.2,0.2,sum(id.hitlc));
            x_hit_hc = 2*ones(1,sum(id.hithc)) + linspace(-0.2,0.2,sum(id.hithc));
            x_miss_lc = 3*ones(1,sum(id.misslc)) + linspace(-0.2,0.2,sum(id.misslc));
            x_miss_hc = 4*ones(1,sum(id.misshc)) + linspace(-0.2,0.2,sum(id.misshc));
            scatter(x_hit_lc, max_hit_lc.(subjses).(regions{reg}), [], conf.color.hit1, 'filled');
            scatter(x_hit_hc, max_hit_hc.(subjses).(regions{reg}), [], conf.color.hit2, 'filled');
            scatter(x_miss_lc, max_miss_lc.(subjses).(regions{reg}), [], conf.color.miss2, 'filled');
            scatter(x_miss_hc, max_miss_hc.(subjses).(regions{reg}), [], conf.color.miss1, 'filled');
            errorbar([1,2], [mean(max_hit_lc.(subjses).(regions{reg})) mean(max_hit_hc.(subjses).(regions{reg}))],....
                [ci(max_hit_lc.(subjses).(regions{reg})) ci(max_hit_hc.(subjses).(regions{reg}))], 'Color', conf.color.hit, 'LineWidth',2);
            errorbar([3,4], [mean(max_miss_lc.(subjses).(regions{reg})) mean(max_miss_hc.(subjses).(regions{reg}))],...
                [ci(max_miss_lc.(subjses).(regions{reg})) ci(max_miss_hc.(subjses).(regions{reg}))], 'Color', conf.color.miss, 'LineWidth',2);
            ylabel('max(w*hg)')
            xlabel("Stimulus type")
            set(gca,'xtick',[1, 2, 3, 4])
            set(gca,'xticklabel',{'Hit low conf', 'Hit high conf', 'Miss low conf', 'Miss high conf'})
            offset = (max(max_hit_hc.(subjses).(regions{reg}))-min(max_miss_lc.(subjses).(regions{reg})))/10;
            axis([0.7 4.3 min(max_miss_lc.(subjses).(regions{reg}))-offset max(max_hit_hc.(subjses).(regions{reg}))+offset])
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
            if p_hit_conf < 0.05
                line([1 2], [max(max_hit_hc.(subjses).(regions{reg}))-offset max(max_hit_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([1 1], [max(max_hit_hc.(subjses).(regions{reg}))-offset*2 max(max_hit_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([2 2], [max(max_hit_hc.(subjses).(regions{reg}))-offset*2 max(max_hit_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
            end
            if p_miss_conf < 0.05
                line([3 4], [max(max_miss_hc.(subjses).(regions{reg}))-offset max(max_miss_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([4 4], [max(max_miss_hc.(subjses).(regions{reg}))-offset*2 max(max_miss_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                line([3 3], [max(max_miss_hc.(subjses).(regions{reg}))-offset*2 max(max_miss_hc.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
            end
            if p_hit_conf < 0.001
                text(1.4, max(max_hit_hc.(subjses).(regions{reg})), '***', 'Fontsize', 16);
            elseif p_hit_conf < 0.01
                text(1.44, max(max_hit_hc.(subjses).(regions{reg})), '**', 'Fontsize', 16);
            elseif p_hit_conf < 0.05
                text(1.48, max(max_hit_hc.(subjses).(regions{reg})), '*', 'Fontsize', 16);
            end
            if p_miss_conf < 0.001
                text(3.4, max(max_miss_hc.(subjses).(regions{reg})), '***', 'Fontsize', 16);
            elseif p_miss_conf < 0.01
                text(3.44, max(max_miss_hc.(subjses).(regions{reg})), '**', 'Fontsize', 16);
            elseif p_miss_conf < 0.05
                text(3.48, max(max_miss_hc.(subjses).(regions{reg})), '*', 'Fontsize', 16);
            end
            sgtitle(['Max decoding by confidence response for ' regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/max_conf_' regions{reg}], '-djpeg');
            
            % Max in hit trials plotted vs face onset
            figure(449+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            scatter(ev.face_onset(id.hit), time_axis_fl(ind_hit.(subjses).(regions{reg})), [], conf.color.hit, 'filled');
            xlabel("Time of face onset")
            ylabel("Time of max decoding")
            axis([0.75 2.1 conf.win_fix_locked(1)-0.1 conf.win_fix_locked(2)+0.1])
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs, ' num2str(sum(id.hit)) ' trials)'])
            text(1.8, 2, ['Rho = ' num2str(rho_hit)], 'Fontsize', 16);
            if p_corr_hit < 0.001
                text(1.4, conf.win_fix_locked(2), '***', 'Fontsize', 16);
            elseif p_corr_hit < 0.01
                text(1.4, conf.win_fix_locked(2), '**', 'Fontsize', 16);
            elseif p_corr_hit < 0.05
                text(1.4, conf.win_fix_locked(2), '*', 'Fontsize', 16);
            end
            sgtitle(["Max decoding as a function of face onset for hit trials in " regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/corr_hit_' regions{reg}], '-djpeg');
            
            % Max in miss trials plotted vs face onset
            figure(478+reg)
            set(gcf, 'Position', get(0, 'Screensize'));
            nexttile(); hold on
            scatter(ev.face_onset(id.miss), time_axis_fl(ind_miss.(subjses).(regions{reg})), [], conf.color.miss, 'filled');
            xlabel("Time of face onset")
            ylabel("Time of max decoding")
            axis([0.75 2.1 conf.win_fix_locked(1)-0.1 conf.win_fix_locked(2)+0.1])
            title([subj  num2str(ises) ' (' num2str(nel) ' elecs, ' num2str(sum(id.miss)) ' trials)'])
            text(1.8, 2, ['Rho = ' num2str(rho_miss)], 'Fontsize', 16);
            if p_corr_miss < 0.001
                text(1.4, conf.win_fix_locked(2), '***', 'Fontsize', 16);
            elseif p_corr_miss < 0.01
                text(1.4, conf.win_fix_locked(2), '**', 'Fontsize', 16);
            elseif p_corr_miss < 0.05
                text(1.4, conf.win_fix_locked(2), '*', 'Fontsize', 16);
            end
            sgtitle(["Max decoding as a function of face onset for miss trials in " regions{reg}])
            print([conf.dir.bidsproc 'figs/decoding/corr_miss_' regions{reg}], '-djpeg');
            
            % Max decoding for fa vs crej
            if sum(label.resp_ta) > 0
                figure(761+reg)
                set(gcf, 'Position', get(0, 'Screensize'));
                nexttile(); hold on
                x_fa = ones(1,sum(id.fa)) + linspace(-0.2,0.2,sum(id.fa));
                x_crej = 2*ones(1,sum(id.crej)) + linspace(-0.2,0.2,sum(id.crej));
                scatter(x_fa, max_fa.(subjses).(regions{reg}), [], conf.color.fa, 'filled');
                scatter(x_crej, max_crej.(subjses).(regions{reg}), [], conf.color.crej, 'filled');
                errorbar([1,2], [mean(max_fa.(subjses).(regions{reg})) mean(max_crej.(subjses).(regions{reg}))],...
                    [ci(max_fa.(subjses).(regions{reg})) ci(max_crej.(subjses).(regions{reg}))], 'Color', 'k', 'LineWidth',2);
                ylabel('max(w*hg)')
                xlabel("Stimulus type")
                set(gca,'xtick',[1, 2])
                set(gca,'xticklabel',{'FA', 'Crej'})
                offset = (max([max_crej.(subjses).(regions{reg}); max_fa.(subjses).(regions{reg})])-...
                    min([max_crej.(subjses).(regions{reg}); max_fa.(subjses).(regions{reg})]))/10;
                axis([0.7 2.3 min([max_crej.(subjses).(regions{reg}); max_fa.(subjses).(regions{reg})])-offset...
                    max([max_crej.(subjses).(regions{reg}); max_fa.(subjses).(regions{reg})])+offset])
                title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
                if p_ta < 0.05
                    line([1 2], [max(max_fa.(subjses).(regions{reg}))-offset max(max_fa.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                    line([1 1], [max(max_fa.(subjses).(regions{reg}))-offset*2 max(max_fa.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                    line([2 2], [max(max_fa.(subjses).(regions{reg}))-offset*2 max(max_fa.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                end
                if p_ta < 0.001
                    text(1.4, max(max_fa.(subjses).(regions{reg})), '***', 'Fontsize', 16);
                elseif p_ta < 0.01
                    text(1.44, max(max_fa.(subjses).(regions{reg})), '**', 'Fontsize', 16);
                elseif p_ta < 0.05
                    text(1.48, max(max_fa.(subjses).(regions{reg})), '*', 'Fontsize', 16);
                end
                sgtitle(['Max decoding in target-absent trials for ' regions{reg}])
                print([conf.dir.bidsproc 'figs/decoding/max_ta_' regions{reg}], '-djpeg');
            end
            
            % Max decoding for low vs high intensity in no report task
            if isfield(data, 'Str')
                figure(993+reg)
                set(gcf, 'Position', get(0, 'Screensize'));
                nexttile(); hold on
                x_low = ones(1,sum(id.strlow)) + linspace(-0.2,0.2,sum(id.strlow));
                x_high = 2*ones(1,sum(id.strhigh)) + linspace(-0.2,0.2,sum(id.strhigh));
                scatter(x_low, max_strlow.(subjses).(regions{reg}), [], conf.color.crej1, 'filled');
                scatter(x_high, max_strhigh.(subjses).(regions{reg}), [], conf.color.crej2, 'filled');
                errorbar([1,2], [mean(max_strlow.(subjses).(regions{reg})) mean(max_strhigh.(subjses).(regions{reg}))],...
                    [ci(max_strlow.(subjses).(regions{reg})) ci(max_strhigh.(subjses).(regions{reg}))], 'Color', 'k', 'LineWidth',2);
                ylabel('max(w*hg)')
                xlabel("Stimulus intensity")
                set(gca,'xtick',[1, 2])
                set(gca,'xticklabel',{'Low', 'High'})
                offset = (max([max_strhigh.(subjses).(regions{reg}); max_strlow.(subjses).(regions{reg})])-...
                    min([max_strhigh.(subjses).(regions{reg}); max_strlow.(subjses).(regions{reg})]))/10;
                axis([0.7 2.3 min([max_strhigh.(subjses).(regions{reg}); max_strlow.(subjses).(regions{reg})])-offset...
                    max([max_strhigh.(subjses).(regions{reg}); max_strlow.(subjses).(regions{reg})])+offset])
                title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
                if p_str < 0.05
                    line([1 2], [max(max_strhigh.(subjses).(regions{reg}))-offset max(max_strhigh.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                    line([1 1], [max(max_strhigh.(subjses).(regions{reg}))-offset*2 max(max_strhigh.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                    line([2 2], [max(max_strhigh.(subjses).(regions{reg}))-offset*2 max(max_strhigh.(subjses).(regions{reg}))-offset],'Color','black','LineWidth',0.8);
                end
                if p_str < 0.001
                    text(1.4, max(max_strhigh.(subjses).(regions{reg})), '***', 'Fontsize', 16);
                elseif p_str < 0.01
                    text(1.44, max(max_strhigh.(subjses).(regions{reg})), '**', 'Fontsize', 16);
                elseif p_str < 0.05
                    text(1.48, max(max_strhigh.(subjses).(regions{reg})), '*', 'Fontsize', 16);
                end
                sgtitle(['Max decoding in no report trials for ' regions{reg}])
                print([conf.dir.bidsproc 'figs/decoding/max_str_' regions{reg}], '-djpeg');
            end
            
            %% Decode fa vs crej for permutations and plot results compared 
            % to experimental mean.
            
            if sum(label.resp_ta) > 0
                n_perm = size(lab.psy_ta,2);
                for p = 1:n_perm
                    [perm_max_fa.(subjses).(regions{reg})(:,p), shuf_ind_fa.(subjses).(regions{reg})(:,p)] = max(gen_ta_from_max(shuf_id.fa(:,p),:),[],2);
                    [perm_max_crej.(subjses).(regions{reg})(:,p), shuf_ind_crej.(subjses).(regions{reg})(:,p)] = max(gen_ta_from_max(shuf_id.crej(:,p),:),[],2);                
                end
                
                exp_avg.ta = mean(max_fa.(subjses).(regions{reg}))- mean(max_crej.(subjses).(regions{reg}));
                perm_avg.ta = mean(perm_max_fa.(subjses).(regions{reg}))- mean(perm_max_crej.(subjses).(regions{reg}));
                perm_p.ta = sum(perm_avg.ta > exp_avg.ta)/n_perm; % Careful that this is a directional test!!!
                
                figure(590+reg); 
                set(gcf, 'Position', get(0, 'Screensize'));
                nexttile(); hold on
                histogram(perm_avg.ta);
                plot([exp_avg.ta exp_avg.ta], (ylim()), 'k--', 'Linewidth', 2);
                xtextpos = xlim(); ytextpos = ylim();
                text(xtextpos(1), ytextpos(2)*0.9,['P = ' num2str(perm_p.ta)]);
                ylabel("Permutation counts")
                xlabel("mean(fa) - mean(crej) for w*hg")
                title([subj  num2str(ises) ' (' num2str(nel) ' elecs)'])
                sgtitle(['Difference in max decoding in target-absent trials for ' regions{reg}])
                print([conf.dir.bidsproc 'figs/decoding/perm_max_ta_' regions{reg}], '-djpeg');
            end
        end
    end
end

% Save values from max decoding
fprintf(' [ analysis ] saving in %s\n',[conf.dir.bidsproc 'decode/decode_max.mat']);
save([conf.dir.bidsproc 'decode/decode_max.mat'], 'ind_crej', 'ind_fa', 'ind_hit',...
    'ind_hit_hc', 'ind_hit_lc', 'ind_miss', 'ind_miss_hc', 'ind_miss_lc', 'max_crej',...
    'max_fa', 'max_hit', 'max_hit_hc', 'max_hit_lc', 'max_miss', 'max_miss_hc',...
    'max_miss_lc', 'max_strlow', 'max_strhigh');


% %% Plot average region
% 
% for reg = 1:length(regions) 
%     imgs = gen_adapt.(regions{reg})(~cellfun(@isempty, gen_adapt.(regions{reg})));
%     tmp = cat(3,imgs{:});
%     img_avg = mean(tmp,3);
%     figure()
%     gen_matrix(conf, img_avg, ['Average in ' regions{reg} ' for adapt']); 
% end
