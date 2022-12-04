%========================================================
%
%   Classification for decoding
%   
%========================================================

clear

addpath ../
addpath ../functions/
conf = getconfig();

band = 'smnhg';
regions = {'fusiform', 'spc', 'ifg', 'dlpfc'};

% Decide to train the decoder on stim-locked or unlocked data;
epoch_type = "locked"; % otherwise: "unlocked"

% Decide to run on actual data or permutations
permut = 1; % actual data if 0, permutations if 1

% Threshold for the number of electrodes necessary to run models
thresh = 5;

nsub = length(conf.subjects);

for s = 6%1:nsub

    subj = conf.subjects{s};
    
    % Skip subjects that are not in conf.subjectspreproc
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
    
    for ises = 1%ses

        for reg = 3%1:length(regions) 
            
            % Load data
            dir = [conf.dir.bidsproc '/analysis/' subj '/ses-' num2str(ises)];
            filename = [dir '/data_by_condition_resampled'];
            fprintf('Loading resampled data by condition for %s from %s\n',subj,dir);
            load(filename);

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
            
            % Load shuffled labels
            if permut == 1
                fprintf('Loading shuffled labels for permutations for %s .\n',subj);
                load([conf.dir.bidsproc 'permutations/shuffled_labels/'  subj num2str(ises)], 'lab');
            end
            
            % get behavior and data excluding faceface and 1+ trials
            if epoch_type == "locked"
                tr = string(data.PsyH_and_M.events.trialtype) ~= "faceface" & data.PsyH_and_M.events.resp ~= 2;
                ev = data.PsyH_and_M.events(tr,:);
                hg = data.PsyH_and_M.trial(sel,tr,:);
                tstim = data.PsyH_and_M.time{1};
            elseif epoch_type == "unlocked"
                tr = string(data.PsyH_and_M_unl.events.trialtype) ~= "faceface" & data.PsyH_and_M_unl.events.resp ~= 2;
                ev = data.PsyH_and_M_unl.events(tr,:);
                hg = data.PsyH_and_M_unl.trial(sel,tr,1:end-1);
                tstim = data.PsyH_and_M_unl.time{1}(1:end-1);
            else
                warning("Unknown epoch type")
                break
            end
            
            % Get max over time of hg for each trial (over electrodes)
            max_hg = squeeze(max(hg, [], 3));

            % Get data size and label
            [nel,ntr,nt] = size(hg);
            if permut == 0
                label = ev.resp;
                nperm = 1;
            else
                label = lab.psy_tp;
                nperm = size(lab.psy_tp,2);
            end

            %% Train a decoder on each time point
            
            rng('default');
            
            if permut == 0
                % prealoc
                clas = [];
                clas.acc = nan(1,nt);
                clas.pred = nan(ntr,nt);
                clas.output = nan(ntr,nt);
                clas.cv = cvpartition(ntr,'LeaveOut');
                % clas.cv = cvpartition(ntr,'Kfold',conf.classif.nfoldext);
                fprintf('          ');

                tic
                for ti = 1:nt

                    % select features, exclude zero firing rate
                    feat = hg(:,:,ti)';

                    % for each cross-validation fold (LEAVE ONE OUT procedure)
                    for xv=1:ntr % conf.classif.nfoldext
                        % train decoder
                        idtrain = clas.cv.training(xv);
                        idtest = clas.cv.test(xv);
                        clas.cl_xval{ti,xv} = fitcdiscr(feat(idtrain,:),label(idtrain), ...
                            'FillCoeffs','on','Prior','uniform','ScoreTransform','logit', ...
                            'Gamma',conf.classif.gamma,'Delta',conf.classif.delta,'DiscrimType',conf.classif.discrim, ...
                            'HyperparameterOptimizationOptions',...
                            struct('ShowPlots',0,'Verbose',0,'UseParallel',1,'MaxObjectiveEvaluations',30,...
                            'KFold',conf.classif.nfoldint, 'AcquisitionFunctionName','expected-improvement-plus'), ...
                            'OptimizeHyperparameters',conf.classif.optimize);
                        % save weights and bias
                        w = clas.cl_xval{ti,xv}.Coeffs(2,1).Linear;
                        b = clas.cl_xval{ti,xv}.Coeffs(2,1).Const;

                        % predict test data
                        out = feat(idtest,:)*w+b;
                        pred = out > 0;
                        clas.output(clas.cv.test(xv),ti) = out;
                        clas.pred(clas.cv.test(xv),ti) = pred;

                    end

                    clas.acc(ti) = mean(clas.pred(:,ti) == label); %%% ATTENTION CLASSE NON-CONTREBALANCEE

                    fprintf('\b\b\b\b\b\b\b\b\b\b[%3d - %2.0f]',ti,100*clas.acc(ti));

                end
                fprintf('\n');
                toc
            end
            
            %% Train a decoder on maximum

            % prealoc
            clas_max = [];
            if permut == 0
                clas_max.acc = nan;
                clas_max.pred = nan(ntr,1);
                clas_max.output = nan(ntr,1);
            else
                clas_max.acc = nan(nperm,1);
                clas_max.pred = nan(ntr,nperm);
                clas_max.output = nan(ntr,nperm);
            end
            clas_max.cv = cvpartition(ntr,'LeaveOut');
            fprintf("Running decoder on maximum HGA value for each trial.\n")

            % select features, exclude zero firing rate
            feat = max_hg';

            % for each cross-validation fold (LEAVE ONE OUT procedure)
            tic
            for p = 1:nperm
                for xv=1:ntr % conf.classif.nfoldext
                    % train decoder
                    idtrain = clas_max.cv.training(xv);
                    idtest = clas_max.cv.test(xv);
                    clas_max.cl_xval{xv,p} = fitcdiscr(feat(idtrain,:),label(idtrain, p), ...
                        'FillCoeffs','on','Prior','uniform','ScoreTransform','logit', ...
                        'Gamma',conf.classif.gamma,'Delta',conf.classif.delta,'DiscrimType',conf.classif.discrim, ...
                        'HyperparameterOptimizationOptions',...
                        struct('ShowPlots',0,'Verbose',0,'UseParallel',1,'MaxObjectiveEvaluations',30,...
                        'KFold',conf.classif.nfoldint, 'AcquisitionFunctionName','expected-improvement-plus'), ...
                        'OptimizeHyperparameters',conf.classif.optimize);
                    % save weights and bias
                    w = clas_max.cl_xval{xv,p}.Coeffs(2,1).Linear;
                    b = clas_max.cl_xval{xv,p}.Coeffs(2,1).Const;

                    % predict test data
                    out = feat(idtest,:)*w+b;
                    pred = out > 0;
                    clas_max.output(clas_max.cv.test(xv),p) = out;
                    clas_max.pred(clas_max.cv.test(xv),p) = pred;

                    clas_max.acc(p) = mean(clas_max.pred(:,p) == label(:,p)); %%% ATTENTION CLASSE NON-CONTREBALANCEE

                end
            
                % additional classifier trained on all trials (for
                % generalization to target absent and no report trials)
                clas_max_all.cl_xval{p} = fitcdiscr(feat,label(:,p), ...
                    'FillCoeffs','on','Prior','uniform','ScoreTransform','logit', ...
                    'Gamma',conf.classif.gamma,'Delta',conf.classif.delta,'DiscrimType',conf.classif.discrim, ...
                    'HyperparameterOptimizationOptions',...
                    struct('ShowPlots',0,'Verbose',0,'UseParallel',1,'MaxObjectiveEvaluations',30,...
                    'KFold',conf.classif.nfoldint, 'AcquisitionFunctionName','expected-improvement-plus'), ...
                    'OptimizeHyperparameters',conf.classif.optimize);

                if permut == 1
                    fprintf("Permutation %i out of %i.\n", p, nperm)
                end
            
            end
            toc
            
            %% Save results
            
            confclassif = conf;
            if epoch_type == "locked"
                if permut == 0
                    fprintf(' [ analysis ] saving in %s\n',[conf.dir.bidsproc 'decode/decoder_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
                    save([conf.dir.bidsproc 'decode/decoder_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim','clas', 'clas_max', 'clas_max_all');%,'reghit','regmiss');
                else
                    fprintf(' [ analysis ] saving permutations in %s\n',[conf.dir.bidsproc 'permutations/decode/decoder_permut_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
                    save([conf.dir.bidsproc 'permutations/decode/decoder_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim', 'clas_max', 'clas_max_all','-v7.3');
                end
            elseif epoch_type == "unlocked"
                fprintf(' [ analysis ] saving in %s\n',[conf.dir.bidsproc 'decode/decoder_unlocked_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
                save([conf.dir.bidsproc 'decode/decoder_unlocked_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim','clas', 'clas_max', 'clas_max_all');
            end
        end
    end
end
