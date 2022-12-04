%========================================================
%
%   Decoding analyses with the permutation results
%   
%========================================================

% Need to run decode_prep.m and decode_analysis.m before.

%% Set up

clear
close all

addpath ../
addpath ../functions/
addpath plots/
conf = getconfig();

band = 'smnhg';
regions = {'ifg'};

% Compute generalization on permutations or load previously ran computations
compute_gen = 0;

% nsub = length(conf.subjects);
s=6;
subj = conf.subjects{s};
ises=1;
reg=1;

%% Load everything we need and run generalizations

if compute_gen == 1
    % Load HGA
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

    % Load generalizations on the actual data
    fprintf(' [ analysis ] Loading generalization from %s\n',[conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
    load([conf.dir.bidsproc 'decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
        'genaroc_tl', 'genaroc_fl', 'genaroc_adapt_tl', 'genaroc_high_tl', 'generalization_tl',...
        'generalization_fl', 'gen_fl_from_max',...
        'gen_ta_from_max', 'gen_str_from_max', 'gen_tl_from_max');

    % Load decoders trained on permutations
    fprintf(' [ analysis ] Loading permutations from %s\n',[conf.dir.bidsproc 'permutations/decode/decoder_permut_' subj '_ses' num2str(ises) '_' regions{reg} '.mat']);
    load([conf.dir.bidsproc 'permutations/decode/decoder_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],'confclassif','tstim', 'clas_max', 'clas_max_all');

    % Load permuted labels
    fprintf('Loading shuffled labels for permutations for %s .\n',subj);
    fprintf("This takes a long time as the file is big.\n")
    tic
    load([conf.dir.bidsproc 'permutations/shuffled_labels/'  subj num2str(ises)], 'lab');
    toc

    % Get actual labels
    tr = string(data.PsyH_and_M.events.trialtype) ~= "faceface" & data.PsyH_and_M.events.resp ~= 2;
    label.resp_tp = data.PsyH_and_M.events(tr,:).resp;
    label.str = data.Str.events.trialtype == "high";

    % Get required data
    hg.tp_tl = data.PsyH_and_M.trial(sel,tr,:);
    hg.str = data.Str.trial(sel,:,:);

    % Data size
    [nel,ntr_tp,nt_tl] = size(data.PsyH_and_M.trial(sel,tr,:));
    ntr_str = size(data.Str.trial(sel,:,:),2);
    nperm = size(lab.psy_tp,2);

    % Check that the number of trials in permutations and data match
    if ntr_tp ~= length(label.resp_tp) || ntr_str ~= length(label.str)
        warning("There's something fishy with the number of trials!!!")
    end
    
    % Compute AUROC on actual data (based on max)
    fprintf(' [ analysis ] Computing decoder generalization AUROC (stim-locked) based on trial max.\n');
    genaroc_max_del = nan(nt_tl,1);
    genaroc_max_str = nan(nt_tl,1);
    for k1 = 1:nt_tl
        [~,~,~,genaroc_max_del(k1)] = perfcurve(label.resp_tp,gen_tl_from_max(:,k1),1);
        [~,~,~,genaroc_max_str(k1)] = perfcurve(label.str,gen_str_from_max(:,k1),1);
    end

    % Run the generalization models on permuted data and compute AUROC
    tic
    gen_tl_from_max_permut = nan(ntr_tp,nt_tl);
    gen_str_from_max = nan(ntr_str,nt_tl);
    genaroc_max_del_permut = nan(nt_tl,nperm);
    genaroc_max_str_permut = nan(nt_tl,nperm);
    for p = 1:nperm
        for xv = 1:clas_max.cv.NumTestSets
            idtest = clas_max.cv.test(xv);
            w_max = clas_max.cl_xval{xv,p}.Coeffs(2,1).Linear;
            b_max = clas_max.cl_xval{xv,p}.Coeffs(2,1).Const;
            f_all = squeeze(hg.tp_tl(:,idtest,:));
            gen_tl_from_max_permut(idtest,:) = reshape(w_max.'*f_all+b_max,[nt_tl,sum(idtest)]).';
        end
        w_max = clas_max_all.cl_xval{p}.Coeffs(2,1).Linear;
        b_max = clas_max_all.cl_xval{p}.Coeffs(2,1).Const;
        for idtest = 1:ntr_str
            f = squeeze(hg.str(:,idtest,:));
            gen_str_from_max_permut(idtest,:) = reshape(w_max.'*f+b_max,[nt_tl,1]).';
        end
        for k1 = 1:nt_tl
            [~,~,~,genaroc_max_del_permut(k1,p)] = perfcurve(lab.psy_tp(:,p),gen_tl_from_max_permut(:,k1),1);
            [~,~,~,genaroc_max_str_permut(k1,p)] = perfcurve(lab.str(:,p),gen_str_from_max_permut(:,k1),1);
        end
        fprintf("Computed AUROC for permutation %i out of %i.\n", p, nperm)
    end
    toc
    
    save([conf.dir.bidsproc 'permutations/decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
        'genaroc_max_del_permut', 'genaroc_max_str_permut', 'genaroc_max_del',...
        'genaroc_max_str');
    
else
    load([conf.dir.bidsproc 'permutations/decode/genaroc_' subj '_ses' num2str(ises) '_' regions{reg} '.mat'],...
        'genaroc_max_del_permut', 'genaroc_max_str_permut', 'genaroc_max_del',...
        'genaroc_max_str');
end

%% Plots

% Define time axis
time_axis = conf.win_disp(1):1/conf.resample_freq:conf.win_disp(2);
time_axis = time_axis(1:end-1);

% Find thresholds and area to shade
p_thresh_del = quantile(genaroc_max_del_permut, 0.95,2);
shade_del = genaroc_max_del-p_thresh_del;
shade_del(shade_del<0) = 0;
p_thresh_str = quantile(genaroc_max_str_permut, 0.95,2);
shade_str = genaroc_max_str-p_thresh_str;
shade_str(shade_str<0) = 0;
y_scale = [min(min(genaroc_max_del), min(genaroc_max_str))-0.01,...
    max(max(genaroc_max_del), max(genaroc_max_str))+0.01];

% Figures
figure(1); clf
set(gcf,'Position',[100 100 650 500])
nexttile();  hold on
% plot(time_axis, genaroc_max_del_permut, 'k')
% plot(time_axis, quantile(genaroc_max_del_permut, 0.95,2), 'k')
plot(time_axis, genaroc_max_del, 'Color', conf.color.hit, 'LineWidth', 2)
a = area(time_axis, [p_thresh_del shade_del], 'EdgeAlpha', 0);
a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
a(2).FaceColor = conf.color.hit; a(2).FaceAlpha = 0.1;
line(xlim(), [0.5 0.5], 'Color', 'k', 'Linestyle', '-');
line([0 0], y_scale, 'Color', 'k', 'Linestyle', '--');
ylim(y_scale)
xlim([time_axis(1) 1.2])
ylabel('Area under ROC', 'Fontsize', 24)
xlabel('Time from face onset (sec)', 'Fontsize', 24)
ax = gca;
ax.FontSize = 18;
grid
print([conf.dir.bidsproc 'figs/decoding/' subj '_' regions{reg} '_AUROC_permut_del'], '-djpeg');
print([conf.dir.bidsproc 'figs/decoding/' subj '_' regions{reg} '_AUROC_permut_del'], '-dpdf');

figure(2); clf
set(gcf,'Position',[100 100 650 500])
nexttile(); hold on
% plot(genaroc_max_str_permut, 'k')
% plot(time_axis, quantile(genaroc_max_str_permut, 0.95,2), 'k')
plot(time_axis, genaroc_max_str, 'Color', conf.color.crej2, 'LineWidth', 2)
a = area(time_axis, [p_thresh_str shade_str], 'EdgeAlpha', 0);
a(1).FaceColor = 'w'; a(1).FaceAlpha = 0; 
a(2).FaceColor = conf.color.crej2; a(2).FaceAlpha = 0.1;
line(xlim(), [0.5 0.5], 'Color', 'k', 'Linestyle', '-');
line([0 0], y_scale, 'Color', 'k', 'Linestyle', '--');
ylim(y_scale)
xlim([time_axis(1) 1.2])
ylabel('Area under ROC', 'Fontsize', 24)
xlabel('Time from face onset (sec)', 'Fontsize', 24)
ax = gca;
ax.FontSize = 18;
grid
print([conf.dir.bidsproc 'figs/decoding/' subj '_' regions{reg} '_AUROC_permut_str'], '-djpeg');
print([conf.dir.bidsproc 'figs/decoding/' subj '_' regions{reg} '_AUROC_permut_str'], '-dpdf');
