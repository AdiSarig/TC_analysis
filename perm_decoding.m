function [stat] = perm_decoding(celldat)
% Run the permutaion test on decoded data

% Define permutation parameters
cfg                  = [];
cfg.parameter        = 'trial';
% cfg.channel          = 'all';
cfg.avgovertime      = 'no';
% cfg.avgoverchan      = 'yes';
cfg.method           = 'montecarlo';
cfg.statistic        = 'mvpa';
cfg.alpha            = 0.05/12; % divided by the number of hypotheses tested
cfg.numrandomization = 360; % 1000
% cfg.alpha            = 0.5; % divided by the number of hypotheses tested
% cfg.numrandomization = 3; % 1000
% cfg.minnbchan        = 2;
% cfg.layout           = 'standard_1020.elc';
cfg.correctm         = 'no';

%-------------------------------------------------------
% for using cluster based permutation
% load('biosemi64_neighb.mat');
% cfg.neighbours=neighbours;
cfg.neighbours = [];
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % thershold for point-by-point t-test
cfg.clusterstatistic = 'wcm';

cfg.clusterthreshold = 'nonparametric_individual'; 
% the statistics routines also allow you to use a
% nonparametric threshold for cluster-member candidates, based on the
% generated distribution of the test statistic under the null
% hypothesis. To use this, simply specify cfg.clusterthreshold =
% 'nonparametric_individual' or cfg.clusterthreshold =
% 'nonparametric_common'. The difference between the two is that the
% former computes a threshold per voxel, and the latter uses the same
% threshold for all voxels. Which one is appropriate for you I don't
% know. (Good reasons for using 'nonparametric_individual' might be a
% strong variation of your test statistic with frequency. I know for a
% fact this is the case with certain quantifications of phase-amplitude
% coupling; these show much higher values in the low frequencies even
% when computed on noise.)
%-----------------------------------------------------------------

cfg.correcttail      = 'prob';
cfg.tail             = 1;
cfg.uvar             = 2;
cfg.ivar             = 1;

sub_num         = size(celldat,2);
design1 = [];
design2 = [];
for sub = 1:sub_num
    n1 = length(celldat{1, sub}.trial);
    n2 = length(celldat{2, sub}.trial);
    current_design = [ones(1,n1); sub*ones(1, n1)];
    design1 = [design1 current_design];
    current_design = [2*ones(1,n2); sub*ones(1, n2)];
    design2 = [design2 current_design];
end
cfg.design = [design1, design2];

cfg.latency          = [-0.2, 0.8];
cfg.features         = 'chan';
cfg.mvpa.classifier  = 'svm';
cfg.mvpa.metric      = 'auc';
cfg.mvpa.k           = 5;
cfg.mvpa.repeat      = 1;


[stat] = ft_timelockstatistics(cfg,celldat{1,:},celldat{2,:});
% stat is a structure that contains the bootstrapping results.



end

