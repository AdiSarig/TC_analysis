function [stat] = permutation_erp(celldat)
% Run the permutaion test

% Define permutation parameters
cfg                  = [];
cfg.parameter        = 'avg';
cfg.channel          = 'all';
cfg.latency          = [0.14, 0.6]; % changed from 'all'
cfg.avgovertime      = 'no';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.clusteralpha     = 0.05; % thershold for point-by-point t-test
cfg.alpha            = 0.05; % divided by the number of hypotheses tested (2 main effects and one interaction)
cfg.numrandomization = 2000;
cfg.minnbchan        = 2;
cfg.clusterstatistic = 'maxsum';
cfg.layout           = 'standard_1020.elc';
cfg.correctm         = 'cluster';
cfg.correcttail      = 'prob';
cfg.tail             = 0;
cfg.uvar             = 1;
cfg.ivar             = 2;

sub_num         = length(celldat);
cfg.design=[1:sub_num 1:sub_num; ones(1,sub_num),2*ones(1,sub_num)];

load('biosemi64_neighb.mat');
cfg.neighbours=neighbours;

[stat] = ft_timelockstatistics(cfg,celldat{1,:},celldat{2,:});
% stat is a structure that contains the bootstrapping results.


end

