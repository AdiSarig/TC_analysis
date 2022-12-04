function stat = decoding(cond1,cond2)

% Configuration
n1 = numel(cond1.trial);
n2 = numel(cond2.trial);

cfg = [] ;  
cfg.method           = 'mvpa';
cfg.latency          = [-0.2, 0.8];
cfg.features         = 'chan';
cfg.mvpa.classifier  = 'lda';
cfg.mvpa.metric      = 'auc';
cfg.mvpa.k           = 10;
cfg.mvpa.repeat      = 5;
cfg.design           = [ones(n1,1); 2*ones(n2,1)];

stat = ft_timelockstatistics(cfg, cond1, cond2);

% mv_plot_result(stat.mvpa, stat.time)

end

