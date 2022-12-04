%========================================================
%
%   Random permutations for decoding
%   
%   Pereira et al., 2021, 
%   Evidence accumulation determines perceptual 
%       consciousness and monitoring
%      
%   Michael Pereira <michael.pereira@univ-grenoble-alpes.ch>
%   17 Dec 2020
%========================================================

% 
clear

addpath ../conf
addpath ../functions/
conf.exp = 2;
conf.session = 1;
config_spikes;

frate = 0.1;
t = linspace(-0.5,1.5,2001);

selstat = find(t>=conf.tstat(1) & t<conf.tstat(2));
selbaseline = find(t>=conf.tbaseline(1) & t<conf.tbaseline(2));

d = dir(['../conf/neurons/e2*_gauss' num2str(frate) '_stim.mat']);
[cumfeat,feat,neuron] = getfeatures(d,conf);

norm = mean(cumfeat(:,selbaseline,:),2);
cumfeat = bsxfun(@minus,cumfeat,norm);

id = getid(neuron.behav,conf.badtrials);

[nn,nt,ntr] = size(cumfeat);

nfold = 10;
cv = cvpartition(ntr,'Kfold',nfold);
rng('default');
%%
fprintf('     ');

label = id.idhit(id.idstim);
confidence = id.trueconf(id.idstim);
nt = length(conf.selt);

rndclas = [];
rndclas.acc = nan(1,nt,conf.nperm);
rndclas.pred = nan(length(label),nt,conf.nperm);
rndclas.output = nan(length(label),nt,conf.nperm);

saveconf = conf;

ntr = length(label);

rndclas.cv = cvpartition(length(label),'Kfold',conf.classif.nfoldext);

% loop across permuations
for rnd = 1:conf.nperm
    
    fprintf('\b\b\b\b\b[%3d]',rnd);
    
    % precompute permutation 
    perm = randperm(ntr);
    
    % loop across time points
    for ti=1:nt
        %%
        f = squeeze(cumfeat(:,conf.selt(ti),:));
        featsel = find(any(f.'~=0));
        f = (f(featsel,id.idstim)).';
        f = f(perm,:);
        % loop across cross-validation folds
        for xv=1:conf.classif.nfoldext
            % train classifier
            idtrain = rndclas.cv.training(xv);
            idtest = rndclas.cv.test(xv);
            rndcl_xval = fitcdiscr(f(idtrain,:),label(idtrain), ...
                'FillCoeffs','on','Prior','uniform','scoreTransform','logit', ...
                'Gamma',conf.classif.gamma,'Delta',conf.classif.delta,'DiscrimType',conf.classif.discrim, ...
                'HyperparameterOptimizationOptions',...
                struct('ShowPlots',0,'Verbose',0,'UseParallel',1,'MaxObjectiveEvaluations',30,...
                'KFold',conf.classif.nfoldint, 'AcquisitionFunctionName','expected-improvement-plus'), ...
                'OptimizeHyperparameters',conf.classif.optimize);
            % get weights
            w = rndcl_xval.Coeffs(2,1).Linear;
            b = rndcl_xval.Coeffs(2,1).Const;
            % predict test set
            out = f(idtest,:)*w+b;
            pred = out > 0;
            rndclas.output(rndclas.cv.test(xv),ti,rnd) = out;
            rndclas.pred(rndclas.cv.test(xv),ti,rnd) = pred;
        end
        rndclas.acc(ti,rnd) = mean(rndclas.pred(:,ti) == label);
      
        
    end
    if mod(rnd,100) == 0
        save('save/rnddecoder_exp2.mat','rndclas','saveconf');
    end
end
fprintf('\n');

save('save/rnddecoder_exp2.mat','rndclas','saveconf');

