

% examples for generating bootstrap distributions and statistics with stats toolbox
boot = bootstrp(10,@stat_bootfun,[1:EEG.CAT.trials],EEG,EEG.CAT.configs.est_fitMVAR);

% passing data directly through bootstrap fun
boot = bootstrp(1,@(tr,varargin) varargin{1},1:2,EEG.CAT.PConn);

% computing confidence intervals
bootconf = bootci(10,{@stat_bootfun,[1:EEG.CAT.trials],EEG,EEG.CAT.configs.est_fitMVAR},'type','student','nbootstd',5);
