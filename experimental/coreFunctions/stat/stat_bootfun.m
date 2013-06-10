function cmat = stat_bootfun(rsmpIdx,varargin)

EEG = varargin{1};
EEG.CAT.srcdata = EEG.CAT.srcdata(:,:,rsmpIdx);
EEG.CAT.trials = length(rsmpIdx);

args = varargin{2};

M = est_fitMVAR(EEG,args,'verb',0);

cmat = M.AR{1};
