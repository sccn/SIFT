function Model = hlp_sift_emptymodel(varargin)
% create a new (empty) Model datastructure with default fields initialized
% varargin is an optional set of <name,value> pairs indicating fields and 
% associated values to store into the set

if nargin > 1 && mod(length(varargin),2)
    error('SIFT:hlp_sift_emptymodel', ...
        ['Additional arguments must be a list of <name,value> pairs. ' ...
         'At least one name is missing its corresponding value.']);
end

Model.AR = {};
Model.PE = {};
Model.RC = {};
Model.mu = {};
Model.th = {};
Model.winStartTimes = []; %EEG.CAT.times(g.winStartIdx)/1000;
Model.morder        = [];
Model.winstep       = [];
Model.winlen        = [];
Model.algorithm     = '';
Model.modelclass    = '';
Model.timeelapsed   = [];
Model.normalize     = [];
Model.modelapproach = '';
Model.taperFcn      = '';

% append additional arguments
if nargin > 1
    Model = hlp_varargin2struct(varargin,Model);    
end


