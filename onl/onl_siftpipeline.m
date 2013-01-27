function [signal cfg] = onl_siftpipeline(varargin)
% implements a basic SIFT pipeline (Preprocessing, Modeling, Connectivity)

% if ~exp_beginfun('filter'), return; end
% declare_properties('name','SIFT');

EEG = arg_extract(varargin,'EEG',struct([]));

cfg = arg_define([0 Inf],varargin, ...
        arg_norep({'EEG','Signal','signal'}), ...
        arg_sub({'preproc','Preprocessing'},{'EEG',EEG},@pre_prepData,'Pre-processing options'), ...
        arg_subswitch({'modeling','Modeling'},{'Segmentation VAR' 'EEG',EEG}, ...
            hlp_getModelingApproaches, 'Select a modeling approach and define parameters.'), ...
        arg_sub({'connectivity','Connectivity'},{'EEG',EEG,'MODEL',struct([])},@est_mvarConnectivity,'Select connectivity methods'), ...
        arg_subtoggle({'validation','Validation'},{'EEG',EEG},@est_validateMVAR,'Validate Model fit'), ...
        arg({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
        );


% Implement basic SIFT pipeline (Preprocessing, Modeling, Connectivity)
% -------------------------------------------------------------------------

% rereference to average
% cfg.EEG.data = bsxfun(@minus,cfg.EEG.data,mean(cfg.EEG.data)); end

% pre-process data
cfg.EEG = pre_prepData('EEG',cfg.EEG,cfg.preproc,'verb',cfg.verb);

% if strcmp(cfg.preproc.sigtype.arg_selection,'Sources')
%     % if using sources, construct the dipfit matrix
%     cfg.EEG.dipfit = hlp_microcache('dipfit',@hlp_makeDipfitStruct,hmObj.sourceSpace,cfg.EEG.roiVertices);
% end
                
% get the m-file name of the function implementing the modeling approach
modelingFuncName = hlp_getModelingApproaches('mfileNameOnly',cfg.modeling.arg_selection);

% fit model
cfg.EEG.CAT.MODEL = feval(modelingFuncName,'EEG',cfg.EEG,cfg.modeling,'verb',cfg.verb);

% calculate connectivity
cfg.EEG.CAT.Conn = est_mvarConnectivity('EEG',cfg.EEG,'MODEL',cfg.EEG.CAT.MODEL,cfg.connectivity,'verb',cfg.verb);

% perform model validation (optional)
if cfg.validation.arg_selection
    [whitestats PCstats stability] = est_validateMVAR('EEG',cfg.EEG,cfg.validation);

    if isempty(whitestats)
        whitenessCriteria = {};
    else
        whitenessCriteria = cfg.validation.checkWhiteness.whitenessCriteria;
    end
    vis_validation(whitestats,PCstats,stability,whitenessCriteria);
    
    % store results
    cfg.EEG.CAT.validation.whitestats = whitestats;
    cfg.EEG.CAT.validation.PCstats    = PCstats;
    cfg.EEG.CAT.validation.stability  = stability;
end

% return SIFT-augmented dataset
signal = cfg.EEG;

if nargout>1
    cfg = rmfield(cfg,'EEG');
end