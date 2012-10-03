


function [m_norm surrogate_m] = est_PMGC(EEG,varargin)
% Calculate the event-related phase-modulated granger-causality estimates for a range of
% frequency bands


EEG = arg_extract(varargin,'EEG',1);
MyComponentNames = [];
if ~isempty(EEG)
    EEG = EEG(1);
    if isfield(EEG,'CAT') && isfield(EEG.CAT,'curComponentNames') && ~isempty(EEG.CAT.curComponentNames)
        MyComponentNames = EEG.CAT.curComponentNames;
    elseif ~isfield(EEG,'CAT')
        MyComponentNames = 1:size(EEG.icaact,1);
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    else
        MyComponentNames = EEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    clear EEG
end


g = arg_define([0 1], varargin, ...
        arg_norep({'EEG'},mandatory),...
        arg({'ampband','AmplitudePassBand'},[80 150],[],'The [lo hi] pass-band to use for amplitude (Hz)','shape','row'), ...
        arg({'phaseband','PhasePassBand'},[3 7],[],'The [lo hi] pass-band to use for phase (Hz)','shape','row'), ...
        arg({'components','Components'},true,MyComponentNames,'Components to use','type','logical'), ...
        arg({'numsurrogate','NumberOfSurrogateValues'},200,[],'Number of surrogate values for computing statistics') ...
    );




