function [EEGsim EEGtrue fwdModel] = sim_eegdata(varargin)
% simulate (realistic) scalp EEG data using a source-level dynamical
% system and a forward model
%
% This returns an EEG dataset or [nchs x npnts x ntr] matrix of simulated 
% source and scalp data. Source data is stored in EEG.srcpot. Scalp data is
% stored in EEG.data
%
% The function also can return the 'ground truth' in an EEGLAB dataset
%
% See Also: sim_fwdProj()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2013, SCCN/INC, UCSD.
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

HeadModelPath = [hlp_getSiftRoot filesep 'resources' filesep 'headmodels' filesep 'standard-Colin27-385ch.mat'];

arg_define(varargin, ...
    arg_subswitch({'srcdyn','SourceDynamics'},'VAR', ...
        {'VAR',@sim_varmodel},'Modeling approach for generating source dynamics. ''VAR'' simulates a Vector AutoRegressive system','cat','Source Data Generation'), ...
    arg_sub({'fwdproj','ForwardModel'},{'hmObj',HeadModelPath},@sim_fwdProj,'Set up the forward model. This projects source activations (see ''SourceDynamics'' option) through a forward head model to generate EEG channel data','cat','Scalp Data Generation'), ...
    arg({'makedipfit','MakeDipfitStruct'},true,[],'Create dipfit structure. Source centroids and ROI surfaces will be stored here'), ...
    arg_subtoggle({'vismodel','VisualizeModel'},'off',@vis_csd,'Visualize 3D model.'), ...
    arg({'verb','VerbosityLevel'},int32(2),{int32(0) int32(1) int32(2)},'Verbose','type','int32') ...
    );

% Error checking
if length(srcdyn.sim.expr) ...
    ~= max([length(fwdproj.sourceShape.roiAtlasLabels), ...
           length(fwdproj.sourceShape.roiOrdered),      ...
           size(fwdproj.sourceShape.sourceCoords,1)])
     error('You must have an equal number of sources in the dynamical model as source ROIs in the forward model');
end

% Simulate the source data
switch srcdyn.arg_selection
    case 'VAR'
        [EEGsim EEGtrue] = sim_varmodel(srcdyn,'verb',verb);
end

drawnow;

% Simulate scalp data 
[scalpData, srcData, LFM, chaninds, centroids_LFM_idx, roiVertices sourceCoords] = sim_fwdProj('sourceAmps',EEGsim.data,fwdproj);
hmObj = hlp_validateHeadModelObject(fwdproj.hmObj);

% ...prune channel space
hmObj.channelSpace = hmObj.channelSpace(chaninds,:);
hmObj.label        = hmObj.label(chaninds);

% Store data in EEG structures
chanlabels     = hmObj.getChannelLabels();
elocs          = hmObj.channelSpace;
EEGsim.srcpot     = EEGsim.data;
EEGsim.data       = scalpData;
EEGsim.nbchan     = size(scalpData,1);
EEGsim.srcpot_all = srcData;
% build chanlocs
EEGsim = rmfield(EEGsim,'chanlocs');
for k=1:EEGsim.nbchan
    EEGsim.chanlocs(k) = struct(...
        'X',elocs(k,1), ...
        'Y',elocs(k,2), ...
        'Z',elocs(k,3), ...
        'labels',chanlabels{k}, ...
        'type','EEGsim', ...
        'theta',[], ...
        'radius',0.5, ...
        'sph_theta',[],...
        'sph_phi',[],...
        'sph_radius',[],...
        'urchan',k, ...
        'ref','');
end
EEGsim.chanlocs = convertlocs(EEGsim.chanlocs,'cart2all');
% update source names
if isempty(srcdyn.makeEEGset.sourceNames)
    if ~isempty(fwdproj.sourceShape.roiOrdered)
        EEGsim.CAT.curComponentNames = fwdproj.sourceShape.roiOrdered;
    elseif ~isempty(fwdproj.sourceShape.roiAtlasLabels)
        EEGsim.CAT.curComponentNames = fwdproj.sourceShape.roiAtlasLabels;
    end
end
% create the dipfit struct
if makedipfit
    % load the original source space (non-reduced)
    tmp = load(hmObj.surfaces);
    fn  = fieldnames(tmp);
    sourceSpace = tmp.(fn{1})(3);  % dim3 = cortical surface

    brainStructsToRemove = setdiff_bc(unique_bc(hmObj.atlas.label),fwdproj.sourceAtlasLabels);
    reducedSpace = getSourceSpace4PEB(hmObj,brainStructsToRemove);

    EEGsim.dipfit = hlp_makeDipfitStruct(sourceSpace,roiVertices,reducedSpace,[],sourceCoords);
end

% Visualize the model
if vismodel.arg_selection
    gobj=vis_csd('signal',EEGsim,vismodel,'hmObj',hmObj);
end

% Prepare outputs
if nargout > 1
    if ~isempty(EEGtrue)
        tmp     = EEGtrue;
        EEGtrue = EEGsim;
        EEGtrue.CAT = tmp.CAT;
        EEGtrue.CAT.curComponentNames = EEGsim.CAT.curComponentNames;
        EEGtrue.setname   = tmp.setname;
        EEGtrue.condition = tmp.condition;
    end
end
if nargout > 2
    fwdModel.LFM = LFM;
    fwdModel.centroids_LFM_idx = centroids_LFM_idx;
end

