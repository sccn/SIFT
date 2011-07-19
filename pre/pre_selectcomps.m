
function [EEG g] = pre_selectcomps(varargin)
%
% Remove some components from icaact matrix and dipfit structure.
%
% Inputs:
%
%   EEG:        EEG data structure
%
% Optional:     <'Name',Value> pairs
%
%     VerbosityLevel:   Verbosity level. 0 = no output, 1 = text, 2 = graphical                                   
%                       Possible values: 0,1,2                                                                    
%                       Default value  : 0                                                                        
%                       Input Data Type: real number (double)                                                     
% 
%     ComponentsToKeep: Select components to analyze                                                              
%                       This should be a cell array of strings containing the IDs of components you wish to keep  
%                       Input Data Type: boolean    
% Outputs:
%
%   EEG:        processed EEG structure
%   g:          argument specification structure
%
% See Also: pre_prepData()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Section 6.5.1 
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
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


% Extract some defaults and ranges
MyComponentNames = [];
if ~isempty(varargin)
    EEG = arg_extract(varargin,'EEG',1);
    if ~isempty(EEG) && isfield(EEG,'icaact')
        % user called function with EEG argument
        if isfield(EEG.CAT,'curComponentNames') && ~isempty(EEG.CAT.curComponentNames)
            MyComponentNames = EEG.CAT.curComponentNames;
        else
            MyComponentNames = EEG.CAT.curComps;
            MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
        end
        clear EEG
    else
        % user did not provide EEG argument but passed in list of component names
        MyComponentNames = arg_extract(varargin,'MyComponentNames');
    end
end


% setup the argument structure
g = arg_define([0 1],varargin, ...
    arg_norep('EEG',mandatory), ...
    arg_nogui({'verb','VerbosityLevel'},0,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical'), ...
    arg({'ComponentsToKeep'},MyComponentNames,MyComponentNames,'Select components to analyze. This should be a cell array of strings containing the IDs of components you wish to keep','type','logical'), ...
    arg_norep('MyComponentNames',[],[],'A dummy argument needed to set range of ComponentToKeep','shape','row','type','cellstr') ...
    );


nbcomps = size(g.EEG.CAT.srcdata,1);
ComponentIndicesToKeep = find(ismember(MyComponentNames,g.ComponentsToKeep));

if g.verb, fprintf('Selecting components...\n'); end
       
% select a subset of components from g.EEG
rmcomps = setdiff(1:nbcomps,ComponentIndicesToKeep);
% g.EEG.icawinv(:,rmcomps)    = [];
% g.EEG.icaweights(rmcomps,:) = [];
g.EEG.CAT.srcdata(rmcomps,:)     = [];
g.EEG.CAT.srcdata = reshape(g.EEG.CAT.srcdata, size(g.EEG.CAT.srcdata,1), g.EEG.pnts, g.EEG.trials);
% g.EEG.icachansind           = 1:size(g.EEG.icaweights,1);

% recompute the ICA activations
% g.EEG = hlp_icaact(g.EEG);

% remove source structures of excluded components
% if isfield(g.EEG,'dipfit')
%     g.EEG.dipfit.model(rmcomps) = [];
% end

g.EEG.CAT.curComps            = ComponentIndicesToKeep;
g.EEG.CAT.curComponentNames   = g.ComponentsToKeep;

EEG = g.EEG;

if g.verb, disp('Done!'); end
