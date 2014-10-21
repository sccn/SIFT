function [GEEG clustIndsToKeep clustids] = grp_hbmm_selSigClust(varargin)
% Select significant clusters following Multi-view Hierarchical Bayesian
% Learning. 
% 
% ------------------------------------------------------------------------------------------------------------------------------
% Input                    Information                                                                                          
% ------------------------------------------------------------------------------------------------------------------------------
% GroupEEG:                EEG structure returned by grp_hbmm                                                                   
%                          This must contain the final MCMC state in GEEG.CAT.MCMC_LastState                                    
%                          Input Range  : Unrestricted                                                                          
%                          Default value: MANDATORY INPUT                                                                       
%                          Input Data Type: string                                                                              
%                                                  
% ConnMethod:              Which connmethod to prune                                                                            
%                          We can only do one at a time                                                                         
%                          Input Range  : Unrestricted                                                                          
%                          Default value: MANDATORY INPUT                                                                       
%                          Input Data Type: string    
%
% PruneDipoles:            Remove non-significant dipoles (if present)                                                          
%                          Input Range  : Unrestricted                                                                          
%                          Default value: 1                                                                                     
%                          Input Data Type: boolean                                                                             
%                                                                                                                               
% SignificanceThresholds:  Significance thresholds                                                                              
%                          Input Range  : Unrestricted                                                                          
%                          Default value: n/a                                                                                   
%                          Input Data Type: string                                                                                                                                                       
%                                                                                                                               
%     | SubjProbThresh:    Probability threshold for assigning subjects to clusters                                             
%                          A subject is assigned to a cluster if her probability of membership exceeds p_subj                   
%                          Input Range  : [0  1]                                                                                
%                          Default value: 0.5                                                                                   
%                          Input Data Type: real number (double)                                                                
%                                                                                                                               
%     | ClusterProbThresh: Probability threshold for keeping a cluster                                                          
%                          A cluster is "retained" if the proportion of subjects assigned to the cluster (with Prob > p_subj)   
%                          exceeds p_clust                                                                                      
%                          Input Range  : [0  1]                                                                                
%                          Default value: 0.5                                                                                   
%                          Input Data Type: real number (double)                                                                
%                                                                                                                               
%     | VerbosityLevel:    Verbosity level. 0 = no output, 1 = text, 2 = graphical                                              
%                          Possible values: 0,1,2                                                                               
%                          Default value  : 2                                                                                   
%                          Input Data Type: real number (double)    
% 
% ------------------------------------------------------------------------------------------------------------------------
% Output              Information                                                                                          
% ------------------------------------------------------------------------------------------------------------------------
% GEEG              Pruned dataset
% clustIdsToKeep    Indices of significant clusters
% clustids          A {N x 1} cell array of cluster assignments for N subjects where clustids{n}(k) is the cluster index for the kth dipole from the nth subject.
%
% See also: grp_hbmm(), grp_hbmm_mcmc()
%
% Author: Tim Mullen, 2014 SCCN/INC

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

g = arg_define([0 Inf],varargin, ...
        arg_norep({'GEEG','GroupEEG'},mandatory,[],'EEG structure returned by grp_hbmm. This must contain the final MCMC state in GEEG.CAT.MCMC_LastState','type','char'), ...
        arg({'cm','ConnMethod'},mandatory,[],'Which connmethod to prune. We can only do one at a time'), ...
        arg({'pruneDipoles','PruneDipoles'},false,[],'Remove non-significant dipoles (if present). Note that this may require adjustment of EEG.CAT.curComps to reference new indexing'), ...
        arg_sub({'sigthresh','SignificanceThresholds'},{},@grp_hbmm_selclust,'Significance thresholds') ...
        );

if ~isfield(g.GEEG.CAT,'MCMC_LastState') || ~isfield(g.GEEG.CAT.MCMC_LastState.(g.cm),'Z')
    error('SIFT:MissingZMatrix','Missing indicator matrix GEEG.CAT.MCMC_LastState.(%s).Z',g.cm);
end

% get the list of clusters to keep
[clustIndsToKeep clustids] = grp_hbmm_selclust(g.sigthresh,'Z',g.GEEG.CAT.MCMC_LastState.(g.cm).Z);
% prune the EEG dataset
g.GEEG.CAT = hlp_selectConnNodes(g.GEEG.CAT,clustIndsToKeep,g.cm);
% prune dipole array to keep only the significant centroids
if g.pruneDipoles && ~isempty(g.GEEG.dipfit)
    g.GEEG.dipfit.model = g.GEEG.dipfit.model(clustIndsToKeep);
    % FIME: adjust GEEG.CAT.curComps to account for newly pruned dipole indices
end

% return pruned structure
GEEG = g.GEEG;
