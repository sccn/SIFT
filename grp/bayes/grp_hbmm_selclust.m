function [clustIndsToKeep clustids] = grp_hbmm_selclust(varargin)
% Select significant clusters following Multi-view Hierarchical Bayesian
% Learning. 
% 
% ------------------------------------------------------------------------------------------------------------------------
% Input              Information                                                                                          
% ------------------------------------------------------------------------------------------------------------------------
% Z:                 Cluster indicators returned by grp_hbmm_mcmc                                                         
%                    This is typically stored in MCMC_InitState.Z and for M subjects 
%                    is a [1 x M] cell array containing [M_i x M x NumIters] indicator matrices.                                  
%                    Input Range  : Unrestricted                                                                          
%                    Default value: MANDATORY INPUT                                                                       
%                    Input Data Type: string                                                                              
%                                                                                                                         
% SubjProbThresh:    Probability threshold for assigning subjects to clusters                                             
%                    A subject is assigned to a cluster if her probability of membership exceeds p_subj                   
%                    Input Range  : [0  1]                                                                                
%                    Default value: 0.5                                                                                   
%                    Input Data Type: real number (double)                                                                
%                                                                                                                         
% ClusterProbThresh: Probability threshold for keeping a cluster                                                          
%                    A cluster is "retained" if the proportion of subjects assigned to the cluster (with Prob > p_subj)   
%                    exceeds p_clust                                                                                      
%                    Input Range  : [0  1]                                                                                
%                    Default value: 0.5                                                                                   
%                    Input Data Type: real number (double) 
% 
% ------------------------------------------------------------------------------------------------------------------------
% Output              Information                                                                                          
% ------------------------------------------------------------------------------------------------------------------------
% clustIdsToKeep    List of clusters to keep (linear indices)
% 
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

arg_define([0 Inf],varargin, ...
    arg_nogui({'Z'},mandatory,[],'Cluster indicators returned by grp_hbmm_mcmc. This is typically stored in MCMC_InitState.Z and for M subjects (observations) is a [1 x M] cell array containing [M_i x M x NumIters] indicator matrices.'), ...
    arg({'p_subj','SubjProbThresh'},0.5,[0 1],'Probability threshold for assigning subjects to clusters. A subject is assigned to a cluster if her probability of membership exceeds p_subj'), ...
    arg({'p_clust','ClusterProbThresh'},0.5,[0 1],'Probability threshold for keeping a cluster. A cluster is "retained" if the proportion of subjects assigned to the cluster (with Prob > p_subj) exceeds p_clust'), ...
    arg({'verb','VerbosityLevel'},2,{int32(0) int32(1) int32(2)},'Verbosity level. 0 = no output, 1 = text, 2 = graphical') ...
    );

if verb
    fprintf('Selecting clusters for which %0.5g%% of subjects have a greater than %0.5g%% probability of membership\n',p_clust*100,p_subj*100);
end

nsubj  = length(Z);
nclust = size(Z{1},2);

% get the posterior mean of the Z-matrix for each subject
Z_mean=cellfun(@(Zi) mean(Zi,ndims(Zi)), Z,'UniformOutput',false);

% estimate p_M(k,i), the probability that subject i has membership in cluster k
p_M=zeros(nclust,nsubj); 
for i=1:nsubj
    for k=1:nclust
        p_M(k,i)=max(Z_mean{i}(:,k));
    end
end
% determine which clusters to keep
pr_k = sum(p_M > p_subj,2)/nsubj;
clustIndsToKeep = pr_k > p_clust;
% convert to linear indices
clustIndsToKeep = find(clustIndsToKeep);
clustIndsToKeep = clustIndsToKeep(:)';

for i=1:nsubj
    % store the ids of the retained clusters to 
    % which each source is most probably assigned
    Zm = Z_mean{i}(:,clustIndsToKeep);
    [~,clustids{i}] = max(Zm,[],2);
end

if verb
    fprintf('Keeping the following (%d) clusters out of (%d) total: \n\t%s\n',length(clustIndsToKeep),nclust,hlp_tostring(clustIndsToKeep));
end