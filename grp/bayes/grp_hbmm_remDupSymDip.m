function [ALLEEG clustind] = grp_hbmm_remDupSymDip(varargin)
% Remove redundant dipoles (and connectivity) assigned to same cluster and remove one half of
% dual symmetric dipole pairs.
%
% Author: Tim Mullen, 2014, SCCN/INC

arg_define([0 Inf],varargin,...
    arg_norep({'ALLEEG','EEG'},mandatory,[],'Array of EEG datasets','type','expression'), ...
    arg({'clustind','ClusterAssignments'},[],[],'A [N x 1] cell array of cluster assignments. Here clustids{n}(k) is the cluster assignment (index) of the kth dipole from the nth subject. If empty, then ignore the assignments'), ...
    arg({'symDipDirToKeep','SymDipDirToKeep'},'left',{'left','right'},'Which hemisphere to keep for dual symmetric dipoles. This assumes that dipfit.model.posxyz(:,2)<0 specifies a left-hemisphere dipole while posxyz(:,2)>0 specifies right hemis.') ...
    );

% for subject i, for M_i sources
% first identify and strip duplicate dipoles
for si=1:length(ALLEEG)
    symDipRem  = [];
    dipmodel    = ALLEEG(si).dipfit.model;
    if ~isempty(clustind)
        % ensure each subject has no more than one dipole in each cluster 
        % check if more than one dipole are assigned to the same cluster
        % and, if so, keep only the first dipole
        [tmp,dipToKeep] = unique_bc(clustind{si},'first');
        [dipToKeep, ib] = sort(dipToKeep); % preserve ordering
        clustind{si} = tmp(ib);
        % delete excess dipoles
        dipmodel = dipmodel(dipToKeep);
        % store the indices (in original dataset) of removed dipoles
        ALLEEG(si).dipfit.dupDipRem = setdiff(1:length(dipmodel),dipToKeep);
        % remove associated connectivity (if present)
        try ALLEEG(si) = hlp_selectConnNodes(ALLEEG(si),dipToKeep); catch; end
    end

    for k=1:length(dipmodel)
        % if we have more than one dipole...
        if size(dipmodel(k).posxyz,1) > 1
            % ...check if this is a valid duplicate
            if length(dipmodel(k).select) > 1 || ~all(dipmodel(k).posxyz(2,:)==0)
                % if so, update list of duplicates
                symDipRem = [symDipRem k];
            end
            
            if all(dipmodel(k).posxyz(2,:)==0)
                dir = 1;
            else
                % ...keep only one of the dipoles
                switch symDipDirToKeep
                    case 'right'
                        dir = find(dipmodel(k).posxyz(:,2)>=0);
                    case 'left'
                        dir = find(dipmodel(k).posxyz(:,2)<0);
                end
            end
            if isempty(dir), dir = 1; end
            dipmodel(k).posxyz = dipmodel(k).posxyz(dir,:);
            dipmodel(k).momxyz = dipmodel(k).momxyz(dir,:);
            dipmodel(k).select = 1;
        end
    end
    ALLEEG(si).dipfit.model      = dipmodel;
    ALLEEG(si).dipfit.symDipRem = symDipRem;
end

