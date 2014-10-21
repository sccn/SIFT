function result = hlp_lookupAnatLabels(varargin)
% Lookup anatomical labels for source locations, using Talairach [1] or LONI atlas [2].
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Input                Information                                                                                                                                                                                      
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ...| SourcePos:      EEGLAB dipfit structure                                                                                                                                                                          
%                      Can also be an [N x 3] matrix of X,Y,Z dipole locations                                                                                                                                          
%                      Possible values: 'Unrestricted'                                                                                                                                                                  
%                      Default value  : 'n/a'                                                                                                                                                                           
%                      Input Data Type: any evaluable Matlab expression.                                                                                                                                                
%                                                                                                                                                                                                                       
% ...| ConfusionRange: Confusion Radius (mm)                                                                                                                                                                            
%                      Assumed standard deviation of uncertainty in the dipole fits, in millimeters.                                                                                                                    
%                      Possible values: [1 30]                                                                                                                                                                          
%                      Default value  : 4                                                                                                                                                                               
%                      Input Data Type: real number (double)                                                                                                                                                            
%                                                                                                                                                                                                                       
% ...| BrainAtlas:     Brain atlas to use                                                                                                                                                                               
%                      Talairach has a larger repertoire of areas, but is (in this version) non-probabilistic. The LONI LBPA40 atlas is a high-quality probabilistic atlas; however, you have to download it yourself.  
%                      Possible values: {'Talairach', 'LBPA40'}                                                                                                                                                         
%                      Default value  : 'Talairach'                                                                                                                                                                     
%                      Input Data Type: string  
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Output                Information                                                                                                                                                                                      
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% result               dipfit structure with annotated .model fields
%                       .structures
%                       .probabilities 
%                       .labels
% 
% References:
%   [1] Lancaster JL, Woldorff MG, Parsons LM, Liotti M, Freitas CS, Rainey L, Kochunov PV, Nickerson D, Mikiten SA, Fox PT, "Automated Talairach Atlas labels for functional brain mapping". 
%       Human Brain Mapping 10:120-131, 2000
%   [2] Shattuck DW, Mirza M, Adisetiyo V, Hojatkashani C, Salamon G, Narr KL, Poldrack RA, Bilder RM, Toga AW, "Construction of a 3D probabilistic atlas of human cortical structures."
%       Neuroimage. 2008 39(3):1064-80
%
% Author: Tim Mullen, SCCN/INC/UCSD, 2013
%         ** This function is closely derived from the BCILAB function set_fit_dipoles by Christian Kothe
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
        arg({'dipfit','SourcePos'},[],[],'EEGLAB dipfit structure. Can also be an [N x 3] matrix of X,Y,Z dipole locations','type','expression'), ...
        arg({'confusion_range','ConfusionRange'}, 4, [1 30], 'Confusion Radius (mm). Assumed standard deviation of uncertainty in the dipole fits, in millimeters.'),...    
        arg({'brain_atlas','BrainAtlas'}, 'Talairach', {'Talairach','LBPA40'}, 'Brain atlas to use. Talairach has a larger repertoire of areas, but is (in this version) non-probabilistic. The LONI LBPA40 atlas is a high-quality probabilistic atlas; however, you have to download it yourself.'),...
        arg({'xyzorder'},[2 1 3],[],'Order of X,Y,Z coordinates in dipfit.posxyz.') ...
    );

if all(cellfun(@isempty,strfind(javaclasspath,'talairach.jar')))
    javaaddpath(which('talairach.jar'));
end

if ismatrix(g.dipfit)
    hlp_posxyz2dipfit(g.dipfit,'MNI');
end

switch lower(g.brain_atlas)
    
    case 'talairach'
        db = org.talairach.Database;
        db.load(fullfile(hlp_getSiftRoot,'resources','talairach.nii'));
        for k=1:length(g.dipfit.model)
            try
                p = icbm_spm2tal(g.dipfit.model(k).posxyz(g.xyzorder));
                g.dipfit.model(k).labels = cellfun(@(d)char(d),cell(db.search_range(p(1),p(2),p(3),1.5*g.confusion_range)),'UniformOutput',false);
                % and compute structure probabilities within the selected volume
                [structures,x,idxs] = unique_bc(hlp_split(sprintf('%s,',g.dipfit.model(k).labels{:}),',')); %#ok<ASGLU>
                probabilities = mean(bsxfun(@eq,1:max(idxs),idxs'));
                [probabilities,reindex] = sort(probabilities,'descend');
                structures = structures(reindex);
                mask = ~strcmp(structures,'*');
                g.dipfit.model(k).structures = structures(mask);
                g.dipfit.model(k).probabilities = probabilities(mask)*5; % there are 5 partitions
                % figure('Position',[0 0 2560 900]); topoplot(tmp.icawinv(:,k),tmp.chanlocs); title([hlp_tostring(structures(mask)) 10 hlp_tostring(5*probabilities(mask))]); pop_dipplot(tmp,k);
            catch
                g.dipfit.model(k).labels = {};
                g.dipfit.model(k).structures = {};
                g.dipfit.model(k).probabilities = [];
            end
        end
        
    case 'lbpa40'
        mask = find(~cellfun('isempty',{g.dipfit.model.posxyz}));
        % build query
        coords = vertcat(g.dipfit.model.posxyz);
        coords = [coords g.confusion_range * ones(size(coords,1),1)];
        % look up from atlas (note: slow!)
        [probs,labels] = label_dipoles(coords);
        % determine overlapped specialty regions to match up with Talairach's labeling scheme
        L = ~cellfun('isempty',strfind(labels,' L '));
        R = ~cellfun('isempty',strfind(labels,' R '));
        B = ~cellfun('isempty',strfind(labels,'Brainstem'));
        C = ~cellfun('isempty',strfind(labels,'Cerebellum'));
        allprobs = [sum(probs(:,L),2) sum(probs(:,R),2) sum(probs(:,C),2)*[0.5 0.5] sum(probs(:,B),2)];
        allstructs = {'Left Cerebrum' 'Right Cerebrum' 'Left Cerebellum' 'Right Cerebellum' 'Brainstem'};
        % go through all gyri and add up left & right probabilities
        gyri = unique(cellfun(@(l)l(12:end),labels(1:end-2),'UniformOutput',false))';
        allstructs = [allstructs gyri];
        for g=1:length(gyri)
            curgyrus = gyri{g};
            matches = ~cellfun('isempty',strfind(labels,curgyrus));
            allprobs = [allprobs sum(probs(:,matches),2)];
        end
        for k=1:length(mask)
            % retain only those with non-zero probability
            sel = allprobs(k,:) ~= 0;
            probabilities = allprobs(k,sel);
            structures = allstructs(sel);
            % sort probs & associated structs by descending probability
            [probabilities,reindex] = sort(probabilities,'descend');
            structures = structures(reindex);
            % store in the model...
            g.dipfit.model(mask(k)).structures = structures;
            g.dipfit.model(mask(k)).probabilities = probabilities;
            % figure('Position',[0 0 2560 900]);topoplot(tmp.icawinv(:,k),tmp.chanlocs); title([hlp_tostring(structures) 10 hlp_tostring(probabilities)]); pop_dipplot(tmp,k);
        end
end

result = g.dipfit;