function handles = vis_plotModelValidation(whitestats,PCstats,stabilitystats,varargin)
% Visualize the results of model validation as computed by
% pop_est_validateMVAR()
%
% Inputs:
%
%       whitestats:     Cell array of structs each containing results of
%                       tests for whiteness of residuals for a single
%                       dataset/condition
%                       (as output by est_checkMVARWhiteness())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%       PCstats:        Cell array of structs each containing results of
%                       tests for consistency of the model for a single
%                       dataset/condition
%                       (as output by est_checkMVARConsistency())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%       stabilitystats: Cell array of structs each containing results of
%                       tests for stability of the model for a single
%                       dataset/condition
%                       (as output by est_checkMVARStability())
%                       If not available, pass in an empty matrix to skip
%                       plotting this result.
%
%
% Optional:
%
%       conditions:     Cell array of condition labels
%                       e.g. conditions = {ALLEEG.condition}
%
% Outputs:
%
%       Figures will be rendered showing model validation results for each
%       dataset. The handles to the figures are returned.
%
% See Also: pop_est_validateMVAR(),est_checkMVARWhiteness(),
%           est_checkMVARStability(), est_checkMVARConsistency()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen, 2010, SCCN/INC, UCSD.
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


% if ~isempty(varargin) && isstruct(varargin{1})
%     varargin = [hlp_struct2varargin(varargin{1}) varargin];
% end

whitenessCriteria = {};

% find the names of all whiteness tests stored in whitestats
if ~isempty(whitestats{1})
    fn = fieldnames(whitestats{1});
    for i=1:length(fn)
        if isfield(whitestats{1}.(fn{i}),'fullname')
            whitenessCriteria = [whitenessCriteria fn{i}];
        end
    end
end

g = finputcheck(varargin, ...
    {'whitenessCriteria'   'cell'  whitenessCriteria   whitenessCriteria; ...
    'checkWhiteness',     'boolean'   []          true; ...
    'checkConsistency'    'boolean'   []          true; ...
    'checkStability'      'boolean'   []          true; ...
    'conditions'          'cell'      {}          {}; ...
    },'mode','ignore','quiet');


if isempty(whitestats{1}),         g.checkWhiteness    = false;    end
if isempty(PCstats{1}),            g.checkConsistency  = false;    end
if isempty(stabilitystats{1}),     g.checkStability    = false;    end

num_conds = max([length(whitestats),length(PCstats),length(stabilitystats)]);

numrows = sum([g.checkWhiteness g.checkConsistency g.checkStability]);
numcols = 1;

for cond = 1:num_conds
    
    if isempty(g.conditions)
        g.conditions{cond} = sprintf('Condition %d',cond);
    end
    
    % plot results
    handles(cond) = figure('Name',sprintf('%s - Model Validation Results',g.conditions{cond}));
    curplot=1;
    
    if g.checkWhiteness
        
        if ~iscell(whitestats), whitestats = {whitestats}; end
        
        subplot(numrows,numcols,curplot);
        for i = 1:length(g.whitenessCriteria)
            wcstr = lower(hlp_variableize(g.whitenessCriteria{i}));
            wc = whitestats{cond}.(wcstr);
            pvals(i,:) = wc.pval;
            if size(pvals,2)>1
                lgnd{i} = sprintf('%s (%d/%d white)',wc.fullname, sum(wc.w),length(wc.w));
            else
                lgnd{i} = sprintf('%s (%swhite)',wc.fullname, fastif(wc.w,'','not '));
            end
        end
        
        if size(pvals,2)>1
            % more than one window -- make lineplot
            plot(1:length(whitestats{cond}.winStartIdx),pvals');
            xlabel('Window number');
            legend(lgnd);
            set(gca,'xlim',[0 length(whitestats{cond}.winStartIdx)+1],'Ylim',[max(0,min(pvals(:))-0.5), min(1,max(pvals(:))+0.5)]);
        else
            % single window -- create bar plots
            h = bar(pvals);
            ch = get(h,'Children');
            
            set(gca,'xticklabel',lgnd);
            colors = [[1 0 0];[0 0 1];[0 1 0];[0 0 0];[1 0 1];[0 1 1]];
            colors = colors(1:length(pvals),:);
            set(ch,'FaceVertexCData',colors);
%             set(gca,'xlim',[0 length(length(g.whitenessCriteria))+1],'Ylim',[max(0,min(pvals(:))-0.5), min(1,max(pvals(:))+0.5)]);
        end
        
        hl=hline(whitestats{cond}.alpha);
        set(hl,'linestyle','--','linewidth',2);
        ylabel({'Whiteness p-value'});
        axcopy(gca);
        
        curplot=curplot+1;
        
        if ismember('acf',g.whitenessCriteria)
            hl=hline(1-whitestats{cond}.alpha,'r');
            set(hl,'linestyle','--','linewidth',2);
        end
    end
    
    
    
    if g.checkConsistency
        
        if ~iscell(PCstats), PCstats = {PCstats}; end
        
        ax=subplot(numrows,numcols,curplot);
        if length(PCstats{cond}.PC)>2
            % more than one window -- make lineplot
            plot(1:length(PCstats{cond}.winStartIdx),PCstats{cond}.PC);
            axes(ax);
            text(0.98,0.9,sprintf('Mean PC: %0.2f%%',mean(PCstats{cond}.PC)), ...
                'units','normalized','horizontalalignment','right', ...
                'edgecolor','k','backgroundcolor','w');
            xlabel('Window number');
        else
            % single window -- make barplot
            bar(PCstats{cond}.PC);
            legend(sprintf('(%0.2f%% Consistent)',PCstats{cond}.PC));
        end
        set(gca,'xlim',[0 length(PCstats{cond}.winStartIdx)+1],'ylim',[min(PCstats{cond}.PC)-50 min(max(PCstats{cond}.PC)+50,100)]);
        ylabel('Percent Consistency');
        axcopy(gca);
        curplot = curplot+1;
    end
    
    if g.checkStability
        
        if ~iscell(stabilitystats), stabilitystats = {stabilitystats}; end
        
        % plot stability results
        subplot(numrows,numcols,curplot);
        if length(stabilitystats{cond}.stability)>2
            % more than one window -- make lineplot
            %boxplot(real(lambda)');
            maxlambda = max(real(stabilitystats{cond}.lambda),[],2);
            plot(1:length(stabilitystats{cond}.winStartIdx),maxlambda);
            xlabel('Window number')
        else
            % single window -- make barplot
            maxlambda = max(real(stabilitystats{cond}.lambda(:)));
            bar(maxlambda);
        end
        
        
        %     set(gca,'ylim',[max(0,0.7*min(abs(lambda(:)))) max(1.3,1.3*max(abs(lambda(:))))]);
        set(gca,'xlim',[0 length(stabilitystats{cond}.winStartIdx)+1],'ylim',[1.2*min(maxlambda(:)) max(0.01,1.3*max(maxlambda(:)))]);
        %     axis auto
        hl=hline(0);
        set(hl,'linestyle','--','linewidth',2);
        ylabel({'Stability Index','(should be < 0)'});
        numstable = sum(stabilitystats{cond}.stability);
        legend(sprintf('(%d/%d stable)',numstable,length(stabilitystats{cond}.stability)));
        axcopy(gca);
    end
    
    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
    
end
