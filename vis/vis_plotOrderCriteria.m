function handles = vis_plotOrderCriteria(IC,conditions,icselector,minimizer,prclim)
% Visualize the results of model order selection as computed by
% est_selModelOrder()
% 
% Input: 
%
%       IC:             Cell array of structures containing results of
%                       model order selection for one or more datasets 
%                       (see pop_est_selModelOrder())
%
% Optional: 
%       conditions:     Cell array of strings containing names of
%                       conditions for figure labeling
%                       Default: {}
%       icselector:     Cell array of string denoting which information
%                       criteria to plot (must be fields of IC{.}).
%                       Default: all selectors in IC structure
%       minimizer:      'min' - find model order that minimizes info crit.
%                       'elbow' - find model order corresponding to elbow
%                                 of info crit. plot
%                       Default: 'min'
%       prclim:         upper limit for percentile
%                       Default: 90
%                       
% Output:
%
%       Figures will be rendered showing model selection results for each
%       dataset. The handles to the figures are returned.
%
% See Also: est_selModelOrder(), pop_est_selModelOrder()
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

if nargin<2
    conditions = [];
end

if nargin<3 || isempty(icselector)
    icselector = IC{1}.selector;
end

if nargin<4
    minimizer = 'min';
end

if nargin<5
    prclim = 90;
end

for cond = 1:length(IC)
    
    if isempty(conditions)
        conditions{cond} = sprintf('Condition %d',cond);
    end
    
    numinfocriteria = length(icselector);
    numWins = size(IC{cond}.(icselector{1}).ic,2);


    % plot the results
    % ----------------------
    handles(cond) = figure('name',sprintf('%s - Model Order Selection Results (%s)',conditions{cond},fastif(strcmpi(minimizer,'min'),'min ic','elbow ic')));
    for i=1:numinfocriteria
        allinfocrit(i,:) = mean(IC{cond}.(lower(icselector{i})).ic,2);
    end

    if numWins>1
        numrows = floor(sqrt(numinfocriteria))+1;
        numcols = ceil(numinfocriteria/(numrows-1));
    else
        numrows = 1;
        numcols = 1;
    end

    set(gca,'Position',[0 0 1 1]);

    % plot information criteria
    subplot(numrows,numcols,1:numcols); 
    plot(IC{cond}.pmin:IC{cond}.pmax,allinfocrit,'linewidth',2);
    set(gca,'xtick',IC{cond}.pmin:IC{cond}.pmax);
    axis auto
    xlabel('model order    ','fontsize',12);
    ylabel('Information criteria (bits)    ','fontsize',12);
    set(gca,'xgrid','on');
    title('Mean IC across sampled windows   ','fontsize',12);
    axcopy(gca);

    hold on

    defcolororder = get(0,'defaultaxescolororder');
    
    % make markers at minima
    for i=1:numinfocriteria
        sel = icselector{i};
        
        switch minimizer
            case 'min'
                [minic popt] = min(allinfocrit(i,:)); %ceil(mean(IC{cond}.(lower(sel)).popt));
            case 'elbow'
                [minic popt] = hlp_findElbow(allinfocrit(i,:));
        end
        
%         minic = allinfocrit(i,popt);
        plot(popt+IC{cond}.pmin-1,minic,'rs','markersize',10,'MarkerEdgeColor','k','markerfacecolor',defcolororder(i,:));
        lmin = vline(popt+IC{cond}.pmin-1);
        set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
        legendstr{i} = sprintf('%s (%d)',icselector{i},popt+IC{cond}.pmin-1);
    end
    
%     for i=1:numinfocriteria
%         sel = icselector{i};
%         popt = ceil(mean(IC{cond}.(lower(sel)).popt));
%         minic = allinfocrit(i,popt-IC{cond}.pmin+1);
%         plot(popt,minic,'rs','markersize',10,'MarkerEdgeColor','k','markerfacecolor',defcolororder(i,:));
%         lmin = vline(popt);
%         set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
%         legendstr{i} = sprintf('%s (%d)',icselector{i},popt);
%     end

    xlim([IC{cond}.pmin IC{cond}.pmax]);
    yl=ylim;
    ylim([yl(1)-0.001*(yl(2)-yl(1)) yl(2)]);
    
    legend(legendstr,'linewidth',2);
    
    % popt = popt+IC{cond}.pmin-1;
    % plot(popt',minic);

    % popt = popt(:);
    % ylim = get(gca,'Ylim');
    % ylim = ylim(ones(1,numinfocriteria),:);
    % line([popt popt],ylim);
    
    hold off
    % 


    if numWins>1
        % plot histograms of optimal model order across windows
        for i=1:numinfocriteria
            ax=subplot(numrows,numcols,i+numcols);
            xScale = IC{cond}.pmin:IC{cond}.pmax;
            
            switch minimizer
                case 'min'
                    bar(xScale,histc(IC{cond}.(lower(icselector{i})).popt,xScale),'k');         %  plot histogram
                    popt = ceil(mean(IC{cond}.(lower(icselector{i})).popt));                    %  mean
                    poptstd = ceil(std(IC{cond}.(lower(icselector{i})).popt));                  %  stdev
                    poptprctile = ceil(prctile(IC{cond}.(lower(icselector{i})).popt,prclim));   %  upper 95th prctile
                case 'elbow'
                    bar(xScale,histc(IC{cond}.(lower(icselector{i})).pelbow,xScale),'k');
                    popt = ceil(mean(IC{cond}.(lower(icselector{i})).pelbow));
                    poptstd = ceil(std(IC{cond}.(lower(icselector{i})).pelbow));
                    poptprctile = ceil(prctile(IC{cond}.(lower(icselector{i})).pelbow,prclim));
            end
            axes(ax);
            
            % shade region for stdev
            [patchHandles textHandles] = hlp_vrect([popt-poptstd popt+poptstd], ...
                'patchProperties', ...
                {'FaceColor',defcolororder(i,:), ...
                 'FaceAlpha',0.2,...
                 'EdgeColor',[0.2 0.2 0.2],...
                 'EdgeAlpha',0.2 ...
                 });
             
            % mark mean popt
            lmin = vline(popt,'-',[num2str(popt) '+-' num2str(poptstd)],[-0.01 0.95],gca,defcolororder(i,:));
            set(lmin,'color',defcolororder(i,:),'linestyle','-','linewidth',2);
            
            % mark upper prctile
            lprctile = vline(poptprctile,':',sprintf('%g%%',prclim),[-0.01 0.05],gca,defcolororder(i,:));
            set(lprctile,'color',defcolororder(i,:),'linewidth',2);
            
            title([icselector{i} '  '],'fontsize',12,'fontweight','bold','color',defcolororder(i,:));
            axcopy(gca);
            xlabel('opt. model order   ','fontsize',12);
            ylabel('histogram count   ','fontsize',12);
            xlim([IC{cond}.pmin-1 IC{cond}.pmax+1]);
            yl=ylim;
            ylim([yl(1)-0.001*(yl(2)-yl(1)) yl(2)]);
        end
    %     axis on
        %set(gca,'yaxislocation','right');
    end


    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
    
end