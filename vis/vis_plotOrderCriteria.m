function handles = vis_plotOrderCriteria(IC,conditions,icselector)
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
%       ALLEEG:         Array of EEG structures for all conditions (only
%                       used to obtain additional information about the 
%                       dataset for plotting). 
%                       Default: []
%       icselector:     Cell array of string denoting which information
%                       criteria to plot (must be fields of IC{.}).
%                       Default: all selectors in IC structure
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

if nargin<3
    icselector = IC{1}.selector;
end

for cond = 1:length(IC)
    
    if isempty(conditions)
        conditions{cond} = sprintf('Condition %d',cond);
    end
    
    numinfocriteria = length(icselector);
    numWins = size(IC{cond}.(icselector{1}).ic,2);


    % plot the results
    % ----------------------
    handles(cond) = figure('name',sprintf('%s - Model Order Selection Results',conditions{cond}));
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
    xlabel('model order','fontsize',12);
    ylabel('Information criteria (bits)','fontsize',12);
    set(gca,'xgrid','on');
    title('Mean IC across sampled windows','fontsize',12);
    axcopy(gca);

    hold on

    defcolororder = get(0,'defaultaxescolororder');
    
    % make markers at minima
    for i=1:numinfocriteria
        sel = icselector{i};
        [minic popt] = min(allinfocrit(i,:)); %ceil(mean(IC{cond}.(lower(sel)).popt));
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

    legend(legendstr);
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
            bar(xScale,histc(IC{cond}.(lower(icselector{i})).popt,xScale),'k');
            popt = ceil(mean(IC{cond}.(lower(icselector{i})).popt));
            poptstd = ceil(std(IC{cond}.(lower(icselector{i})).popt));
            axes(ax);
            lmin = vline(popt,'--',[num2str(popt) '+-' num2str(poptstd)],1.05);
            set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
            lstd = vline([popt-poptstd popt+poptstd]);
            set(lstd,'color',defcolororder(i,:),'linestyle',':','linewidth',1);
            title(icselector{i},'fontsize',12,'fontweight','bold');
            axcopy(gca);
            xlabel('opt. model order','fontsize',12);
            ylabel('histogram count','fontsize',12);
            xlim([IC{cond}.pmin IC{cond}.pmax]);
        end
    %     axis on
        %set(gca,'yaxislocation','right');
    end


    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;
    
end