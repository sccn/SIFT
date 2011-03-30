function IC = pop_est_selModelOrder(ALLEEG,typeproc,varargin)
%
% Fit a series of MVAR models up to a specified model order and compute the
% model order selection (information) criteria. If the model fitting
% algorithm is 'vierra-morf' then we fit a single model up to maximum order
% and downdate the prediction errors. This function generates a figure
% containing the results of the model fitting.
% If nargin == 2, an input GIU will be generated. Otherwise
% est_selModelOrder is called directly.
%
% Input:
%
%   EEG                Preprocessed EEG structure. Must contain .CAT
%   typeproc           Reserved for future use. Use 0.
%
% Optional:            arguments to est_selModelOrder()
%
%     'icselector'         cell array of strings denoting which model order 
%                          selection criteria to estimate
%     'downdate'           [true, false] whether or not to use downdated noise
%                          covariance matrices
%                          ('aic','sbc','fpe','hq')
%     'algorithm',         string denoting which algorithm to use for model
%                           fitting ('vierra-morf','arfit')
%     'winStartIdx'        vector of sample points (start of windows) at which to estimate windowed VAR model
%     'morder',            [min max] VAR model order to fit
%     'winlen',            window length (sec)
%     'winstep',           window step size (sec)
%     'epochTimeLims',     time range to analyze (sec) where 0 = start of the epoch
%     'prctWinToSample',   percent of time windows to randomly select  [0 100]
%     'verb',              verbosity level (0=no output, 1=text, 2=gui)
%
% Output:
%
%   IC                 a structure containing results of model order selection
%                      IC.selector     - the chosen information criteria
%                      IC.pmin         - the minimum model order tested
%                      IC.pmax         - the maximum model order tested
%                      IC.('sel') contains results for a selector 'sel'.
%                      This consists of subfields
%                           .ic         - [P numwins] matrix of information
%                                         critera for all P model orders tested
%                                         P = morder(2)-morder(1)+1 is the
%                                         number of model orders tested
%                           .minic      - the minimum of ic across model
%                                         orders
%                           .popt       - the model order that minimizes ic
% 
% See Also: est_selModelOrder()
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%   Springer.
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



IC = {};

% display help if not enough arguments
% ------------------------------------
if nargin < 2 || typeproc~=0
	help pop_est_selModelOrder;
	return;
end;	
lastcom = [];

var = hlp_mergeVarargin(varargin{:});
g = finputcheck(var, hlp_getDefaultArglist('est'), 'pop_est_selModelOrder','ignore');
if ischar(g), error(g); end
if ~isfield(g,'icselector'), g.icselector = {'sbc','aic'}; end
    
popup = 1; %nargin<3;

% possible information criteria
orderCriteria = {'AIC','FPE','SBC','HQ'};

if ~isempty(g.morder) && length(g.morder) == 2
    pmin = g.morder(1);
    pmax = g.morder(2);
else
    pmin = 1;
    pmax = 30;
end

% pop up window
% -------------
if popup
% 	[txt vars] = gethelpvar('pop_est_selModelOrder.m');
	
	geomhoriz = {1 1 1 1 [3 1 0.5 1] [3 2.5] };
    uilist = { ...
               { 'Style', 'text', 'string', 'Select order criteria to estimate' }...
               { 'Style', 'text', 'string', '(hold Ctrl to select multiple)' }...
               { 'Style', 'listbox', 'string', orderCriteria, 'tag', 'lstOrderCriteria','Value',1,'Min',1,'Max',20} ...
               { 'Style', 'checkbox', 'string', 'Downdate model', 'value', fastif(any(ismember({'arfit','vieira-morf-cpp'},g.algorithm)),false,true),'tag', 'chkDowndate','enable',fastif(any(ismember({'arfit','vieira-morf-cpp'},g.algorithm)),'off','on')} ...
               { 'Style', 'text', 'string', 'model order range: ' }...
               { 'Style', 'edit', 'string', pmin, 'tag','edtMin'}...
               { 'Style', 'text', 'string', '-'} ...
               { 'Style', 'edit', 'string', pmax, 'tag', 'edtMax' }...
               { 'Style', 'text', 'string', '% windows to sample'} ...
               { 'Style', 'edit', 'string', 20, 'tag','prctWinToSample' } ...
			 };

	[ tmp1 tmp2 strhalt result ] = inputgui( 'geometry', geomhoriz, 'geomvert',[1 1 3.5 1 1 1], ...
                        'uilist',uilist, 'helpcom','pophelp(''pop_est_selModelOrder'');', ...
					    'title','Plot Information Criteria');
	if isempty( tmp1 ), return; end;
    
    if ~isempty(result.prctWinToSample)
        g.prctWinToSample = str2num(result.prctWinToSample);
    else
        g.prctWinToSample = 100;
    end
    
    if isempty(result.lstOrderCriteria)
        errordlg2('You must select at least one order criterion','Model order selection');
        return;
    end
    
    g.downdate   = result.chkDowndate;
    g.icselector = lower(orderCriteria(result.lstOrderCriteria));
    g.morder  = [str2num(result.edtMin) str2num(result.edtMax)];
    
    if length(g.morder)~=2
        errordlg2('Model order range must be specified','Model order selection');
        return;
    end
    
    g.verb = 2;
end


for cond=1:length(ALLEEG)
    % calculate the information criteria
    IC{cond} = est_selModelOrder(ALLEEG(cond),g);

    numinfocriteria = length(g.icselector);
    numWins = size(IC{cond}.(g.icselector{1}).ic,2);


    % plot the results
    % ----------------------
    figure('name',sprintf('%s - Model Order Selection Results',ALLEEG(cond).condition));
    for i=1:numinfocriteria
        allinfocrit(i,:) = mean(IC{cond}.(lower(g.icselector{i})).ic,2);
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
        sel = g.icselector{i};
        popt = round(mean(IC{cond}.(lower(sel)).popt));
        minic = allinfocrit(i,popt-IC{cond}.pmin+1);
        plot(popt,minic,'rs','markersize',10,'MarkerEdgeColor','k','markerfacecolor',defcolororder(i,:));
        lmin = vline(popt);
        set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
        legendstr{i} = sprintf('%s (%d)',g.icselector{i},popt);
    end
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
            bar(xScale,histc(IC{cond}.(lower(g.icselector{i})).popt,xScale),'k');
            popt = round(mean(IC{cond}.(lower(g.icselector{i})).popt));
            axes(ax);
            lmin = vline(popt);
            set(lmin,'color',defcolororder(i,:),'linestyle','--','linewidth',2);
            title(g.icselector{i},'fontsize',12,'fontweight','bold');
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

end

   

