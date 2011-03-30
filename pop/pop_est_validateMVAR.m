function varargout = pop_est_validateMVAR(ALLEEG,typeproc,varargin)
%
% Validate a fitted VAR model. With two inputs, this function generates a
% GUI where the validation scheme can be specified. Validation consists of
% statistical tests for "whiteness" of fitted VAR model residuals [1,2], 
% consistency of the fitted model [1,3] and stability of fitted model [1-3]
% 
% Input:
%   
%   ALLEEG:     Array of EEGLAB data structures containing fitted MODEL
%   typeproc:   reserved for future use. Use 0
%
% Optional:
%
%     'whitenessCriteria':    Cell array containing names of residual whiteness 
%                             test criteria to evaluate. See [1, 2] for details.
%                             Possible Values: {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'}
%                             Default Value  : all
%                             Data Input Type: cell array
%
% 
%     'checkWhiteness':       Whether or not to check whiteness.
%                             Default Value  : true
%                             Data Input Type: boolean
% 
%     'checkConsistency':     Whether or not to check consistency. See [1,3]
%                             for details.
%                             Default Value  : true
%                             Data Input Type: boolean
% 
%     'checkStability':       Whether or not to check stability. See
%                             [1-3] for details.
%                             Default Value  : true
%                             Data Input Type: boolean
% 
%     'alpha':                significance level for determining whiteness
%                             Data Input Range: [0 1]
%                             Default Value   : 0.05
%                             Data Input Type : real number (double)
%
%     'prctWinToSample':      percent of time windows to randomly select  
%                             Data Input Range: [0 100]
%                             Default Value   : 100
%                             Data Input Type : real number (double)
%
%     'verb':                 verbosity level (0=no output, 1=text, 2=gui)
%   
%
% Output:
%   
%     lastcom:                Reserved for future use
%
%     whitestats:             Structure containing whiteness statistics.
%                             See est_checkMVARWhiteness() for details on
%                             structure format
%
%     PC:                     Vector of percent consistency estimates for 
%                             each window.
%
%     stability:              Vector of stability estimates for each window
%     
% See Also: est_checkMVARWhiteness(), est_checkMVARStability(),
%           est_checkMVARConsistency, pop_est_fitMVAR()
%
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 3.6 and 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% [2] Lutkepohl, H. (2007) New Introduction to Time Series Analysis.
%   Springer.
% 
% [3] Ding M, Bressler SL, Yang W, Liang H (2000) Short-window spectral 
%   analysis of cortical event-related potentials by adaptive multivariate 
%   autoregressive modeling: data preprocessing, model validation, and 
%   variability assessment. Biol. Cybern. 83:35-45 
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



% input new model structure
for i=1:nargout
    varargout{i} = [];
end

% display help if not enough arguments
% ------------------------------------
if nargin < 2
	help pop_est_validateMVAR;
	return;
end;	
lastcom = [];
whitestats = [];
PC = [];
stability = [];

if nargin < 4
    popup = 1;
else
    popup = 0;
end

for i=1:length(ALLEEG)
    if ~isfield(ALLEEG(i).CAT,'MODEL')
        error('ALLEEG.CAT.MODEL must be present in all datasets');
    else
        MODEL(i) = ALLEEG(i).CAT.MODEL;
    end
end

whitenessCriteria = {'Ljung-Box','ACF','Box-Pierce','Li-McLeod'};  %{'ACF'}; %{'Ljung-Box','ACF','Durbin-Watson','Portmanteau','LMP'};
numWins = length(MODEL(1).AR);

% pop up window
% -------------
if popup
% 	[txt vars] = gethelpvar('pop_est_selModelOrder.m');
	
	geomhoriz = {1 1 [3 1] 1 1 [3 1] };
    uilist = { ...
               { 'Style', 'checkbox', 'string', 'Check Whiteness of Residuals' ,'tag', 'chkWhiteness' ,'Value',1}...
               { 'Style', 'listbox', 'string', whitenessCriteria, 'tag', 'lstWhitenessCriteria','Value',1,'Min',1,'Max',20} ...
               { 'Style', 'text', 'string', 'significance level: '} ...
               { 'Style', 'edit', 'string', 0.05, 'tag','txtAlpha' } ...
               { 'Style', 'checkbox', 'string', 'check percent consistency ' ,'tag', 'chkConsistency'}...
               { 'Style', 'checkbox', 'string', 'check model stability ' ,'tag', 'chkStability' , 'enable','on'}...
               { 'Style', 'text', 'string', '% windows to sample'} ...
               { 'Style', 'edit', 'string', 100, 'tag','prctWinToSample' } ...
			 };

	[ tmp1 tmp2 strhalt result ] = inputgui( 'geometry', geomhoriz, 'geomvert',[1 3.5 1 1 1 1], 'uilist',uilist, 'helpcom','pophelp(''pop_est_validateMVAR'');', ...
					   'title','Select Model Validation Methods');
	if isempty( tmp1 ), return; end;
    
    if ~isempty(result.prctWinToSample)
        g.prctWinToSample = str2num(result.prctWinToSample);
    else
        g.prctWinToSample = 100;
    end
    
    g.checkConsistency = result.chkConsistency;
    g.checkWhiteness = result.chkWhiteness;
    g.checkStability = result.chkStability;
    g.whitenessCriteria = lower(hlp_variableize(whitenessCriteria(result.lstWhitenessCriteria)));
    if ischar(g.whitenessCriteria), g.whitenessCriteria={g.whitenessCriteria}; end
    g.verb = 2;
    g.winStartIdx = [];
    g.alpha = str2num(result.txtAlpha);
else  
    var = hlp_mergeVarargin(varargin{:});
    myargs = {'whitenessCriteria'   'cell'  whitenessCriteria   whitenessCriteria; ...
              'checkWhiteness',     'boolean'   []      true; ...
              'checkConsistency'    'boolean'   []      true; ...
              'checkStability'      'boolean'   []      true; ...
              'alpha'               'real'      [0 1]   0.05; ...    
              };
    g = finputcheck(var, [myargs ; hlp_getDefaultArglist('est')], 'est_checkWhiteness','ignore');
    if ischar(g), error(g); end
    g.alpha = 0.06;
end

numrows = sum([g.checkWhiteness g.checkConsistency g.checkStability]);
numcols = 1;

% determine which windows to use
if isempty(g.winStartIdx)
    % starting point of each window (points)
    g.winStartIdx  = floor(MODEL(1).winStartTimes*ALLEEG(1).srate)+1;    
end
if g.prctWinToSample<100 
    % randomly select percentage of windows to work with
    randwin = randperm(length(g.winStartIdx));
    randwin = sort(randwin(1:ceil(length(g.winStartIdx)*g.prctWinToSample/100)));
    g.winStartIdx = g.winStartIdx(randwin);
    for cnd = 1:length(MODEL);
        ALLEEG(cnd).CAT.MODEL.AR = ALLEEG(cnd).CAT.MODEL.AR(randwin);
        ALLEEG(cnd).CAT.MODEL.PE = ALLEEG(cnd).CAT.MODEL.PE(randwin);
        ALLEEG(cnd).CAT.MODEL.winStartTimes = ALLEEG(cnd).CAT.MODEL.winStartTimes(randwin);
    end
    
    g.prctWinToSample = 100;
end

g.winlen = ALLEEG(1).CAT.MODEL.winlen;
g.winstep = ALLEEG(1).CAT.MODEL.winstep;
g.morder = ALLEEG(1).CAT.MODEL.morder;

for cond=1:length(ALLEEG)
    if g.verb
        fprintf('checking condition %s...\n',ALLEEG(cond).condition);
    end
    
    if g.checkWhiteness
        [whitestats] = est_checkMVARWhiteness(ALLEEG(cond),MODEL(cond),typeproc,g);
    end
    
    if g.checkConsistency
        PC = est_checkMVARConsistency(ALLEEG(cond),MODEL(cond),typeproc,g);
    end
    
    if g.checkStability
        [stability lambda] = est_checkMVARStability(ALLEEG(cond),MODEL(cond),typeproc,g);
    end
    


    % plot results
    figure('Name',sprintf('%s - Model Validation Results',ALLEEG(cond).condition));
    curplot=1;
    if g.checkWhiteness

        subplot(numrows,numcols,curplot);
        for i = 1:length(g.whitenessCriteria)
            wcstr = lower(hlp_variableize(g.whitenessCriteria{i}));
            wc = whitestats.(wcstr);
            pvals(i,:) = wc.pval;
            lgnd{i} = sprintf('%s (%d/%d white)',wc.fullname, sum(wc.w),length(wc.w));
        end

        if size(pvals,2)>1
            % more than one window -- make lineplot
            plot(1:length(g.winStartIdx),pvals);
            xlabel('Window number');
            legend(lgnd);
        else
            % single window -- make barplot
            h = bar(pvals);
            ch = get(h,'Children');

            set(gca,'xticklabel',lgnd);
            colors = [[1 0 0];[0 0 1];[0 1 0];[0 0 0];[1 0 1];[0 1 1]];
            colors = colors(1:length(pvals),:);
            set(ch,'FaceVertexCData',colors);

        end

        hl=hline(g.alpha); 
        set(gca,'xlim',[0 length(g.winStartIdx)+1],'Ylim',[max(0,min(pvals(:))-0.5), min(1,max(pvals(:))+0.5)]);
        set(hl,'linestyle','--','linewidth',2);
        ylabel({'Significance of whiteness','(larger is better)'});
        axcopy(gca);

        curplot=curplot+1;
    end

    if g.checkConsistency
        ax=subplot(numrows,numcols,curplot);
        if length(PC)>2
            % more than one window -- make lineplot       
            plot(1:length(g.winStartIdx),PC);
            axes(ax);
            text(0.98,0.9,sprintf('Mean PC: %0.2f%%',mean(PC)), ...
                'units','normalized','horizontalalignment','right', ...
                'edgecolor','k','backgroundcolor','w');
            xlabel('Window number');
        else
            % single window -- make barplot
            bar(PC);
            legend(sprintf('(%0.2f%% Consistent)',PC));
        end
        set(gca,'xlim',[0 length(g.winStartIdx)+1],'ylim',[min(PC)-50 min(max(PC)+50,100)]);
        ylabel('Percent Consistency');
        axcopy(gca);
        curplot = curplot+1;
    end

    if g.checkStability
        % plot stability results
        subplot(numrows,numcols,curplot);
        if length(stability)>2
            % more than one window -- make lineplot    
            %boxplot(real(lambda)');
            maxlambda = max(real(lambda),[],2);
            plot(1:length(g.winStartIdx),maxlambda);
            xlabel('Window number')
        else
            % single window -- make barplot
            maxlambda = max(real(lambda(:)));
            bar(maxlambda);
        end


    %     set(gca,'ylim',[max(0,0.7*min(abs(lambda(:)))) max(1.3,1.3*max(abs(lambda(:))))]);
        set(gca,'xlim',[0 length(g.winStartIdx)+1],'ylim',[1.2*min(maxlambda(:)) max(0.01,1.3*max(maxlambda(:)))]);
    %     axis auto
        hl=hline(0);
        set(hl,'linestyle','--','linewidth',2);
        ylabel({'Stability Index','(should be < 0)'});
        numstable = sum(stability);
        legend(sprintf('(%d/%d stable)',numstable,length(stability)));
        axcopy(gca);
    end

    try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;

    varargout{1}{cond} = lastcom;
    varargout{2}{cond} = whitestats;
    varargout{3}{cond} = PC;
    varargout{4}{cond} = stability;
end

