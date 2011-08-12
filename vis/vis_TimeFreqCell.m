
function g = vis_TimeFreqCell(varargin)
%
% Plot bidirectional information flow between two sources. This is a
% detailed expansion of one cell of a TimeFrequencyGrid.
%
% Inputs:
% 
%       ConnMatrix:     [N x N x T x F] connectivity matrix, where
%                       N=numvars, T=numtimes, F=numfreqs
%       alltimes:       Times. Vector of timepoints (ordinate). See
%                       Conn.erWinCenterTimes
%       allfreqs:       'Frequencies. Vector of frequencies (abscissa). See
%                       Conn.freqs')
%       StatsMatrix:    Matrix of statistics
%
% Optional:
%
%     topovec:                    [2 x nchs] matrix of topoplots                                                                        
%                                 Input Data Type: real number (double)                                                                 
% 
%     dipfitstruct:               EEG.dipfit structure containing only models for [ch_j ch_i]                                           
%                                 Where ch_j -> ch_i                                                                                    
%                                 Input Data Type: real number (double)                                                                 
% 
%     chanlocs:                   Chanlocs structure                                                                                    
%                                 Input Data Type: string                                                                               
% 
%     chaninfo:                   Chaninfo structure                                                                                    
%                                 Input Data Type: string                                                                               
% 
%     connmethod:                 Connectivity method name                                                                              
%                                 Input Data Type: real number (double)                                                                 
% 
%     NodeLabels:                 Labels for the two nodes                                                                              
%                                 Input Data Type: real number (double)                                                                 
% 
%     TimeRange:                  Time Range to plot                                                                                    
%                                 Input Data Type: real number (double)                                                                 
% 
%     Frequencies:                Frequencies to plot                                                                                   
%                                 Input Data Type: real number (double)                                                                 
% 
%     Baseline:                   Time range of baseline [Min Max] (sec)                                                                
%                                 Will subtract baseline from each point. Leave blank for no baseline.                                  
%                                 Input Data Type: real number (double)                                                                 
% 
%     FrequencyScale:             Frequency Scale                                                                                       
%                                 Make the y-scale logarithmic or linear                                                                
%                                 Possible values: {'log','linear'}                                                                     
%                                 Default value  : 'linear'                                                                             
%                                 Input Data Type: string                                                                               
% 
%     LineWidth:                  Linewidth for marginals                                                                               
%                                 Input Range  : [1  Inf]                                                                               
%                                 Default value: 2                                                                                      
%                                 Input Data Type: real number (double)                                                                 
% 
%     EventMarkers:               Event marker time and style                                                                           
%                                 Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. {    
%                                 { 0.2 'y' ':' 2} { 1.5 'r' ':' 2}} will render two dotted-line event makers, yellow at 200 ms and     
%                                 red at 1500 ms                                                                                        
%                                 Input Data Type: any evaluable Matlab expression.                                                     
% 
%     Bidirectional:              Plot both directions                                                                                  
%                                 Input Data Type: boolean                                                                              
% 
%     SourceMarginPlot:           Source location plotting                                                                              
%                                 Options: 'Topoplot': plot source scalp projection. 'Dipole': plot dipole                              
%                                 Possible values: {'none','topoplot','dipole'}                                                         
%                                 Default value  : 'dipole'                                                                             
%                                 Input Data Type: string                                                                               
% 
%     DipolePlottingOptions:      Options for dipole plotting                                                                           
%                                 Input Data Type: string                                                                               
%     ----------------------                                                                                                            
% 
%         mri:                    Dipplot MRI structure                                                                                 
%                                 Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a    
%                                 path to a Matlab file containing MRI structure. Default uses MNI brain.                               
%                                 Input Data Type: string                                                                               
% 
%         DipoleCoordinateFormat: Coordinate format for dipplot                                                                         
%                                 Possible values: {'spherical','mni'}                                                                  
%                                 Default value  : 'mni'                                                                                
%                                 Input Data Type: string                                                                               
% 
%         DipplotOptions:         Additional dipplot options                                                                            
%                                 Cell array of <'name',value> pairs of additional options for dipplot (see 'doc dipplot')              
%                                 Input Data Type: any evaluable Matlab expression.                                                     
% 
%     TitleString:                Figure time string                                                                                    
%                                 Input Data Type: string                                                                               
% 
%     TitleFontSize:              Title Font Size                                                                                       
%                                 Input Data Type: real number (double)                                                                 
% 
%     AxesFontSize:               Axes Font Size                                                                                        
%                                 Input Data Type: real number (double)                                                                 
% 
%     TextColor:                  Text color                                                                                            
%                                 See 'doc ColorSpec'.                                                                                  
%                                 Input Data Type: any evaluable Matlab expression.                                                     
% 
%     BackgroundColor:            Background Color                                                                                      
%                                 See 'doc ColorSpec'.                                                                                  
%                                 Input Data Type: any evaluable Matlab expression.                                                     
% 
%     ColorLimits:                Color scaling limits                                                                                  
%                                 If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                 scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                 prctile(abs(Conn),scalar)                                                                             
%                                 Input Data Type: real number (double)                                                                 
% 
%     Thresholding:               Thresholding options                                                                                  
%                                 You can choose to use statistics (passed in as 'stats' structure), or simple percentile or absolute   
%                                 thresholds.                                                                                           
%                                 Possible values: {'None','Statistics','Simple'}                                                       
%                                 Default value  : 'None'                                                                               
%                                 Input Data Type: string                                                                               
%     -------------                                                                                                                     
% 
%         AlphaSignificance:      P-value threshold for significance. e.g., 0.05 for p<0.05                                             
%                                 Input Range  : [0  1]                                                                                 
%                                 Default value: 0.05                                                                                   
%                                 Input Data Type: real number (double)                                                                 
% 
%         PercentileThreshold:    Percentile threshold                                                                                  
%                                 If of form [percentile, dimension], percentile is applied elementwise across the specified            
%                                 dimension.                                                                                            
%                                 Input Data Type: real number (double)                                                                 
% 
%         AbsoluteThreshold:      Exact threshold                                                                                       
%                                 Input Data Type: real number (double)     
%
% 
% See Also: pop_vis_TimeFreqGrid(), vis_TimeFreqGrid(), 
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



if ~isempty(varargin)
    % determine whether statistics are present
    if isfield(varargin{1},'icaact') && length(varargin)==2
        stats = [];
    elseif isfield(varargin{1},'icaact')
        stats = arg_extract(varargin(3:end),'stats');
    else
        stats = arg_extract(varargin,'stats');
    end
    
    if isempty(stats)
        usestatsdef = [];  % false
    else
        usestatsdef = {};  % true
    end
    clear stats;
end


g = arg_define([0 Inf],varargin, ...
    arg_norep({'ConnMatrix','ConnectivityMatrix'},mandatory), ...
    arg_norep({'alltimes','AllTimes'},mandatory,[],'Times. Vector of timepoints (ordinate). See Conn.erWinCenterTimes'), ...
    arg_norep({'allfreqs','AllFrequencies'},mandatory,[],'Frequencies. Vector of frequencies (abscissa). See Conn.freqs'), ...
    arg_norep({'StatsMatrix','StatisticsMatrix'},[],[],'Matrix of statistics'), ...
    arg_nogui({'topovec'},[],[],'[2 x nchs] matrix of topoplots.'), ...
    arg_nogui({'dipfitstruct'},[],[],'EEG.dipfit structure containing only models for [ch_j ch_i]. Where ch_j -> ch_i'), ...
    arg_nogui({'elocs','chanlocs'},mandatory,[],'Chanlocs structure'), ...
    arg_nogui({'chaninfo'},mandatory,[],'Chaninfo structure'), ...
    arg_nogui({'connmethod'},[],[],'Connectivity method name'), ...
    arg_nogui({'nodelabels','NodeLabels'},[],[],'Labels for the two nodes'), ...
    arg({'timeRange','TimeRange'},[],[],'Time Range to plot','cat','DisplayProperties'), ...
    arg({'freqValues','Frequencies'},[],[],'Frequencies to plot','cat','DisplayProperties'), ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point. Leave blank for no baseline.','cat','DataProcessing'), ...
    arg({'freqscale','FrequencyScale'},'linear',{'log','linear'},'Frequency Scale. Make the y-scale logarithmic or linear','cat','DisplayProperties'), ...
    arg({'linewidth','LineWidth'},2,[1 Inf],'Linewidth for marginals','cat','DisplayProperties'), ...
    arg({'events','EventMarkers'},[],[],'Event marker time and style. Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. { { 0.2 ''y'' '':'' 2} { 1.5 ''r'' '':'' 2}} will render two dotted-line event makers, yellow at 200 ms and red at 1500 ms','type','expression','shape','row','cat','DisplayProperties'), ...
    arg_norep({'colormap','Colormap'},'jet(300)',[],'Colormap. Matlab expression denoting colormap to use (e.g., ''jet(64)''). See ''help colormap''.','type','expression','cat','DisplayProperties'), ...
    arg({'bidir','Bidirectional'},true,[],'Plot both directions','cat','DisplayProperties'), ...
    arg({'topoplot','SourceMarginPlot'},'dipole',{'none','topoplot','dipole'},'Source location plotting. Options: ''Topoplot'': plot source scalp projection. ''Dipole'': plot dipole','cat','DisplayProperties'), ...
    arg_sub({'dipplot','DipolePlottingOptions'},[], ...
    { ...
        arg_nogui({'mri'},'',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
        arg_nogui({'coordformat','DipoleCoordinateFormat'},'mni',{'spherical','mni'},'Coordinate format for dipplot','type','char','shape','row'), ...
        arg_nogui({'dipplotopt','DipplotOptions'},'{}','','Additional dipplot options. Cell array of <''name'',value> pairs of additional options for dipplot (see ''doc dipplot'')','type','expression','shape','row') ...
    },'Options for dipole plotting'), ...
    arg({'titleString','TitleString'},'','','Figure time string','type','char','cat','TextAndFont'), ...
    arg({'titleFontSize','TitleFontSize'},12,[],'Title Font Size','cat','TextAndFont'), ...
    arg({'axesFontSize','AxesFontSize'},10,[],'Axes Font Size','cat','TextAndFont'), ...
    arg({'textColor','TextColor'},[0 0 0],[],'Text color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg({'backgroundColor','BackgroundColor'},[0 0 0],[],'Background Color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg({'clim','ColorLimits'},[0 100],[],'Color scaling limits. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','shape','row','cat','DisplayProperties'), ...
    arg_subswitch({'thresholding','Thresholding'},'None', ...
    {'None' { ...
    arg_norep({'dummy1'},[],[],'dummy') ...
    }, ...
    'Statistics' {...
    arg({'alpha','AlphaSignificance'},0.05,[0 1],'P-value threshold for significance. e.g., 0.05 for p<0.05') ...
    }, ...
    'Simple' {...
    arg({'prcthresh','PercentileThreshold'},100,[],'Percentile threshold. If of form [percentile, dimension], percentile is applied elementwise across the specified dimension.','shape','row','cat','Thresholding'), ...
    arg({'absthresh','AbsoluteThreshold'},[],[],'Exact threshold.','cat','Thresholding') ...
    } ...
    }, 'Thresholding options. You can choose to use statistics (passed in as ''stats'' structure), or simple percentile or absolute thresholds.','cat','Thresholding') ...
    );

% check inputs and handle defaults
% ---------------------------------------------------
g.applyThreshold = ~strcmpi(g.thresholding.arg_selection,'none');

if ndims(g.ConnMatrix) == 2
    % make ConnMatrix a [1 x N x M] matrix
    tmp          = [];
    tmp(1,:,:)   = g.ConnMatrix;
    g.ConnMatrix = tmp;
end

if g.applyThreshold && strcmpi(g.thresholding.arg_selection,'Statistics') && isempty(g.StatsMatrix)
    error('SIFT:StatsMatrixRequired', ...
        'You must provide a Statistics matrix to use statistical thresholding');
end

if g.bidir && size(g.ConnMatrix,1) == 1
    warning('SIFT:InvalidNumberOfDirections', ...
        'Missing both feedforward and feedback matrices in g.ConnMatrix. Plotting only the first direction');
    g.bidir = 1;
end

if strcmpi(g.topoplot,'dipplot') && isempty(g.dipfitstruct)
    error('SIFT:DipfitmodelRequired',...
        'Argument ''dipfitstruct'' must be an valid dipfit structure containing 2 EEG dipfit.model structures');
end

if strcmpi(g.topoplot,'topoplot') && isempty(g.topovec)
    error('SIFT:TopovecRequired',...
        'Argument ''topovec'' must contain a valid matrix of topographic maps');
end

if isempty(g.freqValues)
    g.freqValues    = g.allfreqs;   end
if isempty(g.timeRange)
    g.timeRange     = g.alltimes([1 end]);   end

% get indices of selected time range
timeIndices     = getindex(g.alltimes,g.timeRange);
timeIndices     = timeIndices(1):timeIndices(2);
g.alltimes      = g.alltimes(timeIndices);

% get indices of selected freq range
freqIndices = getindex(g.allfreqs,g.freqValues);
g.allfreqs  = g.allfreqs(freqIndices);

% select times and frequencies
g.ConnMatrix                = g.ConnMatrix(:,freqIndices,timeIndices);
if ~isempty(g.StatsMatrix)
    switch ndims(g.StatsMatrix)
        case 3
            g.StatsMatrix    = g.StatsMatrix(:,freqIndices,timeIndices); 
        case 4
            g.StatsMatrix    = g.StatsMatrix(:,freqIndices,timeIndices,:);
    end
end

cmaplen = length(g.colormap);

% Rangle = angle(ConnOrig);

if g.bidir
    ordinate(1)  = 0.67;
    ordinate(2)  = 0.1;
    height       = 0.33;
    numdirs = 2;
else
    ordinate(1)  = 0.1;
    height       = 0.9;
    numdirs = 1;
end

% compute angles
% --------------
% Rangle = angle(ConnOrig);
% if ~isreal(ConnOrig)
%     ConnOrig = abs(ConnOrig);
%     Rraw =ConnOrig; % raw coherence values
%     setylim = 1;
%     if ~isnan(g.baseline)
%      	ConnOrig = ConnOrig - repmat(baseline',[1 g.timesout]); % remove baseline mean
%     end;
% else
%     Rraw = ConnOrig;
%     setylim = 0;
% end;

%
%
% if g.bidir
%     setylim = 1;
% else
%     setylim = 0;
% end



% create figure
figure('name',g.titleString,'DefaultAxesFontSize',g.axesFontSize); %,'color',g.backgroundColor
colormap(g.colormap);

% COLORMAP
%    map=hsv(300); % install circular color map - green=0, yellow, orng, red, violet = max
%    %                                         cyan, blue, violet = min
map = g.colormap;
% % cmapmid =  floor(cmaplen/2)+1;
% map = flipud([map(251:end,:);map(1:250,:)]);
% map(151,:) = map(151,:)*0.9; % tone down the (0=) green!

pos = get(gca,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];
axis('off')

[leftMargAxLim botMargAxLim] = deal(zeros(numdirs,4));

for dir=1:numdirs
    %
    % Image the coherence [% perturbations]
    %
    
    Conn = squeeze(g.ConnMatrix(dir,:,:,:,:));
    
    if all(isnan(Conn(:)))
        continue;
    end
    
    ConnOrig = Conn;
    
    %
    %     if numdirs>1
    %        % bidirectional
    %        Conn = squeeze(g.ConnMatrix(dir,:,:,:,:));
    %        ConnOrig = Conn;
    %     else
    %        Conn = g.ConnMatrix;     % for thresholding
    %        ConnOrig = Conn;         % unthresholded
    %     end
    %
    
    % remove baseline
    if ~isempty(g.baseline)
        % baseline mean coherence at each frequency
        baseidx  = getindex(g.alltimes,g.baseline);
        baseline = ConnOrig(:,baseidx(1):baseidx(2));
        baselineMean = mean(baseline,ndims(baseline));
        
        
        Conn = Conn-repmat(baselineMean,1,size(Conn,2));  %hlp_rmbaseline(Conn,g.baseline,g.alltimes);
        ConnOrig = Conn;
    else
        baselineMean = Conn;
    end
    
    % perform statistical thresholding
    % zero out (and 'green out') nonsignif. ConnOrig values
    if g.applyThreshold
        Stat = squeeze(g.StatsMatrix(dir,:,:,:,:));
        
        %         if numdirs>1
        %             Stat = squeeze(g.StatsMatrix(dir,:,:,:,:));
        %         else
        %             Stat = g.StatsMatrix;
        %         end
        
        switch dims(Stat)
            case 3, Conn  (Conn > Stat(:,:,1) & (Conn < Stat(:,:,2))) = 0;
                %                    Rraw(Conn > Stat(:,:,1) & (Conn < Stat(:,:,2))) = 0;
            case 2, Conn  (Conn < Stat) = 0;
                %                    Rraw(Conn < Stat) = 0;
            case 1, Conn  (Conn < repmat(Stat(:),[1 size(Conn,2)])) = 0;
                %                    Rraw(Conn < repmat(Stat(:),[1 size(Rraw,2)])) = 0;
        end;
    end
    
    
    % plot the information flow
    h(6) = axes('Units','Normalized', 'Position',[.1 ordinate(dir) .8 height].*s+q);
    colormap(map);
    
    if ~strcmpi(g.freqscale, 'log')
        try, imagesc(g.alltimes,g.allfreqs,Conn,max(Conn(:))*[-1 1]);
        catch, imagesc(g.alltimes,g.allfreqs,Conn,[-1 1]); end;
    else
        try, imagesclogy(g.alltimes,g.allfreqs,Conn,max(Conn(:))*[-1 1]);
        catch, imagesclogy(g.alltimes,g.allfreqs,Conn,[-1 1]); end;
    end
    set(gca,'ydir','norm');
    
    if ~isempty(g.clim)
        caxis([g.clim(1) g.clim(2)]);
    end
    tmpscale = caxis;
    
    % draw event markers
    hold on
    if ~isempty(g.events)
        for i=1:length(g.events)
            events = g.events{i};
            
            % set defaults
            if length(events) < 4
                events{4} = 2;      end
            if length(events) < 3
                events{3} = ':';    end
            if length(events) < 2
                events{2} = 'r';     end
            
            vl = vline(events{1});
            set(vl,'color',events{2},'linestyle',events{3},'linewidth',events{4});
        end
    end
    hold off
    
    
    % create legend
    if dir==1,
        ch1=1; ch2=2;
    else
        ch1=2; ch2=1;
    end
    
    if ~isempty(str2num(g.nodelabels{ch1}))
        ltext = sprintf('IC%s -> IC%s',g.nodelabels{ch1},g.nodelabels{ch2});
    else
        ltext = sprintf('%s -> %s',g.nodelabels{ch1},g.nodelabels{ch2});
    end
    text(0.98,0.92,ltext,'parent',gca,'units','normalized','fontweight','bold','fontsize',g.axesFontSize,'horizontalalignment','right','backgroundcolor',[1 1 1],'edgecolor',[0.3 0.3 0.3]);
    
    
    set(h(6),'YTickLabel',[],'YTick',[])
    set(h(6),'XTickLabel',[],'XTick',[])
    
    % create colorbar
    h(8) = axes('Position',[.95 ordinate(dir) .05 height].*s+q);
    
    cbar(h(8),1:length(map),[tmpscale(1) tmpscale(2)]);
    %
    %    if setylim
    %        cbar(h(8),1:length(map),[0 tmpscale(2)]);
    % %         cbar(h(8),151:300, [0 tmpscale(2)]); % use only positive colors (gyorv)
    %    else
    %        cbar(h(8),1:length(map)  , [-tmpscale(2) tmpscale(2)]);
    % %        cbar(h(8),1:300  , [-tmpscale(2) tmpscale(2)]); % use only positive colors (gyorv)
    %    end;
    %
    
    
    
  
    % Plot bottom marginal
    % ----------------------------
    
    h(10) = axes('Units','Normalized','Position',[.1 ordinate(dir)-0.1 .8 .1].*s+q);
    Emax = max(ConnOrig,[],1); % envelope of infoflow at each time point
    Emin = min(ConnOrig,[],1); % envelope infoflow at each time point
    
    %    Emax = mean(ConnOrig,1);  % max of data collapsed across time
    %    Emin = Emax;
    
    plotdat = [];
    [dummy idx] = max(abs(ConnOrig));
    for i=1:length(idx)   % CRAPPY CODE!!
        plotdat(i) = ConnOrig(idx(i),i);
    end
    
    %    plot(g.alltimes,Emin, g.alltimes, Emax, 'LineWidth',g.linewidth); hold on;
    plot(g.alltimes,plotdat,'LineWidth',g.linewidth); hold on;
    plot([g.alltimes(1) g.alltimes(length(g.alltimes))],[0 0],'LineWidth',0.7);
    
    % plot bootstrap significance limits (base mean +/-)
    if ~isempty(g.StatsMatrix) && dims(Stat) > 1
        switch dims(Stat)
            case 2, plot(g.alltimes,mean(Stat(:,:),1),'g' ,'LineWidth',g.linewidth);
                plot(g.alltimes,mean(Stat(:,:),1),'k:','LineWidth',g.linewidth);
            case 3, plot(g.alltimes,mean(Stat(:,:,1),1),'g' ,'LineWidth',g.linewidth);
                plot(g.alltimes,mean(Stat(:,:,1),1),'k:','LineWidth',g.linewidth);
                plot(g.alltimes,mean(Stat(:,:,2),1),'g' ,'LineWidth',g.linewidth);
                plot(g.alltimes,mean(Stat(:,:,2),1),'k:','LineWidth',g.linewidth);
        end;
        axis([min(g.alltimes) max(g.alltimes) min([Emin(:)' Stat(:)'])*1.2-10^-5 max([Emax(:)' Stat(:)'])*1.2+10^-5]);
    else
        if ~all(isnan(Emin(:)))
            axis([min(g.alltimes)-eps max(g.alltimes)+eps min(Emin)*1.2-10^-5 max(Emax)*1.2+10^-5]);
        end
    end;
    
    
    tick = get(h(10),'YTick');
    set(h(10),'YTick',[tick(1) ; tick(length(tick))])
    set(h(10),'YAxisLocation','right')
    xlabel('Time (sec)')
    set(gca,'tag','botmarginal');
    
    % Create left marginal
    % ---------------------------
    
    h(11) = axes('Units','Normalized','Position',[0 ordinate(dir) .1 height].*s+q);
    
    % DETERMINE WHAT TO PLOT ON THE LEFT SIDE OF THE IMAGE
    Emax = max(ConnOrig,[],2);  % max of data collapsed across time
    Emin = min(ConnOrig,[],2);
    
    %    Emax = mean(ConnOrig,2);  % max of data collapsed across time
    %    Emin = Emax;
    
    % plot the maximal +- values
    plotdat = [];
    [dummy idx] = max(abs(ConnOrig),[],2);
    for i=1:length(idx)   % CRAPPY CODE!!
        plotdat(i) = ConnOrig(i,idx(i));
    end
    
    if ~strcmpi(g.freqscale, 'log')
        plot(g.allfreqs,plotdat,'b','LineWidth',g.linewidth); % plot baseline
    else
        semilogx(g.allfreqs,Emax,'b','LineWidth',g.linewidth); % plot baseline
        set(h(11),'View',[90 90])
        divs = linspace(log(g.allfreqs(1)), log(g.allfreqs(end)), 10);
        set(gca, 'xtickmode', 'manual');
        divs = ceil(exp(divs)); divs = unique(divs); % ceil is critical here, round might misalign
        % out-of border label with within border ticks
        set(gca, 'xtick', divs);
    end;
    if g.applyThreshold % plot bootstrap significance limits (base max +/-)
        hold on
        if ~strcmpi(g.freqscale, 'log')
            switch dims(Stat)
                case 1, plot(g.allfreqs,Stat(:),'g' ,'LineWidth',g.linewidth);
                    plot(g.allfreqs,Stat(:),'k:','LineWidth',g.linewidth);
                case 2, plot(g.allfreqs,mean(Stat(:,:),2),'g' ,'LineWidth',g.linewidth);
                    plot(g.allfreqs,mean(Stat(:,:),2),'k:','LineWidth',g.linewidth);
                case 3, plot(g.allfreqs,mean(Stat(:,:,1),2),'g' ,'LineWidth',g.linewidth);
                    plot(g.allfreqs,mean(Stat(:,:,1),2),'k:','LineWidth',g.linewidth);
                    plot(g.allfreqs,mean(Stat(:,:,2),2),'g' ,'LineWidth',g.linewidth);
                    plot(g.allfreqs,mean(Stat(:,:,2),2),'k:','LineWidth',g.linewidth);
            end;
        else
            switch dims(Stat)
                case 1, semilogy(g.allfreqs,Stat(:),'g' ,'LineWidth',g.linewidth);
                    semilogy(g.allfreqs,Stat(:),'k:','LineWidth',g.linewidth);
                case 2, semilogy(g.allfreqs,mean(Stat(:,:),2),'g' ,'LineWidth',g.linewidth);
                    semilogy(g.allfreqs,mean(Stat(:,:),2),'k:','LineWidth',g.linewidth);
                case 3, semilogy(g.allfreqs,mean(Stat(:,:,1),2),'g' ,'LineWidth',g.linewidth);
                    semilogy(g.allfreqs,mean(Stat(:,:,1),2),'k:','LineWidth',g.linewidth);
                    semilogy(g.allfreqs,mean(Stat(:,:,2),2),'g' ,'LineWidth',g.linewidth);
                    semilogy(g.allfreqs,mean(Stat(:,:,2),2),'k:','LineWidth',g.linewidth);
            end;
        end;
        if ~isnan(max(Emax))
            axis([g.allfreqs(1) g.allfreqs(end) min([Emin(:)' Stat(:)'])*1.2-10^-5 max([Emax(:)'  Stat(:)'])*1.2+10^-5]);
        end;
    else             % plot marginal mean coherence only
        if ~isnan(max(Emax))
            axis([g.allfreqs(1)-eps g.allfreqs(end)+eps min(Emin)*1.2-10^-5 max(Emax)*1.2+10^-5]);
        end;
    end
    
    
    tick = get(h(11),'YTick');
    set(h(11),'YTick',[tick(1) ; tick(length(tick))]); % crashes for log
    set(h(11),'View',[90 90])
    xlabel('Freq. (Hz)')
    set(gca,'xdir','rev');
    set(gca,'tag','leftmarginal')
end;

leftMargAx = findobj(gcf,'tag','leftmarginal');
botMargAx  = findobj(gcf,'tag','botmarginal');

if length(leftMargAx) > 1
    leftMargAxLim = cell2mat(axis(leftMargAx));
    set(leftMargAx,'Ylim',[min(leftMargAxLim(:,3)) max(leftMargAxLim(:,4))]);
end

if length(botMargAx) > 1
    botMargAxLim = cell2mat(axis(botMargAx));
    set(botMargAx,'Ylim',[min(botMargAxLim(:,3)) max(botMargAxLim(:,4))]);
end

% draw shaded region for baselines
for i=1:length(botMargAx)
    axes(botMargAx(i));
    v = axis(gca);
    if ~isempty(g.baseline)
        % shade in the baseline region
        base = g.alltimes(baseidx);
        patch([base(1) base(1) base(2) base(2)],[v(3) v(4) v(4) v(3)],[0.7 0.7 1],'FaceAlpha',0.5,'EdgeColor','none');
    end

    % plot event markers
    if ~isempty(g.events)
        for ev=1:length(g.events)
            vl = vline(g.events{ev}{1});
            set(vl,'color',g.events{ev}{2},'linestyle',g.events{ev}{3},'linewidth',g.events{ev}{4});
        end
    end
    
end
    
%%%%%%%%%%%%%%% plot topoplot() and/or dipoles %%%%%%%%%%%%%%%%%%%%%%%

switch lower(g.topoplot)
    case 'topoplot'
        h(15) = subplot('Position',[-.1 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(1),g.elocs,'electrodes','off', ...
                'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
        else
            topoplot(g.topovec(1,:),g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
        end;
        dippos = get(h(15),'position');
        lpos = dippos([1 2])+[dippos(3)/2 2*dippos(4)];
        if ~isempty(str2num(g.nodelabels{ch2}))
            ltext = sprintf('IC%s',g.nodelabels{ch2});
        else
            ltext = sprintf('%s',g.nodelabels{ch2});
        end
        text(1.2,0.5,ltext,'horizontalalignment','left','units','normalized','parent',h(15),'fontweight','bold','fontsize',g.axesFontSize);
        axis('square')
        %         axcopy(gca);
        
        h(16) = subplot('Position',[.9 .43 .2 .14].*s+q);
        if size(g.topovec,2) <= 2
            topoplot(g.topovec(2),g.elocs,'electrodes','off', ...
                'style', 'blank', 'emarkersize1chan', 10, 'chaninfo', g.chaninfo);
        else
            topoplot(g.topovec(2,:),g.elocs,'electrodes','off', 'chaninfo', g.chaninfo);
        end;
        dippos = get(h(16),'position');
        lpos = dippos([1 2])+[dippos(3)/2 2*dippos(4)];
        
        if ~isempty(str2num(g.nodelabels{ch1}))
            ltext = sprintf('IC%s',g.nodelabels{ch1});
        else
            ltext = sprintf('%s',g.nodelabels{ch1});
        end
        text(-0.2,0.5,ltext,'horizontalalignment','right','units','normalized','parent',h(16),'fontweight','bold','fontsize',g.axesFontSize);
        axis('square')
        %         axcopy(gca);
    case 'dipole'
        h(15) = subplot('Position',[-.1 .43 .2 .14].*s+q);
        %         dipplot(g.dipfitstruct(1),'color',{'r'},'verbose','off', ...
        %             'dipolelength',0.01,'dipolesize',20,'view',[1 0 0], ...
        %             'projimg', 'off', 'projlines', 'on', 'axistight',  ...
        %             'on', 'cornermri', 'on', 'normlen', 'on','gui','off');
        
        % view [1 0 0] % saggital
        % view [0 -0.99 0.01] for zeynep model
        pop_dipplot(struct('dipfit',g.dipfitstruct),1,'color',{'r'},'verbose','off','dipolelength',0.01,...
            'dipolesize',20,'view',[1 0 0],'projimg', 'off',  ...
            'projlines', 'on', 'axistight', 'on',            ...
            'cornermri', 'on', 'normlen', 'on','gui','off','mri',g.dipplot.mri,'coordformat',g.dipplot.coordformat);
        
        if ~isempty(str2num(g.nodelabels{1}))
            ltext = sprintf('IC%s',g.nodelabels{1});
        else
            ltext = g.nodelabels{1};
        end
        dippos = get(h(15),'position');
        lpos = dippos([1 2])+[dippos(3)/2 2*dippos(4)];
        text(1.2,0.5,ltext,'horizontalalignment','left','units','normalized','parent',h(15),'fontweight','bold','fontsize',g.axesFontSize);
        
        %         ylabel(g.nodelabels{1});
        %         axcopy(gca);
        
        h(16) = subplot('Position',[.9 .43 .2 .14].*s+q);
        %         dipplot(g.dipfitstruct(2),'color',{'r'},'verbose','off', ...
        %             'dipolelength',0.01,'dipolesize',20,'view',[1 0 0], ...
        %             'projimg', 'off', 'projlines', 'on', 'axistight',  ...
        %             'on', 'cornermri', 'on', 'normlen', 'on','gui','off');
        pop_dipplot(struct('dipfit',g.dipfitstruct),2,'color',{'r'},'verbose','off','dipolelength',0.01,...
            'dipolesize',20,'view',[1 0 0],'projimg', 'off',  ...
            'projlines', 'on', 'axistight', 'on',            ...
            'cornermri', 'on', 'normlen', 'on','gui','off','mri',g.dipplot.mri,'coordformat',g.dipplot.coordformat);
        
        if ~isempty(str2num(g.nodelabels{2}))
            ltext = sprintf('IC%s',g.nodelabels{2});
        else
            ltext = g.nodelabels{2};
        end
        %         lpos = [1 .73].*s(1:2)+q(1:2);
        dippos = get(h(16),'position');
        lpos = dippos([1 2])+[dippos(3)/2 2*dippos(4)];
        text(-0.2,0.5,ltext,'horizontalalignment','right','units','normalized','parent',h(16),'fontweight','bold','fontsize',g.axesFontSize);
        %         axcopy(gca);
end


try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;



if (length(g.titleString) > 0) % plot titleString
    axes('Position',pos,'Visible','Off');
    h(13) = text(-.05,1.01,[g.titleString '   ']);
    set(h(13),'VerticalAlignment','bottom')
    set(h(13),'HorizontalAlignment','left')
    set(h(13),'FontSize',g.titleFontSize)
end

try, axcopy(gcf); catch, end;



% HELPER FUNCTIONS

function res = dims(array)
res = min(ndims(array), max(size(array,2),size(array,3)));
