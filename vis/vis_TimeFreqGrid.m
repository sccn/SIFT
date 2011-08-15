
function [figureHandles g] = vis_TimeFreqGrid(varargin)
%
% Create a Time-Frequency Grid from a connectivity matrix. For details on
% the Interactive Time-Frequency Grid see [1].
%
% Inputs:
% 
%       ALLEEG:     Array of EEGLAB datasets
%       Conn:       SIFT Connectivity Structure
%
% Optional:
%
%     Stats:                          A structure containing statistics.
%                                     Input Data Type: structure
%
%     VisualizationMode:              Visualization Modes                                                                                   
%                                     Create Time-Frequency imageplots, Causality x Frequency plots (collapsing across time), Causality x   
%                                     Time plots (collapsing across frequency)                                                              
%                                     Possible values: {'TimeXFrequency','TimeXCausality','FrequencyXCausality'}                            
%                                     Default value  : 'TimeXFrequency'                                                                     
%                                     Input Data Type: string                                                                               
% 
%     MatrixLayout:                   Select the measure and layout                                                                         
%                                     Possible values: {'Full','Partial'}                                                                   
%                                     Default value  : 'Full'                                                                               
%                                     Input Data Type: string                                                                               
%     -------------                                                                                                                         
% 
%         UpperTriangle:              Estimator to render on upper triangle                                                                 
%                                     Possible values: {'none',''}                                                                          
%                                     Default value  : 'n/a'                                                                                
%                                     Input Data Type: string                                                                               
% 
%         LowerTriangle:              Estimator to render on upper triangle                                                                 
%                                     Possible values: {'none',''}                                                                          
%                                     Default value  : 'n/a'                                                                                
%                                     Input Data Type: string                                                                               
% 
%         Diagonal:                   Estimator to render on diagonal                                                                       
%                                     Possible values: {'none',''}                                                                          
%                                     Default value  : 'n/a'                                                                                
%                                     Input Data Type: string                                                                               
% 
%         Estimator:                  Estimator to visualize                                                                                
%                                     Possible values: {''}                                                                                 
%                                     Default value  : 'n/a'                                                                                
%                                     Input Data Type: string                                                                               
% 
%     ColorLimits:                    Color/Y-axis scaling limits                                                                           
%                                     If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If        
%                                     scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is                           
%                                     prctile(abs(Conn),scalar)                                                                             
%                                     Input Data Type: real number (double)                                                                 
% 
%     TimesToPlot:                    [Min Max] Time range to image (sec)                                                                   
%                                     Leave blank to use all timewindows                                                                    
%                                     Input Data Type: real number (double)                                                                 
% 
%     FrequenciesToPlot:              Vector of frequencies (Hz) to image                                                                   
%                                     Leave blank to use all frequencies                                                                    
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     TimeWindowsToPlot:              Time window centers (sec)                                                                             
%                                     If a vector of times, will plot a separate curve for each specified time                              
%                                     Input Data Type: real number (double)                                                                 
% 
%     LineColor:                      Color of line for single-window plots                                                                 
%                                     Input Data Type: real number (double)                                                                 
% 
%     PlotConfidenceIntervals:        Plot confidence intervals (if available)                                                              
%                                     Does not apply to for time-frequency images.                                                          
%                                     Input Data Type: boolean                                                                              
% 
%     PlotContour:                    Plot contours around significant regions                                                              
%                                     Input Data Type: boolean                                                                              
%     ------------                                                                                                                          
% 
%         ContourColor:               Contour Color                                                                                         
%                                     Can use any allowable Matlab color specification (see 'help ColorSpec').                              
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     Thresholding:                   Thresholding options                                                                                  
%                                     You can choose to use statistics (passed in as 'stats' structure), or simple percentile or absolute   
%                                     thresholds.                                                                                           
%                                     Possible values: {'None','Statistics','Simple'}                                                       
%                                     Default value  : 'None'                                                                               
%                                     Input Data Type: string                                                                               
%     -------------                                                                                                                         
% 
%         AlphaSignificance:          P-value threshold for significance. e.g., 0.05 for p<0.05                                             
%                                     Input Range  : [0  1]                                                                                 
%                                     Default value: 0.05                                                                                   
%                                     Input Data Type: real number (double)                                                                 
% 
%         PercentileThreshold:        Percentile threshold                                                                                  
%                                     If of form [percentile, dimension], percentile is applied elementwise across the specified            
%                                     dimension.                                                                                            
%                                     Input Data Type: real number (double)                                                                 
% 
%         AbsoluteThreshold:          Exact threshold                                                                                       
%                                     Input Data Type: real number (double)                                                                 
% 
%     Baseline:                       Time range of baseline [Min Max] (sec)                                                                
%                                     Will subtract baseline from each point. Leave blank for no baseline.                                  
%                                     Input Data Type: real number (double)                                                                 
% 
%     FigureHandles:                  Vector of figure handles to superimpose new graph onto                                                
%                                     New figures and grid will *not* be created. Old grid will be used and new subplots overlaid           
%                                     Input Data Type: real number (double)                                                                 
% 
%     Smooth2D:                       Smooth time-freq image                                                                                
%                                     This will apply nearest-neighbor interpolation.                                                       
%                                     Input Data Type: boolean                                                                                                                                           
% 
%     XTickLabels:                    Labels for X-Tickmarks                                                                                
%                                     Must equal number of time windows                                                                     
%                                     Input Data Type: real number (double)                                                                 
% 
%     YTickLabels:                    Labels for Y-Tickmarks                                                                                
%                                     Must equal number of time windows                                                                     
%                                     Input Data Type: real number (double)                                                                 
% 
%     PlottingOrder:                  Specify index order                                                                                   
%                                     Subset of [1:nbchan] in which to arrange columns/rows. Useful for grouping channels.                  
%                                     Input Data Type: real number (double)                                                                 
% 
%     SourceMarginPlot:               What to plot on margins                                                                               
%                                     Options: 'Topoplot': plot source scalp projection. 'Dipole': plot dipole                              
%                                     Possible values: {'none','topoplot','dipole'}                                                         
%                                     Default value  : 'dipole'                                                                             
%                                     Input Data Type: string                                                                               
% 
%     DipolePlottingOptions:          Options for dipole plotting                                                                           
%                                     Input Data Type: string                                                                               
%     ----------------------                                                                                                                
% 
%         mri:                        Dipplot MRI structure                                                                                 
%                                     Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a    
%                                     path to a Matlab file containing MRI structure. Default uses MNI brain.                               
%                                     Input Data Type: string                                                                               
% 
%         DipoleCoordinateFormat:     Coordinate format for dipplot                                                                         
%                                     Possible values: {'spherical','mni'}                                                                  
%                                     Default value  : 'mni'                                                                                
%                                     Input Data Type: string                                                                               
% 
%         DipplotOptions:             Additional dipplot options                                                                            
%                                     Cell array of <'name',value> pairs of additional options for dipplot (see 'doc dipplot')              
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     NodeLabels:                     List of labels for each node. e.g., {'Node1','Node2',...}                                             
%                                     Leave blank to use defaults.                                                                          
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     FrequencyMarkers:               Vector of frequencies (Hz) at which to draw horizontal lines                                          
%                                     Input Data Type: real number (double)                                                                 
% 
%     FrequencyMarkerColor:           Coloring for frequency markers                                                                        
%                                     If an [1 x 3] array of RBG values, then color all lines using this color. If an [N x 3] matrix of     
%                                     RBG values, then color the kth line with the colorspec from the kth row. If empty then cycle          
%                                     through colorlist                                                                                     
%                                     Input Data Type: real number (double)                                                                 
% 
%     ClusterMaps:                    Cell matrix of mean cluster maps to topoplot                                                          
%                                     Input Data Type: real number (double)                                                                 
% 
%     EventMarkers:                   Event marker time and style                                                                           
%                                     Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. {    
%                                     { 0.2 'y' ':' 2} { 1.5 'r' ':' 2}} will render two dotted-line event makers, yellow at 200 ms and     
%                                     red at 1500 ms                                                                                        
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     FrequencyScale:                 Make the y-scale logarithmic or linear                                                                
%                                     Possible values: {'linear'}                                                                           
%                                     Default value  : 'linear'                                                                             
%                                     Input Data Type: string                                                                               
% 
%     Transform:                      transform the data (logarithmically or other)                                                         
%                                     Possible values: {'log','linear',''}                                                                  
%                                     Default value  : 'n/a'                                                                                
%                                     Input Data Type: string                                                                               
% 
%     TitleString:                    Figure title string                                                                                   
%                                     Input Data Type: string                                                                               
% 
%     TitleFontSize:                  Title Font Size                                                                                       
%                                     Input Data Type: real number (double)                                                                 
% 
%     AxesFontSize:                   Axes Font Size                                                                                        
%                                     Input Data Type: real number (double)                                                                 
% 
%     TextColor:                      Text color                                                                                            
%                                     See 'doc ColorSpec'.                                                                                  
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     Colormap:                       Colormap                                                                                              
%                                     Matlab expression denoting colormap to use (e.g., 'jet(64)'). See 'help colormap'.                    
%                                     Input Data Type: any evaluable Matlab expression.                                                     
% 
%     BackgroundColor:                Background Color                                                                                      
%                                     See 'doc ColorSpec'.                                                                                  
%                                     Input Data Type: any evaluable Matlab
%                                     expression. 
% 
% Outputs:
% 
%       figureHandles:                Handles to figures.      
%
%       g:                            Argument specification. Can be
%                                     supplied in lieu of <name, value>
%                                     argument pairs to exactly reconstruct
%                                     TF-grid.
%
%
% See Also: pop_vis_TimeFreqGrid() 
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


figureHandles = [];

% extract some stuff from inputs for arg defaults
Conn = arg_extract(varargin,'Conn',2);
if ~isempty(Conn)
    Conn = Conn(1);
    ConnNames   = hlp_getConnMethodNames(Conn);
    conndef     = ConnNames{1};
    freqrange   = [Conn.freqs(1) Conn.freqs(end)];
    freqdef     = ['[' num2str(freqrange(1)) ':' num2str(freqrange(end)) ']'];
    timerange   = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    timedef     = [timerange(1) timerange(end)];
    clear Conn;
else
    ConnNames = {''};
    conndef = '';
    [freqrange, freqdef, timerange, timedef] = deal([]);
end


% get some defaults from ALLEEG
ALLEEG = arg_extract(varargin,'ALLEEG',1);
[MyComponentNames MyChannelNames] = deal([]);
if ~isempty(ALLEEG)
    if isfield(ALLEEG(1).CAT,'curComponentNames') && ~isempty(ALLEEG(1).CAT.curComponentNames)
        MyComponentNames = ALLEEG(1).CAT.curComponentNames;
    else
        MyComponentNames = ALLEEG(1).CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    if isfield(ALLEEG(1),'chanlocs') && ~isempty(ALLEEG(1).chanlocs)
        MyChannelNames = {ALLEEG(1).chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:ALLEEG(1).nbchan)'))');
    end
    
    % set default range for marginal plots
    sourceMarginOptions = {'none'};
    if isfield(ALLEEG(1),'chanlocs') && isfield(ALLEEG(1).chanlocs,'X')
        sourceMarginOptions = [sourceMarginOptions 'topoplot'];
    end
    if isfield(ALLEEG(1),'dipfit') && ~isempty(ALLEEG(1).dipfit)
        sourceMarginOptions = [sourceMarginOptions 'dipole'];
    end
    
    clear ALLEEG
else
    sourceMarginOptions = {'none','topoplot','dipole'};
end
    

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
%     methods = fieldnames(stats);
%     StatThreshMethods = fieldnames(stats.methods{1});
end
StatThreshMethods = {'pval'};

clear stats;

% setup the argument list
% -----------------------------------------------------

g = arg_define([0 2],varargin, ...
    arg_norep({'ALLEEG'},mandatory),...
    arg_norep({'Conn'},mandatory),...
    arg_norep({'stats','Stats'},[],[],'A structure containing statistics.'), ...
    arg_nogui({'vismode','VisualizationMode'},'TimeXFrequency',{'TimeXFrequency','TimeXCausality','FrequencyXCausality'},'Visualization Modes. Create Time-Frequency imageplots, Causality x Frequency plots (collapsing across time), Causality x Time plots (collapsing across frequency)'), ...
    arg_norep({'connmethods'},ConnNames{1},ConnNames,'Connectivity Measures to visualize','shape','row'), ...
    arg_norep({'msubset'},'all',{'tril','triu','diag','nodiag','all'},'Subset of the full matrix to keep. Lower/upper triangle (''tril''/''triu''), diagonals (''diag''), everything except diagonal (''nodiag''), everything (''all'').'), ...
    arg_subswitch({'MatrixLayout'},'Full', ...
    {'Full', ...
        { ...
            arg({'estimator','Estimator'},ConnNames{1},ConnNames,'Estimator to visualize','shape','row') ...
        } ...
     'Partial', ...
        {...
            arg({'triu','UpperTriangle'},ConnNames{1},['none' ConnNames],'Estimator to render on upper triangle.','shape','row'), ...
            arg({'tril','LowerTriangle'},ConnNames{1},['none' ConnNames],'Estimator to render on upper triangle.','shape','row'), ...
            arg({'diag','Diagonal'},ConnNames{1},['none' ConnNames],'Estimator to render on diagonal.','shape','row') ...
         }, ...
    },'Select the measure and layout','cat','DisplayProperties'), ...
    arg({'clim','ColorLimits'},100,[],'Color/Y-axis scaling limits. If [min max], scale by [min max]. If scalar, and all(Conn>0), limits are set to [0 maxprc]. If scalar, and any(Conn<0), limits are set to [-maxprc maxprc] where maxprc is prctile(abs(Conn),scalar)','type','denserealdouble','shape','row','cat','DisplayProperties'), ...
    arg({'timeRange','TimesToPlot'},timedef,[],'[Min Max] Time range to image (sec). Leave blank to use all timewindows','shape','row','type','denserealdouble','cat','DisplayProperties'), ...
    arg({'freqValues','FrequenciesToPlot'},freqdef,[],'Vector of frequencies (Hz) to image. Leave blank to use all frequencies','type','expression','shape','row','cat','DisplayProperties'), ...
    arg_nogui({'windows','TimeWindowsToPlot'},[],[],'Time window centers (sec). If a vector of times, will plot a separate curve for each specified time','shape','row','cat','DisplayProperties'),...
    arg_nogui({'linecolor','LineColor'},[0 0 0],[],'Color of line for single-window plots','shape','row','type','denserealdouble'), ...
    arg_nogui({'plotci','PlotConfidenceIntervals'},false,[],'Plot confidence intervals (if available). Does not apply to for time-frequency images.'), ...
    arg_subtoggle({'pcontour','PlotContour'},[], ...
        {...
        arg({'contourcolor','ContourColor'},[0 0 0],[],'Contour Color. Can use any allowable Matlab color specification (see ''help ColorSpec'').','shape','row','type','expression','cat','DisplayProperties') ...
        }, 'Plot contours around significant regions','cat','DisplayProperties'), ...
    arg_subswitch({'thresholding','Thresholding'},'None', ...
        {'None' { ...
            arg_norep({'dummy1'},[],[],'dummy') ...
            }, ...
         'Statistics' {...
            arg({'sigthreshmethod','ThresholdingMethod'},StatThreshMethods{1},StatThreshMethods,'Method to use for significance masking') ...
            arg({'alpha','AlphaSignificance'},0.05,[0 1],'P-value threshold for significance. e.g., 0.05 for p<0.05') ...
            }, ...
         'Simple' {...
            arg({'prcthresh','PercentileThreshold'},100,[],'Percentile threshold. If of form [percentile, dimension], percentile is applied elementwise across the specified dimension.','type','denserealdouble','shape','row','cat','Thresholding'), ...
            arg({'absthresh','AbsoluteThreshold'},[],[],'Exact threshold.','cat','Thresholding') ...
            } ...
       }, 'Thresholding options. You can choose to use statistics (passed in as ''stats'' structure), or simple percentile or absolute thresholds.','cat','Thresholding'), ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point. Leave blank for no baseline.','shape','row','type','denserealdouble','cat','DataProcessing'), ...
    arg_nogui({'fighandles','FigureHandles'},[],[],'Vector of figure handles to superimpose new graph onto. New figures and grid will *not* be created. Old grid will be used and new subplots overlaid'), ...
    arg({'smooth','Smooth2D'},false,[],'Smooth time-freq image. This will apply nearest-neighbor interpolation.','cat','DataProcessing'), ...
    arg_nogui({'xord','XTickLabels'},[],[],'Labels for X-Tickmarks. Must equal number of time windows','cat','DisplayProperties'), ...
    arg_nogui({'yord','YTickLabels'},[],[],'Labels for Y-Tickmarks. Must equal number of time windows','cat','DisplayProperties'), ...
    arg_norep({'channels','VariablesToKeep'},[],[],'List of indices of channels to keep. Can be [vector], a subset of [1:nbchan]'), ...
    arg({'plotorder','PlottingOrder'},[],[],'Specify index order. Subset of [1:nbchan] in which to arrange columns/rows. Useful for grouping channels.','cat','DisplayProperties'), ...
    arg({'topoplot','SourceMarginPlot'},sourceMarginOptions{end},sourceMarginOptions,'What to plot on margins. Options: ''Topoplot'': plot source scalp projection. ''Dipole'': plot dipole','cat','DisplayProperties'), ...
    arg_sub({'dipplot','DipolePlottingOptions'},[], ...
    { ...
        arg_nogui({'mri'},'',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
        arg_nogui({'coordformat','DipoleCoordinateFormat'},'mni',{'spherical','mni'},'Coordinate format for dipplot','type','char','shape','row'), ...
        arg_nogui({'dipplotopt','DipplotOptions'},'{}','','Additional dipplot options. Cell array of <''name'',value> pairs of additional options for dipplot (see ''doc dipplot'')','type','expression','shape','row') ...
    },'Options for dipole plotting'), ...
    arg({'nodelabels','NodeLabels'},MyComponentNames,{},'List of labels for each node. e.g., {''Node1'',''Node2'',...}. Leave blank to use defaults.','shape','row','type','expression','cat','DisplayProperties'),...
    arg({'foilines','FrequencyMarkers'},[],[],'Vector of frequencies (Hz) at which to draw horizontal lines','cat','FrequencyMarkers'), ...
    arg({'foilinecolor','FrequencyMarkerColor'},[],[],'Coloring for frequency markers. If an [1 x 3] array of RBG values, then color all lines using this color. If an [N x 3] matrix of RBG values, then color the kth line with the colorspec from the kth row. If empty then cycle through colorlist','shape','matrix','cat','FrequencyMarkers'), ...
    arg_nogui({'clustmaps','ClusterMaps'},[],[],'Cell matrix of mean cluster maps to topoplot','shape','row','cat','DisplayProperties'), ...
    arg({'events','EventMarkers'},{{0 'r' ':' 2}},[],'Event marker time and style. Specify event markers with a cell array of {time linecolor linestyle linewidth} cell arrays. Ex. { { 0.2 ''y'' '':'' 2} { 1.5 ''r'' '':'' 2}} will render two dotted-line event makers, yellow at 200 ms and red at 1500 ms','type','expression','shape','row','cat','DisplayProperties'), ...
    arg({'freqscale','FrequencyScale'},'linear',{'linear'},'Make the y-scale logarithmic or linear','cat','DisplayProperties'), ...
    arg_nogui({'transform','Transform'},'',{'log','linear',''},'transform the data (logarithmically or other)'), ...
    arg({'titleString','TitleString'},[],[],'Figure title string','type','char','shape','row','cat','TextAndFont'), ...
    arg({'titleFontSize','TitleFontSize'},12,[],'Title Font Size','cat','TextAndFont'), ...
    arg({'axesFontSize','AxesFontSize'},10,[],'Axes Font Size','cat','TextAndFont'), ...
    arg({'textColor','TextColor'},[1 1 1],[],'Text color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg({'colormap','Colormap'},'jet(300)',[],'Colormap. Matlab expression denoting colormap to use (e.g., ''jet(64)''). See ''help colormap''.','type','expression','cat','DisplayProperties'), ...
    arg({'backgroundColor','BackgroundColor'},[0 0 0],[],'Background Color. See ''doc ColorSpec''.','type','expression','shape','row','cat','TextAndFont'), ...
    arg_norep({'report_args'},[],[],'Need this to allow recursive calls...') ...
    );

%     arg_sub({'subplotargs','SubplotExpansionProperties'},[],@vis_TimeFreqCell,'Additional arguments for subplot callback function.','cat','SubplotExpansion'), ...

    
% Commit ALLEEG and Conn variables to workspace   
[data g] = hlp_splitstruct(g,{'ALLEEG','Conn'});        
arg_toworkspace(data);
clear data;


    
% handle plotting multiple estimators on the grid
switch lower(g.MatrixLayout.arg_selection)
    case 'full'
        g.connmethods = g.MatrixLayout.estimator;
        g.msubset     = 'all';
    case 'partial'
        layouts = {g.MatrixLayout.triu g.MatrixLayout.tril g.MatrixLayout.diag};
        if sum(~strcmpi(layouts,'none')) > 1
            % if there is more than one layout to plot...
            if all(strcmpi(layouts,layouts{1}))
                % ... and all layouts use the same conn. estimator
                % then just plot the full TimeFreqGrid
                figureHandles = vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',Conn,g,'MatrixLayout',{'Full','estimator',g.MatrixLayout.triu});
                return;
            elseif ~strcmpi(g.MatrixLayout.diag,'none') ...
                    && isequal(g.MatrixLayout.triu,g.MatrixLayout.tril)
                % same estimator on upper/lower triangles, with different
                % estimator on diagonal
                % ... first plot diagonals...
                g.arg_direct = 0;
                figureHandles = vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',Conn,g, ...
                        'FigureHandles', [g.fighandles], ...
                        'MatrixLayout',{'Partial','tril','none','diag',g.MatrixLayout.diag,'triu','none'});  
                % then continue with off-diagonals
                g.fighandles = [g.fighandles figureHandles];
                g.connmethods = g.MatrixLayout.tril;
                g.msubset = 'nodiag';
                
            else
                % different connectivity estimators for upper, lower, diag
                % so plot each separately via recursive calls to
                % vis_TimeFreqGrid
                if ~strcmpi(g.MatrixLayout.triu,'none')
                    % plot upper triangle
                    figureHandles = vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',Conn,g, ...
                        'FigureHandles', [g.fighandles figureHandles], ...
                        'MatrixLayout',{'Partial','tril','none','diag','none','triu',g.MatrixLayout.triu});
                end
                if ~strcmpi(g.MatrixLayout.tril,'none')
                    % plot upper triangle
                    figureHandles = vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',Conn,g, ...
                        'FigureHandles', [g.fighandles figureHandles], ...
                        'MatrixLayout',{'Partial','triu','none','diag','none','tril',g.MatrixLayout.tril});
                end
                if ~strcmpi(g.MatrixLayout.diag,'none')
                    % plot upper triangle
                    figureHandles = vis_TimeFreqGrid('ALLEEG',ALLEEG,'Conn',Conn,g, ...
                        'FigureHandles', [g.fighandles figureHandles], ...
                        'MatrixLayout',{'Partial','tril','none','triu','none','diag',g.MatrixLayout.diag});
                end

                return;
            end
        elseif sum(~strcmpi(layouts,'none')) == 1
            % there is just one layout to plot
            g.connmethods = layouts{~strcmpi(layouts,'none')};
            allowedLayouts = {'triu','tril','diag'};
            g.msubset = allowedLayouts{~strcmpi(layouts,'none')};
        else
            return;
        end
end
    


if isempty(g.connmethods)
    % if no connectivity methods specified, select all of them
    g.connmethods   = hlp_getConnMethodNames(Conn(1));          end
if isempty(g.freqValues)
    g.freqValues    = freqdef;                                  end
if isempty(g.timeRange)
    g.timeRange     = timedef;                                  end


CEstimator = g.connmethods;


% do some error checking
if ~isfield(Conn(1),'erWinCenterTimes') || isempty(Conn(1).erWinCenterTimes)
    error('Conn.erWinCenterTimes not found!'); end

if isempty(g.clim)
    g.clim = 100; end    

if isempty(g.channels)
    g.channels = 1:ALLEEG(1).CAT.nbchan; end

if ~isempty(g.plotorder)
    for cnd=1:length(Conn)
        sz = size(Conn(cnd).(CEstimator));
        if sz(1)~=length(unique(g.plotorder)) || any(g.plotorder>sz(1))
            fprintf('error: each channel index in connectivity matrix must appear in ''plotorder'' exactly once\n');
            figureHandles = [];
            return;
        end
    end
else
    g.plotorder = 1:size(Conn(1).(CEstimator),1);
end

if ~isfield(Conn(1),'freqs') || isempty(Conn(1).freqs)
    error('''freqs'' must be a field of Conn');
end

% load the MRI file for dipole plotting, if necessary
if strcmpi(g.topoplot,'dipole')
    if isempty(g.dipplot.mri) 
        g.dipplot.mri = ALLEEG(1).dipfit.mrifile;
        tmp = load(g.dipplot.mri);
        fn = fieldnames(tmp);
        g.dipplot.mri = tmp.(fn{1});
    elseif evalin('base',['exist(''' g.dipplot.mri ''',''var'')'])==1
        % MRI variable is in workspace, so copy it
        g.dipplot.mri = evalin('base',g.dipplot.mri);
        if ~isfield(g.dipplot.mri,'anatomy')
            error('MRI structure is invalid format (see ''dipplot'' for more info)');
        end
    elseif isdir(fileparts(g.dipplot.mri))
        % User specified path to MRI file, so load it up
        tmp = load(g.dipplot.mri);
        fn = fieldnames(tmp);
        g.dipplot.mri = tmp.(fn{1});
    else
        % User specified an invalid path to MRI file
        error('Unable to load MRI file for dipplot');
    end
end


% INITIALIZE DEFAULT VARS
%--------------------------------------------------------------------------
gridmargin  = [0.05 0.05];                      % margin (normalized units) around grid of subplots [horiz vert]
pmargin     = 0.005;                            % margin between subplots
OFFSET      = 0.05;
colorlist   = {'k','g','b','c','m','y','r'};    % list of colors for sequential overlapping plots of different time windows
StatsMatrix = [];
TwoSidedThresholding = false;

% setup some labels
if length(Conn)>1
    if length(ALLEEG)<2
        error('To get between-condition labels, ALLEEG must contain datasets for both conditions');
    end
    condstring = sprintf('(%s) - (%s)', ...
    ALLEEG(1).condition, ...
    ALLEEG(2).condition);
else
    condstring = sprintf('(%s)', ...
        ALLEEG(1).condition);
end

% set up the node labels
% ---------------------------------
if isempty(g.nodelabels)
    if isfield(ALLEEG(1).CAT,'curComponentNames') && ~isempty(ALLEEG(1).CAT.curComponentNames)
        g.nodelabels = ALLEEG(1).CAT.curComponentNames;
    else
        g.nodelabels = cell(size(ALLEEG(1).CAT.curComps));
        for i=1:length(g.nodelabels)
            g.nodelabels{i} = fastif(ALLEEG(1).chanlocs(i).labels, ...
                ALLEEG(1).chanlocs(i).labels,num2str(ALLEEG(1).chanlocs(i).urchan)); 
        end
    end
end

% determine whether or not we will plot sources on row-col margins
g.PlotSourceOnMargin = ~strcmpi(g.topoplot,'none');

% extract some variables for convenience
nch   = ALLEEG(1).CAT.nbchan;
erWinCenterTimes = Conn(1).erWinCenterTimes;


% get indices of selected time range
timeIndices                 = getindex(erWinCenterTimes,g.timeRange);
timeIndices                 = timeIndices(1):timeIndices(2);
erWinCenterTimes            = erWinCenterTimes(timeIndices);
g.erWinCenterTimesSelected  = erWinCenterTimes;


% get indices of selected freq range
freqIndices = getindex(Conn(1).freqs,g.freqValues);


% select time and frequency range from Connectivity ...
for i=1:length(Conn)
    Conn(i).(CEstimator) = Conn(i).(CEstimator)(:,:,freqIndices,timeIndices);
end

% ... and from Statistics (if it exists)
if ~isempty(g.stats)
    if ~isequal(size(g.stats.(CEstimator).(g.sigthreshmethod)),size(Conn(1).(CEstimator)))
        error('Stats matrix must be same size as Connectivity matrix');
    end

    g.stats.(CEstimator).(g.sigthreshmethod) ...
        = g.stats.(CEstimator).(g.sigthreshmethod)(:,:,freqIndices,timeIndices);
    
    % select window for confidence interval
    if isfield(g.stats.(CEstimator),'ci')
        g.stats.(CEstimator).ci = g.stats.(CEstimator).ci(:,:,:,freqIndices,timeIndices);
    end
end


% do logarithmic transform if desired
if strcmpi(g.transform,'log')
    for i=1:length(Conn)
        Conn(i).(CEstimator) = log(Conn(i).(CEstimator));
    end
end


% specify new x- and y-axes (TODO: remove this)
if ~isempty(g.xord), erWinCenterTimes = g.xord; end
if ~isempty(g.yord), freqValues = g.yord; end
    


% if there's more than one window, compute the step size
% if length(erWinCenterTimes)>1
%     g.winstep = abs(erWinCenterTimes(2)-erWinCenterTimes(1));
% else
%     g.winstep = 1;
% end


% create the Time Frequency Grid [nch x nch] for current connectivity measure
%--------------------------------------------------------------------------

sz = size(Conn(1).(CEstimator));

if length(sz)==3 && sz(1)~=sz(2)
    % we have a univariate measure (ERSP, mCOH, etc) 
    % so expand univariate onto diagonals
    if ~isempty(g.stats)
        error('stats not currently compatible with univariate connectivity measures');
    end

    for cnd=1:length(Conn)
        tmp(cnd).(CEstimator) = zeros(sz(1),sz(1),sz(2),sz(3),'single');
        for ch=1:sz(1)
            tmp(cnd).(CEstimator)(ch,ch,:,:) = Conn(cnd).(CEstimator)(ch,:,:);
        end
    end
    Conn = tmp; clear tmp;
end


% generate figure and subplot array
numSubplotRows  = nch + g.PlotSourceOnMargin;
numSubplotCols  = nch + g.PlotSourceOnMargin;
g.titleString = sprintf('Subj %s. Cond %s. %s', ...
                        ALLEEG(1).subject, ...
                        condstring,g.titleString);
if ~isempty(g.fighandles)
    % set focus to the selected figure
    figureHandles(end+1)  = figure(g.fighandles);
else
    % create a new figure
    figureHandles(end+1)  = figure('units','normalized','visible','off'); 
    % initialize subplot array
    subplot1(numSubplotRows,numSubplotCols, ...
            'Min',gridmargin,'Max',1-gridmargin,...
            'Gap',[pmargin pmargin], ...
            'YTickL','RightMargin');
    set(figureHandles(end),'name', g.titleString);
end



set(figureHandles(end),'color',g.backgroundColor);
colormap(g.colormap);

if length(Conn)>1
    % we are creating difference plot
    ConnMatrix  = Conn(1).(CEstimator) - Conn(2).(CEstimator);
    TwoSidedThresholding = true;
else
    ConnMatrix  = Conn.(CEstimator);
end

if ~isempty(g.baseline)
    % subtract baseline from connectivity matrix
    ConnMatrix = hlp_rmbaseline(ConnMatrix,g.baseline,erWinCenterTimes);
    TwoSidedThresholding = true;
end


if ~strcmpi(g.msubset,'all')
    % only plot a subset of the full TF Grid matrix
    switch lower(g.msubset)
        case 'nodiag' 
            % NaN diagonal elements of nch x nch grid
            Dd = single(~eye(nch));
            Dd(Dd==0)=nan;
        case 'diag'
            % NaN off-diagonal elements of nch x nch grid
            Dd = eye(nch);
            Dd(Dd==0)=nan;
        case 'tril'
            % NaN everything except lower triangle of nch x nch grid
            Dd = tril(ones(nch),-1);
            Dd(Dd==0)=nan;
        case 'triu'
            % NaN everything except upper triangle of nch x nch grid
            Dd = triu(ones(nch),1);
            Dd(Dd==0)=nan;
    end
    ConnMatrix=ConnMatrix.*(repmat(Dd,[1,1,size(ConnMatrix,3),size(ConnMatrix,4)]));
end


% Apply statistics and thresholding
%--------------------------------------------------------------------
if any(ConnMatrix(:)<0), TwoSidedThresholding = true; end

switch lower(g.thresholding.arg_selection)
case 'simple'
    if ~isempty(g.thresholding.prcthresh)
        % percentile Thresholding

        if length(g.thresholding.prcthresh)>1
            % get percentiles across a specified dimension
            dim=g.thresholding.prcthresh(2);
            sz=size(ConnMatrix); nd=length(sz);
            odims = setdiff(1:nd,dim);
            Cntmp = permute(ConnMatrix,[odims(1) dim odims(2:end)]);    % put dim in second dim (for reshape)
            Cntmp = reshape(Cntmp,[prod(sz(odims)) sz(dim)]);           % we'll take percentiles for ea. col

            if TwoSidedThresholding
                % apply two-sided thresholding
                StatsMatrix(2,:)= prctile(Cntmp,g.thresholding.prcthresh(1),1);             % upper limit
                StatsMatrix(1,:)= prctile(Cntmp,100-g.thresholding.prcthresh(1),1);         % lower limit

                % expand StatsMatrix to same size as ConnMatrix
                SMtmp(2,:,:,:,:,:) = repmat(StatsMatrix(2,:),[sz(odims(1)) 1 sz(odims(2:end))]);
                SMtmp(1,:,:,:,:,:) = repmat(StatsMatrix(1,:),[sz(odims(1)) 1 sz(odims(2:end))]);
                StatsMatrix = ipermute(SMtmp,[1 1+[odims(1) dim odims(2:end)]]);   
                clear SMtmp
            else
                % apply single-sided thresholding
                StatsMatrix = prctile(Cntmp,g.thresholding.prcthresh(1),1);

                % expand StatsMatrix to same size as ConnMatrix
                StatsMatrix = repmat(StatsMatrix,[sz(odims(1)) 1 sz(odims(2:end))]);
                StatsMatrix = ipermute(StatsMatrix,[odims(1) dim odims(2:end)]);
            end

            clear Cntmp

        else
            % get percentiles of complete data matrix
            if TwoSidedThresholding
                StatsMatrix(2,:)= prctile(ConnMatrix(:),g.thresholding.prcthresh(1),1);     % upper limit
                StatsMatrix(1,:)= prctile(ConnMatrix(:),100-g.thresholding.prcthresh(1),1); % lower limit
            else
                StatsMatrix = prctile(ConnMatrix(:),g.thresholding.prcthresh);
            end
        end
    elseif ~isempty(g.absthresh)
        % use scalar threshold
        StatsMatrix = g.absthresh(size(ConnMatrix));    
    end

case 'statistics'

    % Use Statistics Structure for thresholding
    if nargin>3 && isfield(g.stats,CEstimator)
        if length(g.stats)>1
            fprintf('WARNING: Stats contains more than one structure. Taking the first one...\n');
            g.stats = g.stats(1);
        end
        if isstruct(g.stats.(CEstimator))
            if ~isfield(g.stats.(CEstimator),g.sigthreshmethod)
                fprintf('ERROR: %s is not a field of g.stats.%s\n',g.sigthreshmethod,CEstimator);
                return;
            end
            StatsMatrix = g.stats.(CEstimator).(g.sigthreshmethod);
        else
            StatsMatrix = g.stats.(CEstimator);
        end
    else
        StatsMatrix = [];
    end

end

    
% Apply significance mask
% ---------------------------------------------------------------

% preserve the original connectivity matrix
OrigConnMatrix = ConnMatrix;

if ~isempty(StatsMatrix) && isempty(g.windows) && ~g.pcontour.arg_selection

    if TwoSidedThresholding % two-sided thresholds (x < lothresh | x > hithresh = nan)
        if isscalar(StatsMatrix) % && ~g.pcontour.arg_selection
            % uniform two-sided threshold
            ConnMatrix(abs(ConnMatrix) < abs(StatsMatrix)) = 0;
        elseif length(StatsMatrix)==2  %&& ~g.pcontour.arg_selection
            % uniform two-sided threshold
            ConnMatrix(ConnMatrix > StatsMatrix(1) & ConnMatrix < StatsMatrix(2)) = 0;
        elseif isequal(size(StatsMatrix),size(ConnMatrix)) && islogical(StatsMatrix)
            % logical thresholding (e.g., p-value)
            ConnMatrix(~StatsMatrix) = 0;    
        elseif isequal(size(StatsMatrix),[2 size(ConnMatrix)])
            % two-sided numeric thresholding
            ConnMatrix(ConnMatrix > squeeze(StatsMatrix(1,:,:,:,:))  ...
                     & ConnMatrix < squeeze(StatsMatrix(2,:,:,:,:))) = 0;
        else
            error('unknown statistical thresholding paradigm');
        end

    else % single-sided thresholding (x < thresh = 0)
        if isscalar(StatsMatrix) % && ~g.pcontour.arg_selection
            ConnMatrix(ConnMatrix < StatsMatrix) = 0;
        elseif isvector(StatsMatrix) && ndims(ConnMatrix)>3
%                 sz=size(StatsMatrix);
%                 % expand vector thresh to dims of ConnMatrix
%                 StatsMatrix = repmat(StatsMatrix,fastif(sz(1)==1,[length(freqs) 1],[1 length(erWinCenterTimes)]));
        elseif isequal(size(StatsMatrix),size(ConnMatrix))
            if islogical(StatsMatrix)
                % logical thresholding (e.g., p-value)
                ConnMatrix(~StatsMatrix) = 0;
            else
                % single-sided numeric thresholding
                ConnMatrix(ConnMatrix < StatsMatrix) = 0;
            end
        else
            error('unknown statistical thresholding paradigm');
        end
    end

end
% ---------------------------------------------------------------


if ~isempty(g.windows)
    % Instead of Time-Frequency plots, user wishes to plot causal
    % spectra for individual selected window(s)
    g.windows = unique(g.windows);
    windowIndex = getindex(erWinCenterTimes,g.windows);
    ConnMatrix=ConnMatrix(:,:,:,windowIndex);
    if ~isempty(StatsMatrix) && ~isscalar(StatsMatrix)
        StatsMatrix = StatsMatrix(:,:,:,windowIndex);
    end

    OrigConnMatrix = OrigConnMatrix(:,:,:,windowIndex);
    
    % select window for confidence interval
    if ~isempty(g.stats) && isfield(g.stats.(CEstimator),'ci') && g.plotci
        g.stats.(CEstimator).ci = g.stats.(CEstimator).ci(:,:,:,:,windowIndex);
    end
    
end
% ---------------------------------------------------------------------


[nch nch nfreqs ntime] = size(ConnMatrix);


% set up the color limits
% ---------------------------------
if ~isempty(g.clim)
    if isscalar(g.clim)
        % use percentile colorlimits
        if ~TwoSidedThresholding %all(ConnMatrix(~isnan(ConnMatrix))>=0)
            g.clim=[0 prctile(ConnMatrix(:),g.clim)];
        else
            maxprc=prctile(abs(ConnMatrix(:)),g.clim);
            g.clim=[-maxprc maxprc];
        end
    end
end

if any(isnan(g.clim)) || diff(g.clim)<=0
    g.clim = [0 1];  % clims are nan or non-increasing
end

%     nodelabels = cell(size(ALLEEG(1).CAT.params.newchans));
%     for i=1:length(ALLEEG(1).chanlocs), nodelabels{i} = ALLEEG(1).chanlocs(i).labels; end


% plot the source locations or topoplots on marginal
% ---------------------------------
if ~strcmpi(g.topoplot,'none')

    subplot1(sub2ind([numSubplotRows,numSubplotCols],1,1));
    
    set(gca,'visible','off');
    chidx=0;
    for ch=g.plotorder
        chidx = chidx+1;
        
        % row marginals
        % --------------
        subplot1(sub2ind([numSubplotRows,numSubplotCols],1,chidx+1));
        plotmarginal(ALLEEG,ch,g,'view',[1 0 0]); % [0 -1 0] for zeynep
        pos = get(gca,'position');
        th=annotation('textbox',[pos(1)-0.01 pos(2) 0.01 pos(4)]);
        set(th,'string',g.nodelabels(ch),'color',g.textColor,  ...
            'horizontalalignment','center','fontsize',g.axesFontSize, ...
            'verticalalignment','middle','edgecolor','none');
        
        % column marginals
        % ----------------
        subplot1(sub2ind([numSubplotRows,numSubplotCols],chidx+1,1));
        plotmarginal(ALLEEG,ch,g,'view',[0 0 1])  % [0 0 1] for zeynep
        title(g.nodelabels(ch),'color',g.textColor,'fontsize',g.axesFontSize);
        set(gco,'tag','coltitle');
    end
end

% Plot the causality: 
% column ch_j --> row ch_i
for ch_i=1:nch
    for ch_j=1:nch

        % plot in specified order
        i = g.plotorder(ch_i);
        j = g.plotorder(ch_j);

        % Index the appropriate subplot. 
        % subplot1 counts linearly by columns, so we flip ch_i,ch_j
        subplot1(sub2ind([numSubplotRows,numSubplotCols],ch_j+numSubplotRows-nch,ch_i+numSubplotCols-nch));
        set(gca,'tag','causalplot');

        % create column/row titles, if necessary
        if ch_i==1 && strcmpi(g.topoplot,'none')
            ht=title(g.nodelabels(j));
            set(ht,'tag','title','color',g.textColor, ...
                   'fontsize',g.axesFontSize);
        end
        if ch_j==1 && strcmpi(g.topoplot,'none')
            pos = get(gca,'position');
            th=annotation('textbox',[pos(1)-0.01 pos(2) 0.01 pos(4)]);
            set(th,'string',g.nodelabels(i),'color',g.textColor,  ...
                'horizontalalignment','center','fontsize',g.axesFontSize, ...
                'verticalalignment','middle','edgecolor','none');
        end
        
        % check if we want to image this cell
        if (ch_i==ch_j    &&  strcmpi(g.msubset,'nodiag'))    || ...
           (ch_i~=ch_j    &&  strcmpi(g.msubset,'diag'))      || ...
           (ch_i>=ch_j    &&  strcmpi(g.msubset,'triu'))      || ...
           (ch_i<=ch_j    &&  strcmpi(g.msubset,'tril'))         ...
            
           % if we get here, then we don't want to actually image this cell
           set(gca,'color',get(gcf,'color'));
           continue; 
        end


        % MODE1: Time x Frequency Images
        % -----------------------------------------------------
        if ntime > 1 && nfreqs > 1 && isempty(g.windows) 

            C=squeeze(ConnMatrix(i,j,:,:));

            if g.smooth
                h=pcolor(erWinCenterTimes,g.freqValues,double(C));
                shading interp
            else
                if strcmpi(g.freqscale,'log')
%                     options = { 'chanlocs', struct([]), 'electrodes', 'off', 'cbar', 'off', ...
%                     'cmode', 'separate', 'logfreq','on' };       
%         
%                     tftopo( C, erWinCenterTimes, g.freqValues,options{:}); 
% %                     imagesclogy(erWinCenterTimes,g.freqValues,C);
                    h=imagesc(erWinCenterTimes,g.freqValues,C);
%                     h=gca;
                else                        
                    h=imagesc(erWinCenterTimes,g.freqValues,C);
                end
            end
            
            set(gca,'Clim',g.clim,'YDir','normal');
            
            % extract the stats matrix for this pair
            if isequal(size(StatsMatrix),size(ConnMatrix))
                Sji = squeeze(StatsMatrix(i,j,:,:));
                Sij = squeeze(StatsMatrix(j,i,:,:));
            elseif size(StatsMatrix,1)==2 && ndims(StatsMatrix)==5
                Sji = permute(squeeze(squeeze(StatsMatrix(:,i,j,:,:))),[2 3 1]);
                Sij = permute(squeeze(squeeze(StatsMatrix(:,j,i,:,:))),[2 3 1]);
            else
                Sji = StatsMatrix;
                Sij = StatsMatrix;
            end
                
            % Plot contour
            if g.pcontour.arg_selection 
                if isscalar(Sji) && any(C(:)-C(1))   
                    % use contour for constant threshold
                    hold on; 
                    contour(erWinCenterTimes,g.freqValues,C,[Sji Sji], ...
                        'color',g.pcontour.contourcolor); 
                    hold off 
                elseif length(Sji)==2 && any(C(:)-C(1)) 
                    hold on; 
                    contour(erWinCenterTimes,g.freqValues,C,Sji, ...
                        'color',g.pcontour.contourcolor); 
                    hold off 
                end
            end
         
        
            
            % prepare the subplot function
            subargs.topovec     = squeeze(ALLEEG(1).icawinv(:,ALLEEG(1).CAT.curComps([j i])))';
            if ~isempty(ALLEEG(1).dipfit)
                subargs.dipfitstruct = ALLEEG(1).dipfit;
                subargs.dipfitstruct.model = subargs.dipfitstruct.model(ALLEEG(1).CAT.curComps([j i]));
            else
                subargs.dipfitstruct = [];
            end
            subargs.elocs       = ALLEEG(1).chanlocs;
            subargs.chaninfo    = ALLEEG(1).chaninfo;
            subargs.alltimes    = erWinCenterTimes;
            subargs.allfreqs    = g.freqValues;
            
            if ~isempty(Sji)
                    subargs.StatsMatrix(1,:,:,:,:) = Sji;
                    subargs.StatsMatrix(2,:,:,:,:) = Sij;
            else
                subargs.StatsMatrix = [];
            end
            
            subargs.ConnMatrix(1,:,:)  = squeeze(OrigConnMatrix(i,j,:,:));
            subargs.ConnMatrix(2,:,:)  = squeeze(OrigConnMatrix(j,i,:,:));
            subargs.baseline    = g.baseline;
            subargs.freqscale   = g.freqscale;
            subargs.events      = g.events;
            subargs.topoplot    = g.topoplot;
            subargs.titleString = g.titleString;
            subargs.titleFontSize   = g.titleFontSize;
            subargs.axesFontSize    = g.axesFontSize;
            subargs.textColor       = g.textColor;
            subargs.backgroundColor = g.backgroundColor;
            subargs.clim            = g.clim;
            subargs.thresholding    = g.thresholding;
            subargs.bidir           = true;
            subargs.connmethod      = CEstimator;
            subargs.nodelabels      = g.nodelabels([j i]);
            subargs.dipplot         = g.dipplot;
            
%             subplotargs = hlp_varargin2struct({g.subplotargs subargs});
           
            set(gca,'userdata',subargs)
            set([gca h],'buttondownfcn','vis_TimeFreqCell(get(gca,''UserData''));');
%             set([gca h],'tooltip',sprintf('%s --> %s. Click to expand',g.nodelabels{j},g.nodelabels{i}));
            
            set(gca,'Xlim',[erWinCenterTimes(1)+OFFSET erWinCenterTimes(end)-OFFSET]);  
            % [erWinCenterTimes(1)-winlen/(2*ALLEEG(1).srate) erWinCenterTimes(end)+winlen/(2*ALLEEG(1).srate)]
            set(gca,'Ylim',g.freqValues([1 end]));
            
            set(gca,'XColor',g.textColor,'YColor',g.textColor);
            set(gca,'fontsize',g.axesFontSize);
            

            % draw event markers
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

            % draw horizontal lines at frequencies of interest
            if ~isempty(g.foilines)
                for ln=1:length(g.foilines)
                    hl = hline(g.foilines(ln));
                    if isempty(g.foilinecolor)
                        color = colorlist{mod(ln-1,length(colorlist))+1};
                    elseif size(g.foilinecolor,1) > 1
                        color = g.foilinecolor(ln,:);
                    elseif size(g.foilinecolor,1) == 1
                        color = g.foilinecolor;
                    end
                    
                    
                    set(hl,'color',color,'linestyle','-','linewidth',1);
                    set(hl,'tag','foilines');
                end
            end
            
            % create a red border around diagonal plots
            if ch_i==ch_j
%                 subplot1(sub2ind([numSubplotRows,numSubplotCols],ch_j+numSubplotRows-nch,ch_i+numSubplotCols-nch));
                pos = get(gca,'position'); 
                hborder = annotation('rectangle',pos,'edgecolor',[1 0 0],'linewidth',2);
                set(hborder,'userdata',gca);
                set(hborder,'buttondownfcn','vis_TimeFreqCell(get(get(gco,''UserData''),''UserData''));');
            end
            
         
        % MODE 2: Causality x Frequency lineplots
        % -----------------------------------------------------
        elseif nfreqs > 1   

            if isequal(size(StatsMatrix),size(ConnMatrix))
                S = squeeze(StatsMatrix(i,j,:,:));
            else
                S = StatsMatrix; 
            end

            hold on
            for tt=1:ntime
                
                % plot confidence intervals
                if ~isempty(g.stats) && isfield(g.stats.(CEstimator),'ci') && g.plotci
                   ci = g.stats.(CEstimator).ci;
                   if ndims(ci)>=4 && size(ci,1)==2
%                            if j==3 && i==2, keyboard; end
%                            plot(g.freqValues,squeeze(ci(1,i,j,:,tt)),'r');
%                            plot(g.freqValues,squeeze(ci(2,i,j,:,tt)),'b');
                       ciplot(squeeze(ci(1,i,j,:,tt)),squeeze(ci(2,i,j,:,tt)),g.freqValues,[0.7 0.7 0.7],0,'Ylim',g.clim,'FaceAlpha',0.5,'EdgeColor',[0.2 0.2 0.2]);
                   else
%                            ciplot(-squeeze(ci(i,j,:,tt)),squeeze(ci(i,j,:,tt)),g.freqValues,[0.7 0.7 0.7],1,'FaceAlpha',0.5,'EdgeColor',[0.2 0.2 0.2]);
                   end
                end
                
                % plot causality trace
                if ntime==1
                    h=plot(g.freqValues,squeeze(ConnMatrix(i,j,:,:)),'color',g.linecolor);
                    set(gca,'color',1-g.linecolor);
%                         h=area(g.freqValues,squeeze(ConnMatrix(i,j,:,:)));
                else
                    h=plot(g.freqValues,squeeze(ConnMatrix(i,j,:,tt)),'color',colorlist{mod(tt-1,length(colorlist))+1});
                    set(gca,'color',1-g.linecolor);
                end


            end


            if ~isempty(S)
                if isscalar(S)
                    plot(g.freqValues,S(ones(1,length(g.freqValues))),'r:');
                elseif isvector(S)
                    plot(g.freqValues,S,'r:');
                end
            end
            hold off

            set(gca,'Ylim',g.clim);
            set(gca,'Xlim',[g.freqValues(1) g.freqValues(end)]);
            set(gca,'tag','lineplot');
%             set(gca,'userdata',struct('dat',squeeze(ConnMatrix(i,j,:,:)),'freqs',g.freqValues,'from',g.nodelabels(j),'to',g.nodelabels(i)));
%             try, axcopy(gca); catch; end
%                 set([gca h],'buttondownfcn',cbfun);

        % MODE 3: Causality x Time (for single freq) lineplots
        % -----------------------------------------------------
        elseif ntime > 1   
%                 if isequal(size(StatsMatrix),size(ConnMatrix)), S = squeeze(StatsMatrix(i,j,:,:));
%                 else S = StatsMatrix; end
            S = [];

            hold on
            if nfreqs==1
                h=plot(erWinCenterTimes,squeeze(ConnMatrix(i,j,:,:)));
            else
                %h=plot(erWinCenterTimes,squeeze(ConnMatrix(i,j,:,:)),colorlist{mod(tt-1,length(colorlist))+1});
            end

%                 % plot threshold (if available)
%                 if ~isempty(S)
%                     if isscalar(S)
%                         plot(erWinCenterTimes,S(ones(1,length(erWinCenterTimes))),'r:');
%                     elseif isvector(S)
%                         plot(erWinCenterTimes,S,'r:');
%                     end
%                 end
%                 hold off
%                 
%                 set(gca,'Ylim',clim);
%                 set(gca,'Xlim',[erWinCenterTimes(1) erWinCenterTimes(end)]);
% %                 set(gca,'userdata',struct('dat',squeeze(ConnMatrix(i,j,g.freqValues,:)),'freqs',g.freqValues,'from',g.nodelabels(j),'to',g.nodelabels(i)));
                set(gca,'tag','lineplot');
%                 try, axcopy(gca); catch; end
        end
      
    end
end


set(gcf,'color',g.backgroundColor);

for ch_i=1:nch
    subplot1(sub2ind([numSubplotRows,numSubplotCols],ch_i+numSubplotRows-nch,ch_i+numSubplotCols-nch));
    % give diagonals a red border
    set(gca,'color','r');
    set(gca,'xcolor','r');
    set(gca,'ycolor','r');
end
           

% axcopy all lineplots
lineplots = findobj(gcf,'tag','lineplot');
for i=1:length(lineplots)
    axcopy(lineplots(i));
end

% axcopy all topoplots
topos = findobj(gcf,'tag','topo');
for i=1:length(topos)
    axcopy(topos(i));
end
    
% create legend
subplot1(1);
% axes(gca);
if isempty(findobj(gcf,'tag','legendborder'))
    % create border around legend
    pos = get(gca,'position'); 
    annotation('rectangle',pos,'edgecolor',[1 0 0],'tag','legendborder');
end

legendsettings = {'horizontalalignment','left','units','normalized','parent',gca,'edgecolor','none','color',g.textColor,'fontsize',g.axesFontSize};
if strcmpi(g.MatrixLayout.arg_selection,'full')
    text(0.1,0.5, CEstimator,legendsettings{:})
else
    upperstr = g.MatrixLayout.triu;
    diagstr  = g.MatrixLayout.diag;
    lowerstr = g.MatrixLayout.tril;
    leghandle = [];
    if ~strcmpi(upperstr,'none')
        leghandle = [leghandle text(0.1,0.85, sprintf('upper: %s',upperstr),legendsettings{:})];
    end
    if ~strcmpi(diagstr,'none')
        leghandle = [leghandle text(0.1,0.5, sprintf('diag: %s',diagstr),legendsettings{:})];
    end
    if ~strcmpi(lowerstr,'none')
        leghandle = [leghandle text(0.1,0.15, sprintf('lower: %s',lowerstr),legendsettings{:})];
    end
end
            
% text(0.1,0.9,sprintf('upper: dDTF08\ndiag: S\nlower: dDTF08'),'horizontalalignment','left','units','normalized','parent',gca,'edgecolor','w','color','w');



function h = plotmarginal(ALLEEG,curch,g,varargin)
% render topoplot or dipplot in the current axis
% varargin can be a list of <'name',value> argument pairs for pop_dipplot

dipplotdefs ={'color',{'r'},'verbose','off','dipolelength',0.01,...
            'dipolesize',20,'view',[1 0 0],'projimg', 'off',  ...
            'projlines', 'on', 'axistight', 'on',            ...
            'cornermri', 'on', 'normlen', 'on','gui','off','mri',g.dipplot.mri,'coordformat',g.dipplot.coordformat};
        
dipplotargs = [hlp_struct2varargin(hlp_varargin2struct(varargin,dipplotdefs{:})) g.dipplot.dipplotopt];
    
if strcmpi(g.topoplot,'dipole') && ~strcmpi(g.msubset,'diag')
    % plot dipole locations
    cbstr=pop_dipplot(ALLEEG(1),ALLEEG(1).CAT.curComps(curch), ...
        dipplotargs{:});
    set(gca,'buttondownfcn',['figure;' cbstr]);
end

if strcmpi(g.topoplot,'topoplot') && ~strcmpi(g.msubset,'diag')
    % plot topoplots
    color = get(gcf,'color');
    topoplot(squeeze(ALLEEG(1).icawinv(:,ALLEEG(1).CAT.curComps(curch))), ...
             ALLEEG(1).chanlocs,'electrodes','off');
    set(gcf,'color',color)
    set(findobj(gca,'type','patch'),'facecolor',color);
%                 ylabel(g.nodelabels(curch),'color','w');
%     cbstr = 'topoplot(squeeze(ALLEEG(1).icawinv(:,ALLEEG(1).CAT.curComps(curch))),fastif(isfield(ALLEEG(1),''chanlocs''),ALLEEG(1).chanlocs,ALLEEG(1).chanlocs),''electrodes'',''off'');';
end

if strcmpi(g.topoplot,'clustmaps') && ~strcmpi(g.msubset,'diag')
   % plot mean cluster topoplots
   toporeplot(g.clustmaps{curch},'plotrad',.75,'intrad',.75);
end

set(gca,'tag','topo');

        