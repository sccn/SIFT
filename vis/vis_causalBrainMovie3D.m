
function [cfg handles BMout] = vis_causalBrainMovie3D(varargin)
%
% Create an interactive 3D BrainMovie from a connectivity matrix. See [1]
% for more details on the Interactive BrainMovie3D.
%
% Inputs:
% 
%       ALLEEG:     Array of EEGLAB datasets
%       Conn:       SIFT Connectivity Structure
%
% Optional:
%
%     Stats:                        Name of variable in base containing statistics                                                                    
%                                   Input Data Type: string                                                                                           
% 
%     ConnectivityMethod:           Connectivity Measure to visualize                                                                                 
%                                   Possible values: {'Determined_by_data'}                                                                           
%                                   Default value  : 'Determined_by_data'                                                                             
%                                   Input Data Type: string                                                                                           
% 
%     MovieTimeRange:               Time Range for Movie [Min Max] (sec)                                                                              
%                                   Specified w.r.t to event time (e.g., [-1 2]). Leave blank for complete epoch.                                     
%                                   Input Data Type: real number (double)                                                                             
% 
%     FrequenciesToCollapse:        Frequencies over which to collapse (Hz)                                                                           
%                                   E.g., [1:50] for 1-50Hz. Leave blank for all frequencies                                                          
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%     FreqCollapseMethod:           Method for collapsing frequency dimension                                                                         
%                                   Possible values: {'mean','max','peak','integrate'}                                                                
%                                   Default value  : 'mean'                                                                                           
%                                   Input Data Type: string                                                                                           
% 
%     TimeResamplingFactor:         Time resampling factor                                                                                            
%                                   If 0, don't resample. If < 1, downsample timecourse by this factor. If > 1, upsample by this                      
%                                   factor. Uses resample() from Sigproc Toolbox                                                                      
%                                   Input Range  : [0  20]                                                                                            
%                                   Default value: 0                                                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%     SubtractConditions:           Subtract conditions                                                                                               
%                                   If true, then plot difference between conditions. If false, then render two brainmovies                           
%                                   side-by-side.                                                                                                     
%                                   Input Data Type: boolean                                                                                          
% 
%     NodeLabels:                   List of labels for each node. e.g., {'Node1','Node2',...}                                                         
%                                   Leave blank to omit labels.                                                                                       
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%     NodesToExclude:               Exclude these sources from Brainmovie                                                                             
%                                   Specify using the Name/ID of the source to exclude.                                                               
%                                   Input Data Type: boolean                                                                                          
% 
%     EdgeColorMapping:             Specify mapping for edge color                                                                                    
%                                   This determines how we index into the colormap. If 'None', edge color is not modulated. If                        
%                                   'Connectivity', use connectivity strength. If 'PeakFreq', use index of peak frequency                             
%                                   Possible values: {'None','Connectivity','PeakFreq','Directionality'}                                              
%                                   Default value  : 'Connectivity'                                                                                   
%                                   Input Data Type: string                                                                                           
% 
%     EdgeSizeMapping:              Specify mapping for edge size                                                                                     
%                                   If 'None', edges are not rendered. If 'Connectivity', use connectivity strength. If                               
%                                   'ConnMagnitude', use connectivity magnitude (absval). If 'PeakFreq', use index of peak frequency.                 
%                                   If 'Directionality', map directionality to the lower and upper extremes of the colormap (e.g.,                    
%                                   i->j: blue, j->i: red)                                                                                            
%                                   Possible values: {'None','ConnMagnitude','Connectivity'}                                                          
%                                   Default value  : 'ConnMagnitude'                                                                                  
%                                   Input Data Type: string                                                                                           
% 
%     NodeColorMapping:             Specify mapping for node color                                                                                    
%                                   This determines how we index into the colormap. Options are as follows. None: node color is not                   
%                                   modulated. Outflow: sum connectivity strengths over outgoing edges. Inflow: sum connectivity                      
%                                   strengths over incoming edges. CausalFlow: Outflow-Inflow. Asymmetry Ratio: node colors are defined               
%                                   by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1                  
%                                   for exclusive outflow, and 0.5 for balanced inflow/outflow                                                        
%                                   Possible values: {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'}  
%                                   Default value  : 'Outflow'                                                                                        
%                                   Input Data Type: string                                                                                           
% 
%     NodeSizeMapping:              Specify mapping for node size. Options are as follows:                                                            
%                                   None: node size is not modulated.                                                                                 
%                                   Outflow: sum connectivity strengths over outgoing edges.                                                          
%                                   Inflow: sum connectivity strengths over incoming edges.                                                           
%                                   CausalFlow: Outflow-Inflow.                                                                                       
%                                   Asymmetry Ratio: node size is defined by the equation C = 0.5*(1 +                                                
%                                   outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5                
%                                   for balanced inflow/outflow                                                                                       
%                                   Possible values: {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'}  
%                                   Default value  : 'Outflow'                                                                                        
%                                   Input Data Type: string                                                                                           
% 
%     Baseline:                     Time range of baseline [Min Max] (sec)                                                                            
%                                   Will subtract baseline from each point. Leave blank for no baseline.                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%     NormalizeConn:                Normalize edge and node values to [0 1]                                                                           
%                                   Values mapped to edge/node width and color are devided by max to put in [0 1] range. Recommended!                 
%                                   Input Data Type: boolean                                                                                          
% 
%     UseStatistics:                Use Statistics                                                                                                    
%                                   Input Data Type: boolean                                                                                          
%     --------------                                                                                                                                  
% 
%         Thresholding:             Type of thresholding for stats                                                                                    
%                                   If 'both' then stats should have upper and lower thresholds                                                       
%                                   Possible values: {'single','both','lessthan'}                                                                     
%                                   Default value  : 'single'                                                                                         
%                                   Input Data Type: string                                                                                           
% 
%         AlphaSignificance:        Significance threshold. e.g., 0.05 for p<0.05                                                                     
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.05                                                                                               
%                                   Input Data Type: real number (double)                                                                             
% 
%     PercentileThreshold:          Percentile threshold                                                                                              
%                                   Fraction of "strongest" connections to display. E.g: PercentileThreshold=0.05 will display only the               
%                                   top 5% of connections                                                                                             
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: n/a                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%     AbsoluteThreshold:            Exact threshold                                                                                                   
%                                   If a single value, then render only connections with strength above this threshold. If [low hi]                   
%                                   then render only connections with strength between [lo hi]. Overrides PercentileThreshold                         
%                                   Input Data Type: real number (double)                                                                             
% 
%     FooterPanelDisplaySpec:       Configure footer panel displayed at the bottom of the figure                                                      
%                                   If 'off', don't render footer. If 'ICA_ERP_Envelope', then display the ERP envelope of                            
%                                   backprojected components. If 'Chan_ERP_Envelope' then display the ERP envelope of selected channels               
%                                   Possible values: {'off','ICA_ERPenvelope','Chan_ERPenvelope'}                                                     
%                                   Default value  : 'off'                                                                                            
%                                   Input Data Type: string                                                                                           
%     -----------------------                                                                                                                         
% 
%         icaenvelopevars:          Select components to use in the display                                                                           
%                                   Input Data Type: boolean                                                                                          
% 
%         backprojectedchans:       List of channels to use in the backprojection                                                                     
%                                   Input Data Type: boolean                                                                                          
% 
%         chanenvelopevars:         Select channels to use in the display                                                                             
%                                   Input Data Type: boolean                                                                                          
% 
%     BrainMovieOptions:            Additonal options for rendering the brainmovie                                                                    
%                                   Input Data Type: string                                                                                           
%     ------------------                                                                                                                              
% 
%         Visibility:               Figure visibility when rendering movie                                                                            
%                                   If 'on,' render frames on screen (slower). If 'off,' keep them hidden (faster).                                   
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         LatenciesToRender:        Subset of latencies to render (sec)                                                                               
%                                   Must be in TimeRange. Can be a vector The time point closest to the latency given are plotted. If                 
%                                   empty, render all latencies in TimeRange.                                                                         
%                                   Input Data Type: real number (double)                                                                             
% 
%         FramesToRender:           Vector of frame indices to compute                                                                                
%                                   E.g. [1:2] only computes the first two frames. If empty, render all frames                                        
%                                   Input Data Type: real number (double)                                                                             
% 
%         FigureHandle:             Handle to a figure to render brainmovie in                                                                        
%                                   Input Data Type: real number (double)                                                                             
% 
%         RotationPath3D:           Specify the rotation path for the BrainMovie                                                                      
%                                   Possible values: {'none','automatic','manual'}                                                                    
%                                   Default value  : 'none'                                                                                           
%                                   Input Data Type: string                                                                                           
%         ---------------                                                                                                                             
% 
%             AngleFactor:          Angle multiplicative factor                                                                                       
%                                   Input Data Type: real number (double)                                                                             
% 
%             PhaseFactor:          Phase multiplicative factor                                                                                       
%                                   Input Data Type: real number (double)                                                                             
% 
%         InitialView:              3D static starting view for the movie                                                                             
%                                   See 'help view'. Default is [50 36].                                                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%         ProjectGraphOnMRI:        Project the 3D graph onto the 2D MRI slices                                                                       
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'off'                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%         RenderCorticalSurface:    Superimpose semi-transparent smoothed cortex on brainmovie                                                        
%                                   UseOpenGL must be set to "on"                                                                                     
%                                   Input Data Type: boolean                                                                                          
%         ----------------------                                                                                                                      
% 
%             VolumeMeshFile:       Filename or path to mesh volume file                                                                              
%                                   Input Data Type: string                                                                                           
% 
%             Transparency:         Transparency of the cortical surface                                                                              
%                                   Real number in [0 1], where 0=opaque, 1=transparent.                                                              
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.7                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             VolumeMeshFile:       Filename or path to mesh volume file                                                                              
%                                   Input Data Type: string                                                                                           
% 
%             Transparency:         Transparency of the cortical surface                                                                              
%                                   Real number in [0 1], where 0=opaque, 1=transparent.                                                              
%                                   Input Range  : [0  1]                                                                                             
%                                   Default value: 0.7                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%         UseOpenGL:                OpenGL usage                                                                                                      
%                                   OpenGL may cause rendering problems with MacOSX                                                                   
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         EventFlashTimes:          Vector of time indices at which the background flashes                                                            
%                                   Specify the color of the flash with a cell array of [1,2] cell arrays. Ex. { { 200 'y' } { 1500 '5'               
%                                   }} will generate two flashes, yellow at 200 ms and red at 1500 ms                                                 
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%         Square:                   ?                                                                                                                 
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         DisplayLegendPanel:       Display legends in BrainMovie                                                                                     
%                                   Possible values: {'on','off'}                                                                                     
%                                   Default value  : 'on'                                                                                             
%                                   Input Data Type: string                                                                                           
% 
%         ShowLatency:              Display latency of current frame                                                                                  
%                                   This will render in lower left corner.                                                                            
%                                   Input Data Type: boolean                                                                                          
% 
%         DisplayRTProbability:     Display reaction time probabilty (if RT available)                                                                
%                                   This will render a small bar the height of which will vary based on the probability of response.                  
%                                   Input Data Type: boolean                                                                                          
% 
%         BackgroundColor:          Background color                                                                                                  
%                                   Can use any allowable Matlab color specification (see 'help ColorSpec').                                          
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%         GraphColorAndScaling:     Graph and Color Scaling                                                                                           
%                                   Options for coloring and scaling components of the directed graph                                                 
%                                   Input Data Type: string                                                                                           
%         ---------------------                                                                                                                       
% 
%             NodeSizeLimits:       [Min Max] limits for node size (pixels)                                                                           
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeColorLimits:      [Min Max] limits for node color (colormap index)                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeSizeLimits:       [Min Max] limits for edge size (pixels)                                                                           
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeColorLimits:      [Min Max] limits for edge color (colormap index)                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeSizeDataRange:    [Min Max] range of node size data                                                                                 
%                                   Input Data Type: real number (double)                                                                             
% 
%             NodeColorDataRange:   [Min Max] range of node color data                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeSizeDataRange:    [Min Max] range for edge size data                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             EdgeColorDataRange:   [Min Max] limits for edge color data                                                                              
%                                   Input Data Type: real number (double)                                                                             
% 
%             CenterDataRange:      Make 0 in the center of the colormap/datarange                                                                    
%                                   Input Data Type: boolean                                                                                          
% 
%             EdgeColorMap:         Expression defining the colormap for edges                                                                        
%                                   E.g., jet(64). See 'help colormap'.                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             NodeColormap:         Expression defining the colormap for nodes                                                                        
%                                   E.g., jet(64). See 'help colormap'.                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             DiskScalingFactor:    Numeric value that scales the size of disks                                                                       
%                                   Input Range  : [0  Inf]                                                                                           
%                                   Default value: 0.3                                                                                                
%                                   Input Data Type: real number (double)                                                                             
% 
%             MagnificationFactor:  Magnification factor for graphics                                                                                 
%                                   Input Range  : [0  Inf]                                                                                           
%                                   Default value: 1                                                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%         OutputFormat:             Options for saving the movie/figures                                                                              
%                                   Input Data Type: string                                                                                           
%         -------------                                                                                                                               
% 
%             ImageOutputDirectory: Output directory to save images                                                                                   
%                                   If 'prompt', then you will be prompted to select the folder from a dialog. If blank, don't save                   
%                                   images                                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%             ImageOutputFormat:    Format for saving images                                                                                          
%                                   Possible values: {'jpg','eps','ppm','tif','ai','bmp','emf','pbm','pcx','pdf','pgm','png','fig'}                   
%                                   Default value  : 'jpg'                                                                                            
%                                   Input Data Type: string                                                                                           
% 
%             MovieOutputFilename:  Movie filename                                                                                                    
%                                   E.g, 'movie.avi'. If 'prompt', then you will be prompted to select the file from a dialog. If                     
%                                   blank, don't save movie.                                                                                          
%                                   Input Data Type: string                                                                                           
% 
%             MovieOpts:            Cell array of movie options for avifile function                                                                  
%                                   See "help avifile".                                                                                               
%                                   Input Data Type: any evaluable Matlab expression.                                                                 
% 
%             ImageSize:            Image size (pixels)                                                                                               
%                                   Input should be [widthcond height]. If more than one condition is being plotted horizontally, then                
%                                   widthcond is the width of each condition subplot                                                                  
%                                   Input Data Type: real number (double)                                                                             
% 
%         mri:                      Dipplot MRI structure                                                                                             
%                                   Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a                
%                                   path to a Matlab file containing MRI structure. Default uses MNI brain.                                           
%                                   Input Data Type: string                                                                                           
% 
%         DipoleCoordinateFormat:   Coordinate format for dipplot                                                                                     
%                                   Possible values: {'spherical','mni'}                                                                              
%                                   Default value  : 'spherical'                                                                                      
%                                   Input Data Type: string                                                                                           
% 
%         DipplotOptions:           Additional dipplot options                                                                                        
%                                   Cell array of <'name',value> pairs of additional options for dipplot (see 'doc dipplot')                          
%                                   Input Data Type: any evaluable Matlab expression.  
%        
% Outputs:
%
%       cfg:        Configuration structure. Any brainmovie can be replicated 
%                   via the command vis_causalBrainMovie3D(ALLEEG,Conn,cfg);
%       handles:    Handles to figure and other objects
%       BMout:      Reserved for future use
%
%
% See Also: pop_vis_causalBrainMovie3D(), brainmovie3d_causal(), 
%
%
% References: 
% 
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual. Chapter 6.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift
%
% Author: Tim Mullen and Arnaud Delorme, 2010, SCCN/INC, UCSD. 
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



[cfg handles BMout] = deal([]);

subconddef = false;

% extract some stuff from inputs for arg defaults
Conn = arg_extract(varargin,'Conn',2);
if ~isempty(Conn)
    Conn = Conn(1);
    connnames   = hlp_getConnMethodNames(Conn);
    conndef     = connnames{1};
    freqrange   = [Conn.freqs(1) Conn.freqs(end)];
    freqdef     = ['[' num2str(freqrange(1)) ':' num2str(freqrange(end)) ']'];
    timerange   = [Conn.erWinCenterTimes(1) Conn.erWinCenterTimes(end)];
    timedef     = [timerange(1) timerange(end)];
    
    if length(Conn)>1
        subconddef = true;
    else
        subconddef = false;
    end
    
    clear Conn;
end

ALLEEG = arg_extract(varargin,'ALLEEG',1);
[MyComponentNames MyChannelNames] = deal([]);
if ~isempty(ALLEEG)
    ALLEEG = ALLEEG(1);
    if isfield(ALLEEG.CAT,'curComponentNames') && ~isempty(ALLEEG.CAT.curComponentNames)
        MyComponentNames = ALLEEG.CAT.curComponentNames;
    else
        MyComponentNames = ALLEEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    if isfield(ALLEEG,'chanlocs')
        MyChannelNames = {ALLEEG.chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:ALLEEG.nbchan)'))');
    end
    clear ALLEEG
end

% extract stats
% if isfield(varargin{1},'icaact') && length(varargin)==2
%     stats = [];
% elseif isfield(varargin{1},'icaact') && ~isempty(varargin(3:end))
%     stats = arg_extract(varargin(3:end),'stats');
% elseif ~isempty(varargin(3:end))
%     stats = arg_extract(varargin,'stats');
% elseif isempty(varargin(3:end))
%     stats = [];
% end
% 
% if isempty(stats)
%     usestatsdef = [];  % false
% else
%     usestatsdef = {};  % true
% end
% clear stats;

usestatsdef = [];  % false

% if length(varargin)>1 && iscell(varargin{1})
%     % arg_guidialog called the function...
%     % extract some arguments
%     varargin{1} = varargin{1}{1};
%     connnames   = hlp_getConnMethodNames(varargin{1});
%     conndef     = connnames{1};
%     freqrange   = [varargin{1}.freqs(1) varargin{1}.freqs(end)];
%     freqdef     = ['[' num2str(freqrange(1)) ':' num2str(freqrange(end)) ']'];
%     timerange   = [varargin{1}.erWinCenterTimes(1) varargin{1}.erWinCenterTimes(end)];
%     timedef     = [timerange(1) timerange(end)];
%     varargin    = varargin(2:end);
% else
%     % user called the function
%     [freqdef timedef conndef, ...
%      freqrange timerange connnames] = deal([]);
% end

% ensure we have row vectors
MyComponentNames = MyComponentNames(:)';
MyChannelNames = MyChannelNames(:)';


NodeColorMappingOpts = {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'};
NodeSizeMappingOpts = {'None','Outflow','Inflow','CausalFlow','Outdegree','Indegree','CausalDegree','AsymmetryRatio'};
EdgeColorMappingOpts = {'None','Connectivity','PeakFreq','Directionality'};
EdgeSizeMappingOpts = {'None','ConnMagnitude','Connectivity'};

g = arg_define([0 2],varargin, ...
    arg_norep({'ALLEEG'},mandatory),...
    arg_norep({'Conn'},mandatory),...
    arg_nogui({'stats','Stats'},'',[],'Name of variable in base containing statistics.','type','char','shape','row'), ...
    arg({'connmethod','ConnectivityMethod'},conndef,connnames,'Connectivity Measure to visualize','shape','row','cat','DataProcessing'), ...
    arg({'timeRange','MovieTimeRange'},timedef,{},'Time Range for Movie [Min Max] (sec). Specified w.r.t to event time (e.g., [-1 2]). Leave blank for complete epoch.','cat','DataProcessing'), ...
    arg({'freqsToCollapse','FrequenciesToCollapse'},freqdef,[],'Frequencies over which to collapse (Hz). E.g., [1:50] for 1-50Hz. Leave blank for all frequencies','type','expression','shape','row','cat','DataProcessing'), ...
    arg({'collapsefun','FreqCollapseMethod'},'mean',{'mean','max','peak','integrate'},'Method for collapsing frequency dimension.','cat','DataProcessing'),...
    arg({'resample','TimeResamplingFactor'},0,[0 20],'Time resampling factor. If 0, don''t resample. If < 1, downsample timecourse by this factor. If > 1, upsample by this factor. Uses resample() from Sigproc Toolbox','cat','DataProcessing'), ...
    arg({'subtractconds','SubtractConditions'},subconddef,[],'Subtract conditions. If true, then plot difference between conditions. If false, then render two brainmovies side-by-side.','cat','DataProcessing'), ...
    arg({'nodelabels','NodeLabels'},MyComponentNames,{},'List of labels for each node. e.g., {''Node1'',''Node2'',...}. Leave blank to omit labels.','shape','row','type','expression','cat','DisplayProperties'),...
    arg({'nodesToExclude','NodesToExclude'},false,MyComponentNames,'Exclude these sources from Brainmovie. Specify using the Name/ID of the source to exclude.','type','logical','cat','DisplayProperties'),...
    arg({'edgeColorMapping','EdgeColorMapping'},'Connectivity',EdgeColorMappingOpts,'Specify mapping for edge color. This determines how we index into the colormap. If ''None'', edge color is not modulated. If ''Connectivity'', use connectivity strength. If ''PeakFreq'', use index of peak frequency','cat','DisplayProperties'),...
    arg({'edgeSizeMapping','EdgeSizeMapping'},'ConnMagnitude',EdgeSizeMappingOpts,'Specify mapping for edge size. If ''None'', edges are not rendered. If ''Connectivity'', use connectivity strength. If ''ConnMagnitude'', use connectivity magnitude (absval). If ''PeakFreq'', use index of peak frequency. If ''Directionality'', map directionality to the lower and upper extremes of the colormap (e.g., i->j: blue, j->i: red)','cat','DisplayProperties'),...
    arg({'nodeColorMapping','NodeColorMapping'},'Outflow',NodeColorMappingOpts,'Specify mapping for node color. This determines how we index into the colormap. Options are as follows. None: node color is not modulated. Outflow: sum connectivity strengths over outgoing edges. Inflow: sum connectivity strengths over incoming edges. CausalFlow: Outflow-Inflow. Asymmetry Ratio: node colors are defined by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5 for balanced inflow/outflow','cat','DisplayProperties'), ...
    arg({'nodeSizeMapping','NodeSizeMapping'},'Outflow',NodeSizeMappingOpts,{'Specify mapping for node size. Options are as follows:' 'None: node size is not modulated.' 'Outflow: sum connectivity strengths over outgoing edges.' 'Inflow: sum connectivity strengths over incoming edges.' 'CausalFlow: Outflow-Inflow.' 'Asymmetry Ratio: node size is defined by the equation C = 0.5*(1 + outflow-inflow/(outflow+inflow)). This is 0 for exclusive inflow, 1 for exclusive outflow, and 0.5 for balanced inflow/outflow'},'cat','DisplayProperties'), ...
    arg({'baseline','Baseline'},[],[],'Time range of baseline [Min Max] (sec). Will subtract baseline from each point. Leave blank for no baseline.','cat','DataProcessing'),...
    arg_nogui({'normalize','NormalizeConn'},true,[],'Normalize edge and node values to [0 1]. Values mapped to edge/node width and color are devided by max to put in [0 1] range. Recommended!','cat','DataProcessing'), ...
    arg_subtoggle({'useStats','UseStatistics'},usestatsdef, ...
            { ... 
            arg({'threshside','Thresholding'},'single',{'single','both','lessthan'},'Type of thresholding for stats. If ''both'' then stats should have upper and lower thresholds'), ...
            arg({'alpha','AlphaSignificance'},0.05,[0 1],'Significance threshold. e.g., 0.05 for p<0.05') ...
            }, 'Use Statistics','cat','Thresholding'), ...
    arg({'prcthresh','PercentileThreshold'},[],[0 1],'Percentile threshold. Fraction of "strongest" connections to display. E.g: PercentileThreshold=0.05 will display only the top 5% of connections','type','denserealdouble','cat','Thresholding'), ...
    arg({'absthresh','AbsoluteThreshold'},[],[],'Exact threshold. If a single value, then render only connections with strength above this threshold. If [low hi] then render only connections with strength between [lo hi]. Overrides PercentileThreshold','type','denserealdouble','shape','row','cat','Thresholding'), ...
    arg_subswitch({'footerPanelSpec','FooterPanelDisplaySpec'},'off', ...
        {'off' { ...
            arg_norep({'dummy1'},[],[],'dummy') ...
            }, ...
         'ICA_ERPenvelope' {...
            arg({'icaenvelopevars'},true,MyComponentNames,'Select components to use in the display','type','logical'), ...
            arg({'backprojectedchans'},true,MyChannelNames,'List of channels to use in the backprojection'), ...
            } ...
          'Chan_ERPenvelope' {...
            arg({'chanenvelopevars'},true,MyChannelNames,'Select channels to use in the display','type','logical'), ...
            } ...
        },'Configure footer panel displayed at the bottom of the figure. If ''off'', don''t render footer. If ''ICA_ERP_Envelope'', then display the ERP envelope of backprojected components. If ''Chan_ERP_Envelope'' then display the ERP envelope of selected channels','cat','DisplayProperties'), ...
    arg_sub({'BMopts','BrainMovieOptions'},[], ...
        { ...
            arg({'visible','Visibility'},'on',{'on','off'},'Figure visibility when rendering movie. If ''on,'' render frames on screen (slower). If ''off,'' keep them hidden (faster).','cat','DisplayProperties'), ...
            arg_nogui({'latency','LatenciesToRender'},[],[],'Subset of latencies to render (sec). Must be in TimeRange. Can be a vector The time point closest to the latency given are plotted. If empty, render all latencies in TimeRange.','cat','DataProcessing'), ...
            arg_nogui({'frames','FramesToRender'},[],[],'Vector of frame indices to compute. E.g. [1:2] only computes the first two frames. If empty, render all frames','cat','DataProcessing'), ...
            arg_nogui({'figurehandle','FigureHandle'},[],[],'Handle to a figure to render brainmovie in'), ...
            arg_subswitch({'rotationpath3d','RotationPath3D'},'none', ...
                {'none' {arg_norep({'junk3'},[],[],'')}, ...
                 'automatic' {arg_norep({'junk2'},[],[],'')}, ...
                 'manual' { ...
                            arg({'AngleFactor'},1,[],'Angle multiplicative factor.'), ...
                            arg({'PhaseFactor'},0.75,[],'Phase multiplicative factor.') ...
                          } ...
                 },'Specify the rotation path for the BrainMovie.','cat','DisplayProperties'), ...
            arg({'view','InitialView'},[50 36],[],'3D static starting view for the movie.  See ''help view''. Default is [50 36].','cat','DisplayProperties'), ...
            arg({'project3d','ProjectGraphOnMRI'},'off',{'on','off'},'Project the 3D graph onto the 2D MRI slices','cat','DisplayProperties'), ...
            arg_subtoggle({'plotCortex','RenderCorticalSurface','plotcortex'},{}, ...
                        { 
                        arg({'cortexVolumeFile','VolumeMeshFile'},'standard_BEM_vol.mat',[],'Filename or path to mesh volume file.','type','char','shape','row'), ...
                        arg({'cortexTransparency','Transparency','cortextransparency'},0.7,[0 1],'Transparency of the cortical surface. Real number in [0 1], where 0=opaque, 1=transparent.') ...
                        },'Superimpose semi-transparent smoothed cortex on brainmovie. UseOpenGL must be set to "on"'), ...
            arg({'opengl','UseOpenGL'},'on',{'on','off'},'OpenGL usage. OpenGL may cause rendering problems with MacOSX','cat','DisplayProperties'), ...
            arg({'flashes','EventFlashTimes'},[],[],'Vector of time indices at which the background flashes.  Specify the color of the flash with a cell array of [1,2] cell arrays. Ex. { { 200 ''y'' } { 1500 ''5'' }} will generate two flashes, yellow at 200 ms and red at 1500 ms','type','expression','shape','row','cat','DisplayProperties'), ...
            arg_nogui({'square','Square'},'on',{'on','off'},'?','cat','Miscellaneous'), ...
            arg({'caption','DisplayLegendPanel'},'on',{'on','off'},'Display legends in BrainMovie','cat','DisplayProperties'), ...
            arg({'showLatency','ShowLatency'},true,[],'Display latency of current frame. This will render in lower left corner.','cat','DisplayProperties'), ...
            arg({'dispRT','DisplayRTProbability'},false,[],'Display reaction time probabilty (if RT available). This will render a small bar the height of which will vary based on the probability of response.','cat','DisplayProperties'), ...
            arg({'backcolor','BackgroundColor'},[0 0 0],[],'Background color.  Can use any allowable Matlab color specification (see ''help ColorSpec'').','shape','row','type','expression','cat','DisplayProperties'), ...
            arg_sub({'graphColorAndScaling','GraphColorAndScaling'},{}, ...
                { ...
                arg({'nodeSizeLimits','NodeSizeLimits'},[0.1 1],[],'[Min Max] limits for node size (pixels).','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColorLimits','NodeColorLimits'},[0 1],[],'[Min Max] limits for node color (colormap index).','shape','row','cat','DisplayProperties'), ...
                arg({'edgeSizeLimits','EdgeSizeLimits'},[0.1 1],[],'[Min Max] limits for edge size (pixels)','shape','row','cat','DisplayProperties'), ...
                arg({'edgeColorLimits','EdgeColorLimits'},[0 1],[],'[Min Max] limits for edge color (colormap index).','shape','row','cat','DisplayProperties'), ...
                arg({'nodeSizeDataRange','NodeSizeDataRange'},[],[],'[Min Max] range of node size data.','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColorDataRange','NodeColorDataRange'},[],[],'[Min Max] range of node color data.','shape','row','cat','DisplayProperties'), ...
                arg({'edgeSizeDataRange','EdgeSizeDataRange'},[],[],'[Min Max] range for edge size data.','shape','row','cat','DisplayProperties'), ...
                arg({'edgeColorDataRange','EdgeColorDataRange'},[],[],'[Min Max] limits for edge color data.','shape','row','cat','DisplayProperties'), ...
                arg({'centerDataRange','CenterDataRange'},false,[],'Make 0 in the center of the colormap/datarange','cat','DisplayProperties'), ...
                arg({'edgeColormap','EdgeColorMap'},'jet(64)',[],'Expression defining the colormap for edges. E.g., jet(64). See ''help colormap''.','type','expression','shape','row','cat','DisplayProperties'), ...
                arg({'nodeColormap','NodeColormap'},'jet(64)',[],'Expression defining the colormap for nodes. E.g., jet(64). See ''help colormap''.','type','expression','shape','row','cat','DisplayProperties'), ...
                arg({'diskscale','DiskScalingFactor'},0.3,[0 Inf],'Numeric value that scales the size of disks.','cat','DisplayProperties'), ...
                arg({'magnify','MagnificationFactor'},1,[0 Inf],'Magnification factor for graphics','cat','DisplayProperties') ...
                },'Graph and Color Scaling. Options for coloring and scaling components of the directed graph'), ...
            arg_sub({'outputFormat','OutputFormat'},{}, ...   
                { ...
                arg({'framefolder','ImageOutputDirectory'},[],[],'Output directory to save images. If ''prompt'', then you will be prompted to select the folder from a dialog. If blank, don''t save images','shape','row','type','char'), ...
                arg({'framesout','ImageOutputFormat'},'jpg',{'jpg','eps','ppm','tif','ai','bmp','emf','pbm','pcx','pdf','pgm','png','fig'},'Format for saving images'), ...
                arg({'moviename','MovieOutputFilename'},[],[],'Movie filename. E.g, ''movie.avi''. If ''prompt'', then you will be prompted to select the file from a dialog. If blank, don''t save movie.','shape','row','type','char'),...
                arg({'movieopts','MovieOpts'},{'videoname',''},[],'Cell array of movie options for avifile function. See "help avifile".','type','expression','shape','row'), ...
                arg({'size','ImageSize'},[600 600],[],'Image size (pixels). Input should be [widthcond height]. If more than one condition is being plotted horizontally, then widthcond is the width of each condition subplot','cat','DisplayProperties') ...
                },'Options for saving the movie/figures','cat','MovieOutput'), ...
            arg({'mri'},'',[],'Dipplot MRI structure. Can be the name of matlab variable (in the base workspace) containing MRI structure. May also be a path to a Matlab file containing MRI structure. Default uses MNI brain.','type','char','shape','row'), ...
            arg({'coordformat','DipoleCoordinateFormat'},'spherical',{'spherical','mni'},'Coordinate format for dipplot','type','char','shape','row'), ...
            arg({'dipplotopt','DipplotOptions'},'{}','','Additional dipplot options. Cell array of <''name'',value> pairs of additional options for dipplot (see ''doc dipplot'')','type','expression','shape','row'), ...
            arg_nogui('renderBrainMovie',true,[],'Special option. Determines whether to actually render the brainmovie (or just return values)') ...
            }, ...
    'Additonal options for rendering the brainmovie','cat','DisplayProperties') ...
    );

       
    
    % insert an arg_sub here that calls @brainmovie3D to populate
    % arguments list

    %  TODO: edit brainmovie3d to add arg specification

    
% Commit data variables to workspace    
% [data g] = hlp_splitstruct(g,{'ALLEEG','Conn','Stats'});        
% arg_toworkspace(data);
% clear data;
    
% copy the stats structure from base to current workspace
if ~isempty(g.stats)
    g.stats = evalin('base',g.stats);
end


% VALIDATE INPUTS
% ---------------------------------------------

if ~isfield(g,'renderBrainMovie')
    g.renderBrainMovie = true;
end

if ~ismember(g.collapsefun,{'max','peak'}) && strcmpi(g.edgeColorMapping,'peakfreq')
    error('To use PeakFreq EdgeColorMapping, you must select ''max'' or ''peak'' as the FreqCollapseMethod');
end

if length(setdiff(MyComponentNames,g.nodesToExclude)) < 2
    error('You must include at least two nodes in the BrainMovie');
end
    
% if no stats present, disable stats
if isempty(g.stats)
    g.useStats.arg_selection = false;  end
    
% check that frequencies are in valid range
if any(g.freqsToCollapse < freqrange(1)) || ...
   any(g.freqsToCollapse > freqrange(end))
    error('FreqsToCollapse contains elements outside valid range [%0.1f %0.1f]', ...
        freqrange(1),freqrange(end));
end
    
% check that times are in valid range
if g.timeRange(1)+1e-5 < timerange(1) || ...
   g.timeRange(2)-1e-5 > timerange(end)
    error('TimeRange contains elements outside valid range [%0.1f %0.1f]', ...
        timerange(1),timerange(end));
end

% check that latencies are in TimeRange
if any(g.BMopts.latency < g.timeRange(1)) || ...
   any(g.BMopts.latency > g.timeRange(end))
    error('Latency (%0.1f s) is outside the specified TimeRange [%0.1f %0.1f]. Please provide a valid latency', ...
            g.BMopts.latency,g.timeRange(1),g.timeRange(1));
end

if ~isempty(g.BMopts.mri) 
    if evalin('base',['exist(''' g.BMopts.mri ''',''var'')'])==1
        % MRI variable is in workspace, so copy it
        g.BMopts.mri = evalin('base',g.BMopts.mri);
        
        if ~isfield(g.BMopts.mri,'anatomy')
            error('MRI structure is invalid format (see ''dipplot'' for more info)');
        end
    elseif isdir(fileparts(g.BMopts.mri))
        % User specified path to MRI file, so load it up
        tmp = load(g.BMopts.mri);
        fn = fieldnames(tmp);
        g.BMopts.mri = tmp.(fn{1});
    else
        % User specified an invalid path to MRI file
        error('Invalid path to MRI matlab file');
    end
end

if g.BMopts.plotCortex.arg_selection
    tmp = load(g.BMopts.plotCortex.cortexVolumeFile);
    fn = fieldnames(tmp);
    g.BMopts.csf.vol = tmp.(fn{1});
end

% ---------------------------------------------

    

% HANDLE DEFAULTS
% ---------------------------------------------

% do some cleanup
g.nodeSizeMapping   = lower(g.nodeSizeMapping);
g.nodeColorMapping  = lower(g.nodeColorMapping);
g.edgeSizeMapping   = lower(g.edgeSizeMapping);
g.edgeColorMapping  = lower(g.edgeColorMapping);


if isempty(MyComponentNames)
    if isfield(g.ALLEEG.CAT,'curComponentNames') && ~isempty(g.ALLEEG.CAT.curComponentNames)
        MyComponentNames = g.ALLEEG.CAT.curComponentNames;
    else
        MyComponentNames = g.ALLEEG.CAT.curComps;
        MyComponentNames = strtrim(cellstr(num2str(MyComponentNames'))');
    end
    
    if isfield(g.ALLEEG,'chanlocs')
        MyChannelNames = {g.ALLEEG.chanlocs.labels};
    else
        MyChannelNames = strtrim(cellstr(num2str((1:g.ALLEEG.nbchan)'))');
    end
end
    
if isempty(g.connmethod)
    % if no connectivity methods specified, select first one
    g.connmethod = hlp_getConnMethodNames(g.Conn(1));
    g.connmethod = g.connmethod{1};                         end
if isempty(g.freqsToCollapse)
    g.freqsToCollapse = g.Conn(1).freqs;                    end
if isempty(g.timeRange)
    g.timeRange = g.Conn(1).erWinCenterTimes([1 end]);      end


% make a copy for convenience
erWinCenterTimes = g.Conn(1).erWinCenterTimes;

% identify the indices of the components to keep
g.BMopts.selected = find(~ismember(MyComponentNames,g.nodesToExclude));

% obtain factors for resampling
if g.resample
    [p q] = rat(g.resample);
    g.resample = [p q];
end

if g.resample
    erWinCenterTimes = resample(erWinCenterTimes,g.resample(1),g.resample(2));
end


% number of sources
N=g.ALLEEG(1).CAT.nbchan;

% indices of nodes we will exclude
nodeIndicesToExclude = find(~ismember(1:N,g.BMopts.selected));

% get the indices of frequencies to collapse
freqIndicesToCollapse = getindex(g.Conn(1).freqs,g.freqsToCollapse);

% get indices of desired time windows ...
timeIndices =  getindex(erWinCenterTimes,g.timeRange(1)):getindex(erWinCenterTimes,g.timeRange(2));

% ... and construct the time vector for the Brainmovie
BrainMovieTimeRangeInMs = 1000*erWinCenterTimes(timeIndices);

% ---------------------------------------------


% TRANSFORM DATA INTO BRAINMOVIE FORMAT
% ---------------------------------------------
for cnd = 1:length(g.Conn)
    % select the connectivity measure
    Conn = g.Conn(cnd).(g.connmethod);
    
    % select the appropriate stats
    if ~isempty(g.stats)
        Stats = g.stats(cnd).(g.connmethod);
    end

    % remove the baseline
    if ~isempty(g.baseline)
        Conn = hlp_rmbaseline(Conn,g.baseline,g.Conn(cnd).winCenterTimes);
    end

    % Apply statistics and thresholding
%     if g.absthresh > 0
%         % absolute threshold
%         Conn(Conn < g.absthresh) = 0;
%     elseif g.prcthresh<1
%         % simple percentile threshold
%     %     prcthresh = repmat(prctile(c.Conn,4),(1-g.prcthresh)*100,3),[1 1 size(Conn,3) size(Conn,4)]);
    if g.useStats.arg_selection
        if isempty(Stats)
            error('IFT:StatsThresholding:NoStats','You must supply a ''stats'' structure. See help');
            return;
        end
        % apply statistical threshold
%         Conn(Stats.p > g.useStats.alpha) = 0;

        switch g.useStats.threshside
            case 'both'
                Conn((Conn > squeeze(Stats.thresh(1,:,:,:,:,:))) & (Conn < squeeze(Stats.thresh(2,:,:,:,:,:)))) = 0;
            case 'single'
                Conn(abs(Conn) < abs(Stats.thresh)) = 0;
            case 'lessthan'
                Conn(Conn < squeeze(Stats.thresh(2,:,:,:,:,:))) = 0;
        end
    end

    
    % Format edge data (size/color)
    % -----------------------------
    [EdgeSize EdgeColor] = deal(zeros(N,N,length(timeIndices),class(Conn(1))));
    for ch1=1:N
        for ch2=1:N

            if ch1==ch2 ...
               || ~isempty(intersect(nodeIndicesToExclude,[ch1 ch2])) ...
               || all(ismember({g.edgeSizeMapping, g.edgeColorMapping},'none'))
           
                % don't compute data for this channel pair
%                 EdgeSize(ch1,ch2,:)=zeros(1,length(timeIndices),class(Conn(1)));
%                 EdgeColor(ch1,ch2,:) = EdgeSize{ch1,ch2};
                continue;
            end
            
            % extract data
            causality = squeeze(Conn(ch1,ch2,:,:));
            
%             % make row vector if necessary
%             if any(size(causality)==1)
%                 causality = causality(:)';
%             end
            
            % collapse matrix across frequencies
            if length(g.freqsToCollapse)>1
                [causality peakidx] = hlp_collapseFrequencies( causality, ...
                  g.collapsefun,freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
            elseif size(causality,2)==1
                % only one frequency
                causality = causality(:,freqIndicesToCollapse);
                peakidx = freqIndicesToCollapse;
            else
                % multiple freqs, but only one frequency selected
                causality = causality(freqIndicesToCollapse,:);
                peakidx = freqIndicesToCollapse;
            end

             
%                 interp1(1:1:size(causality,2), ...
%             causality', ...
%             1:1/INTERPFAC:size(causality,2))'; 
                
            if g.resample
                causality   = resample(causality.',g.resample(1),g.resample(2)).';
                peakidx     = resample(peakidx.',g.resample(1),g.resample(2)).';
            end

            
            switch g.edgeSizeMapping
                case 'connectivity'
                    EdgeSize(ch1,ch2,:) = causality;
                case 'peakfreq'
                    EdgeSize(ch1,ch2,:) = peakidx;
                case 'connmagnitude'
                    EdgeSize(ch1,ch2,:) = abs(causality);    
%                 case 'none'
%                     EdgeSize(ch1,ch2,:) =
%                     zeros(1,length(timeIndices),class(Conn(1)));
            end
            
            if strcmpi(g.edgeSizeMapping,g.edgeColorMapping)
                EdgeColor(ch1,ch2,:) = EdgeSize(ch1,ch2,:);
            else
                switch g.edgeColorMapping
                    case 'connectivity'
                        EdgeColor(ch1,ch2,:) = causality;
                    case 'peakfreq'
                        EdgeColor(ch1,ch2,:) = peakidx;
%                     case 'none'
%                         EdgeColor(ch1,ch2,:) = zeros(1,length(timeIndices),class(Conn(1)));
                    case 'directionality'
                        EdgeColor(ch1,ch2,:) = causality;
                    case 'connmagnitude'
                        EdgeColor(ch1,ch2,:) = abs(causality);
                end
            end
            
            
            
        end % for ch2
    end % for ch1

    
    if ~isempty(g.absthresh)
        if isscalar(g.absthresh)
            edgeSizeThresh = [-Inf g.absthresh];
        else
            edgeSizeThresh = g.absthresh;
        end
        
        EdgeSize(EdgeSize   > edgeSizeThresh(1)     ...
               & EdgeSize   < edgeSizeThresh(2))    = 0;   
    end
    
    if ~isempty(g.prcthresh) && g.prcthresh < 1
        % apply percentile thresholding of the edges
%         if any(EdgeColor(:) < 0)
%             % two-sided thresholds
%             prcthresh = [g.prcthresh/2 1-g.prcthresh/2]*100;
%         else
%             prcthresh = [0 1-g.prcthresh]*100;
%         end
        
        prcthresh = [0 1-g.prcthresh]*100;

        % get thresholds
        edgeSizeThresh  = prctile(abs(EdgeSize(:)) ,prcthresh);
%         edgeColorThresh = prctile(abs(EdgeColor(:)),prcthresh);
        
        EdgeSize(abs(EdgeSize)   > edgeSizeThresh(1)     ...
               & abs(EdgeSize)   < edgeSizeThresh(2))    = 0;   
%         EdgeColor(abs(EdgeColor) > edgeColorThresh(1) 	...
%                 & abs(EdgeColor) < edgeColorThresh(2))   = 0;     
    end
    

    
    
    if ~all(ismember({g.nodeSizeMapping, g.nodeColorMapping},'none'))
        % we will need the connectivity data in matrix format so we can
        % obtained graph statistics for mapping Node Size/Color
        
        if strcmpi(g.edgeSizeMapping,'connectivity')
            causality = squeeze(EdgeSize(:,:,:));
        elseif strcmpi(g.edgeColorMapping,'connectivity')
            causality = squeeze(EdgeColor(:,:,:));
        else % Conn matrices have not been collapsed ...
             %  ... we have to collapse matrix across frequencies
             causality = zeros(N,N,length(timeIndices));
             for ch1=1:N
                for ch2=1:N
                    [causality(ch1,ch2,:)] = hlp_collapseFrequencies( ...
                            squeeze(Conn(ch1,ch2,:,:)), g.collapsefun, ...
                            freqIndicesToCollapse,timeIndices,g.freqsToCollapse(2)-g.freqsToCollapse(1));
                end
             end
        end
    end


    % Format node data (size/color)
    [NodeSize NodeColor] = deal(zeros(N,length(timeIndices),class(Conn(1))));
    for ch1=1:N
        if ismember(ch1,nodeIndicesToExclude) || ...
           all(ismember({g.nodeSizeMapping, g.nodeColorMapping},'none'))
%             NodeSize(ch1,:) = zeros(size(EdgeSize{1,2}),class(Conn(1)));
%             NodeColor(ch1,:) = NodeSize{ch1,1};
            continue;

        end
        
        % compute desired graph measure for this node and map to size
        % we add eps=10^-16 just to ensure that values are never exactly
        % zero (or they will get set to NaN later)
        NodeSize(ch1,:) = hlp_computeGraphMeasure(causality,ch1,g.BMopts.selected,g.nodeSizeMapping)+eps;
        
        if strcmp(g.nodeSizeMapping,g.nodeColorMapping)
            NodeColor(ch1,:) = NodeSize(ch1,:);
        else
            % compute desired graph measure for this node and map to color
            NodeColor(ch1,:) = hlp_computeGraphMeasure(causality,ch1,g.BMopts.selected,g.nodeColorMapping)+eps;
        end
                
        if g.resample
            NodeSize(ch1,:)  = resample(NodeSize(ch1,:).',g.resample(1),g.resample(2)).';
            NodeColor(ch1,:) = resample(NodeColor(ch1,:).',g.resample(1),g.resample(2)).';
        end
    end % loop over channels
    
    
    % check for negative values
%     HaveNegativeValues = (any(EdgeSize(:)<0) ||  any(NodeColor(:)<0));
    
    % set to NaN any nonsignificant values
    NodeSize(NodeSize==0)   = nan;
    NodeColor(NodeColor==0) = nan;
    EdgeSize(EdgeSize==0)   = nan;
    EdgeColor(EdgeColor==0) = nan;
    
    
    if nargout>2
        BMout.NodeSize  = NodeSize;
        BMout.NodeColor = NodeColor;
        BMout.EdgeSize  = EdgeSize;
        BMout.EdgeColor = EdgeColor;
    end
        

    % convert to Brainmovie format
    % ----------------------------
    for ch1=1:N
        for ch2=1:N
            if ch1==ch2, continue; end
                EdgeSizeCell{ch1,ch2,cnd}   = squeeze(EdgeSize(ch1,ch2,:))';
                EdgeColorCell{ch1,ch2,cnd}  = squeeze(EdgeColor(ch1,ch2,:))';
        end
    end
    
    for ch1=1:N
        NodeSizeCell{ch1,cnd}     = squeeze(NodeSize(ch1,:));
        NodeColorCell{ch1,cnd}    = squeeze(NodeColor(ch1,:));
    end
    
    
    % Create footer panel data
    % ---------------------------
    if strcmp(g.footerPanelSpec.arg_selection,'off')
        envdata     = [];
        envtimes    = [];
    else
        switch g.footerPanelSpec.arg_selection
            
            case 'Chan_ERPenvelope'
                % compute ERP of selected channels
                erpvars     = ismember(MyChannelNames, ...                      % select chans
                                       g.footerPanelSpec.chanenvelopevars);
                erpdata     = mean(g.ALLEEG(cnd).data(erpvars,:,:),3)';         % compute ERP
                if size(erpdata,1)>1    % if more than one channel selected
                    erpdata     = env(erpdata);                                 % compute envelope
                else
                    erpdata = erpdata;
                end
            case 'ICA_ERPenvelope'
                % compute ERP of selected backprojected components
                erpvars     = g.ALLEEG(cnd).CAT.curComps( ...
                                ismember(MyComponentNames(g.BMopts.selected),...  % select comps
                                       g.footerPanelSpec.icaenvelopevars));
                erpchans    = ismember(MyChannelNames,g.footerPanelSpec.backprojectedchans);
                erpdata     = mean(g.ALLEEG(cnd).icaact(erpvars,:,:),3)';       % obtain ERP
                erpdata     = g.ALLEEG(cnd).icawinv(erpchans,erpvars)*erpdata';        % backproject ERPs to scalp
                if size(erpdata,1)>1    % if more than one channel selected
                    erpdata     = env(erpdata);                                 % compute envelope
                else
                    erpdata = erpdata;
                end
        end

        erptimes    = g.ALLEEG(cnd).times;
        erptlims    = getindex(erptimes,g.timeRange*1000);
        erpdata     = erpdata(:,erptlims(1):erptlims(2));
        erptimes    = erptimes(erptlims(1):erptlims(2));

        if g.resample
            % resample the envelope timecourse to make consistent
            envdata(:,:,cnd)    = resample(envdata' ,g.resample(1),g.resample(2))';
            envtimes(:,cnd)     = resample(erptimes,g.resample(1),g.resample(2));
        else
            envdata(:,:,cnd)    = erpdata;
            envtimes(:,cnd)     = erptimes;
        end
    end
    
end  % loop over conditions


% subtract multiple conditions
if g.subtractconds && size(EdgeSizeCell,3)==2  
    for ch1=1:N
        for ch2=1:N
            if ch1==ch2, continue; end
                tmpedgesize{ch1,ch2,1}   = EdgeSizeCell{ch1,ch2,1}  - EdgeSizeCell{ch1,ch2,2};
                tmpedgecolor{ch1,ch2,1}  = EdgeColorCell{ch1,ch2,1} - EdgeColorCell{ch1,ch2,2};
        end
    end
    EdgeSizeCell    = tmpedgesize;
    EdgeColorCell   = tmpedgecolor;

    for ch1=1:N
        tmpnodesize{ch1,1}     = NodeSizeCell{ch1,1}  - NodeSizeCell{ch1,2};
        tmpnodecolor{ch1,1}    = NodeColorCell{ch1,1} - NodeColorCell{ch1,2};
    end
    
    NodeSizeCell    = tmpnodesize;
    NodeColorCell   = tmpnodecolor;
    
    g.BMopts.condtitle = sprintf('%s - %s', g.ALLEEG(1).condition,g.ALLEEG(2).condition); 
else
    g.BMopts.condtitle    = {g.ALLEEG.condition};
end



% SETUP THE BRAINMOVIE OPTIONS
% ---------------------------------------------
BMopts              = g.BMopts;
BMopts.causality    = true;
BMopts.nodelabels   = g.nodelabels;
BMopts.title = BMopts.condtitle;
% extract the coordinates of dipoles
% for cnd=1:length(g.Conn)
%     coords{cnd} = {g.ALLEEG(cnd).dipfit.model.posxyz};
%     coords{cnd} = coords{cnd}(g.ALLEEG(cnd).CAT.curComps);
% end
coords = {g.ALLEEG(cnd).dipfit.model.posxyz};
coords = coords(g.ALLEEG(cnd).CAT.curComps);
BMopts.coordinates = coords;

% determine whether to render and/or modulate properties of edges/nodes 
BMopts.crossfphasecolor = fastif(strcmpi('none',g.edgeColorMapping),'off','on');    % do/don't modulate edgecolor
BMopts.crossf           = fastif(strcmpi('none',g.edgeSizeMapping),'off','on');     % do/don't render edges
BMopts.itc              = fastif(strcmpi('none',g.nodeColorMapping),'off','on');    % do/don't modulate nodecolor
BMopts.power            = fastif(strcmpi('none',g.nodeSizeMapping),'off','on');     % do/don't modulate nodesize

SLASH = fastif(isunix,'/','\');

% Select the folder to save frames
if strcmpi(BMopts.outputFormat.framefolder,'prompt')
    framefolder = uigetdir(pwd,'Select directory to save all movie frames');
    if framefolder~=0
        BMopts.outputFormat.framefolder = framefolder;
    else
        BMopts.outputFormat.framefolder = '';
    end
elseif ~isempty(BMopts.outputFormat.framefolder) && ~isdir(BMopts.outputFormat.framefolder)
    resp = questdlg2(sprintf('You asked to save figures here:\n ''%s''\nThis folder does not exist. Should I create it?',BMopts.outputFormat.framefolder),'Yes','No','No');
    if strcmpi(resp,'Yes')
        mkdir(BMopts.outputFormat.framefolder);
    else
        BMopts.outputFormat.framefolder = '';
    end
end
    

% select the movie filename
if strcmpi(BMopts.outputFormat.moviename,'prompt')
    moviename = uiputfile('*.mov, *.avi','Create output movie file');
    if moviename~=0
        BMopts.outputFormat.moviename = moviename;
    else
        BMopts.outputFormat.moviename = '';
    end
end


BMopts.times        = envtimes;
BMopts.envelope     = envdata;

% get latencies for movie
if isempty(g.BMopts.latency) && isempty(g.BMopts.frames)
    % default: use all latencies in time range
    BMopts.latency = BrainMovieTimeRangeInMs;
else
    % convert to ms
    BMopts.latency = 1000*g.BMopts.latency;
end


% setup the rotation path for the movie
switch lower(BMopts.rotationpath3d.arg_selection)
    case 'none'
        BMopts.path3d = 'off';
    case 'automatic'
        BMopts.path3d = 'on';
    case 'manual'
        % polar rotation factors specified 
        BMopts.path3d = ...
        [BMopt.rotationpath3d.AngleFactor, BMopt.rotationpath3d.PhaseFactor];
end
    
if ~BMopts.plotCortex.arg_selection
    % don't superimpose cortex
    BMopts.cortexTransparency = 1;
end

if strcmpi(g.edgeColorMapping,'Directionality')
    BMopts.EdgeColorMappedToDirectionality = true;
end


% check if we have negative values (and need +- colormaps)
for i=1:length(NodeColorCell)
    if ~isempty(NodeColorCell{i})
        if any(NodeColorCell{i} < 0)
            BMopts.nodeColorPolarity = 'posneg';
        end
    end
end
for i=1:length(EdgeColorCell(:))
    if ~isempty(EdgeColorCell{i})
        if any(EdgeColorCell{i} < 0)
            BMopts.edgeColorPolarity = 'posneg';
        end
    end
end
% 
% if any(any(any(cell2mat(NodeColorCell)<0)))
%     BMopts.nodeColorPolarity = 'posneg'; end
% if any(any(any(cell2mat(EdgeColorCell)<0)))
%     BMopts.edgeColorPolarity = 'posneg'; end

BMopts.nodeColorMapping = g.nodeColorMapping;
BMopts.edgeColorMapping = g.edgeColorMapping;
BMopts.nodeSizeMapping  = g.nodeSizeMapping;
BMopts.edgeSizeMapping  = g.edgeSizeMapping;
BMopts.ConnMethod       = g.connmethod;

% convert options structure to ('name',value) arglist
bmargs = hlp_struct2varargin(hlp_flattenStruct(BMopts,'exclude',{'mri','csf'}),'suppress', ...
                             {'arg_selection','arg_direct',...
                             'rotationpath3d','AngleFactor',...
                             'PhaseFactor'});

if g.renderBrainMovie
    % Call BrainMovie3D
    brainmovie3d_causal( ...
        NodeSizeCell,NodeColorCell,EdgeSizeCell,EdgeColorCell, ...
        BrainMovieTimeRangeInMs,1,g.BMopts.selected,bmargs{:});
end

% EOF






        
        

