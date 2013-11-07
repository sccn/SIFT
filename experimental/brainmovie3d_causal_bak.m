% brainmovie3d() - generate a sequence of images showing event-related coherence,
%               event-related spectral perturbations, and inter-trial coherence
%               of localized EEG waveforms. Uses outputs of timef() and cross().
% Usage:
%   >> brainmovie3d(ersps,itcs,crossfs_amp,crossfs_phase,TIMES,freqs,SELECTED,...
%                         'keyword1',value1,...); % creates files image0001.eps, etc.
%
% Inputs:
% ersps         - Cell array (components,conditions) of ERSP arrays
% (freqs,TIMES)
%                 ERSP = event-related spectral perturbation; returned by
%                 timef()
% itcs          - Cell array (components,conditions) of ITC arrays
% (freqs,TIMES)
%                 ITC = inter-trial coherence; returned by timef()
% crossfs_amp   - Cell array (components,components,conditions) of crossf()
%
%                 amplitude output arrays of size (freqs,TIMES).
% crossfs_phase - Cell array (components,components,conditions) of crossf() phase
%                 output arrays of size (freqs,TIMES). (Only the upper diagonal part
%                 of the matrix is taken into account).
% TIMES         - Array of TIMES returned by timef() or crossf()
% freqs         - Indices into the array of freqs returned by timef() or
% crossf()
%                 (v.henv.g., [1:2] means plot the mean of the first two frequencies).
%                 These indexes determine for which freqs plotting will be
%                 performed.
% SELECTED      - Component indices to plot (default all)
%
% Optional 'keyword' parameters:
% 'latency'   - plot only a subset of latencies. The time point closest to the
%               latency given are plotted. Default = empty, all latencies.
% 'frames'    - vector of frame indices to compute. [1:2] only computes the
%               first two frames.
% 'envelope'  - (2,points,conditions) envelopes of the average data (ERP)
%               in each condition
%               (envelope =  min and max traces of each ERP across all
%               channels and TIMES)
% 'rt'        - cell array of vector containing reaction TIMES of the subject in
%               each conditions. This will plot a small bar which height will vary
%               based on the probability of response (default {} -> ignored)
% 'flashes'   - vector of time indices at which the background flashes.  Specify the color
%               of the flash with a cell array of [1,2] cell arrays.
%               Ex. { { 200 'y' } { 1500 '5' }} will generate two flashes,
%               yellow at 200 ms and red at 1500 ms
%
% Movie ITC, Power and Crossf options:
% 'power'     - ['on'|'off'] vary the size of the component disks according to spectral power
%                                                           {default: on}
% 'itc'       - ['on'|'off'] vary component disk colors according to inter-trial coherence
%							    {default: on}
% 'crossf'    - ['on'|'off'] plot | do not plot coherence   {default: on}
% 'crossfcoh' - ['on'|'off'] vary the width of the connecting arc
%                               according to cross-coherence magnitude {def: on}
% 'crossfphasecolor'   -['on'|'off'] vary the arc color according to coherence {default: on}
% 'crossfphasespeed'   - ['on'|'off'] vary the arc speed according to
%                                      cross-coherence phase {def: off}
% 'crossfphaseunit'    - ['degree'|'radian']. Coherence phase angle unit
% {Default is degree}.
% 'edgeColormap'       - colormap array for arcs {default: hsv(64) with green as 0}
% 'nodeColormap'       - colormap array for disks {default: hot(64)}
% 'nodeSizeDataRange'  - [min max] data range for disk size {default is min/max of data}
% 'nodeSizeLimits'     - [min max] limits for node size {default 0.1 to 1}
% 'nodeColorDataRange' - [min max] data range for disk color {default is min/max of data}
% 'nodeColorLimits'    - [min max] limits for node color {default full colormap 0 to 1}
% 'edgeSizeDataRange'  - [min max] data range for edge size {default is min/max of data}
% 'edgeSizeLimits'     - [min max] limits for edge size {default 0.1 to 1}
% 'edgeColorDataRange' - [min max] data range for edge color {default is
% min/max of data}
% 'edgeColorLimits'    - [min max] limits for edge color {default full colormap 0 to 1}
% 'polarity'  - ['pos'|'posneg'] polarity for ITC and crossf. 'pos' = only positive values
%               'posneg' = positive and negative values.
%
% Movie coordinates and axis options:
% 'magnify'   - integer magnification factor for graphics. Default is 1.
% 'diskscale'   - numeric value that scales the size of disks {default: [1.0]}
% 'xlimaxes'    - x-axis limits axis for the component locations {default: [-1 1]}
% 'ylimaxes'    - y-axis limits axis for the component locations {default: [-1 to 1]}
% 'coordinates' - 3-column array of [x y z] coordinates of the SELECTED components
%                 {default: spaced evenly around the head circle boundary}
% 'square'    - ['on'|'off'] re-square all coordinates (so X and Y width is the same)
%               default is 'on';
% 'project3d' - ['on'|'off'] project disks on each 3-D axis. Default is 'off'.
% 'circfactor'  - (ncomps,ncomps) array of arc curvatures (0=straight; 1=half-round,
%                 positive or negative values give the sense of rotation)
%                 {def: 0s}
% 'envylabel'   - ordinate label for envelope. {Default 'Potential \muV'}
% 'envvert'     - cell array of time indices at which to draw vertical lines.
%                 Can also be a cell array of cell to specify line aspect. For instance
%                 { { 0 'color' 'b' 'linewidth' 2 } {1000 'color' 'r' }} would draw two
%                 lines, one blue thick line at latency 0 and one thin red line at latency 1000.
% 'rthistloc' - location and size of rt histograms in individual axes.
%               [abscissa ordinate width maxheight].
% 'title'       - (string) main movie title
% 'condtitle'   - (string array) condition titles (one condition title per row)
% 'condtitleformat' - list of title properties. Ex: { 'fontize', 12,
% 'fontweight', 'bold' }
% 'plotorder'   - [integer vector] component plot order from 1 to the number of SELECTED
%                 components.
% 'backcolor' - [float array] background color. Default is [1 1 1] (white).
%
% Picture and movie output options:
% 'moviename'  - ['string'] Movie file name. Default is "output.avi".
% 'movieopts'  - [cell] Movie options for avifile function. See "help
% avifile".
% 'framesout'  - ['eps'|'ppm'|'fig'|'tiff'|'none'] Default format for saving frames on disk.
%                Default is 'tiff'.
% 'framefolder' - [string] frames output folder. Default uses current directory.
%               the directory is created if it does not exist.
% 'visible'    - ['on'|'off'] show the images on the screen or keep them hidden {default 'on'}
% 'size'      - [widthcond height] output image size {default [400,400]}
%               widthcond is the width of a single condition plot (in pixels)
% 'view'      - 3D static starting view.  See help view. Default is [1 0 0].
% 'path3d'    - ['on'|'off'|[thetafact phifact]] 'on' activate automatic rotation in 3-D. Use
%               [exttheta extphi] to specify theta and phi multiplicative factor (default is
%               [1 0.75]. Use parameter 'view' to specify starting view point. Default is
% 'stereo'    - [Real] Create a stereo movie. The figure should contain a [left right]
%               display of two identical 3-D plots. The left plot view will follow the
%               given 'path' (see above). The right plot axis will be 3-D rotated by an
%               additional horizontal disparity angle specified by the 'stereo' argument:
%               6 (degrees) suggested. Default is [] = mono display.
%               'off'.
%Outputs to disk:
% imageX      - brainmovie3d() saves an output.avi movie (see 'moviename' option above)
%               and a sequence of image files to disk (image0001.eps, as define in the
%               'framesout' option).
%Example:
%
% % Given ICA activations in array icaact (size ncomps,nframes,ntrials), animate (here)
% % activity at/between two components at 176 points per epoch (from -100 ms to 600 ms
% % re stimulus onset) assuming a 250-Hz sampling rate and 100 output
% frames
%
% >> [ersps{1,1},itcs{1,1},powbase,TIMES,freqs] = ...                          % timef for
%                timef(icaact(1,:),176,[-100 600],'Component
%                1',250,1,32,100); %     1st comp
% >> [ersps{2,1},itcs{2,1},powbase,TIMES,freqs] = ...                          % timef for
%                timef(icaact(2,:),176,[-100 600],'Component 2',250,1,32,100); %     2nd comp
% >> [crossfs_amp{1,2},mcoh,TIMES,freqs,cohboot,crossfs_phase{1,2}] = ...      % crossf for
%      crossf_(icaact(1,:),icaact(2,:),176,[-100 600],'Crossf 1 and
%      2',250,1,32,100); % both
%
% >> brainmovie3d( ersps, itcs, crossfs_amp, crossfs_phase, TIMES, [1:2] );
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 30 Mai 2003
%
% Modified: Tim Mullen, SCCN/INC, UCSD 2010
%
% Note: Better resolution movies can be generated by .eps -> .ppm -> .avi,
%       (or, under a planned upgrade to brainmovie3d, from Matlab6 to .avi directly).
% >> !/usr/local/bin/convert images*.eps movie.mpg % ImageMagic 'convert' may be
%                                                  % used to generate the movie.

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2003

% This program is free software; you can redistribute it and/or
% modify it.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% $Log: brainmovie3d.m,v $
% Revision 1.24  2008/11/25 23:19:28  arno
% test for cirfactor
%
% Revision 1.23  2008/08/20 22:42:15  arno
% fixed last changes
%
% Revision 1.22  2008/08/20 22:28:37  arno
% Reworded documentation, fixed various problems
%
% Revision 1.21  2008/05/09 23:30:05  arno
% new coordinate format for dipoles
%
% Revision 1.20  2007/03/09 18:58:33  arno
% axes size and coherence thickness
%
% Revision 1.19  2007/03/06 18:56:01  arno
% *** empty log message ***
%
% Revision 1.18  2007/03/02 22:42:56  arno
% trying to fix the function
%
% Revision 1.17  2005/07/26 21:24:47  arno
% add option 'project3d' to function header
%
% Revision 1.16  2004/06/30 01:08:31  arno
% test if width>1280
%
% Revision 1.15  2004/06/29 23:31:19  arno
% recompute size
%
% Revision 1.14  2004/04/13 18:08:43  arno
% dipole projections
%
% Revision 1.13  2003/10/21 21:55:19  arno
% + for dB
%
% Revision 1.12  2003/10/21 19:24:30  arno
% new look for caption ...
%
% Revision 1.11  2003/10/09 01:01:34  arno
% fixing coordinate problem for multiple conditions
%
% Revision 1.10  2003/10/08 23:59:03  arno
% removing s to framesfolder
%
% Revision 1.9  2003/10/08 23:55:12  arno
% fix framefolder
%
% Revision 1.8  2003/10/08 23:21:15  arno
% 3d problem
%
% Revision 1.7  2003/10/08 22:54:16  arno
% updating for 3 conditions
%
% Revision 1.6  2003/09/18 22:47:06  arno
% adding frames output
%
% Revision 1.5  2003/08/06 14:42:19  arno
% *** empty log message ***
%
% Revision 1.4  2003/07/03 16:50:20  arno
% header typo
%
% Revision 1.3  2003/07/03 16:49:44  arno
% test for nan itc
%
% Revision 1.2  2003/07/03 16:07:22  arno
% disable axis redrawing
%
% Revision 1.1  2003/07/03 00:23:01  arno
% Initial revision
%

function [alltimepoints mov g] = brainmovie3d_causal(ALLERSP,ALLITC,ALLCROSSF,ALLCROSSFANGLE,TIMES,FREQS,SELECTED,varargin)

if nargin < 6
    help brainmovie3d_causal;
    return;
end;

% create structure for option if necessary
%-----------------------------------------
if ~isempty( varargin ),
    for index=1:length(varargin)
        if iscell(varargin{index})
            varargin{index} = { varargin{index}};
        end;
    end;
    g=struct(varargin{:});
else
    g= [];
end;

if nargin < 7
    SELECTED = 1:size(ALLERSP, 1);
end;

nbconditions = size(ALLERSP,2);
nbcomponents = size(ALLERSP,1);


try, g.dipplotopt;      catch, g.dipplotopt = {}; end
try, g.mri;             catch, g.mri = ''; end
try, g.figurehandle;    catch, g.figurehandle = []; end
try, g.opengl;          catch, g.opengl='on'; end
try, g.head; 			catch, g.head=''; end;
try, g.visible; 		catch, g.visible='on'; end;
try, g.square; 		    catch, g.square='on'; end;
try, g.moviename; 	    catch, g.moviename='output.avi'; end;
try, g.movieopts; 	    catch, g.movieopts={}; end;
try, g.rt; 	        	catch, g.rt={}; end;
try, g.power; 	    	catch, g.power='on'; end;
try, g.latency; 	   	catch, g.latency=[]; end;
try, g.itc; 	    	catch, g.itc='on'; end;
try, g.magnify; 	    catch, g.magnify=1; end;
try, g.crossf; 			catch, g.crossf='on'; end;
try, g.crossfcoh; 		catch, g.crossfcoh='on'; end;
try, g.size; 			catch, g.size=[400 400]; end;
try, g.crossfphasecolor;catch, g.crossfphasecolor='on'; end;
try, g.crossfphasespeed;catch, g.crossfphasespeed='off'; end;
try, g.crossfphaseunit; catch, g.crossfphaseunit='degree'; end;

try, g.nodeSizeDataRange;  catch, g.nodeSizeDataRange  = []; end;
try, g.nodeSizeLimits;     catch, g.nodeSizeLimits     = [0.1 1]; end;
try, g.nodeColorDataRange; catch, g.nodeColorDataRange = []; end;
try, g.nodeColorLimits;    catch, g.nodeColorLimits    = [0 1]; end;

try, g.edgeSizeDataRange;  catch, g.edgeSizeDataRange  = []; end;
try, g.edgeSizeLimits;     catch, g.edgeSizeLimits     = [0.1 1]; end;
try, g.edgeColorDataRange; catch, g.edgeColorDataRange = []; end;
try, g.edgeColorLimits;    catch, g.edgeColorLimits    = [0 1]; end;

try, g.diskscale;       catch, g.diskscale = 1; end;
try, g.framefolder;     catch, g.framefolder = ''; end;
try, g.envelope;        catch, g.envelope = []; end;
try, g.caption;			catch, g.caption = 'on'; end;
try, g.frames;			catch, g.frames = []; end;
try, g.envvert;			catch, g.envvert = {}; end;
try, g.flashes;			catch, g.flashes = []; end;
try, g.polarity;		catch, g.polarity = 'pos'; end;
try, g.framesout;	    catch, g.framesout = 'tiff'; end;
try, g.condtitle;		catch, g.condtitle = []; end;
try, g.condtitleformat;	catch, g.condtitleformat = {'fontsize', 14', 'fontweight', 'bold' }; end;
try, g.title;			catch, g.title = []; end;
try, g.envylabel;		catch, g.envylabel = '\muV'; end;
try, g.plotorder;       catch, g.plotorder = SELECTED; end;
try, g.coordformat;     catch, g.coordformat = 'spherical'; end;
try, g.stereo;          catch, g.stereo = []; end;
try, g.backcolor;       catch, g.backcolor = [0 0 0]; end;
try, g.path3d;          catch, g.path3d = 'off'; end;
try, g.project3d;       catch, g.project3d = 'off'; end;
try, g.view;            catch, g.view = [43.6650 30.4420]; end;
try, g.LONImesh;        catch, g.LONImesh = []; end                         %% TM
try, g.LONITransparency;catch, g.LONITransparency = 1; end                  %% TM
try, g.causality;       catch, g.causality = 0; end                         %% TM: added
try, g.nodelabels;      catch, g.nodelabels = {}; end                       %% TM: added
try, g.collapsefun;     catch, g.collapsefun = 'mean'; end                  %% TM: added  (can be 'integrate','mean','max','absmax','peak')
try, g.TIMES;           catch, g.TIMES = []; end                            %% TM: added (for envelope plot)
try, g.showLatency;        catch, g.showLatency = 0; end                    %% TM: added
try, g.verb;            catch, g.verb = 0;  end                             %% TM: added (display progress bar for making movie)
try, g.csf;             catch, g.csf = []; end                              %% TM
try, g.cortexTransparency;  catch, g.cortexTransparency = 1; end            %% TM: transparency of superimposed cortex (1 = don't plot)
try, g.speedy;          catch, g.speedy = true; end                         %% TM: for fast rendering -- some features disabled
try, g.mode;        catch, g.mode = 'init_and_render'; end                  %% TM: initialization mode. if 'init', then initialize only and current state of BM will be returned in g.vars, pass g in as subsequent input. If 'render', then skip init and use g.vars for rendering. if 'init_and_render' then do both (overwrites g.vars)
try, g.vars;            catch, g.vars = []; end                             %% TM: a structure to hold initialization variables.
try, g.EdgeColorMappedToDirectionality; catch, g.EdgeColorMappedToDirectionality = false; end

% some parameters for captions
try, g.nodeColorMapping; catch, g.nodeColorMapping = ''; end
try, g.EdgeColorMapping; catch, g.EdgeColorMapping = ''; end
try, g.NodeSizeMapping; catch, g.NodeSizeMapping = ''; end
try, g.EdgeSizeMapping; catch, g.EdgeSizeMapping = ''; end
try, g.ConnMethod; catch, g.ConnMethod = 'Connectivity'; end
try, g.nodeColorPolarity; catch, g.nodeColorPolarity = 'pos'; end
try, g.edgeColorPolarity; catch, g.edgeColorPolarity = 'pos'; end
try, g.centerDataRange;   catch, g.centerDataRange = true; end
try, g.colorshadow;     catch, g.colorshadow = 1; end
try, g.nodeColormap;       catch,
    colormtmp = hot(64);
    colormtmp(end,3) = (colormtmp(end,3)+colormtmp(end-1,3))/2; % white does not come out when the
    g.nodeColormap = colormtmp;                                    % the figure is printed to ppm
    g.nodeColormap(:,1) =  colormtmp(:,2);
    g.nodeColormap(:,2) =  colormtmp(:,3);
    g.nodeColormap(:,3) =  colormtmp(:,1);
    g.nodeColormap = [ g.nodeColormap; colormtmp(end:-1:1,:)];
    g.nodeColormap = jet(64);
end;
try, g.edgeColormap; catch,
    g.edgeColormap = jet(64);
    %g.edgeColormap = hsv(64);
    %g.edgeColormap = [ g.edgeColormap(55:end,:);
    %g.edgeColormap(1:54,:)]; g.edgeColormap = g.edgeColormap(linspace(64, 1, 64),:); % reorganize the colormap
    %g.edgeColormap = hsv(64);
    %g.edgeColormap = [g.edgeColormap(16:end,:); g.edgeColormap(1:5,:)];
end;
try, g.xlimaxes; 		catch, g.xlimaxes = [-1 1]; end;
try, g.ylimaxes; 		catch, g.ylimaxes = [-1 1]; end;
try, g.rthistloc; 	    catch, g.rthistloc(1) = (g.xlimaxes(2)-g.xlimaxes(1))*0.74 + g.xlimaxes(1); % abscicia
    g.rthistloc(3) = (g.xlimaxes(2)-g.xlimaxes(1))*0.1; % width
    g.rthistloc(2) = (g.ylimaxes(2)-g.ylimaxes(1))*0.34 + g.ylimaxes(1); % ordinate
    g.rthistloc(4) = (g.ylimaxes(2)-g.ylimaxes(1))*0.1; % max height
end;
try, g.coordinates; catch,
    % coordinates around a circle
    g.coordinates = zeros( nbcomponents, 3 );
    count = 0;
    for index = SELECTED
        if length(SELECTED) > 1
            g.coordinates( index,:) = [ cos(count/length(SELECTED)*2*pi) sin(count/length(SELECTED)*2*pi) 0 ] * 0.7;
        else	g.coordinates(index,:) = [ 0.01 0.01 0 ];
        end;
        count = count + 1;
    end;
end;
try, g.circfactor; catch, g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
if isempty(g.circfactor), g.circfactor = ones( nbcomponents, nbcomponents )*0.01; end;
if isstr(g.path3d)
    switch g.path3d
        case 'on', g.path3d = [ 1 0.75];
        case 'off', g.path3d = [];
    end;
else
    if length(g.path3d) ~= 2, error('path3d length have to be a string or a 2 element vector'); end;
end;

% add defaults ranges
%--------------------
if isempty(g.nodeSizeDataRange)
    g.nodeSizeDataRange = [Inf -Inf];
    for i=1:length(ALLERSP)
        g.nodeSizeDataRange(1) = min(g.nodeSizeDataRange(1), min(ALLERSP{i}(:)));
        g.nodeSizeDataRange(2) = max(g.nodeSizeDataRange(2), max(ALLERSP{i}(:)));
    end;

    % make 0 in the center of colormap
    if g.centerDataRange
        g.nodeSizeDataRange(1) = -max(g.nodeSizeDataRange);
        g.nodeSizeDataRange(2) =  max(g.nodeSizeDataRange);
    end
    if g.verb, fprintf('Node size data range automatically set to %1.2f to %1.2f\n', g.nodeSizeDataRange(1), g.nodeSizeDataRange(2)); end
end
if isempty(g.nodeColorDataRange)
    g.nodeColorDataRange = [Inf -Inf];
    for i=1:length(ALLITC)
        g.nodeColorDataRange(1) = min(g.nodeColorDataRange(1), min(ALLITC{i}(:)));
        g.nodeColorDataRange(2) = max(g.nodeColorDataRange(2), max(ALLITC{i}(:)));
    end;
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.nodeColorDataRange(1) = -max(g.nodeColorDataRange);
        g.nodeColorDataRange(2) =  max(g.nodeColorDataRange);
    end
    
    if g.verb, fprintf('Node color data range automatically set to %1.2f to %1.2f\n', g.nodeColorDataRange(1), g.nodeColorDataRange(2)); end
end;
if isempty(g.edgeSizeDataRange)
    g.edgeSizeDataRange = [Inf -Inf];
    for i=1:length(ALLCROSSF(:))
        if ~isempty(ALLCROSSF{i})
            g.edgeSizeDataRange(1) = min(g.edgeSizeDataRange(1), min(ALLCROSSF{i}(:)));
            g.edgeSizeDataRange(2) = max(g.edgeSizeDataRange(2), max(ALLCROSSF{i}(:)));
        end;
    end;
    
    % make 0 in the center of colormap
    if g.centerDataRange
        g.edgeSizeDataRange(1) = -max(g.edgeSizeDataRange);
        g.edgeSizeDataRange(2) =  max(g.edgeSizeDataRange);
    end
    if g.verb, fprintf('Edge size data range automatically set to %1.2f to %1.2f\n', g.edgeSizeDataRange(1), g.edgeSizeDataRange(2)); end
end;
if isempty(g.edgeColorDataRange)
    g.edgeColorDataRange = [Inf -Inf];
    for i=1:length(ALLCROSSFANGLE(:))
        if ~isempty(ALLCROSSFANGLE{i})
            g.edgeColorDataRange(1) = min(g.edgeColorDataRange(1), min(ALLCROSSFANGLE{i}(:)));
            g.edgeColorDataRange(2) = max(g.edgeColorDataRange(2), max(ALLCROSSFANGLE{i}(:)));
        end;
    end;
    
     % make 0 in the center of colormap
    if g.centerDataRange
        g.edgeColorDataRange(1) = -max(g.edgeColorDataRange);
        g.edgeColorDataRange(2) =  max(g.edgeColorDataRange);
    end
    
    if g.verb, fprintf('Edge color data range automatically set to %1.2f to %1.2f\n', g.edgeColorDataRange(1), g.edgeColorDataRange(2)); end
end;

% check size of inputs
% --------------------
try
    if ~all(size(ALLERSP) == size(ALLITC))
        disp('Error: ERSP and ITC cells array must be the same size'); return;
    end;
    if ~isempty(ALLCROSSF)
        if ~all(size(ALLCROSSF) == size(ALLCROSSFANGLE))
            disp('Error: Crossf amplitude and Crossf angle cells array must be the same size'); return;
        end;
        if ~(size(ALLCROSSF,2) == size(ALLERSP,1))
            disp('Error: number of components different in ERSP and Crossf arrays'); return;
        end;
        if ~(size(ALLCROSSF,3) == size(ALLERSP,2))
            disp('Error: number of conditions different in ERSP and Crossf arrays'); return;
        end;
        if ~(size(ALLCROSSF{1,2,1},1) == size(ALLERSP{1,1},1))
            disp('Error: number of frequencies (rows) different in ERSP and Crossf arrays'); return;
        end;
        if ~(size(ALLCROSSFANGLE{1,2,1},2) == size(ALLITC{1,1},2))
            disp('Error: number of time points (columns) different in ERSP and Crossf arrays'); return;
        end;
        if ~(size(ALLCROSSF{1,2,1},2) == length(TIMES))
            disp('Error: number of time points (columns) different in TIMES and Crossf arrays'); return;
        end;
    end;
    try tmp = ALLERSP{1,1}; tmp(FREQS,:); catch, disp('Error: unable to access the defined frequencies in ERSPs (out of bounds) '); return; end;
    try ALLERSP{SELECTED,1}; catch, disp('Error: unable to access the defined components in ERSPs (out of bounds)'); return; end;
catch
    disp('Error accessing one of the variable. Remember: Except for SELECTED, freqs, TIMES and circfactor, all vars are cell arrays. Check also: dimensions and content.'); return;
end;

% check structure content
% -----------------------
if ~isempty(g.rt)
    if length(g.rt) ~= nbconditions
        disp('Error: Rt must be either an array of the size of the number of conditions (might be 0 for some conditions)'); return;
    end;
end;
switch lower(g.visible)
    case {'on', 'off'} ;
    otherwise disp('Error: Visibility must be either ''on'' or ''off'''); return;
end;
switch lower(g.square)
    case {'on', 'off'} ;
    otherwise disp('Error: Square must be either ''on'' or ''off'''); return;
end;
switch lower(g.power)
    case {'on', 'off'} ;
    otherwise disp('Error: Power must be either ''on'' or ''off'''); return;
end;
switch lower(g.itc)
    case {'on', 'off'} ;
    otherwise disp('Error: Itc must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossf)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossf must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossfcoh)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfcoh must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossfphasecolor)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfphasecolor must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossfphasespeed)
    case {'on', 'off'} ;
    otherwise disp('Error: Crossfphasespeed must be either ''on'' or ''off'''); return;
end;
switch lower(g.crossfphaseunit)
    case {'degree', 'radian'} ;
    otherwise disp('Error: Crossfphaseunit must be either ''degree'' or ''radian'''); return;
end;
switch lower(g.caption)
    case {'on', 'off'} ;
    otherwise disp('Error: Caption must be either ''on'' or ''off'''); return;
end;
switch lower(g.polarity)
    case {'pos', 'posneg'} ;
    otherwise disp('Error: Polarity must be either ''pos'' or ''posneg'''); return;
end;
if ~isempty(g.envvert),
    if ~iscell(g.envvert) && ~( isstruct(g.envvert{1}) || isnumeric(g.envvert{1}) )
        disp('Error: Invalid type for Envvert.'); return;
    end
end
if ~isempty(g.latency) && ~isnumeric(g.latency)
    disp('Error: Latency must be a vector'); return;
end;
if length(g.nodeSizeDataRange) ~= 2
    disp('Error: Scalepower must be a 2-element array'); return;
end;
if length(g.edgeSizeDataRange) ~= 2
    disp('Error: Scalecoher must be a 2-element array'); return;
end;
if (length(g.diskscale) ~= 1 || g.diskscale < 0)
    disp('Error: Diskscale must be a scalar value >= 0.'); return;
end
if size(g.nodeColormap,2) ~= 3
    disp('Error: Colmapcoh must be a colormap (3 columns)'); return;
end;
if size(g.edgeColormap,2) ~= 3
    disp('Error: Colmapcrossf must be a colormap (3 columns)'); return;
end;
if size(g.circfactor,1) ~= size(g.circfactor,2)
    disp('Error: Circfactor must be a square matrix'); return;
end;
if ~iscell(g.coordinates) && ~isempty(g.circfactor)
    if size(g.circfactor,1) ~= size(g.coordinates,1)
        disp('Error: Circfactor must have the same number of rows as the number of rows of coordinates'); return;
    end;
    if nbcomponents ~= size(g.coordinates,1)
        disp('Error: The array of SELECTED components must have length nrows of the array coordinates'); return;
    end;
end;
if ~isstr(g.envylabel)
    disp('Error: envelope label must be a string'); return;
end;
if ~isempty(g.envelope) && isempty(g.TIMES)
    if (size( g.envelope,1 ) ~=2) || (size( g.envelope,2 ) ~= length(TIMES)) || (size( g.envelope,3 ) ~= nbconditions)
        fprintf('Error: Envelope array does not have the right size (%s), i.v.henv. 2 x %d (number of time points) x %d (number of conditions)\n', int2str(size( g.envelope)), length(TIMES), nbconditions); return;
    end;
end;
if ~isempty(g.condtitle)
    if iscell(g.condtitle), g.condtitle = strvcat(g.condtitle{:}); end;
    if ~isempty(g.condtitle) && size( g.condtitle,1 ) ~= nbconditions
        fprintf('Error: The number of rows in the title array(%d) must match the number of conditions (%d)\n', size(g.condtitle,1), nbconditions); return;
    end;
end;
if length(g.plotorder) ~= length(SELECTED)
    error([ 'Error: ''plotorder'' must be the same size as the number of SELECTED components:' int2str(length(SELECTED)) ]);
end;
if max(g.plotorder) > max(SELECTED)
    error([ 'Error: ''plotorder'' must be below the number of SELECTED components:' int2str(max(SELECTED)) ]);
end;

SLASH = fastif(isunix,'/','\');

if ~isempty(g.framefolder) && ~isdir(g.framefolder)
    error('%s does not exist',g.framefolder);
%     [tmp1 tmp2] = mkdir(SLASH, g.framefolder(2:end) );
end;

% append a slash to framefolder, if necessary
if ~isempty(g.framefolder) && g.framefolder(end) ~= SLASH 
    g.framefolder = [g.framefolder SLASH]; 
end;
   
if g.verb==2
    g.vars.hwaitbar=waitbar(0,'Initializing Brainmovie...');
end

mov = [];
alltimepoints = [];

% other variables
% ---------------
%limits: power -6 to 6
%limits: ITC 0-1
%limits: coherence 0-1
%limits: coherence angle -180 to 180
g.factproj = [-88 88 -71]; %[-71 88 -71];
g.projcolor = [0.35 0.35 0.35];
g.rthistcolor  = [1 1 1];
g.resmult = 1;
    
if ismember_bc(lower(g.mode),{'init','init_and_render'})

    % create movie
    % ------------
    if ~isempty(g.moviename)
    %     fprintf('A movie is being saved under %s (movie parameters shown below):',g.moviename);
        mov = avifile(g.moviename, g.movieopts{:});
    end


    currentphase   = zeros( length(SELECTED), length(SELECTED), nbconditions);
    tmp = ALLERSP{1,1};
    nwin = size(tmp,2);

    %for index = 1:64
    %	circle(1+index,1, 0.5, g.colormaphsv(index, :));
    %end;

    % optional resqure of all coordinates
    % -----------------------------------
    g.magnify = g.magnify/4;

    % compute RT distribution
    % -----------------------
    if ~isempty(g.rt)
        RTdist = zeros(nbconditions,nwin);
        for index = 1:nbconditions
            if ~isempty(g.rt{index})
                timestep = (TIMES(2)-TIMES(1))/2;
                for indeximage = 1:nwin
                    RTdist(index, indeximage) = length( intersect_bc( find( g.rt{index} > TIMES(indeximage)-timestep ) , ...
                        find(  g.rt{index} <= TIMES(indeximage)+timestep ) ) );
                end;
                RTdist(index,:) = RTdist(index,:)/max(RTdist(index,:));
            end;
        end;
        RTdist = RTdist/max(RTdist(:));
    end;

    if ~isempty(g.figurehandle)
        try 
            figure(g.figurehandle);
            clf
        catch
            g.figurehandle = figure( 'position', [100, 100, ceil(nbconditions*g.size(1)/4)*4, ceil(g.size(2)/4)*4], ...
                'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible,'tag','BrainMovieFigure'); %'paperorientation', 'landscape' );
        end
    else
       g.figurehandle = figure( 'position', [100, 100, ceil(nbconditions*g.size(1)/4)*4, ceil(g.size(2)/4)*4], ...
            'PaperPositionMode', 'auto', 'papertype', 'A1', 'visible',g.visible,'tag','BrainMovieFigure'); %'paperorientation', 'landscape' );
    end

    axis off
    if 	strcmpi(g.framesout, 'ppm')
        r = 0.8465;
        pos = get(gcf,'position');
        if floor(pos(3)/r)> 1280
            fact = 1280/(pos(3)/r);
            set(gcf, 'position', [ 0 0 1280  floor(pos(4)/r*fact) ]);
        else
            set(gcf, 'position', [ 0 0 floor(pos(3)/r), floor(pos(4)/r) ]);
        end;
    end;
    pos = get(gca,'position');
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)];

    % compute SELECTED latency point
    % ------------------------------
    if ~isempty(g.latency)
        g.vars.alltimepoints = [];
        for index = 1:length(g.latency)
            [tmp tmptimepoint] = min(abs(g.latency(index)-TIMES));
            g.vars.alltimepoints = [ g.vars.alltimepoints tmptimepoint];
        end;
    else
        if isempty(g.frames)
            g.vars.alltimepoints = 1:nwin;
        else
            g.vars.alltimepoints = g.frames;
        end;
    end;

    % make black patch behind figure
    % ------------------------------
    g.vars.hback = axes('position' , [0 0 1 1], 'xtick', [], 'ytick', [], 'box', 'off','tag','figureBackground');
    g.vars.hpatch = patch([0 1 1 0], [0 0 1 1], g.backcolor); xlim([0 1]); ylim([0 1]);
    set(g.vars.hpatch, 'facecolor' , g.backcolor, 'edgecolor', 'none');

    % compute flashes latency
    % -----------------------
    if ~isempty(g.flashes)
        if iscell(g.flashes)
            for index = 1:length(g.flashes)
                flasheslat(index) = g.flashes{index}{1};
                g.vars.flashescol{index} = g.flashes{index}{2};
            end;
        else
            flasheslat = g.flashes;
            for index = 1:length(g.flashes)
                g.vars.flashescol{index} = [0.5 0.5 0.5];
            end;
        end;
        g.vars.allflashes = [];
        for index = 1:length(g.flashes)
            [tmp tmptimepoint] = min(abs(flasheslat(index)-TIMES));
            g.vars.allflashes = [ g.vars.allflashes tmptimepoint];
        end;
        %g.vars.hpatch = patch([ 0.02 .11 .11 0.02], [0.05 0.05 0.925 0.925], [0.5 0.5 0.5]); lateral
        %g.vars.hpatch = patch([ 0 1 1 0], [0 0 1 1], [0.5 0.5 0.5]); full
        %g.vars.hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.92 0.92 1 1], [0.5 0.5 0.5]); %up
        g.vars.hpatch = patch([ 0.13 0.84 0.84 0.13 ], [0.8 0.8 0.93 0.93], [0.5 0.5 0.5]);
        set(g.vars.hpatch, 'facecolor', 'w', 'edgecolor', 'none');
        xlim([0 1]); ylim([0 1]);
        posf = 0; % used as a counter to preserve color
    end;

    % draw axes and display images
    % ----------------------------
    ordinate = 0.2;
    max_ordinate = 1-1.4*ordinate;   % makes space at top for figure title
    maxcoordx    = 1.1-1/nbconditions/4;
    coords = g.coordinates;
    g.coordinates = cell(nbconditions);
    for tmpcond=1:nbconditions

        % plot 3d head (*0.9 added for Nick - Arno).
        % ------------
        g.vars.hh(tmpcond) = axes('position', [0+maxcoordx/nbconditions*(tmpcond-1), ordinate, maxcoordx/nbconditions*0.9, max_ordinate].*s+q );
        gr = [ 0.3 0.3 0.3 ];
        g.dipplotopt = [{ 'coordformat' g.coordformat 'gui', 'off', 'cornermri', 'on', 'color', { gr gr gr gr gr gr gr gr gr } } g.dipplotopt];

        if ~isempty(g.mri)
            g.dipplotopt = [g.dipplotopt 'mri' g.mri];
        end

        if iscell(coords)
            for index = 1:length(coords)
                if size(coords{index},1) == 2 && all(coords{index}(2,:) == 0), coords{index}(2,:) = []; end;
                dipstruct(index).posxyz = coords{index};
                dipstruct(index).momxyz = [1 1 1];
                if size(dipstruct(index).posxyz,1) == 2, dipstruct(index).momxyz(2,:) = [1 1 1]; end;
                dipstruct(index).component = index;
                dipstruct(index).rv = 0.1;
            end;
        else
            for index = 1:size(coords, 1);
                dipstruct(index).posxyz = coords(index,:);
                dipstruct(index).momxyz = [0 0 0];
                dipstruct(index).component = index;
                dipstruct(index).rv = 0.1;
            end;
        end;

        dipplot( dipstruct, 'view', g.view, g.dipplotopt{:}); axis off;

        %     g.maxc = 100;
        %     surface([-2 -2; -2 -2]*g.maxc, [-20 20; -20 20]*g.maxc,[-20 -20; 20 20]*g.maxc, repmat(reshape([0 0 0], 1, 1, 3), [2 2 1]), 'facelighting', 'none');
        %     surface([-20 20; -20 20]*g.maxc,[2 2; 2 2]*g.maxc, [-20 -20; 20 20]*g.maxc,     repmat(reshape([0 0 0], 1, 1, 3), [2 2 1]), 'facelighting', 'none');

        %camproj('perspective');
        brainaxis(tmpcond) = gca;
        set(gca, 'cameraviewanglemode', 'manual'); % disable change size
        set(gca,'tag',['brain' num2str(tmpcond)]);
        axis vis3d % same as above (for security)
        camlight left
        camlight right
        view(g.view)
        %camzoom(1.2)

        g.vars.coordinates = g.coordinates;
        
        for index = 1:length(dipstruct)
            htmp = findobj(gca, 'tag', [ 'dipole' int2str(index) ]);
            for dipindex = 1:length(htmp)
                tmpstruct = get(htmp(dipindex), 'userdata');
                if isstruct(tmpstruct) % look for dipole location % THIS DOES NOT WORK
                    if isfield(tmpstruct, 'pos3d') && ~all(tmpstruct.pos3d == 0)
                        if length(g.vars.coordinates{tmpcond}) >= index && ~isempty(g.vars.coordinates{tmpcond}{index}) && ~isequal(g.vars.coordinates{tmpcond}{index}, tmpstruct.pos3d)
                            g.vars.coordinates{tmpcond}{index}(2,:) = tmpstruct.pos3d;
                        else
                            g.vars.coordinates{tmpcond}{index} = tmpstruct.pos3d;
                        end;
                    elseif isfield(tmpstruct, 'eleccoord') && ~all(tmpstruct.eleccoord == 0)
                        if length(g.vars.coordinates{tmpcond}) >= index && ~isempty(g.vars.coordinates{tmpcond}{index}) && ~isequal(g.vars.coordinates{tmpcond}{index}, tmpstruct.eleccoord)
                            g.vars.coordinates{tmpcond}{index}(2,:) = tmpstruct.eleccoord;
                        else
                            g.vars.coordinates{tmpcond}{index} = tmpstruct.eleccoord;
                        end;
                    else
                        tmpstruct
                        error('Field not found in tmpstruct');
                    end;
                end;
            end;
            delete(htmp);
        end;

        %h = plot3(g.vars.coordinates{tmpcond}(:, 1),  g.vars.coordinates{tmpcond}(:, 2),  g.vars.coordinates{tmpcond}(:, 3), 'r.', 'markersize', 30);
        %dsaf
        xltmp = xlim;
        yltmp = ylim;
        g.vars.dimratio = (xltmp(2) - xltmp(1)) / (yltmp(2) - yltmp(1));

        axis off;
        if ~isempty(g.condtitle)
            h = title(g.condtitle(tmpcond,:));
            if ~isempty(g.condtitleformat)
                set(h, g.condtitleformat{:} );
            end;
        end;

        % this axis is used for the enveloppe but
        % also used to print current time (which is why it is always created
        if ~isempty( g.envelope )
            g.vars.henv(tmpcond) = axes('position', [0.1/nbconditions+maxcoordx/nbconditions*(tmpcond-1), 0, ...
                maxcoordx/nbconditions-0.1/nbconditions, ordinate-0.1].*s+q,'visible', g.visible);
        end;
    end;


    % adjust the colormap as needed
    if strcmpi(g.nodeColorPolarity, 'pos')
        g.nodeColormap = g.nodeColormap(length(g.nodeColormap)/2:end,:); end
    if strcmpi(g.edgeColorPolarity, 'pos')
        g.edgeColormap = g.edgeColormap(length(g.edgeColormap)/2:end,:); end


    % draw captions if necessary
    % --------------------------
    countl = 1;
    switch lower(g.caption)
        case 'on' ,
            xlimnorm = (1.1-maxcoordx)/(maxcoordx/nbconditions) * g.xlimaxes;
            ylimnorm = 0.45/(1-ordinate) * g.ylimaxes;
            switch g.power, case 'on',
                % create 3 spheres to show node size variation
                g.vars.hlgnd(countl) = axes('position', [maxcoordx, -0.1,    (1.1-maxcoordx), 0.45].*s+q, 'xlim', xlimnorm, ...
                    'ylim', ylimnorm,'visible', g.visible, 'color', 'w' );
                % draw 3 spheres
                [xstmp ystmp zs] = sphere(15);
                l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
                normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
                tmpsize = 0.5; xs1 = tmpsize*ystmp; ys1 = tmpsize*xstmp; zs1 = tmpsize*zs;
                tmpsize = 0.9; xs2 = tmpsize*ystmp; ys2 = tmpsize*xstmp; zs2 = tmpsize*zs + 2;
                tmpsize = 0.1; xs3 = tmpsize*ystmp; ys3 = tmpsize*xstmp; zs3 = tmpsize*zs - 1.5;
                colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zs,1) size(zs,2) 1]);
                handles = surf(xs1, ys1, zs1, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                    'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3); hold on;
                handles = surf(xs2, ys2, zs2, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                    'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
                handles = surf(xs3, ys3, zs3, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                    'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
                axis off;
                camlight left
                camlight right
                view([1 0 0])
                lightangle(45,0);
                lighting phong;
                material shiny;
                axis equal;
                set(gca, 'zlim', [-2 4]);
                set(gca, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
                text(0.1, 0.42, {'',g.nodeSizeMapping}, 'units','normalized','horizontalalignment','center','rotation',90,'fontsize', 11, 'fontweight', 'bold');
                text(0, 1.3, 2, num2str(g.nodeSizeDataRange(2),'%0.2f'), 'fontweight', 'bold');    % upper limit
                text(0, 1, 0, num2str((g.nodeSizeDataRange(2)+g.nodeSizeDataRange(1))/2,'%0.2f'), 'fontweight', 'bold');   % midrange
                text(0, 0.5, -1.5, num2str(g.nodeSizeDataRange(1),'%0.2f'), 'fontweight', 'bold'); % lower limit
                %       scalepower(mean(xlimnorm), min(ylimnorm)+0.2, g); % see function at the end
                %       axis off;
                countl = countl + 1;
            end;
            switch g.itc, case 'on',
                % create ball colormap for nodecolor legend
                g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.29 , (1.1-maxcoordx)/2, 0.14].*s+q, ...
                    'visible', g.visible, 'color', 'none' );
                countl = countl + 1;

                if strcmpi(g.nodeColorPolarity, 'posneg') % negative ITCs (difference only) ?
                    cbar( [-1 1], [-1 1], g.nodeColormap, 'vert', 'circle', g);
                    %           ylabel(g.nodeColorMapping, 'fontweight', 'bold');
                    %           set(gca, 'ytick', [-1 0 1], 'yticklabel', [-1 0 1], 'xticklabel', [], 'box', 'off');
                else
                    cbar( [0 1], [0 1], g.nodeColormap, 'vert', 'circle', g);
                    %           cbar( [0 1], [0 1], g.nodeColormap(length(g.nodeColormap)/2:end,:), 'vert', 'circle', g);
                    %           ylabel({'',g.nodeColorMapping},'rotation',90, 'fontweight', 'normal');
                    %           set(gca, 'ytick', [0 1], 'yticklabel', [0 1], 'xticklabel', [], 'box', 'off');
                    %           set(gca, 'ytick', [0], 'yticklabel', [0], 'xticklabel', [], 'box', 'off');
                end;
                axis off;
                set(gca, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
                % yticks
                text(-0.2, 1, num2str(g.nodeColorDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized'); % upper limit
                text(-0.2, 0.5, num2str((g.nodeColorDataRange(2)+g.nodeColorDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized'); % midrange
                text(-0.2, 0, num2str(g.nodeColorDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized'); % lower limit
                text(-0.8, 0.55, {'',g.nodeColorMapping}, 'horizontalalignment','center','rotation',90,'fontsize', 11, 'fontweight', 'bold','units','normalized');
            end;

            switch g.crossf, case 'on',
                % create bars showing variability in edge size
                g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.69,(1.1-maxcoordx)/2, 0.25 ].*s+q, ...
                    'visible', g.visible );
                countl = countl + 1;
                scalecoher([0.02 1], [0.04 0.96], 5, g); % see function at the end
                text(-0.2, 1, num2str(g.edgeSizeDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized'); % upper limit
                text(-0.2, 0.5, num2str((g.edgeSizeDataRange(2)+g.edgeSizeDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized'); % midrange
                text(-0.2, 0, num2str(g.edgeSizeDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized'); % lower limit
                text(-0.8, 0.55, {'',g.edgeSizeMapping}, 'horizontalalignment','center','rotation',90,'fontsize', 11, 'fontweight', 'bold','units','normalized');
                set(gca, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');

                switch g.crossfphasecolor, case 'on',
                    % create colorbar for edgecolor legend
                    g.vars.hlgnd(countl) = axes('position', [maxcoordx+(1.1-maxcoordx)/2, 0.47 , (1.1-maxcoordx)/4, 0.14].*s+q, ...
                        'visible', g.visible, 'color', 'none' );
                    countl = countl + 1;
                    if strcmpi(g.edgeColorPolarity, 'posneg') % negative ITCs (difference only) ?
                        cbar( [-1 1], [-1 1], g.edgeColormap, 'vert', '', g);
                        %           ylabel({'', g.ConnMethod }, 'fontweight', 'bold');
                        %           set(gca, 'ytick', [-1 0 1], 'yticklabel', [g.edgeSizeDataRange(1) 0 g.edgeSizeDataRange(2)], 'xticklabel', []);
                        %         set(gca, 'ytick', [-1 0 1], 'yticklabel', [g.edgeSizeDataRange(1) 0 g.edgeSizeDataRange(2)], 'xticklabel', []);
                    else
                        cbar( [0 1], [0 1], g.edgeColormap, 'vert', '', g);
                        %               cbar( [0 1], [0 1], g.edgeColormap(length(g.edgeColormap)/2:end,:), 'vert', '', g);
                        %           ylabel({'', g.ConnMethod }, 'fontweight', 'bold');
                        %           set(gca, 'ytick', [0 1], 'yticklabel', [g.edgeSizeDataRange(1) g.edgeSizeDataRange(2)], 'xticklabel', []);
                        %         set(gca, 'ytick', [0], 'yticklabel', [0], 'xticklabel', []);
                    end;
                    text(-0.5, 1, num2str(g.edgeColorDataRange(2),'%0.2f'), 'fontsize', 10,'units','normalized'); % upper limit
                    text(-0.5, 0.5, num2str((g.edgeColorDataRange(2)+g.edgeColorDataRange(1))/2, '%0.2f'), 'fontsize', 10,'units','normalized'); % midrange
                    text(-0.5, 0, num2str(g.edgeColorDataRange(1),'%0.2f'), 'fontsize', 10,'units','normalized'); % lower limit
                    text(-1.5, 0.55, {'',g.edgeColorMapping}, 'horizontalalignment','center','rotation',90,'fontsize', 11, 'fontweight', 'bold','units','normalized');

                    set(gca, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'off');
                end;
            end;
        case 'off', maxcoordx = 1;
    end;

    % specify the connectivity method used
    text(0.8927, 0.9459, g.ConnMethod,'fontsize',12,'units','normalized','fontweight','bold','parent',findobj(gcf,'tag','figureBackground'));


    % draw white axis on envelop if flashes DOES NOT WORK WHEN PRINTING IN EPS
    % -------------------------------------
    %if ~isempty(g.flashes)
    %	if ~isempty(g.) % draw axis for the envelope
    %		eflash = axes('position', [0 0 maxcoordx-0.1 ordinate].*s+q, ...
    %					  'xtick', [], 'ytick', [], 'box', 'off', 'visible', g.visible, 'color', 'none');
    %		hpatch2 = patch([ 0 1 1 0], [0 0 1 1], [0.5 0.5 0.5]); set(hpatch2, 'facecolor', 'w', 'edgecolor', 'none');
    %	end;
    %end;

    % scan time windows
    % -----------------
    if strcmpi(g.opengl,'on')
        set(gcf, 'renderer', 'opengl');
    end
    
end % initialization block

alltimepoints = g.vars.alltimepoints;

% return now if we are only initializing the brainmovie
if strcmpi(g.mode,'init')
    return;
end


% Main loop, draw frames
% ----------------------------
for indeximage = g.vars.alltimepoints
    
    switch g.verb
        case 1
            fprintf('Processing image %d\n', indeximage);
        case 2
            waitbar(indeximage/length(g.vars.alltimepoints),g.vars.hwaitbar,sprintf('Rendering timepoint (%d/%d)...',indeximage,length(g.vars.alltimepoints)));
    end
    
    
    % invert background if necessary
    % ------------------------------
    if ~isempty(g.flashes)
        %axes(g.vars.hback); set (gcf, 'visible', g.visible);
        if ~isempty(find(indeximage == g.vars.allflashes, 1))
            posf = find(indeximage == g.vars.allflashes);
            set(g.vars.hpatch, 'facecolor', g.vars.flashescol{posf});
        elseif posf == 0 % allow the color to stay 2 images
            set(g.vars.hpatch, 'facecolor', 'k');
        else
            posf = 0;
        end;
    end;
    
    for tmpcond=1:nbconditions
        g.vars.axhandle = g.vars.hh(tmpcond);
        
%         axes(g.vars.hh(tmpcond)); 
        if ~g.speedy
            set (gcf, 'visible', g.visible);
        end
        
        
        % clean images and update view
        % ----------------------------
        if ~isempty(g.path3d)
            angle = (indeximage-1)/length(g.vars.alltimepoints)*360;
            camorbit( cos(angle/180*pi)*g.path3d(1), sin(angle/180*pi)*g.path3d(2) );
        end;
        
        % get handle to old objects (we'll delete them later)
        oldobjs = findobj( g.vars.hh(tmpcond), 'tag', 'tmpmov');
%         set (gcf, 'visible', g.visible);
        
        % draw correlations
        % -----------------
        switch lower(g.crossf), case 'on',
            for index1 = SELECTED
                for index2 = SELECTED
                    
                    if index1==index2 || (~g.causality && index2 < index1)
                        continue;     % need only do one arc for symmetric connectivity (g.vars.henv.g., cohere)
                    end
                    
                    %% DEBUG
                    %                       if ismember_bc(index2,[3 5 9]), continue; end
                    
                    
                    tmpcrossfpow = ALLCROSSF     	 { index1, index2, tmpcond };
                    tmpcrossfang = ALLCROSSFANGLE    { index1, index2, tmpcond };
                    
                    % collapse data across frequency
                    [tmppower tmppower_peakidx] = hlp_collapseFrequencies(tmpcrossfpow,g.collapsefun,FREQS,indeximage);
                    
                    if g.EdgeColorMappedToDirectionality
                        tmpangle = sign(index1-index2);
                        if ~tmpangle, tmpangle=1; end
                    else
                        % collapse data across frequency
                        [tmpangle tmpangle_peakidx] = hlp_collapseFrequencies(tmpcrossfang,g.collapsefun,FREQS,indeximage);
                        
                    end
                    
                    if strcmp(lower(g.crossfphaseunit), 'radian'), tmpangle = tmpangle/pi*180; end;
                    %fprintf('%d-%d -> power %1.1f\n', index1, index2, tmppower);b
                    drawconnections( g.vars.coordinates{tmpcond}{ index1 }, g.vars.coordinates{tmpcond}{ index2 }, ...
                        tmppower, tmpangle, g.circfactor(index1, index2), g);
                end;
            end;
        end;
        
        % draw circles
        % ------------
        for index1 = g.plotorder(:)'
            tmptimef = ALLERSP{ index1, tmpcond};
            
            % node size is power
            [tmppow tmppow_peakidx] = hlp_collapseFrequencies(tmptimef,g.collapsefun,FREQS,indeximage);
            
            tmptimef = ALLITC{ index1, tmpcond};
            
            % node color is ITC
            [tmpitc tmpitc_peakidx] = hlp_collapseFrequencies(tmptimef,g.collapsefun,FREQS,indeximage);
            
            %index1, tmpitc, tmppow,
            drawcircle( g.vars.coordinates{tmpcond}{ index1 }, tmppow, tmpitc, g);
            if ~isempty(g.nodelabels)
                tmpcoord = g.vars.coordinates{tmpcond}{index1};
                for i=1:size(tmpcoord,1)
                    text(tmpcoord(i,1), ...
                        tmpcoord(i,2), ...
                        tmpcoord(i,3)+5, ...
                        g.nodelabels{index1},'color','k');
                end
            end
        end;
        
        % delete old objects
        delete( oldobjs );
        
        % draw a bar for time probability
        % -------------------------------
        if ~isempty(g.rt)
            if ~isempty(g.rt{tmpcond})
                ll = line([g.rthistloc(1)-g.rthistloc(3)/2 g.rthistloc(1)+g.rthistloc(3)/2], [g.rthistloc(2) g.rthistloc(2)]);
                set(ll, 'linewidth', 2*g.resmult, 'color', 'k');
                barheight = RTdist(tmpcond, indeximage)*g.rthistloc(4);
                x1 = g.rthistloc(1)-0.65*g.rthistloc(3)/2;
                x2 = g.rthistloc(1)+0.65*g.rthistloc(3)/2;
                y1 = g.rthistloc(2);
                y2 = g.rthistloc(2)-barheight;
                ll = patch([x1 x1 x2 x2], [y1 y2 y2 y1], g.rthistcolor, 'linewidth', 2*g.resmult);
            end;
        end;
        
            
        % render each layer
        layers = fieldnames(g.Layers);
        
        for layeridx = 1 : length(layers)
        
            curlayer = layers{layeridx};
        
            if ~isstruct(g.Layers.(curlayer)), continue; end
            
            if g.Layers.(curlayer).transparency < 1 ...
                    && isempty(findobj(g.figurehandle,'tag',sprintf('%smesh',curlayer))) ...

                children = get(g.figurehandle,'children');

                hax = findobj(g.figurehandle,'tag',['brain' num2str(tmpcond)]);

                campos = get(hax,'CameraPosition');

                %         set(gcf,'visible','off')
                hold on;
                if strcmpi(curlayer,'loni')
                    Handle = hlp_plotAtlas(g.Layers.(curlayer).mesh,hax,g.Layers.(curlayer).coloratlas,g.Layers.(curlayer).color);
                else
                    Handle = hlp_plotmesh(g.Layers.(curlayer).mesh.faces, g.Layers.(curlayer).mesh.vertices,[],false,hax,g.Layers.(curlayer).color);
                end
                
                set(Handle, 'facealpha',1-g.Layers.(curlayer).transparency,'tag',sprintf('%smesh',curlayer));
                
                if ~strcmpi(get(g.figurehandle,'Renderer'),'opengl')
                    set(g.figurehandle,'Renderer' ,'opengl');
                end

                %         set(gcf,'visible',g.visible)

                set(g.figurehandle,'children',children);

                set(hax,'CameraPosition',campos);

            end;
        end
        
%             reducepatch(csfHandle, 1/8); % 1/6 crashed the computer
            % each call to dipplot adds two new lights, so we need to remove them all,
            % then add two lights
            %         delete(findobj(gcf, 'type','light'));
            
            camlight left;
            camlight right;
            
            
            
%         % render LONI probabilistic map on cortex
%         if g.LONITransparency < 1 && isempty(findobj(g.figurehandle,'tag','LONImesh'))
%             
%             children = get(g.figurehandle,'children');
%             
%             hax = findobj(g.figurehandle,'tag',['brain' num2str(tmpcond)]);
%             campos = get(hax,'CameraPosition');
%             
% %             csf = load(g.csf_volume); % located in visualization folder
% %             csfVertices = g.csf.vol.bnd(3).pnt;
%             
%             %         set(gcf,'visible','off')
%             hold on;
%             
% %             g.LONImesh = reducepatch(g.LONImesh,1/30);
%             LONIHandle = hlp_plotAtlas(g.LONImesh,hax);
%             
% %             LONIHandle = hlp_plotmesh(g.LONImesh.faces, g.LONImesh.vertices, [], false);
%             set(LONIHandle, 'facealpha',1-g.LONITransparency,'tag','LONImesh');
% %             set(g.figurehandle,'Renderer' ,'opengl');
% 
% %             reducepatch(LONIHandle, 0.5); % 1/6 crashed the computer
% 
%             % each call to dipplot adds two new lights, so we need to remove them all,
%             % then add two lights
%             %         delete(findobj(gcf, 'type','light'));
%             
% %             camlight left;
% %             camlight right;
%             
%             %         set(gcf,'visible',g.visible)
%             
%             set(g.figurehandle,'children',children);
%             
%             set(hax,'CameraPosition',campos);
%             
%             % since each dipplot has a camzoom, we have to zoom out
%             %camzoom(1/(1.2^(2*(length(clusterNumbers)-1))));
%             
%             %         set(gca, 'CameraViewAngleMode', 'manual');
%             
%         end
            
            
        if ~isempty(g.title) & tmpcond == 1
            if iscell(g.title), g.title = g.title{1}; end
            set(gcf,'Name',['BrainMovie3D ' g.title]);
        end;

%         
%         if ~isempty(g.title) & tmpcond == 1
%             t = textsc(g.title{1},'title');
%             set(t,'VerticalAlignment','top', 'fontsize', 15);
%         end;
        
%         % bring brain to front
%         axes(brainaxis(tmpcond));
%         axis off; set (gcf, 'visible', g.visible);
        
    end;  % LOOP OVER CONDITIONS
    
    % put the time in the right bottom corner
    % --------------------------------------
    if g.showLatency
        if ~isfield(g.vars,'hlatency'),
            g.vars.hlatency = text(0.92, 0.05, sprintf('%d ms', round(TIMES(indeximage))), 'unit', 'normalized','parent',g.vars.hback);
            set(g.vars.hlatency, 'fontsize', 12*g.resmult, 'horizontalalignment', 'right', 'tag', 'tmpmov', 'color', 'w');
        else
            set(g.vars.hlatency,'String', sprintf('%d ms', round(TIMES(indeximage))));
        end  
%         uistack(hlatency,'up',10);
    end
    
    
    % draw the envelope of the signal if necessary
    % ---------------------------------------------
%     axes(g.vars.henv(tmpcond)); cla; axis off; set (gcf, 'visible', g.visible);
    if ~isempty( g.envelope )
        minordinate = mxin(min(min(g.envelope)));
        maxordinate = max(max(max(g.envelope)));
        for tmpcond = 1:nbconditions
            cla(g.vars.henv(tmpcond)); axis on; % set (gcf, 'visible', g.visible);
            plot(g.vars.henv(tmpcond),g.TIMES(:,1), g.envelope(:,:,tmpcond), 'k', 'linewidth', 2*g.resmult); hold on;
            set(g.vars.henv(tmpcond), 'ylim', [minordinate maxordinate]);
            set(g.vars.henv(tmpcond), 'xlim', [g.TIMES(1,1) g.TIMES(end,1)]);
            plot(g.vars.henv(tmpcond),[g.TIMES(getindex(g.TIMES(:,1),TIMES(indeximage)),1) g.TIMES(getindex(g.TIMES(:,1),TIMES(indeximage)),1)], [minordinate maxordinate], 'b', 'linewidth', 2*g.resmult);
            line([0 0],[g.TIMES(1,1) g.TIMES(end,1)],'linestyle',':','color','b','linewidth',1,'parent',g.vars.henv(tmpcond));
            line([g.TIMES(1,1) g.TIMES(end,1)],[0 0],'color','k','parent',g.vars.henv(tmpcond));
            xlabel(g.vars.henv(tmpcond),'Time (ms)', 'fontweight', 'bold', 'fontsize', 12*g.resmult); set(g.vars.henv(tmpcond), 'box', 'off');
            set(g.vars.henv(tmpcond), 'fontsize', 10*g.resmult);
            if tmpcond == 1
                ylabel(g.vars.henv(tmpcond),g.envylabel, 'fontweight', 'bold', 'fontsize', 12*g.resmult);
            end;
            
            % draw vertical lines if needed
            % -----------------------------
            if ~isempty(g.envvert)
                drawvert(g.envvert, tmpcond,  [minordinate maxordinate],g.vars.henv(tmpcond));
            end;
        end;
        %
        %           figure(gcf);
        %           h=text(0.8492,0.1508,num2str(
        
        % put the time on the ERP axis
        % ----------------------------
        %coordx1 = (g.xlimaxes(2)-g.xlimaxes(1))*0.1 + g.xlimaxes(1);
        %coordy1 = (g.ylimaxes(2)-g.ylimaxes(1))*0.87 + g.ylimaxes(1);
    end;
    
    
    % last 3-D settings
    % -----------------
    lighting phong;
    material shiny;
    
    % reset text fontcolor to white
    set(findobj(g.vars.hlgnd,'type','text'),'color',[0.99 0.99 0.99]);    
    
%     if ~g.speedy
%         setfont(gcf, 'color', [0.99 0.99 0.99]); % warning, for some reasons white does not print
%         if isfield(g.vars,'hlgnd')
%             for index = 1:length(g.vars.hlgnd)
%                 axes(g.vars.hlgnd(index)); % bring back legend to front
%                 set(gcf,'visible',g.visible);  % TM - restore visibility
%             end;
%         end
%     end
    
    drawnow
    
    if ~isempty(g.moviename)
        % save the file for a movie
        % -------------------------
        movframes = getframe(gcf);
        mov = addframe(mov,movframes);
    end
    

%     figure(findobj('tag','BrainMovieFigure'));
    

    
    if ~isempty(g.framefolder) 
        fname = sprintf('%simage%4.4d.%s', g.framefolder, indeximage,g.framesout);
        saveas(gcf,fname,g.framesout);
        pause(0.2);
%         
%         if strcmpi(g.framesout, 'tiff')
%             command2 = sprintf('print -dtiff %simage%4.4d.tiff', g.framefolder, indeximage);
%             eval(command2);
%         elseif strcmpi(g.framesout, 'eps')
%             command2 = sprintf('print -depsc -loose %simage%4.4d.eps', g.framefolder, indeximage);
%             eval(command2);
%         elseif 	strcmpi(g.framesout, 'ppm')
%             command2 = sprintf('print -dppm -loose %simage%4.4d.ppm', g.framefolder, indeximage);
%             eval(command2);
%         elseif 	strcmpi(g.framesout, 'jpg')
%             command2 = sprintf('print -djpeg -loose %simage%4.4d.jpg', g.framefolder, indeximage);
%             eval(command2);
%         else % fig format
%             hgsave(sprintf('%simage%4.4d.fig', g.framefolder, indeximage));
%             if strcmp(g.visible, 'on')
%                 drawnow;
%             end;
%         end;
    end
end;

if ~isempty(g.moviename)
    mov = close(mov);
end

if g.verb==2
    close(g.vars.hwaitbar);
end

% EOF












% function to draw circles
% ------------------------
function [tmpsize, tmpcolor, handles] = drawcircle( tmpcoord, tmpersp, tmpitc, g);
% tmpcoord         coordinate of the circle
% tmpersp          erps power -> radius
% tmpitc           itc -> color
% g                preference

% deal with dual dipoles
% ----------------------
handles = [];
if size(tmpcoord,1) > 1
    [tmpsize1, tmpcolor1, handles1] = drawcircle( tmpcoord(1,:), tmpersp, tmpitc, g);
    [tmpsize2, tmpcolor2, handles2] = drawcircle( tmpcoord(2,:), tmpersp, tmpitc, g);
    tmpsize  = [ tmpsize1 tmpsize2];
    tmpcolor = [ tmpcolor1; tmpcolor2];
    handles  = [ handles1 handles2];
    return;
end;

% Node size
tmpsize = scalevalue(tmpersp, g.nodeSizeDataRange, g.nodeSizeLimits);
if isnan(tmpersp),          tmpsize = 0; end;
if strcmpi(g.power, 'off') || tmpersp==eps, tmpsize = g.nodeSizeLimits(1); end;
tmpsize = 0.05*tmpsize*10;
tmpsize = g.diskscale*tmpsize*100;

% Node color
tmpcolor = scalevalue(tmpitc, g.nodeColorDataRange, g.nodeColorLimits);
indexcolor = floor(max(tmpcolor-0.001,0)*length(g.nodeColormap))+1;
tmpcolor   = g.nodeColormap( indexcolor,: );
if isnan(tmpitc) || strcmpi(g.itc, 'off'), tmpcolor = [1 1 1]; end; %g.nodeColorLimits(1);

if tmpersp == 0, dashed = 1;
else             dashed = 0;
end;

if tmpsize > 0
    if length(tmpcoord) > 2
        [xstmp ystmp zs] = sphere(15);
        l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
        normals = reshape([xstmp./l ystmp./l zs./l],[16 16 3]);
        xs = tmpcoord(1) + tmpsize*ystmp*g.vars.dimratio;
        ys = tmpcoord(2) + tmpsize*xstmp;
        zs = tmpcoord(3) + tmpsize*zs;
        colorarray = repmat(reshape(tmpcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
        %if tmpitc ~= 0
        %else colorarray = repmat(reshape([1 1 1],  1,1,3), [size(zs,1) size(zs,2) 1]);
        %end;
        %figure;gca; hold on;
        handles = surf(g.vars.axhandle,xs, ys, zs, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
            'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
        %axis off; axis equal; lighting phong; camlight left; rotate3d
        if strcmpi(g.project3d, 'on')
            if g.colorshadow
                colorarray=colorarray.*repmat(reshape(g.projcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
            else
                colorarray = repmat(reshape(g.projcolor, 1,1,3), [size(zs,1) size(zs,2) 1]);
            end
            surf(g.vars.axhandle,xs, ys, g.factproj(3)*ones(size(zs)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
            surf(g.vars.axhandle,xs, g.factproj(2)*ones(size(ys)), zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
            surf(g.vars.axhandle,g.factproj(1)*ones(size(xs)), ys, zs, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
        end;
    else
        circle( tmpcoord(1), tmpcoord(2), tmpsize, tmpcolor, 'k', 0, 360, dashed, fastif(dashed, 2, 1), 50, g.vars.axhandle);
    end;
end;
return;

% function to draw the lines
% --------------------------
function handles = drawconnections( pos1, pos2, crossfpower, crossfangle, circfact, g);
% pos1, pos2		position of the points
% crossfpower       coherence power for width of the line
% crossfangle       coherence angle for color and speed of the line
% cirfact           curvature of the line
% g                 preference
% arrow should point from pos2 to pos1

% deal with dual dipoles
% ----------------------
handles = [];
if crossfpower == 0, return; end;
if size(pos1,1) > 1
    handles1 = drawconnections( pos1(1,:), pos2, crossfpower, crossfangle, circfact, g);
    handles2 = drawconnections( pos1(2,:), pos2, crossfpower, crossfangle, circfact, g);
    handles = [ handles1 handles2];
    return;
end;
if size(pos2,1) > 1
    handles1 = drawconnections( pos1, pos2(1,:), crossfpower, crossfangle, circfact, g);
    handles2 = drawconnections( pos1, pos2(2,:), crossfpower, crossfangle, circfact, g);
    handles = [ handles1 handles2];
    return;
end;

% if the two circle are too close and do not draw the line
% --------------------------------------------------------
distance = sqrt(sum((pos1-pos2).^2));
if distance < 0.05*(g.ylimaxes(2) - g.ylimaxes(1)), return;
end;

% Edge size
% ---------
tmpthick = scalevalue( crossfpower, g.edgeSizeDataRange, g.edgeSizeLimits );
if strcmpi(lower(g.crossfcoh), 'off') || isnan(crossfpower), tmpthick  = 0; end;
tmpthick = 30*tmpthick; % adjusted for Nick
if crossfpower == 0, tmpthick = 0; end;

% Edge color
% ----------
tmpcolor = scalevalue( crossfangle, g.edgeColorDataRange, g.edgeColorLimits );
if strcmpi(lower(g.crossfphasecolor), 'off'), tmpcolor  = 0; end;
indexcolor = floor(max(tmpcolor-0.001,0)*length(g.edgeColormap))+1;
tmpcolor   = g.edgeColormap( indexcolor,: );
if isnan(crossfangle), tmpcolor = [1 1 1]; end;

if tmpthick > 0
    if g.causality
        %             cylwidth=linspace(g.resmult*tmpthick/300*100,0.1,50);
        cylwidth = [0.1 g.resmult*tmpthick/300*100];    % tapered cylinder pointing from pos2->pos1
    else
        cylwidth=g.resmult*tmpthick/300*100;
    end
    [xc yc zc] = cylinder( cylwidth, 10 );   % create unit-length cylinder
    colorarray = repmat(reshape(tmpcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
    handles = surf(g.vars.axhandle,xc, yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', ...
        'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', 'phong', 'ambientstrength', 0.3,'buttondownfcn','msgbox(''foo!'');');
    [xc yc zc] = adjustcylinder2( handles, [pos1(1) pos1(2) pos1(3)], [pos2(1) pos2(2) pos2(3)] );  % stretch and rotate cylinder to match start-end pnts
    
    
    % compute cylinder normals (have to bias normal closer to sphere
    % to get a specular point
    % TM: this is the part that creates the "spotlight" effect on each
    % arc
    cx = mean(xc,2); cx = [(3*cx(1)+cx(2))/4; (cx(1)+3*cx(2))/4];
    cy = mean(yc,2); cy = [(3*cy(1)+cy(2))/4; (cy(1)+3*cy(2))/4];
    cz = mean(zc,2); cz = [(3*cz(1)+cz(2))/4; (cz(1)+3*cz(2))/4];
    tmpx = xc - repmat(cx, [1 11]);
    tmpy = yc - repmat(cy, [1 11]);
    tmpz = zc - repmat(cz, [1 11]);
    l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
    normals = reshape([tmpx./l tmpy./l tmpz./l],[2 11 3]);
    set( handles, 'vertexnormals', normals);
    
    %         if crossfangle>0,k=2; else k=1; end
    %         cx = mean(xc,2); cx = cx(k)-3; %[(cx(1)+cx(2))/4; (cx(1)+cx(2))/4];
    %         cy = mean(yc,2); cy = cy(k); %[(cy(1)+cy(2))/4; (cy(1)+cy(2))/4];
    %         cz = mean(zc,2); cz = cz(k); %[(cz(1)+cz(2))/4; (cz(1)+cz(2))/4];
    %         tmpx = xc - repmat(cx, [2 11]);
    %         tmpy = yc - repmat(cy, [2 11]);
    %         tmpz = zc - repmat(cz, [2 11]);
    %         l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
    %         normals = reshape([tmpx./l tmpy./l tmpz./l],[2 11 3]);
    %         set( handles, 'vertexnormals', normals);
    
    %figure
    %axis off; axis equal; lighting phong; camlight left; rotate3d
    if strcmpi(g.project3d, 'on')
        if g.colorshadow
            colorarray=colorarray.*repmat(reshape(g.projcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
        else
            colorarray = repmat(reshape(g.projcolor, 1,1,3), [size(zc,1) size(zc,2) 1]);
        end
        surf(g.vars.axhandle,xc, yc, g.factproj(3)*ones(size(zc)), colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
        surf(g.vars.axhandle,xc, g.factproj(2)*ones(size(yc)), zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
        surf(g.vars.axhandle,g.factproj(1)*ones(size(xc)), yc, zc, colorarray, 'tag', 'tmpmov', 'edgecolor', 'none', 'facelighting', 'none');
    end;
    %if round(tmpthick) == 7, asdf; end;
end;
return;

% ***************************************************************************************
%                              Caption and tests
% ***************************************************************************************

% function to draw circles at all power
% -------------------------------------
function scalepower(posx, posy, g);

NBCIRCLE = 3;
coordy = posy;
powerscale = [ ceil( g.nodeSizeDataRange(1) ) 0 floor( g.nodeSizeDataRange(2) ) ];
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');

for i=1:NBCIRCLE
    [tmpsize] = drawcircle( [posx coordy], powerscale(i), 0, g);
    if i == 1, tmpsizeori = tmpsize; end;
    
    if i == NBCIRCLE
        hlatency = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g dB', powerscale(i)));
    else hlatency = text( 1.4*(xlim(2) - xlim(1))+xlim(1), coordy , sprintf('%2.1g', powerscale(i)));
    end;
    set(hlatency, 'fontsize', 10*g.resmult, 'horizontalalignment', 'left', 'fontweight', 'bold');
    coordy = coordy + tmpsize + 0.2*(ylim(2)-ylim(1));
    
    %command2 = sprintf('print -depsc -loose scale%d.eps', i);
    %eval(command2);
    %cla;
end;
set(gca, 'xlim', xlim, 'ylim', ylim-tmpsizeori, 'clipping', 'off', 'fontsize', 10*g.resmult);
return;

% function to draw lines at all coherence
% ---------------------------------------
function scalecoher(posx, posy, thickness,g);
compter = -5;
for i=linspace( posy(1), posy(2), 11)
    superline( [ posx(1) posx(2) ], [ i i ], 'b', thickness*g.resmult, mod(compter/10, 1));
    compter = compter + 1;
end;
%ylabel('Phase-Coh', 'fontweight', 'bold', 'fontsize', 12*g.resmult);
axis off;
set(gca, 'ytick', [], 'yticklabel', [], 'xtick',[],'xticklabel', [], 'box', 'on');
% 	set(gca, 'box', 'on', 'ylim', [0 1], 'ytick', [0 0.5 1], ...
% 			 'yticklabel', strvcat('-180�','0�','180�'), 'xlim', [0 1], 'xtick', [], 'xticklabel', [], 'fontsize', 10*g.resmult);
%hold on; ff = fill([0 0.02 0.02 0], [0 0 1 1], 'w'); set(ff, 'edgecolor', 'w');
%hold on; ff = fill([0 0 1 1], [0 0.02 0.02 0], 'w'); set(ff, 'edgecolor', 'w');
return;

% colorbar special
% ----------------
function cbar( X, Y, colors, orientation, style, g )
% colors = colors to plot
% orientation = 'vert' or 'horiz'
% style = shape of the colorbar, 'circle' = circle, bar otherwise

NSEGMENTS = size(colors,1)-1;
compter = 0;
switch lower(orientation)
    case 'horiz'
        inc = (X(2)-X(1))/NSEGMENTS;
        for i=linspace(X(1),X(2)-inc,NSEGMENTS);
            compter = compter + 1;
            hold on;
            h = fill( [i i i+inc i+inc], [Y(1) Y(2) Y(2) Y(1)], colors(size(colors,1)+1-compter, :));
            set(h, 'edgecolor', 'none');
        end;
    case 'vert'
        inc = (X(2)-X(1))/NSEGMENTS;
        for i=linspace(Y(1),Y(2)-(Y(2)-Y(1))/NSEGMENTS,NSEGMENTS);
            hold on;
            switch style
                case 'circle',
                    mid     = (X(2)-X(1))/2;
                    angle   = acos( compter/NSEGMENTS*2-1);
                    angle1  = acos( (compter+1)/NSEGMENTS*2-1);
                    coordx1 = mid - sin( angle )*mid;
                    coordx2 = mid + sin( angle )*mid;
                    coordx3 = mid + sin( angle1 )*mid;
                    coordx4 = mid - sin( angle1 )*mid;
                    coordx = real([coordx1 coordx2 coordx3 coordx4]);
                otherwise,	coordx = [X(1) X(2) X(2) X(1)];
            end;
            compter = compter + 1;
            h = fill( coordx, [i i i+inc i+inc], colors(compter, :));
            set(h, 'edgecolor', 'none');
        end;
    otherwise
        disp('Orientation has to be ''vert'' or ''horiz''');
end;
set(gca, 'fontsize', 10*g.resmult);
if strcmp(style, 'circle'), axis square; end;
return;

% draw vertical lines
% -------------------
function drawvert(tmpev, tmpcond, coords,h)

if nargin<4
    h = gca;
end

if isstruct(tmpev) || isstruct(tmpev{1})
    
    % cooper envert
    %--------------
    if length(tmpev) > 1,
        verts = tmpev{ tmpcond };
    else   verts = tmpev{1};
    end
    
    for v=verts,
        if isstruct(v), ev = v;
        else,           ev.time = v;  ev.color = 'k'; ev.style = '-';
        end
        
        phandle = plot(h,[ev.time ev.time], coords, ev.style, 'linewidth', 1);
        set(phandle,'color',ev.color);
    end;
else
    % standard envvert
    % ----------------
    for index = 1:length(tmpev)
        if ~iscell(tmpev{index}),
            plot(h,[tmpev{index} tmpev{index}], coords, 'k');
        else
            phandle = plot(h,[tmpev{index}{1} tmpev{index}{1}], coords, 'k');
            if length(tmpev{index}) > 2
                set(phandle,tmpev{index}{2:end});
            end;
        end;
    end;
end;
return;

% check the flux
% --------------
for indeximage = 1:nwin-7
    index1 = 1;
    index2 = 2;
    % determine color = coherence phase
    tmpcrossf = ALLCROSSFANGLE     { index1, index2, 1 };
    tmpvalue  = mean(tmpcrossf( 1:2, indeximage));
    tmpcolor  = colormaphsv( ceil((tmpvalue+180)/360*63) + 1, : );    % index for color
    
    % absolute value to 90 degree determine speed
    speed = 1 - abs(90 - abs(tmpvalue))/90; % speed from 1 to 0
    currentphase(index1, index2) = currentphase(index1, index2) + sign(tmpvalue)*speed/3; % 1 cycle in 5 images at max speed
    
    superline( [ 2 1] , [ 1+indeximage 0.8+indeximage], 5, tmpcolor, mod(currentphase(index1, index2),1));
end;
return;

% scale values for node, connector and color
% ------------------------------------------
function value = scalevalue( value, datarange, limits)

value = (value-datarange(1))/(datarange(2)-datarange(1));
value = max(min(value,1),0);
value = value * (limits(2)-limits(1)) + limits(1);










% circle() - draw a circle in the current Matlab axes
%
% Usage:
%   >> [linehandles fillhandle] = circle(X,Y,radius,colorfill,coloredge,oriangle,...
%                                         endangle,dashed,thickness,segments);
% Inputs:
%   X, Y       - circle center coordinates 
%   radius     - circle radius. Can be a vector of 2 values, one
%                for the x dimension and one the y dimension.
%   colorfill  - circle fill color (default:0=none)
%   coloredge  - circle edge color (default:black; 0:no edge)  
%   oriangle   - starting angle (default:0 degrees)
%   endangle   - ending angle (default:360 degrees)
%   dashed     - 0=no, 1=yes (default: no)
%   thickness  - thickness of boarder (default: line default)
%   segments   - number of line segments (default:50)
%
% Outputs:
%   linehandles - handle to the lines of the circle object
%   fillhandle  - handle to the interior of the circle object

% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function [h, h2] = circle( X, Y, radius, colorfill, coloredge, oriangle, endangle, dashed, thickness, segments, h);

if nargin < 3
	error('Not enough arguments. Type help circle');
	return;
end;
if nargin < 4
	colorfill = 0;
end;
if nargin < 5
	coloredge = 'k';		
end;
if nargin < 6
	oriangle = 0;		
end;
if nargin < 7
	endangle = 360;		
end;
if nargin < 8
	dashed = 0;
end;
if nargin < 9
	thickness = 0;
end;
if nargin < 10
	segments = 50;
end;
if nargin < 11
    h = gca;
end
if any(radius <= 0)
    return;
end;

A = linspace(oriangle/180*pi, endangle/180*pi, segments);

% draw surface
% ------------
if any(colorfill)
	A = linspace(oriangle/180*pi, endangle/180*pi, segments);
	h2 = patch( X + cos(A)*radius(1), Y + sin(A)*radius(end), zeros(1,segments), colorfill, 'parent', h);
	set(h2, 'FaceColor', colorfill); 
	set(h2, 'EdgeColor', 'none'); 
end;		

% draw lines
% ----------
if dashed
	compt=0;	
	for i=1:3:segments-2
		compt = compt+1;
		h(compt) = line( X + cos(A(i:i+1))*radius(1), Y + sin(A(i:i+1))*radius(end));
	end;
else
	h = line( X + cos(A)*radius(1), Y + sin(A)*radius(end), 'parent', h);
end;
set(h, 'Color', coloredge); 

if thickness ~= 0
	set(h, 'Linewidth', thickness); 
end;
return;
