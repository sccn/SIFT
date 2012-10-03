% function handles = pop_vis_TimeFreqGrid(ALLEEG,Conn,varargin)
% INPUTS:
%   ALLEEG  -   vector of 1 or 2 EEG datasets. If 2 sets, this function
%               plots difference (and Conn must also be 1x2)
%   Conn    -   connectivity structure with field matrices [nch nch nfreq
%   ntime]
%   
%
% Copyright (C) 2008-2010  Tim Mullen
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
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
%


% TODO: USE INPUTDLG3.m TO CREATE GUI
%       *OR* Use christian's ARG specification to make it work...
function [handles cfg] = pop_vis_TimeFreqGrid(ALLEEG,Conn,varargin)
 

%Spec = arg(Names,Default,Range,Help,Options...)
handles = [];

if length(Conn.erWinCenterTimes) == 1
    % remove some options from the GUI
    subset = {-1,'times','windows'};
else
    subset = {};
end

% cfg = arg_guidialog(@vis_TimeFreqGrid_test,'params',{'ALLEEG',ALLEEG,'Conn',Conn,varargin{:}},'subset',subset);


% Example: rendering buttons
% buttons = {1 {'Style','pushbutton','string','SET TEXTBOX VALUE','callback',@(varargin)(set(findobj('tag','prcthresh'),'string','50')),'tooltipstring',''}, ...
%            3 {'Style','pushbutton','string','LAUCH GUI','callback', @(varargin)pop_vis_TimeFreqGrid(ALLEEG,Conn),'tooltipstring',''}};
% cfg = arg_guidialog(@vis_TimeFreqGrid_test,'params',{'ALLEEG',ALLEEG,'Conn',Conn,varargin{:}},'subset',subset,'buttons',buttons);    %'subset',{'prcthresh',{'Style','commandbutton','callback','msgbox(''foo'')'},'absthresh'}


% Example: using a panel
h = arg_guipanel('Function',@vis_TimeFreqGrid,'params',{'ALLEEG',ALLEEG,'Conn',Conn,varargin{:}});
uiwait(h.Parent);
% determine whether the user selected OK or Cancel

ps = h.GetPropertySpecification;
cfg = arg_tovals(ps,false);

if isempty(cfg)
    return;
end

% execute the function
handles = vis_TimeFreqGrid(ALLEEG,Conn,cfg);






% cfg = arg_define(varargin, ...
%     arg_nogui({'ALLEEG'},mandatory),...
%     arg_nogui({'Conn'},mandatory),...
%     arg({'prcthresh','PercentileThreshold'},100,[],'Percentile threshold. Can be scalar or of form [thresh dim]. If latter, prctile is applied separately for each element of dimension=dim across all remaining dimensions'), ...
%     arg({'absthresh','AbsoluteThreshold'},100,[],'Exact (hard) threshold to apply'), ...
%     arg({'connmethods','ConnectivityMethod'},names{1},names,'Connectivity Measure to visualize'));
    



% 
% cfg = arg_define(varargin, ...
%     ... % core arguments for the paradigm framework (passed by the framework)
%     arg_norep('op',[],{'preprocess','learn','predict'},'Operation to execute on the data. Preprocess the raw data, learn a predictive model, or predict outputs given a model.'), ...
%     arg_norep('data',[],[],'Data to be processed by the paradigm.'), ...
%     arg_norep('model',[],[],'Model according to which to predict.'), ...
%     ... % signal processing arguments (sourced from flt_pipeline)
%     arg_sub({'flt','SignalProcessing'},flt_defaults,@flt_pipeline,'Signal processing stages. These can be enabled, disabled and configured for the given paradigm. The result of this stage flows into the feature extraction stage','cat','Signal Processing'), ...
%     ... % arguments for the feature-extraction plugins (passed by the user paradigms)
%     arg_sub({'fex','FeatureExtraction'},{},fex_declaration,'Parameters for the feature-extraction stage.','cat','Feature Extraction'), ...
%     ... % feature-extraction plugin definitions (passed by the user paradigms)
%     arg_sub({'plugs','PluginFunctions'},fex_defaults,{ ...
%         arg({'adapt','FeatureAdaptor'},@default_feature_adapt,[],'The adaption function of the feature extraction. Function_handle, receives the preprocessed data, an options struct (with feature-extraction), and returns a model (which may just re-represent options).'),...
%         arg({'extract','FeatureExtractor'},@default_feature_extract,[],'The feature extraction function. Function_handle, receives the preprocessed data and the model from the featureadapt stage, and returns a NxF array of N feature vectors (for F features).'), ...
%         arg({'vote','FeatureAdaptorNeedsVoting'},false,[],'Feature-adaption function requires voting. Only relevant if the data contains three or more classes.') ...
%         },'The feature-extraction functions','cat','Feature Extraction'), ...
%     ... % machine learning arguments (sourced from ml_train)
%     arg_sub({'ml','MachineLearning'},ml_defaults,@ml_train,'Machine learning stage of the paradigm. Operates on the feature vectors that are produced by the feature-extraction stage.','cat','Machine Learning'), ...
%     ... % configuration dialog layout
%     arg({'arg_dialogsel','ConfigLayout'},dialog_default,[],'Parameters displayed in the config dialog. Cell array of parameter names to display (dot-notation allowed); blanks are translated into empty rows in the dialog. Referring to a structure argument lists all parameters of that struture, except if it is a switchable structure - in this case, a pulldown menu with switch options is displayed.','type','cellstr','shape','row'));
% 






% 
% 
% if nargin<3
%     popup =1;
% end
% 
% 
% if popup
%     % this is a test
% end
%     
% 
% allconnmethods = fieldnames(Conn(1));
% allconnmethods = setdiff(allconnmethods,{'winCenterTimes','erWinCenterTimes','freqs'});
% connstr = '';
% for i=1:length(allconnmethods)
%     connstr = [connstr allconnmethods{i} '|'];
% end
% connstr(end) = [];
% 
% prompt =   {'Percentile Threshold', ...
%             'Absolute Threshold', ...
%             'Color Limits', ...
%             'Window Centers', ...
%             'Frequencies', ...
%             'Connectivity Measure', ...
%             connstr};
%             
%             
% style = {'edit','edit','edit','edit','edit','text','popupmenu' };      
%             
% default = {100 '' 100 '' ['1:' num2str(freqs(end))] 0 1};
% 
% [result, userdat, strhalt, resstruct] = inputdlg3('prompt',prompt,'style',style,'default',default);
%  
% 
% 
% g = finputcheck(varargin, ...
%     { 'prcthresh'           'real'          [0 100]             []; ...             % percentile threshold (scalar or [thresh dim]. If latter, prctile is applied separately for each element of dimension=dim across all remaining dimensions)
%       'absthresh'           'real'          []                  []; ...             % exact (hard) threshold
%       'clim'                'real'          []                  100;                % [min max] color (or y) limits. If scalar, limits are set to [0 clim-percentile] if all(Conn>0) or [-maxprc maxprc] if any(Conn<0) where maxprc is prctile(abs(Conn),clim-percentile)
%       'evlstyle'            'cell'          {}                  {'k' ':' 2}; ...    % {linecolor linestyle linewidth} of event marker
%       'windows'             'real'          []                  []; ...             % (approximate) centers (seconds) of time windows to plot (using plot)
%       'freqs'               'real'          []                  []; ...             % vector of frequencies (Hz) to examine
%       'connmethods'         'cell'          []                  {}; ...             % which connectivity methods to plot
%       'stats'               ''              []                  []; ...             % Stats structure obtained by cat_calcStats.  If it contains stats.(method).ci, this is the [lo hi] or [-ci ci] confidence intervals to plot (for specific-window plots)
%       'threshside'          'string'        {'single','both','lessthan'}   'single'            % one-sided or two-sided thresholding for stats (if 'both' then stats should have upper and lower thresholds)
%       'pcontour'            'boolean'       [0 1]               0; ...              % plot contour for significance or not
%       'contourcolor'        ''              []                  'k'; ...
%       'baseline'            'real'          []                  []; ...             % ER baseline (seconds)
%       'fighandles'          'real'          []                  []; ...             % vector of figure handles to superimpose new graph onto (new figures and grid will *not* be created). Old grid will be used and new subplots overlaid
%       'interpolate'         'boolean'       [0 1]               0; ...              % whether or not to use interpolation
%       'subplotargs'         'cell'          {}                  {}; ...             % list of <'name',value> arguments to pass to plotconn when a subplot is clicked (these will overwrite defaults)
%       'xord'                'real'          []                  []; ...             % vector of x-ordinates (e.g., window centers) -- length must equal approriate dim of Conn.XXX
%       'yord'                'real'          []                  []; ...             % vector of y-ordinates (e.g., freqs) -- length must equal approriate dim of Conn.XXX
%       'msubset'             'string'        {'tril','triu',...
%                                             'diag','nodiag',''}  ''; ...            % determines which subset of the full matrix to keep/plot: lower/upper triangle ('tril'/'triu'), diagonals ('diag'), everything except diagonal ('nodiag'), everything (''). Data that is not plotted is set to nan.
%       'channels'            'real'          []                  []; ...             % specific list of channel indices to keep ([vector], a subset of [1:nbchan])
%       'plotorder'           'real'          []                  []; ...             % specify index order (subset of 1:nbchan) in which to arrange columns/rows. Useful for grouping channels. 
%       'namestr'             'string'        []                  '';...              % additional stuff to add to figure name
%       'linecolor'           ''              []                  'k'; ...            % linecolor for single-window plots
%       'plotci'              'boolean'       [0 1]               1; ...              % (not for time-frequency images) whether or not to plot confidence intervals (if available)
%       'sigthreshmethod'     'string'        {'p','h','thresh','hfdr',''}  ''; ...   % which field of g.stats.(method) to use for significance thresholding (if '', assumes a matrix in g.stats.(method))
%       'topoplot'            'real'          [0 3]              0; ...
%       'channames'           'cell'          {}                  {}; ...             % list of channel names (overrides default usage of chanlocs names)
%       'foilines'            'real'          []                  []; ...             % vector of frequencies (Hz) at which to draw lines                  
%       'clustmaps'           'cell'          {}                  {}; ...             % cell matrix of mean cluster maps to torepoplot
%   });



