
function axisHandle = vis_dataStreamViewer(varargin)
% A simple datastream viewer
%
% Example:
%
% axisHandle = vis_dataStreamViewer('closeRequestFcn','assignin(''base'',''isclosed'',true); delete(gcbf);');
% isclosed = false;
% while ~isclosed
%     vis_dataStreamViewer('axisHandle',axisHandle);
%     drawnow;
% end

try
    
arg_define([0 3], varargin, ...
    arg({'streamname','MatlabStream'},'mystream','','Matlab stream name'), ...
    arg({'basewsDatasetName'},[],[],'Optional name of EEG data structure in base workspace. Data will be pulled from this instead of streamname'), ...
    arg({'spacing','VerticalSpacing'},30,[0 10^5],'Spacing between lines. Units are in data measurement units (e.g. uV)'), ...
    arg({'forceSpacingUpdate'},false,[],'If true, then spacing is set to the value above. Otherwise, the value of spacing stored in the axis UserData is used'), ...
    arg({'timerange','TimeRange','NumSecondsToDisplay'},10,[0 Inf],'Number of seconds of data to display'), ...
    arg({'forceTimeRangeUpdate'},false,[],'If true, then timerange is set to the value above. Otherwise, the value of timerange stored in the axis UserData is used'), ...
    arg({'channelsToDisplay'},[],[],'Indices of channels to display'), ...
    arg({'axisHandle'},[],[],'Handle to current streamviewer axis. If empty, a new StreamViewer figure is created'), ...
    arg({'closeRequestFcn'},'delete(gcbf)','','Figure close request function. This could delete a timer, for instance'), ...
    arg({'signalType'},'Electrodes',{'ICs','Electrodes'},'Type of signal to display'), ...
    arg('refreshAxis',false,[],'If true, redraw the plot instead of updating xdata,ydata'), ...
    arg('pipeLine','',[],'A pipeline to extract data from instead of onl_peek()'), ...
    arg('draw',false,[],'If true, draw the figure (redraw())') ...
    );

    
if ~isempty(axisHandle) && ishandle(axisHandle)
    ud = get(axisHandle,'UserData');
    
    if ~isempty(ud)
        
        % store new properties
        if forceTimeRangeUpdate
            ud.timerange = timerange;
        end
        
        if forceSpacingUpdate
            ud.spacing = spacing;
        end
        
        set(axisHandle,'UserData',ud);
        
        % ... we always operate on the
        % viewer properties stored in the
        % axis UserData
        timerange   = ud.timerange;
        spacing     = ud.spacing;
        
    end
end

% get a chunk of data
if isempty(basewsDatasetName)
    tmp=onl_peek(streamname,timerange);
else
    tmp = evalin('base',basewsDatasetName);
    if isempty(tmp)
        error('Dataset %s does not exist in base workspace',basewsDatasetName);
    end
end

if isempty(channelsToDisplay)
    channelsToDisplay = 1:tmp.nbchan;
end

% extract the data
data = tmp.data(channelsToDisplay,:)';

[npnts nchs] = size(data);

if isempty(axisHandle) || ~ishandle(axisHandle)
    % construct new figure
    figureHandle = figure('Name',sprintf('Data Stream Viewer: Stream ''%s''',streamname),'CloseRequestFcn',closeRequestFcn, ...
        'KeyPressFcn',@(varargin) updateStreamViewerDisplay(varargin{2}.Key,varargin{1}));
    
    axisHandle = axes('parent',figureHandle);
    
    refreshAxis = true;
    
end

if isempty(axisHandle) || ~ishandle(axisHandle) || refreshAxis
    % draw a new image
    
    h = plot(axisHandle,linspace(tmp.xmin,tmp.xmax,npnts),bsxfun(@plus,data,(0:nchs-1)*spacing-mean(data,1)));
    set(axisHandle,'UserData',struct('spacing',spacing,'timerange',timerange,'lineHandles',h));
    xlabel(axisHandle,'Time (sec)','FontSize',12);
    ylabel(axisHandle,'Voltage (\mu V^2)','FontSize',12);
    set(axisHandle,'tag','StreamViewer');
    
else
    % just update the xdata and ydata
    
    % update the data
    newdata = bsxfun(@plus,data,(0:nchs-1)*spacing-mean(data,1));
    
    for k=1:length(ud.lineHandles)
        set(ud.lineHandles(k),'Ydata',newdata(:,k));
        set(ud.lineHandles(k),'Xdata',linspace(tmp.xmin,tmp.xmax,npnts));
    end
    
end

% update the axis limit and tickmarks
axis(axisHandle,[tmp.xmin tmp.xmax -spacing nchs*spacing+spacing]);
set(axisHandle,'YTick',(0:nchs-1)*spacing,'YTickLabel',{tmp.chanlocs(channelsToDisplay(end:-1:1)).labels});

if draw
    drawnow;
end


catch err
    env_handleerror(err);
%     keyboard;
end



function updateStreamViewerDisplay(Key,figureHandle)
% update the display properties of the streamviewer based on keypress

axHandle = findobj(get(figureHandle,'children'),'tag','StreamViewer');

ud=get(axHandle,'UserData');

switch lower(Key)
    
    case 'uparrow'
        % increase spacing
        ud.spacing = ud.spacing*1.1;
        
    case 'downarrow'
        % decrease spacing
        ud.spacing = ud.spacing*0.9;
        
    case 'rightarrow'
        % increase timerange
        ud.timerange = ud.timerange*1.1;
        
    case 'leftarrow'
        % decrease timerange
        ud.timerange = ud.timerange*0.9;
end

set(axHandle,'UserData',ud);


