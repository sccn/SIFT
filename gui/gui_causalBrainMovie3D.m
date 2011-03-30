function varargout = gui_causalBrainMovie3D(varargin)
%
% GUI_CAUSALBRAINMOVIE3D M-file for gui_causalBrainMovie3D.fig
%      GUI_CAUSALBRAINMOVIE3D, by itself, creates a new GUI_CAUSALBRAINMOVIE3D or raises the existing
%      singleton*.
%
%      H = GUI_CAUSALBRAINMOVIE3D returns the handle to a new GUI_CAUSALBRAINMOVIE3D or the handle to
%      the existing singleton*.
%
%      GUI_CAUSALBRAINMOVIE3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_CAUSALBRAINMOVIE3D.M with the given input arguments.
%
%      GUI_CAUSALBRAINMOVIE3D('Property','Value',...) creates a new GUI_CAUSALBRAINMOVIE3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_causalBrainMovie3D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_causalBrainMovie3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_causalBrainMovie3D

% Last Modified by GUIDE v2.5 17-Nov-2010 14:48:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_causalBrainMovie3D_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_causalBrainMovie3D_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_causalBrainMovie3D is made visible.
function gui_causalBrainMovie3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_causalBrainMovie3D (see VARARGIN)

% Choose default command line output for gui_causalBrainMovie3D
handles.output = hObject;
set(hObject,'Name','BrainMovie3D Control Panel');
% set default termination behavior
handles.ExitButtonClicked = 'Cancel';

% extract some data from command-line input
if isempty(varargin)
    error('You must pass ALLEEG and Conn to gui_CausalBrainMovie3D');
end

% Extract input parameters/data and store
handles.ud.ALLEEG  = varargin{1};  
handles.ud.Conn    = varargin{2};  
varargin([1 2]) = [];

% set the default properties of the Real Time Visualization Panel
if length(handles.ud.Conn(1).erWinCenterTimes)==1
    set(handles.pnlRealTimeVis,'enable','false');
else
    erWinCenterTimes = handles.ud.Conn(1).erWinCenterTimes;
    set(handles.txtMinTime,'String',erWinCenterTimes(1));
    set(handles.txtMaxTime,'String',erWinCenterTimes(end));
    set(handles.txtCurTime,'String',erWinCenterTimes(1));
    set(handles.slideCurTime,'value',1);
    set(handles.slideCurTime,'min',1,'max',length(erWinCenterTimes));
    set(handles.slideCurTime,'sliderstep',[1/length(erWinCenterTimes) 4/length(erWinCenterTimes)]);
end

% set default EEGLAB background and text colors
%-----------------------------------------------
try, icadefs;
catch,
	GUIBACKCOLOR        =  [.8 .8 .8];     
	GUIPOPBUTTONCOLOR   = [.8 .8 .8];    
	GUITEXTCOLOR        = [0 0 0];
end;

allhandlers = hObject;

hh = findobj(allhandlers,'style', 'text');
%set(hh, 'BackgroundColor', get(hObject, 'color'), 'horizontalalignment', 'left');
set(hh, 'Backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hObject, 'color',GUIBACKCOLOR );
% set(hh, 'horizontalalignment', g.horizontalalignment);

hh = findobj(allhandlers, 'style', 'edit');
set(hh, 'BackgroundColor', [1 1 1]); %, 'horizontalalignment', 'right');

hh =findobj(allhandlers, 'parent', hObject, 'style', 'pushbutton');
if ~strcmpi(computer, 'MAC') && ~strcmpi(computer, 'MACI') % this puts the wrong background on macs
    set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
    set(hh, 'foregroundcolor', GUITEXTCOLOR);
end;
hh =findobj(allhandlers, 'parent', hObject, 'style', 'popupmenu');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'checkbox');
set(hh, 'backgroundcolor', GUIBACKCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'listbox');
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hh, 'foregroundcolor', GUITEXTCOLOR);
hh =findobj(allhandlers, 'parent', hObject, 'style', 'radio');
set(hh, 'foregroundcolor', GUITEXTCOLOR);
set(hh, 'backgroundcolor', GUIPOPBUTTONCOLOR);
set(hObject, 'visible', 'on');

set(handles.pnlPropertyGrid,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlPropertyGrid,'foregroundcolor', GUITEXTCOLOR);
set(handles.pnlRealTimeVis,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlRealTimeVis,'foregroundcolor', GUITEXTCOLOR);
%-----------------------------------------------

drawnow

% render the PropertyGrid in the correct panel
handles.PropertyGridHandle = arg_guipanel( ...
                 handles.pnlPropertyGrid, ...
                'Function',@vis_causalBrainMovie3D, ...
                'params',{'ALLEEG',handles.ud.ALLEEG, 'Conn',handles.ud.Conn, varargin{:}});

handles.BrainMovieFigure = [];
            
% Update handles structure
guidata(hObject, handles);

% Wait for user to click OK, Cancel or close figure
uiwait(handles.gui_BrainMovie3D);


% UIWAIT makes gui_causalBrainMovie3D wait for user response (see UIRESUME)
% uiwait(handles.gui_BrainMovie3D);





% --- Outputs from this function are returned to the command line.
function varargout = gui_causalBrainMovie3D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles)
    % user closed the figure
    varargout = {[] hObject};
elseif strcmpi(handles.ExitButtonClicked,'OK')
    % user clicked OK
    % get PropertySpecification
    varargout = {handles.PropertyGridHandle handles.output};
else
    % user clicked cancel
    varargout = {[] handles.output};
end

try, close(hObject);
catch; end


% --- Executes on button press in cmdCancel.
function cmdCancel_Callback(hObject, eventdata, handles)

handles.ExitButtonClicked = 'Cancel';
guidata(hObject,handles);
uiresume(handles.gui_BrainMovie3D)




% --- Executes on button press in cmdOK.
function cmdOK_Callback(hObject, eventdata, handles)

handles.ExitButtonClicked ='OK';
guidata(hObject,handles);
uiresume(handles.gui_BrainMovie3D);




% --- Executes on slider movement.
function slideCurTime_Callback(hObject, eventdata, handles)

% get the PropertyGrid object
% ud = get(handles.gui_BrainMovie3D,'UserData');
% if isempty(ud)
%     error('Could not obtain handle to PropertyGrid! Did you pass it to the UserData of the figure?');
% end

% get the property specification
pg = handles.PropertyGridHandle;
ps = pg.GetPropertySpecification;
cfg = arg_tovals(ps,false);

if isempty(cfg)
    error('could not retrieve property specification from propertygrid');
end

% get the index of the frame the user wants to render
range = get(hObject,'max')-get(hObject,'min');
curframe = round(get(hObject,'Value'));

% update the curTime textbox
set(handles.txtCurTime,'string',handles.ud.Conn(1).erWinCenterTimes(curframe));

minTime = str2num(get(handles.txtMinTime,'String'));
maxTime = str2num(get(handles.txtMaxTime,'String'));

figh = handles.BrainMovieFigure; %findobj('tag','BrainMovieFigure');
try figure(figh)
catch
    figh = [];
end
curpos = get(findobj(figh,'tag','brain'),'CameraPosition');
% execute the brainmovie function to render this frame
vis_causalBrainMovie3D('ALLEEG',handles.ud.ALLEEG,'Conn',handles.ud.Conn,cfg,'MovieTimeRange',[minTime maxTime],'BrainMovieOptions',{cfg.BMopts,'visible','off','FigureHandle',figh,'LatenciesToRender',[],'FramesToRender',curframe,'outputFormat',{'framefolder','','moviename',''}});
figh = gcf;
set(findobj(figh,'tag','brain'),'CameraPosition',curpos);
set(gcf,'visible','on')
handles.BrainMovieFigure = gcf;
guidata(hObject,handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes during object creation, after setting all properties.
function slideCurTime_CreateFcn(hObject, eventdata, handles)


% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object deletion, before destroying properties.
function cmdOK_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to cmdOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdHelp.
function cmdHelp_Callback(hObject, eventdata, handles)
% hObject    handle to cmdHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warndlg2('Coming soon!');

