function varargout = gui_metaControlPanel(varargin)
%
% GUI_METACONTROLPANEL M-file for gui_metaControlPanel.fig
%      GUI_METACONTROLPANEL, by itself, creates a new GUI_METACONTROLPANEL or raises the existing
%      singleton*.
%
%      H = GUI_METACONTROLPANEL returns the handle to a new GUI_METACONTROLPANEL or the handle to
%      the existing singleton*.
%
%      GUI_METACONTROLPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_METACONTROLPANEL.M with the given input arguments.
%
%      GUI_METACONTROLPANEL('Property','Value',...) creates a new GUI_METACONTROLPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_metaControlPanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_metaControlPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_metaControlPanel

% Last Modified by GUIDE v2.5 22-Aug-2012 21:31:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_metaControlPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_metaControlPanel_OutputFcn, ...
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


% --- Executes just before gui_metaControlPanel is made visible.
function gui_metaControlPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_metaControlPanel (see VARARGIN)

% Choose default command line output for gui_metaControlPanel
handles.output = hObject;
handles.gcf    = hObject;
set(hObject,'Name','Meta Control Panel');
% set default termination behavior
handles.ExitButtonClicked = 'Cancel';

% extract some data from command-line input
if isempty(varargin)
    error('First input to gui_metaControlPanel must be a valid bcilab dataset');
end
if length(varargin)<2
    error('Second input to gui_metaControlPanel must be an opts structure');
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

set(handles.pnlPropertyGridFltPip,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlPropertyGridFltPip,'foregroundcolor', GUITEXTCOLOR);

set(handles.pnlPropertyGridSiftPip,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlPropertyGridSiftPip,'foregroundcolor', GUITEXTCOLOR);

set(handles.pnlMisc,'backgroundcolor', GUIBACKCOLOR);
set(handles.pnlMisc,'foregroundcolor', GUITEXTCOLOR);

%-----------------------------------------------

% render the PropertyGrids in the correct panels
handles.SiftPipPropertyGridHandle = arg_guipanel( ...
                handles.pnlPropertyGridSiftPip, ...
                'Function',@onl_siftpipeline, ...
                'Parameters',{'EEG',varargin{1}});
    
handles.FltPipPropertyGridHandle = arg_guipanel( ...
                 handles.pnlPropertyGridFltPip, ...
                'Function',@flt_pipeline, ...
                'Parameters',{'Signal',varargin{1}, 'DataCleaning', {'RetainPhases' true 'HaveBursts' false 'DataSetting' 'drycap'}});

% store opts structure
handles.opts = varargin{2};

% initialize some control contents from opts input
try set(handles.txtWindowLength,'String',num2str(handles.opts.winLenSec));
catch; end
try set(handles.chkDispBenchmarking,'Value',handles.opts.doBenchmark);
catch; end
try set(handles.chkDispBrainMovie,'Value',handles.opts.doBrainMovie);
catch; end
try set(handles.chkDispSpectrum,'Value',handles.opts.dispSpectrum);
catch; end

drawnow


% Update handles structure
guidata(hObject, handles);


% Wait for user to click OK, Cancel or close figure
% uiwait(handles.gui_BrainMovie3D);


% UIWAIT makes gui_metaControlPanel wait for user response (see UIRESUME)
% uiwait(handles.gui_BrainMovie3D);





% --- Outputs from this function are returned to the command line.
function varargout = gui_metaControlPanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout = {hObject};


% --- Executes on button press in cmdCancel.
function cmdCancel_Callback(hObject, eventdata, handles)

if isempty(handles)
    % user closed the figure
    varargout = {[] hObject};
% elseif strcmpi(handles.ExitButtonClicked,'OK')
%     % user clicked OK
%     % get PropertySpecification
%     varargout = {handles.PropertyGridHandle handles.output};
else
    % user clicked cancel
    varargout = {[] handles.output};
end

handles.opts.exitPipeline = true;
assignin('base','opts',handles.opts);
guidata(hObject,handles);

try, close(handles.figurehandle);
catch; end



% --- Executes on button press in cmdUpdate.
function cmdUpdate_Callback(hObject, eventdata, handles)

% handles.ExitButtonClicked ='OK';


% get the property specifications
handles.opts.fltPipCfg  = arg_tovals(handles.FltPipPropertyGridHandle.GetPropertySpecification,true);
handles.opts.siftPipCfg = arg_tovals(handles.SiftPipPropertyGridHandle.GetPropertySpecification,false);

% save pipeline configurations in base workspace and set new data flags
assignin('base','opts',handles.opts);
assignin('base','newPipeline',true);
set(handles.gcf,'UserData','initialized');
% evalin('base','opts.fltpipeline  = evalin(''caller'',''fltcfg'')');
% evalin('base','opts.siftpipeline = evalin(''caller'',''siftcfg'')');

handles.opts.winLenSec = str2num(get(handles.txtWindowLength,'String'));

guidata(hObject,handles);
% guidata(hObject,handles);
% uiresume(handles.gui_BrainMovie3D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% % --- Executes during object creation, after setting all properties.
% function slideCurTime_CreateFcn(hObject, eventdata, handles)
% 
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


% --- Executes during object deletion, before destroying properties.
function cmdUpdate_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to cmdUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
% % --- Executes on button press in cmdHelp.
% function cmdHelp_Callback(hObject, eventdata, handles)
% % hObject    handle to cmdHelp (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% warndlg2('Coming soon!');


% --- Executes on button press in cmdPause.
function cmdPause_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% pause the brainmovie
if strcmpi(get(hObject,'String'),'pause');
    handles.opts.holdPipeline = false;
    set(hObject,'String','Unpause');
else
    handles.opts.holdPipeline = true;
    set(hObject,'String','Pause');
end

assignin('base','opts',handles.opts);
guidata(hObject,handles);

% --- Executes on button press in chkDispBenchmarking.
function chkDispBenchmarking_Callback(hObject, eventdata, handles)
% hObject    handle to chkDispBenchmarking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkDispBenchmarking

handles.opts.doBenchmark = get(hObject,'Value');
assignin('base','opts',handles.opts);
guidata(hObject,handles);

% --- Executes during object deletion, before destroying properties.
function chkDispBenchmarking_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to chkDispBenchmarking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdVisStream.
function cmdVisStream_Callback(hObject, eventdata, handles)
% hObject    handle to cmdVisStream (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% open stream visualizer
vis_filtered;


% --- Executes when pnlMisc is resized.
function pnlMisc_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to pnlMisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chkDispSpectrum.
function chkDispSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to chkDispSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkDispSpectrum

handles.opts.dispSpectrum = get(hObject,'Value');
assignin('base','opts',handles.opts);
guidata(hObject,handles);


% --- Executes on button press in chkDispBrainMovie.
function chkDispBrainMovie_Callback(hObject, eventdata, handles)
% hObject    handle to chkDispBrainMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkDispBrainMovie


handles.opts.doBrainMovie = get(hObject,'Value');
assignin('base','opts',handles.opts);
guidata(hObject,handles);



function txtWindowLength_Callback(hObject, eventdata, handles)
% hObject    handle to txtWindowLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtWindowLength as text
%        str2double(get(hObject,'String')) returns contents of txtWindowLength as a double


% --- Executes during object creation, after setting all properties.
function txtWindowLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtWindowLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
