function varargout = headModelViewer(varargin)
% HEADMODELVIEWER MATLAB code for headModelViewer.fig
%      HEADMODELVIEWER, by itself, creates a new HEADMODELVIEWER or raises the existing
%      singleton*.
%
%      H = HEADMODELVIEWER returns the handle to a new HEADMODELVIEWER or the handle to
%      the existing singleton*.
%
%      HEADMODELVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HEADMODELVIEWER.M with the given input arguments.
%
%      HEADMODELVIEWER('Property','Value',...) creates a new HEADMODELVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before headModelViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to headModelViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help headModelViewer

% Last Modified by GUIDE v2.5 20-Jun-2012 18:00:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @headModelViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @headModelViewer_OutputFcn, ...
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


% --- Executes just before headModelViewer is made visible.
function headModelViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to headModelViewer (see VARARGIN)

% Choose default command line output for headModelViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes headModelViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = headModelViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hLabels,'Visible','off');
else 
    set(obj.hLabels,'Visible','on');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hSensors,'Visible','off');
else 
    set(obj.hSensors,'Visible','on');
end

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hCortex,'Visible','off');
else 
    set(obj.hCortex,'Visible','on');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hSkull,'Visible','off');
else 
    set(obj.hSkull,'Visible','on');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hScalp,'Visible','off');
else 
    set(obj.hScalp,'Visible','on');
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hFiducials,'Visible','off');
else 
    set(obj.hFiducials,'Visible','on');
end


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
obj = get(handles.figure1,'userData');
if get(hObject,'Value')
    set(obj.hCortex,'FaceVertexCData',[ 0 1 0],'LineStyle','-');
else 
    set(obj.hCortex,'FaceVertexCData',obj.streamObj.atlas.color,'LineStyle','none');
end
