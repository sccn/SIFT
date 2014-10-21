function varargout = gui_splashscreen(varargin)
%
% GUI_SPLASHSCREEN M-file for gui_splashscreen.fig
%      GUI_SPLASHSCREEN, by itself, creates a new GUI_SPLASHSCREEN or raises the existing
%      singleton*.
%
%      H = GUI_SPLASHSCREEN returns the handle to a new GUI_SPLASHSCREEN or the handle to
%      the existing singleton*.
%
%      GUI_SPLASHSCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPLASHSCREEN.M with the given input arguments.
%
%      GUI_SPLASHSCREEN('Property','Value',...) creates a new GUI_SPLASHSCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_splashscreen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_splashscreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_splashscreen

% Last Modified by GUIDE v2.5 28-Jun-2012 15:16:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_splashscreen_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_splashscreen_OutputFcn, ...
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


% --- Executes just before gui_splashscreen is made visible.
function gui_splashscreen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_splashscreen (see VARARGIN)

% Choose default command line output for gui_splashscreen
handles.output = hObject;

axes(handles.axLogo);
path = which('splashlogo.jpg');
im = imread(path);
image(im);
axis image
axis off

siftRoot = hlp_getSiftRoot;
splashText = importdata(fullfile(siftRoot,'Readme.txt'),'\n');
vers = importdata(fullfile(hlp_getSiftRoot,'resources','version.txt'));
splashText = strrep(splashText,'[$Version$]',vers);
set(handles.txtSplashText,'String',splashText);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_splashscreen wait for user response (see UIRESUME)
uiwait(handles.gui_splash);


% --- Outputs from this function are returned to the command line.
function varargout = gui_splashscreen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
try,
varargout{1} = handles.output;
close(hObject);
catch
varargout{1} = [];
end


function gui_splash_CreateFcn(hObject, eventdata, handles)


function cmdOK_Callback(hObject, eventdata, handles)

uiresume(handles.gui_splash);


% --- Executes during object creation, after setting all properties.
function txtSplashText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSplashText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function txtSplashText_Callback(hObject, eventdata, handles)
% hObject    handle to txtSplashText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtSplashText as text
%        str2double(get(hObject,'String')) returns contents of txtSplashText as a double


% --- Executes during object creation, after setting all properties.
function txtSplashScreen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtSplashText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
