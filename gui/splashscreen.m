function varargout = splashscreen(varargin)
%
% SPLASHSCREEN M-file for splashscreen.fig
%      SPLASHSCREEN, by itself, creates a new SPLASHSCREEN or raises the existing
%      singleton*.
%
%      H = SPLASHSCREEN returns the handle to a new SPLASHSCREEN or the handle to
%      the existing singleton*.
%
%      SPLASHSCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPLASHSCREEN.M with the given input arguments.
%
%      SPLASHSCREEN('Property','Value',...) creates a new SPLASHSCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before splashscreen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to splashscreen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help splashscreen

% Last Modified by GUIDE v2.5 20-Nov-2010 04:31:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @splashscreen_OpeningFcn, ...
                   'gui_OutputFcn',  @splashscreen_OutputFcn, ...
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


% --- Executes just before splashscreen is made visible.
function splashscreen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to splashscreen (see VARARGIN)

% Choose default command line output for splashscreen
handles.output = hObject;

axes(handles.axLogo);
path = which('splashlogo.jpg');
im = imread(path);
image(im);
axis image
axis off

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes splashscreen wait for user response (see UIRESUME)
uiwait(handles.gui_splash);


% --- Outputs from this function are returned to the command line.
function varargout = splashscreen_OutputFcn(hObject, eventdata, handles) 
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

