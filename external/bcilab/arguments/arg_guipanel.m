function h = arg_guipanel(varargin)
% Create a uipanel that displays an argument property inspector for a Function.
% Handle = arg_guipanel(Options ...)
% Handle = arg_guipanel(Parent, Options ...)
%
% The handle supports the method .GetPropertySpecification(), by means of which the edited argument specification can be retrieved.
% This result can be turned into a valid Function argument using arg_tovals(). Additional Parameters may be passed to the Function,
% in order to override some of its defaults.
%
% In:
%   Parent : optional parent widget
%
%   Options : name-value pairs; possible names are:
%              'Function' : the function for which to display arguments
%
%              'Parameters' : cell array of parameters to the function
%
%              'Position' : position of the panel within the parent widget
%
% Out:
%   Handle : handle to the panel; supports .GetPropertySpecification() to obain the edited specification
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-10-24

% separate the parent, if specified
if isscalar(varargin{1}) && ishandle(varargin{1})
    parent = varargin{1};
    varargin = varargin(2:end);
else
    parent = figure('MenuBar','none');
end

% get the options
opts = hlp_varargin2struct(varargin, {'Function','func'},mandatory, {'Parameters','params'},{}, {'Position','pos'},[0 0 1 1]);

% obtain the argument specification for the function
spec = arg_report('rich', opts.Function, opts.Parameters);
% ... and turn it into an array of PropertyGridField's
properties = PropertyGridField.GenerateFrom(spec);

% instantiate the grid
args = hlp_struct2varargin(opts,'suppress',{'Function','Parameters'});
h = PropertyGrid(parent,args{:},'Properties',properties);
