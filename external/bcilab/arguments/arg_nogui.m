function res = arg_nogui(varargin)
% Like arg(), but not displayed by GUIs.
% Spec = arg_nogui(Names,Default,Range,Options...)
%
% In:
%   Names : The name(s) of the argument. At least one must be specified, and if multiple are specified, they must be passed in a cell array.
%           * The first name specified is the argument's "code" name, as it should appear in the function's code (= the name under which arg_define() 
%             returns it to the function).
%           * The second name, if specified, is the "Human-readable" name, which is exposed in the GUIs (otherwise the code name is displayed).
%           * Further specified names are alternative names for the argument (e.g., for backwards compatibility with older function syntaxes/parameter names).
%
%   Default : Optionally the default value of the argument; can be any data structure (default: []).
%
%   Range : Optionally a range of admissible values (default: []).
%           * If empty, no range is enforced.
%           * If a cell array, each cell is considered one of the allowed values.
%           * If a 2-element numeric vector, the two values are considered the numeric range of the data (inclusive).
%           * If an anonymous function of one of the forms: @(x)a<x<b, @(x)a<x<=b, @(x)a<=x<b, or @(x)a<=x<=b, 
%             then a and b determine the range (either inclusive or not, according to the relational operators).
%
%   Help : The help text for this argument, optional. (default: []).
%
%   Options... : Optional name-value pairs to denote additional properties:
%                 'type' : Specify the primitive type of the parameter (default: [], indicating that it is auto-discovered from the Default and/or Range).
%                          The primitive type is one of the following strings: 
%                             'logical', 'char', 'int8','uint8', 'int16','uint16','int32','uint32','int64','uint64', ...
%                             'denserealsingle','denserealdouble','densecomplexsingle','densecomplexdouble', ...
%                             'sparserealsingle','sparserealdouble','sparsecomplexsingle','sparsecomplexdouble','cellstr','object'.
%                          If auto-discovery was requested, but fails for some reason, the default type is set to 'denserealdouble'.
%
%                 'shape' : Specify the array shape of the parameter (default: [], indicating that it is auto-discovered from the Default and/or Range).
%                           The array shape is one of the following strings: 
%                              'scalar','row','column','matrix','empty'.
%                           If auto-discovery was requested, but fails for some reason, the default shape is set to 'matrix'.
%
% Out:
%   Spec : A cell array, that, when called as spec{1}(reptype,spec{2}{:}), yields a specification of the argument, for use by arg_define.
%          * Generally, this is a cell array (here: one element) of cells formatted as: {Names,Assigner-Function,Default-Value}.
%          * Names is a cell array of admissible names for this argument. 
%          * Assigner-Function is a function that returns the rich specifier with value assigned, when called as Assigner-Function(Value). 
%          * reptype is either 'rich' or 'lean', where in lean mode, the aternatives field remains empty.
%
% Notes:
%   For MATLAB versions older than 2008a, type and shape checking, as well as auto-discovery, are not necessarily executed.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

% we return a function that, when invoked, yields a non-displayable argument
varargin = [varargin cell(1,4-length(varargin))];
res = arg(varargin{:}, 'displayable',false);
