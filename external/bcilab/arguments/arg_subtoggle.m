function res = arg_subtoggle(varargin)
% Specify an argument of a function which is a struct of sub-arguments that can be disabled.
% Spec = arg_subtoggle(Names,Default,Source,Help,Options...)
%
% Accessible to the function as a struct, and visible in the GUI as a an expandable sub-list of arguments (with a checkbox to toggle).
% The special field 'arg_selection' (true/false) indicates whether the argument is enabled or not. The value assigned to the argument determines whether it
% is turned on or off, as determined by the mapper option.
%
% In:
%   Names : The name(s) of the argument. At least one must be specified, and if multiple are specified, they must be passed in a cell array.
%           * The first name specified is the argument's "code" name, as it should appear in the function's code (= the name under which arg_define()
%             returns it to the function).
%           * The second name, if specified, is the "Human-readable" name, which is exposed in the GUIs (otherwise the code name is displayed).
%           * Further specified names are alternative names for the argument (e.g., for backwards compatibility with older function syntaxes/parameter names).
%
%   Defaults : A cell array of arguments to override defaults for the Source; all syntax accepted by the (selected) Source is allowed here,
%              whereas in the case of positional arguments, the leading arg_norep() arguments of the source are implicitly skipped.
%              Note: Whether the argument is turned on or off is determined via the 'mapper' option.
%                    By default, [] and 'off' are mapped to off, whereas {}, non-empty cell arrays and structs are mapped to on.
%
%   Source : A source of argument specifications, usually a function handle (referring to a function which defines arguments via arg_define()).
%
%            For convenience, a cell array with a list of argument declarations, formatted like the Specification part of an arg_define() clause
%            can be given, instead. In this case, the effect is the same as specifying @some_function, for a function implemented as:
%
%               function some_function(varargin)
%               arg_define(Format,varargin,Source{:});
%
%   Help : The help text for this argument (displayed inside GUIs), optional. (default: []).
%          (Developers: Please do *not* omit this, as it is the key bridge between ease of use and advanced functionality.)
%
%          The first sentence should be the executive summary (max. 60 chars), any further sentences are a detailed explanation (examples, units, considerations).
%          The end of the first sentence is indicated by a '. ' followed by a capital letter (beginning of the next sentence). If ambiguous,
%          the help can also be specified as a cell array of 2 cells.
%
%   Options... : Optional name-value pairs to denote additional properties:
%                 'cat' :  The human-readable category of this argument, helpful to present a list of many parameters in a categorized list, and to separate
%                          "Core Parameters" from "Miscellaneous" arguments. Developers: When choosing names, every bit of consistency with other function in the
%                          toolbox helps the uses find their way (default: []).
%
%                 'fmt' : Optional format specification for the Source (if it is a cell array) (default: []). See arg_define() for a detailed explanation.
%
%                 'mapper' : A function that maps the argument list (e.g., Defaults) to a value in the domain of selectors, and a potentially updated argument list.
%                            The mapper is applied to the argument list prior to any parsing (i.e. it faces the raw argument list) to determine the current selection, and its
%                            its second output (the potentially updated argument list) is forwarded to the Source that was selected, for further parsing.
%
%                            The default mapper maps [] and 'off' to off, whereas 'on', empty or non-empty cell arrays and structs are mapped to on.
%
%                 'merge': Whether a value (cell array of arguments) assigned to this argument should completely replace all arguments of the default,
%                          or whether it should instead the two cell arrays should be concatenated ('merged'), so that defaults are only selectively overridden.
%                          Note that for concatenation to make sense, the cell array of Defaults cannot be some subset of all allowed positional arguments, 
%                          but must instead either be the full set of positional arguments (and possibly some NVPs) or be specified as NVPs in the first place.
%
%
% Out:
%   Spec : A cell array, that, when called as spec{1}(reptype,spec{2}{:}), yields a specification of the argument, for use by arg_define. 
%          Technical note: Upon assignment with a value (via the assigner field), the
%          'children' field of the specifier struct is populated according to how the selected (by the mapper) Source parses the value into arguments.
%          The additional struct field 'arg_selection 'is introduced at this point.
%
% Notes:
%   for MATLAB versions older than 2008a, type and shape checking is not necessarily enforced.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

% we return a function that an be invoked to yield a specification (its output is cached for efficiency)
% packed in a cell array together with the remaining arguments
res = {@invoke_argsubtoggle_cached,varargin};


function spec = invoke_argsubtoggle_cached(varargin)
spec = hlp_microcache('arg',@invoke_argsubtoggle,varargin{:});


% the function that does the actual work of building the argument specifier
function spec = invoke_argsubtoggle(reptype,names,defaults,source,help,varargin)

% start with a base specification
spec = arg_specifier('head',@arg_subtoggle, 'fmt',[], 'type','logical', 'shape','scalar', 'mapper',@map_argsubtoggle);

% override properties
if exist('names','var')
    spec.names = names; end
if exist('help','var')
    spec.help = help; end
for k=1:2:length(varargin)
    if isfield(spec,varargin{k})
        spec.(varargin{k}) = varargin{k+1};
    else
        error('BCILAB:arg:no_new_fields','It is not allowed to introduce fields into a specifier that are not declared in arg_specifier.');
    end
end

% do checking
if ~iscell(spec.names)
    spec.names = {spec.names}; end
if isempty(spec.names) || ~iscellstr(spec.names)
    error('The argument must have a name or cell array of names.'); end
if ~exist('source','var') || isempty(source)
    error('BCILAB:args:no_options','The Source argument for arg_subtoggle() may not be omitted.'); end %#ok<*NODEF>
if nargin(spec.mapper) == 1
    spec.mapper = @(x,y,z) spec.mapper(x); end
% parse the help
if ~isempty(spec.help)
    try
        spec.help = parse_help(spec.help);
    catch
        e = lasterror; %#ok<LERR>
        disp(['Problem with the help text for argument ' spec.names{1} ': ' e.message]);
        spec.help = {};
    end
elseif spec.reportable && spec.displayable
    disp(['Please specify a description for argument ' spec.names{1} ', or specify it via arg_nogui() instead.']);
end

% uniformize Source syntax
if iscell(source)
    % args is a cell array instead of a function: we effectively turn this into a regular arg_define-using function (taking & parsing values)
    source = @(varargin) arg_define(spec.fmt,varargin,source{:});
else
    % args is a function: was a custom format specified?
    if isa(spec.fmt,'function_handle')
        source = @(varargin) source(spec.fmt(varargin));
    elseif ~isempty(spec.fmt)
        error('The only allowed form in which the Format of a Source that is a function may be overridden is as a pre-parser (given as a function handle)');
    end
end
spec = rmfield(spec,'fmt');

% find out what boolean flag and value set the default configuration maps to
% this is relevant for the merging option: in this case, we need to pull up the currect default and merge it with the passed value
[default_sel,default_val] = spec.mapper(defaults);

% set up the regular assigner
spec.assigner = @(spec,value) assign_argsubtoggle(spec,value,reptype,source,default_sel,default_val);

% assign the default
if strcmp(reptype,'rich')
    spec = assign_argsubtoggle(spec,defaults,'build',source,NaN,{});
else
	spec = assign_argsubtoggle(spec,defaults,'lean',source,NaN,{});
end


function spec = assign_argsubtoggle(spec,value,reptype,source,default_sel,default_val)
persistent arg_sel arg_desel;
if isempty(arg_sel) || isempty(arg_sel)
    arg_sel = arg_nogui('arg_selection',true); arg_sel = arg_sel{1}([],arg_sel{2}{:});
    arg_desel = arg_nogui('arg_selection',false); arg_desel = arg_desel{1}([],arg_desel{2}{:});
end

% retrieve the values for the realized switch option...
[selected,value] = spec.mapper(value);
% build the complementary alternative, if requested
if strcmp(reptype,'build')
    if selected
        spec.alternatives{1} = arg_desel;
    else
        spec.alternatives{2} = [arg_report('rich',source,{}) arg_sel];
    end
    reptype = 'rich';
end
% obtain the children
if ~selected
    spec.children = arg_desel;
elseif spec.merge && (default_sel==true)
    spec.children = [arg_report(reptype,source,[default_val value]) arg_sel];
else
    spec.children = [arg_report(reptype,source,value) arg_sel];
end
spec.alternatives{selected+1} = spec.children;
% and set the cell's value
spec.value = selected;




function [selected,args] = map_argsubtoggle(args)
if isequal(args,'on')
    selected = true;
    args = {};
elseif isequal(args,'off') || isequal(args,[])
    selected = false;
    args = [];
elseif length(args) == 1 && isstruct(args) && isfield(args,'arg_selection')
    selected = args.arg_selection;
elseif length(args) == 1 && iscell(args) && isstruct(args{1}) && isfield(args{1},'arg_selection')
    selected = args{1}.arg_selection;
else
    selected = true;
end