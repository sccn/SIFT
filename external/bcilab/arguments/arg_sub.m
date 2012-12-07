function res = arg_sub(varargin)
% Specify an argument of a function which is a structure of sub-arguments.
% Spec = arg_sub(Names,Defaults,Source,Help,Options...)
%
% Delivered to the function as a struct, and visible in the GUI as a an expandable sub-list of arguments.
% A function may have an argument which itself consists of several arguments. For example, a function may be passing the contents 
% of this struct as arguments to another function, or may just collect several arguments into sub-fields of a single struct. 
% Differs from the default arg() function by allowing, instead of the Range, either a Source function which exposes a list of arguments 
% (itself using arg_define), or a cell array with argument specifications, identical in format to the Specification part of an arg_define() clause.
%
% In:
%   Names : The name(s) of the argument. At least one must be specified, and if multiple are specified, they must be passed in a cell array.
%           * The first name specified is the argument's "code" name, as it should appear in the function's code (= the name under which arg_define() 
%             returns it to the function).
%           * The second name, if specified, is the "Human-readable" name, which is exposed in the GUIs (otherwise the code name is displayed).
%           * Further specified names are alternative names for the argument (e.g., for backwards compatibility with older function syntaxes/parameter names).
%
%   Defaults : A cell array of arguments to override defaults for the Source; all syntax accepted by the Source is allowed here,
%              whereas in the case of positional arguments, the leading arg_norep() arguments of the source are implicitly skipped.
%              If empty, the defaults of the Source are unaffected.
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
%                 'merge': Whether a value (cell array of arguments) assigned to this argument should completely replace all arguments of the default,
%                          or whether instead the two cell arrays should be concatenated ('merged'), so that defaults are only selectively overridden.
%                          Note that for concatenation to make sense, the cell array of Defaults cannot be some subset of all allowed positional arguments, 
%                          but must instead either be the full set of positional arguments (and possibly some NVPs) or be specified as NVPs in the first place.
%
% Out:
%   Spec : A cell array, that, when called as spec{1}(reptype,spec{2}{:}), yields a specification of the argument, for use by arg_define. 
%          Technical note: Upon assignment with a value (via the assigner field), the 'children' field of the specifier struct is populated according 
%          to how the Source parses the value into arguments.
%
% Notes:
%   for MATLAB versions older than 2008a, type and shape checking is not necessarily enforced.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

% we return a function that an be invoked to yield a specification (its output is cached for efficiency)
% packed in a cell array together with the remaining arguments
res = {@invoke_argsub_cached,varargin};


function spec = invoke_argsub_cached(varargin)
spec = hlp_microcache('arg',@invoke_argsub,varargin{:});

% the function that does the actual work of building the argument specifier
function spec = invoke_argsub(reptype,names,defaults,source,help,varargin)

% start with a base specification
spec = arg_specifier('head',@arg_sub,'fmt',[], 'value','', 'type','char', 'shape','row');

suppressNames = {};

% override properties
if exist('names','var')
    spec.names = names; end
if exist('help','var')
    spec.help = help; end
for k=1:2:length(varargin)
    if isfield(spec,varargin{k})
        spec.(varargin{k}) = varargin{k+1}; 
    elseif strcmpi(varargin{k},'suppress')
        suppressNames = varargin{k+1};
    else
        error('BCILAB:arg:no_new_fields','It is not allowed to introduce fields into a specifier that are not declared in arg_specifier.');
    end
end

% do checking
if ~iscell(spec.names)
    spec.names = {spec.names}; end
if isempty(spec.names) || ~iscellstr(spec.names)
    error('The argument must have a name or cell array of names.'); end
if ~exist('source','var') || isequal(source,[])
    error('BCILAB:args:no_source','The Source argument for arg_sub() may not be omitted.'); end %#ok<*NODEF>
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

if ~isempty(source)
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
end
spec = rmfield(spec,'fmt');

% assignment to this object does not touch the value, but instead creates a new children structure
spec.assigner = @(spec,value) assign_argsub(spec,value,reptype,source,defaults,suppressNames);

% assign the default
spec = assign_argsub(spec,defaults,reptype,source,[],suppressNames);



% function to do the value assignment
function spec = assign_argsub(spec,value,reptype,source,default,suppressNames)
if ~isempty(source)
    if spec.merge
        spec.children = arg_report(reptype,source,[default,value]);
    else
        spec.children = arg_report(reptype,source,value);
    end
end

% toggle the displayable option for children which should be suppressed
if ~isempty(suppressNames)
    % identify which children we want to suppress display
    hidden = find(cellfun(@any,cellfun(@(x,y) ismember(x,suppressNames),{spec.children.names},'UniformOutput',false)));
    % set display flag to false
    for k=hidden(:)'
        spec.children(k).displayable = false;
    end
end