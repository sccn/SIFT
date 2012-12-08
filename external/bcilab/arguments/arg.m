function res = arg(varargin)
% A (rich) specification of a function argument, for use in arg_define() clauses.
% Spec = arg(Names,Default,Range,Help,Options...)
%
% In:
%   Names : The name(s) of the argument. At least one must be specified, and if multiple are specified, they must be passed in a cell array.
%           * The first name specified is the argument's "code" name, as it should appear in the function's code (= the name under which arg_define()
%             returns it to the function).
%           * The second name, if specified, is the "Human-readable" name, which is exposed in the GUIs (otherwise the code name is displayed).
%           * Further specified names are alternative names for the argument (e.g., for backwards compatibility with older function syntaxes/parameter names).
%
%   Default : Optionally the default value of the argument; can be any data structure (default: []). Special values: 
%              * unassigned: the argument is not listed in the function's workspace nor GUI unless explicitly assigned
%              * mandatory: if no value is specified to this argument by the time the function is called, an error is raised
%             Note: If neither Default nor Range are specified, consider specifying the argument's type via the Options... list.
%
%   Range : Optionally a range of admissible values (default: []).
%           * If empty, no range is enforced.
%           * If a cell array, each cell is considered one of the allowed values.
%           * If a 2-element numeric vector, the two values are considered the numeric range of the data (inclusive).
%           * If an anonymous function of one of the forms: @(x)a<x<b, @(x)a<x<=b, @(x)a<=x<b, or @(x)a<=x<=b,
%             then a and b determine the range (either inclusive or not, according to the relational operators).
%           * Note: if neither Default nor Range are specified, consider specifying the argument's type via the Options... list.
%
%   Help : The help text for this argument (displayed inside GUIs), optional. (default: []).
%          (Developers: Please do *not* omit this, as it is the key bridge between ease of use and advanced functionality.)
%
%          The first sentence should be the executive summary (max. 60 chars), any further sentences are a detailed explanation (examples, units, considerations).
%          The end of the first sentence is indicated by a '. ' followed by a capital letter (beginning of the next sentence). If ambiguous,
%          the help can also be specified as a cell array of 2 cells.
%
%   Options... : Optional name-value pairs to denote additional properties:
%                 'cat' : The human-readable category of this argument, helpful to present a list of many parameters in a categorized list, and to separate
%                         "Core Parameters" from "Miscellaneous" arguments. Developers: When choosing names, every bit of consistency with other function in the
%                         toolbox helps the uses find their way (default: []).
%
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
%          The (internal) structure is as follows:
%          * Generally, this is a cell array (here: one element) of cells formatted as: {Names,Assigner-Function,Default-Value}.
%          * Names is a cell array of admissible names for this argument.
%          * Assigner-Function is a function that returns the rich specifier with value assigned, when called as Assigner-Function(Value).
%          * reptype is either 'rich' or 'lean', where in lean mode, the aternatives field remains empty.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

% map range functions to strings (removing their variable references, if any) for proper cacheability
if length(varargin)>=3 && isa(varargin{3},'function_handle')
    varargin{3} = char(varargin{3}); end

% we return a function that an be invoked to yield a specification (its output is cached for efficiency)
% packed in a cell array together with the remaining arguments
res = {@invoke_arg_cached,varargin};


function spec = invoke_arg_cached(reptype,varargin) %#ok<INUSL>
spec = hlp_microcache('arg',@invoke_arg,varargin{:});


% the function that does the actual work of building the argument specifier
function spec = invoke_arg(names,default,range,help,varargin)

% start with a base specification
spec = arg_specifier('head',@arg);

% override properties
if exist('names','var')
    spec.names = names; end
if exist('default','var')
    spec.value = default; end
if exist('range','var')
    spec.range = range; end
if exist('help','var')
    spec.help = help; end
for k=1:2:length(varargin)
    if isfield(spec,varargin{k})
        spec.(varargin{k}) = varargin{k+1}; 
    else
        error('BCILAB:arg:no_new_fields','It is not allowed to introduce fields into a specifier that are not declared in arg_specifier.');
    end
end

% do fixups & checking
if ~iscell(spec.names)
    spec.names = {spec.names}; end
if isempty(spec.names) || ~iscellstr(spec.names)
    error('The argument must have a name or cell array of names.'); end

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

% do type inference
[spec.type,spec.shape,spec.range] = infer_type(spec.type,spec.shape,spec.range,spec.value);

% do minimal fixups on the fully typed data
if isequal(spec.value,[])    
    if strcmp(spec.type,'cellstr')
        spec.value = {}; 
    else
        try
            pt = PropertyType(spec.type,spec.shape,spec.range);
            spec.value = pt.ConvertFromMatLab(spec.value);
        catch
        end
    end
end

if strcmp(spec.type,'logical') && iscell(spec.range) && isscalar(spec.value) && islogical(spec.value)
    if spec.value
        spec.value = spec.range;
    else
        spec.value = {};
    end
end

% infer the type & range of the argument, based on provided info (note: somewhat messy)
function [type,shape,range] = infer_type(type,shape,range,value)
try
    if isempty(type)
        % try to auto-discover the type (or leave empty, if impossible)
        if ~isempty(value)
            type = PropertyType.AutoDiscoverType(value);
        elseif ~isempty(range)
            if isnumeric(range)
                type = PropertyType.AutoDiscoverType(range);
            elseif iscell(range)
                types = cellfun(@PropertyType.AutoDiscoverType,range,'UniformOutput',false);
                if length(unique(types)) == 1
                    type = types{1}; end
            end
        end
    end
    if isempty(shape)
        % try to auto-discover the shape
        if ~isempty(value)
            shape = PropertyType.AutoDiscoverShape(value);
        elseif ~isempty(range)
            if isnumeric(range)
                shape = 'scalar';
            elseif iscell(range)
                shapes = cellfun(@PropertyType.AutoDiscoverShape,range,'UniformOutput',false);
                if length(unique(shapes)) == 1
                    shape = shapes{1}; end
            end
        end
    end
catch
end

% if in doubt, fall back to denserealdouble and/or matrix
if isempty(type)
    type = 'denserealdouble'; end
if isempty(shape)
    shape = 'matrix'; end

% remap function-style ranges into numeric ranges...
if isa(range,'function_handle') || ischar(range)    
    try
        if ischar(range)
            str = range;
        else
            conts = functions(range);
            str = conts.function;
        end
        if ~strncmp(str,'@(x)',4)
            error; end %#ok<*LTARG>
        str = strsplit(str(5:end),'<');
        if length(str) ~= 3
            error; end
        [a,x,b] = str{:};
        if isvarname(a)
            a = conts.workspace{1}.(a);
        else
            a = eval(a);
        end
        if strcmp(x,'=x')
            lequ = true;
        elseif strcmp(x,'x')
            lequ = false;
        else
            error;
        end
        if b(1) == '='
            b = b(2:end);
            requ = true;
        else
            requ = false;
        end
        if isvarname(b)
            b = conts.workspace{1}.(b);
        else
            b = eval(b);
        end
        if ~lequ
            a = a+eps(a); end
        if ~requ
            b = b-eps(b); end
        if a>=b
            error; end   % alternative, but might crash GUIs: range = [NaN NaN]; end
        range = [a b];
    catch
    end
end


function strings = strsplit(string, splitter)
ix = strfind(string, splitter);
strings = cell(1,numel(ix)+1);
ix = [0 ix numel(string)+1];
for k = 2 : numel(ix)
    strings{k-1} = string(ix(k-1)+1:ix(k)-1); end