function m = hlp_serialize(v, n)
% Obtain an executable string equivalent for a data structure.
%
% This code is adapted from Joger Hansegord's FileExchange submission;
% it is extended with support for serializing functions and lambda functions properly (and adds a fix for empty struct serialization,
% as well as optimizations for vector serialization).
%
%   matcode = SERIALIZE(x) generates matlab code of x
%   matcode = SERIALIZE(x, n) generates matlab code of x retaining n digits
%   of precision
%   SERIALIZE() enters self test mode
%
%   SERIALIZE should be able to create matlab code of the following data types:
%   - matrices, vectors, and scalars of any class and dimension
%   - strings
%   - structs, arrays of structs with up to six dimensions
%   - cell arrays
%   - matlab objects with a copy constructor implemented (Not Java)
%   - empty values of any class
%   - any combinations hereof
%
%   The value of x can be obtained by
%     eval(matcode)
%
%   Examples
%     x = [1 2 3; 3 4 5];
%     serialize(x)
%
%     x = uint8(rand(10)*5);
%     matcode = serialize(x)
%
%     x = {rand(3,3,3,4), 'a string value', {1 2 3; 1 3 3 }}
%     matcode = serialize(x, 30)
%
%   See also mat2str, num2str, int2str, sprintf, class, eval
%
%                       adapted, Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-28

% AUTHOR    : J�ger Hanseg�rd
% $DATE     : 29-Jun-2006 17:37:49 $
% $Revision : 1.00 $
% DEVELOPED : 7.2.0.232 (R2006a)
% FILENAME  : serialize.m

if nargin == 1
    n = 15; end

m = serializevalue(v, n);

function val = serializevalue(v, n)
if isnumeric(v) || islogical(v)
    val = serializematrix(v, n);
elseif ischar(v)
    val = serializestring(v, n);
elseif isstruct(v)
    val = serializestruct(v, n);
elseif iscell(v)
    val = serializecell(v, n);
elseif isobject(v)
    val = serializeobject(v, n);
elseif isa(v,'function_handle')
    val = serializehandle(v, n);
else
    error('Unhandled type %s', class(v));
end

function val = serializestring(v,n)
val = ['sprintf(''' v ''')'];
doConvertToUint8 = false;
try
    dummy = eval(val);
catch
    doConvertToUint8 = true;
end
if doConvertToUint8 || ~isequal(eval(val), v)
    val = ['char(' serializevalue(uint8(v), n) ')'];
end

function val = serializematrix(v, n)
if ndims(v) < 3
    if isa(v, 'double')
        if size(v,1) == 1 && length(v) > 3 && isequal(v,v(1):v(2)-v(1):v(end))
            % special case: colon sequence
            if v(2)-v(1) == 1
                val = ['[' num2str(v(1)) ':' num2str(v(end)) ']'];
            else
                val = ['[' num2str(v(1)) ':' num2str(v(2)-v(1)) ':' num2str(v(end)) ']'];
            end
        elseif size(v,2) == 1 && length(v) > 3 && isequal(v',v(1):v(2)-v(1):v(end))
            % special case: colon sequence
            if v(2)-v(1) == 1
                val = ['[' num2str(v(1)) ':' num2str(v(end)) ']'''];
            else
                val = ['[' num2str(v(1)) ':' num2str(v(2)-v(1)) ':' num2str(v(end)) ']'''];
            end        
        else            
            val = mat2str(v, n);
        end
    else
        val = mat2str(v, n, 'class');
    end
else
    if isa(v, 'double')
        val = mat2str(v(:), n);
    else
        val = mat2str(v(:), n, 'class');
    end
    val = sprintf('reshape(%s, %s)', val, mat2str(size(v)));
end

function val = serializecell(v, n)
if isempty(v)
    val = '{}';
    return
end
cellSep = ', ';
if isvector(v) && size(v,1) > 1
    cellSep = '; ';
end
vstr = cellfun(@(val) [serializevalue(val, n) cellSep], v, 'UniformOutput', false);
vstr{end} = vstr{end}(1:end-2);
val = [ '{' vstr{:} '}'];
if ~isvector(v)
    val = ['reshape('  val sprintf(', %s)', mat2str(size(v)))];
end

function val = serializestruct(v, n)
fieldNames   = fieldnames(v);
fieldValues  = struct2cell(v);
if ndims(fieldValues) > 6
    error('Structures with more than six dimensions are not supported');
end
val = 'struct(';
for fieldNo = 1:numel(fieldNames)
    val = [val serializevalue( fieldNames{fieldNo}, n) ', '];
    val = [val serializevalue( permute(fieldValues(fieldNo, :,:,:,:,:,:), [2:ndims(fieldValues) 1]) , n) ];
    val = [val ', '];
end
% fix: empty struct
if ~isempty(fieldNames)
    val = val(1:end-2); end
val = [val ')'];
    
if ~isvector(v)
    val = sprintf('reshape(%s, %s)', val, mat2str(size(v)));
end

function val = serializeobject(v, n)
val = sprintf('%s(%s)', class(v), serializevalue(struct(v), n));

function val = serializehandle(v, n)
tmp = char(v);
if tmp(1) == '@'
    try        
        % read out the function's workspace
        contents = functions(v);
        if ~isempty(strfind(contents.function,'eval('))
            warning('hlp_serialize:eval_handle',['The anonymous function being serialized may show altered behavior if it evaluates variables\n' ...
                'of the workspace in which it was originally declared.']);
        end
        val = sprintf('hlp_makefunction(%s,%s)',serializestring(contents.function,n),serializestruct(contents.workspace{1},n));
    catch
        % fall back to the default serialization without workspace
        val = sprintf('eval(%s)',serializestring(func2str(v),n));
    end
else
    try        
        % read out the function's workspace
        contents = functions(v);
        if ~strcmp(contents.type,'simple')
            warning('hlp_serialize:eval_nonsimple','The function being serialized can probably not be deserialized/evaluated, as it is nested.'); end
    catch,end
    % serialize a regular function handle... (note: can be problematic if pointing to private or nested functions)
    val = sprintf('str2func(''%s'')', func2str(v));
end
