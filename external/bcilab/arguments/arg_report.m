function res = arg_report(type,func,args)
% Report argument-related information of a certain Type from the given Function.
% Result = arg_report(Type,Function,Arguments)
%
% In:
%   Type : Type of information to report, can be one of the following:
%          'rich' : Report a rich declaration of the function's arguments as a struct array, with fields as in arg_specifier.
%          'lean' : Report a rich declaration of the function's arguments as a struct array, with fields as in arg_specifier, excluding the alternatives field.
%          'vals' : Report the values of the function's arguments as a struct, possibly with sub-structs.
%
%          'handle': Report function handles to scoped functions within the Function (i.e., subfunctions). The named of those functions are listed
%                    as a cell string array in place of Arguments.
%
%   Function : a function handle to a function which defines some arguments (via arg_define)
%
%   Arguments : cell array of parameters to be passed to the function; depending on the function's implementation,
%               this can affect the current value assignment (or structure) of the parameters being returned
%               If this is not a cell, it is automatically wrapped inside one (note: to specify the first positional argument as [] 
%               to the function, always pass it as {[]}; this is only relevant if the first argument's default is non-[]).
%
% Out:
%   Result : the reported data.
%
% Notes: 
%   The Function must use arg_define() to define its arguments.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

global bcilab;

if ~exist('args','var')
    args = {}; end
if isequal(args,[])
    args = {}; end
if ~iscell(args)
    args = {args}; end
    
try
    % call one of the arg_report* functions: these essentially communicate to arg_define (via their presence in the call stack) the extraction mode
    res = feval(['arg_report_' lower(type)],func,args);
catch
    e = lasterror; %#ok<LERR>
    if strcmp(e.identifier,'BCILAB:arg:report_args')
        % we received the arguments
        res = bcilab.temp.report_args;
    else
        % this was an actual error; rethrow
        rethrow(e);
    end
end

function res = arg_report_rich(func,args) %#ok<DEFNU>
if exist('exp_eval','file')
    res = exp_eval(func(args{:}),1);
else
    res = func(args{:});
end

function res = arg_report_lean(func,args) %#ok<DEFNU>
if exist('exp_eval','file')
    res = exp_eval(func(args{:}),1);
else
    res = func(args{:});
end

function res = arg_report_vals(func,args) %#ok<DEFNU>
if exist('exp_eval','file')
    res = exp_eval(func(args{:}),1);
else
    res = func(args{:});
end

function res = arg_report_handle(func,args) %#ok<DEFNU>
if exist('exp_eval','file')
    res = exp_eval(func(args{:}),1);
else
    res = func(args{:});
end
