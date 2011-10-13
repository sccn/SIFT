function res = arg_define(varargin)
% Declare function arguments with optional defaults and built-in GUI support.
% Struct = arg_define(Values, Specification...)
% Struct = arg_define(Format, Values, Specification...)
%
% This is essentially a in improved replacement for the parameter declaration line of a function.
% Assigns Values (a cell array of values, typically the "varargin" of the calling function, henceforth named the "Function") to fields in the output Struct,
% with parsing implemented according to a Specification of argument names and their order (optionally with a custom argument Format description).
%
% By default, values can be a list of a fixed number of positional arguments (i.e., the typical MATLAB calling format), optionally followed by a list of name-value
% pairs (NVPs, e.g., as the format accepted by figure()), in which, furthermore, instead of any given NVP, a struct may be passed as well (thus, one may pass a mix
% of 'name',value,struct,'name',value,'name',value, ... parameters). Alternatively, by default the entire list of positional arguments can instead be be specified 
% as a list of NVPs/structs. Only names that are allowed by the Specification may be used, if positional syntax is allowed by the Format (which is the default).
%
% The special feature over hlp_varargin2struct()-like functionality is that arguments defined via arg_define can be reported to the framework (if triggered by
% arg_report()). The resulting specification can be rendered in a GUI or be processed otherwise.
%
% In:
%   Format : Optional format description (default: [0 Inf]):
%            * If this is a function handle, the function is used to transform the Values prior to any other processing into a new Values cell array. The function
%              may specify a new (numeric) Format as its second output argument (if not specified, this is 0).
%            * If this is a number (say, k), it indicates that the first k arguments are specified in a positional manner, and the following arguments
%              are specified as list of name-value pairs and/or structs.
%            * If this is a vector of two numbers [0 k], it indicates that the first k arguments MAY be specified in a positional manner (the following arguments
%              must be be specified as NVPs/structs) OR alternatively, all arguments can be specified as NVPs / structs. Only names that are listed in
%              the specification may be used as names (in NVPs and structs) in this case.
%
%   Values : A cell array of values passed to the function (usually the calling function's "varargin"). Interpreted according to the Format and the Specification.
%
%   Specification... : The specification of the calling function's arguments; this is a sequence of arg(), arg_norep(), arg_nogui(), arg_sub(),
%                      arg_subswitch(), arg_subtoggle() specifiers. The special keyword mandatory can be used in the declaration of default values, which
%                      declares that this argument must be assigned some value via Values (otherwise, an error is raised before the arg is passed to the Function).
%
% Out:
%   Struct : A struct with values assigned to fields, according to the Specification and Format.
%            If this is not captured by the Function in a variable, the contents of Struct are instead assigned to the Function's workspace (default practice).
%
% See also:
%   arg(), arg_nogui(), arg_norep(), arg_splice(), arg_sub(), arg_subswitch(), arg_subtoggle()
%
% Notes:
%   1) If the Struct output argument is omitted by the user, the arguments are not returned as a struct but instead directly copied into the Function's workspace.
%
%   2) Someone may call the user's Function with the request to deliver the parameter specification, instead of following the normal execution flow.
%      arg_define() automatically handles this task, except if the user lists a special argument named 'report_args' is in the Arg-Specification. In this case,
%      it is the user's task to check whether the field is non-empty, and if so, terminate the normal execution flow to return report_args as the Function
%      output. The Function may perform additional computations in between, as long as a timely delivery of the Struct remains ensured.
%
%   3) If a struct with a field named 'arg_direct' (set to true) is passed, all type checking, specification parsing, fallback to default values and reporting 
%      functionality are skipped. This is essentially a fast path to call a function when all the defaults have previously been obtained from it
%      via arg_report.
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24


global bcilab;


% --- get Format, Values and Specification ---

if ~iscell(varargin{1})
    % a Format specifier was given ...
    if isempty(varargin{1})
        % ... but it was empty: use default behavior
        fmt = [0 Inf];
    else
        % ... and was nonempty: use it
        fmt = varargin{1};
    end
    varargin = varargin(2:end);
else
    % default behavior (fmt = # of args that may be specified positionally)
    fmt = [0 Inf];
end

vals = varargin{1};         % Values
spec = varargin(2:end);     % Specification

if isa(fmt,'function_handle')
    % Format is a function: run it
    if nargout(fmt) == 1
        vals = fmt(vals);
        fmt = 0;
    else
        [vals,fmt] = feval(fmt,vals);
    end
end


% --- find out the reporting type ---

% usually, the reporting type is 'none', except if called (possibly indirectly) by arg_report('type', ...): in this case, the reporting type is 'type'
% reporting is a special way to call arg_define, which requests the argument specification, so that it can be displayed by GUIs, etc.
% * 'none' normal execution: arg_define passes a Struct of Values to the Function or assigns the Struct's fields to the Function's workspace
% * 'rich' arg_define yields a rich specifier list to arg_report(), basically an array of specifier structs (see arg_specifier for the field names)
% * 'lean' arg_define yields a lean specifier list to arg_report(), basically an array of specifier structs but without alternatives for multi-option specifiers
% * 'vals' arg_define yields a struct of values to arg_report(), wich can subsequently be used as the full specification of arguments to pass to the Function
try
    throw; %#ok<LTARG>                      % faster than error()
catch
    ctx = lasterror; %#ok<LERR>             % contains the call stack
    names = {ctx.stack(3:min(6,end)).name}; % function names at the considered levels of indirection...
    matches = find(strncmp(names,'arg_report_',11)); % ... which start with 'arg_report_'
    if isempty(matches)
        reporting_type = 'none';            % no report requested (default case)
    else
        reporting_type = names{matches(end)}(11+1:end); % the reporting type is the suffix of the deepest arg_report_* function in the call stack
    end
end


% --- find out if we are in direct (specificationless and therefore fast) mode ---

direct_mode = false;
% this mode is only applicable when all arguments can be passed as NVPs/structs
if any(fmt == 0) && strcmp(reporting_type,'none')
    structs = find(cellfun('isclass',vals,'struct'));
    for k = structs
        if isfield(vals{k},'arg_direct') && vals{k}.arg_direct
            direct_mode = true;
            break;
        end
    end
end


if direct_mode

    % --- we are in direct mode: quickly collect NVPs from the arguments and produce a result
    
    % obtain flat NVP list
    nvps = flatten_structs(vals);
    
    % get names & values
    names = nvps(1:2:end);
    values = nvps(2:2:end);
    
    % use only the last assignment for each name
    [s,indices] = sort(names(:));
    indices( strcmp(s((1:end-1)'),s((2:end)'))) = [];
    
    % build a struct
    res = cell2struct(values(indices),names(indices),2);

    % there are no reports...
    report_args = [];
    
else    
        
    if strcmp(reporting_type,'handle')
        % very special report type: 'handle'--> this asks for function handles to nested / scoped functions. 
        % unfortunately, these cannot be obtained using standard MATLAB functionality...
        if ~iscellstr(vals)
            error('The arguments passed for handle report must denote function names.'); end
        unresolved = {};
        for f=1:length(vals)
            % resolve each function name in the caller scope
            funcs{f} = evalin('caller',['@' vals{f}]); 
            % check if the function could be retrieved
            tmp = functions(funcs{f});
            if isempty(tmp.file)
                unresolved{f} = vals{f}; end      
        end
        if ~isempty(unresolved)
            % resolve remaining functions: evaluate the spec
            for k=1:length(spec)
                spec{k} = spec{k}{1}('rich',spec{k}{2}{:}); end
            spec = [spec{:}];
            for f=find(~cellfun('isempty',unresolved))
                % search each one in the specification
                funcs{f} = find_function(spec,vals{f}); end
        end
        % report it
        bcilab.temp.report_args = funcs;
        error('BCILAB:arg:report_args','This (internal) exception is destined to be caught by arg_report(); please do not interfere with it.');        
    end

    % --- evaluate the Specification list ---

    % evaluate the specification or retrieve from cache
    [spec,all_names,joint_names,remap] = hlp_microcache('spec',@evaluate_spec,spec,reporting_type);
    
    % --- map vals to a pure list of name-value pairs (NVPs) & check name consistency ---

    if length(fmt) == 2
        if fmt(1) ~= 0 || fmt(2) <= 0
            error('For two-element formats, the first possibility must be 0 and the second possibility must be > 0.'); end

        % there are two possible options: either 0 arguments are positional, or k arguments are positional
        % assuming that 0 arguments are positional, splice substructs into NVP list (structs are allowed instead of individual NVPs)
        if any(cellfun('isclass',vals(1:2:end),'struct'))
            nvps = flatten_structs(vals);
        else
            nvps = vals;
        end

        % check if all the resulting names are in the set of allowed names (a disambiguation requirement in this case)
        if iscellstr(nvps(1:2:end))
            try
                disallowed_nvp = fast_setdiff(nvps(1:2:end),[joint_names {'arg_selection','arg_direct'}]);
            catch
                disallowed_nvp = setdiff(nvps(1:2:end),[joint_names {'arg_selection','arg_direct'}]);
            end
        else
            disallowed_nvp = '(or the sequence of names and values was confused)';
        end

        if isempty(disallowed_nvp)
            % the assumption was correct: 0 arguments are positional
            fmt = 0;
        else
            % k arguments are positional, and we enfore strict naming for the rest (i.e. names must be in the Specification).
            strict_names = true;
            fmt = fmt(2);
        end
    elseif fmt == 0
        % 0 arguments are positional
        nvps = flatten_structs(vals);
    elseif fmt > 0
        % k arguments are positional, the rest are NVPs (we do not enforce strict naming here)
        strict_names = false;
    else
        error('Negative or NaN formats are not allowed.');
    end

    if fmt > 0
        % the first k arguments are positional

        % Find out if we are being called by another arg_define; in this case, this definition appears inside an arg_sub/arg_*, and the values passed to the arg_define
        % are part of the defaults declaration of one of these. If these defaults are specified positionally, the first k arg_norep() arguments in Specification
        % are implicitly skipped.
        if ~strcmp(reporting_type,'none') && any(strcmp('arg_define',{ctx.stack(2:end).name}));
            % we implicitly skip the leading non-reportable arguments in the case of positional assignment (assuming that these are supplied by the outer function),
            % by shifting the name/value assignment by the appropriate number of places
            shift_positionals = min(fmt,find([spec.reportable],1)-1);
        else
            shift_positionals = 0;
        end

        % get the effective number of positional arguments
        fmt = min(fmt,length(vals)+shift_positionals);

        % the NVPs begin only after the k'th argument (defined by the Format)
        nvps = vals(fmt+1-shift_positionals:end);

        % splice in any structs
        if any(cellfun('isclass',nvps(1:2:end),'struct'))
            nvps = flatten_structs(nvps); end

        % do minimal error checking...
        if ~iscellstr(nvps(1:2:end))
            error('Some of the specified arguments that should be names or structs, are not.'); end

        if strict_names
            % enforce strict names
            try
                disallowed_pos = fast_setdiff(nvps(1:2:end),[joint_names {'arg_selection','arg_direct'}]);
            catch
                disallowed_pos = setdiff(nvps(1:2:end),[joint_names {'arg_selection','arg_direct'}]);
            end
            if ~isempty(disallowed_pos)
                error(['Some of the specified arguments do not appear in the argument specification; ' format_cellstr(disallowed_pos) '.']); end
        end

        try
            % remap the positionals (everything up to the k'th argument) into an NVP list, using the code names
            poss = [cellfun(@(x)x{1},all_names(shift_positionals+1:fmt),'UniformOutput',false); vals(1:fmt-shift_positionals)];
        catch
            if strict_names
                % maybe the user intended to pass 0 positionals, but used some disallowed names
                error(['Apparently, some of the used argument names are not known to the function: ' format_cellstr(disallowed_nvp) '.']);
            else
                error(['The first ' fmt ' arguments must be passed by position, and the remaining ones must be passed by name (either in name-value pairs or structs).']);
            end
        end
        % ... and concatenate them with the remaining NVPs into one big NVP list
        nvps = [poss(:)' nvps];
    end


    % --- assign values to names ---

    for k=1:2:length(nvps)
        if isfield(remap,nvps{k})
            idx = remap.(nvps{k});
            spec(idx) = spec(idx).assigner(spec(idx),nvps{k+1});
        else
            % append it to the spec (note: this might need some optimization... it is better if the spec automatically contains the 
            % arg_selection field)
            tmp = arg_nogui(nvps{k},nvps{k+1});
            spec(end+1) = tmp{1}([],tmp{2}{:});
        end
    end
    
    
    % --- yield a report, if requested  ---

    report_args = [];
    
    if ~strcmp(reporting_type,'none')
        % but deliver only the reportable arguments, and only if the values are not unassigned
        tmp = spec([spec.reportable] & ~strcmp(unassigned,{spec.value}));
        if strcmp(reporting_type,'vals')
            tmp = arg_tovals(tmp); end
        % find out how the Function prefers to respond arg_report() (either it leaves it to arg_define() (default), or, by declaring a report_args argument,
        % it declares to handle the reporting manually (by returning report_args as result, if that is non-empty).
        if any(strcmp('return_args',nvps(1:2:length(nvps))))
            report_args = tmp;
        else
            bcilab.temp.report_args = tmp;
            error('BCILAB:arg:report_args','This (internal) exception is destined to be caught by arg_report(); please do not interfere with it.');
        end
    end

    
    % --- post-process the outputs (mandatory) and create a result struct to pass to the Function ---

    % generate errors for local mandatory arguments
    missing_entries = strcmp(mandatory,{spec.value});
    if any(missing_entries)
        missing_names = cellfun(@(x)x{1},{spec(missing_entries).names},'UniformOutput',false);
        error(['The arguments ' format_cellstr(missing_names) ' were unspecified but are mandatory.']);
    end

    % strip non-returned arguments, and convert it all to a struct of values
    res = arg_tovals(spec);

end

% assing the pack to be returned to report_args field
res.report_args = report_args;

% place the arguments in the caller's workspace, if requested
if nargout==0
    for fn=fieldnames(res)'
        assignin('caller',fn{1},res.(fn{1})); end
end




% substitute any structs in place of a name-value pair into the name-value list
function args = flatten_structs(args)
k = 1;
while k <= length(args)
    if isstruct(args{k})
        tmp = [fieldnames(args{k}) struct2cell(args{k})]';
        args = [args(1:k-1) tmp(:)' args(k+1:end)];
        k = k+numel(tmp);
    else
        k = k+2;
    end
end



% format a non-empty cell-string array into a string
function x = format_cellstr(x)
x = ['{' sprintf('%s, ',x{1:end-1}) x{end} '}'];



% evaluate a specification
function [spec,all_names,joint_names,remap] = evaluate_spec(spec,reporting_type)
if strcmp(reporting_type,'rich')
    subreport_type = 'rich';
else
    subreport_type = 'lean';
end
% evaluate the functions to get (possibly arrays of) specifier structs
for k=1:length(spec)
    spec{k} = spec{k}{1}(subreport_type,spec{k}{2}{:}); end

% concatenate the structs to one big struct array
spec = [spec{:}];

% obtain the argument names and the joined names
all_names = {spec.names};
joint_names = [all_names{:}];

% create a name/index remapping table
remap = struct();
for n=1:length(all_names)
    for k=1:length(all_names{n})
        remap.(all_names{n}{k}) = n; end
end

% check for duplicate argument names in the Specification
sorted_names = sort(joint_names);
duplicates = all_names(strcmp(sorted_names(1:end-1),sorted_names(2:end)));
if ~isempty(duplicates)
    error(['The names ' format_cellstr(duplicates) ' refer two multiple arguments.']); end



% recursively find a function handle by name in a specification
function r = find_function(spec,name)
r = [];
for k=1:length(spec)
    if isa(spec(k).value,'function_handle') && strcmp(char(spec(k).value),name)
        r = spec(k).value; 
        return;
    elseif ~isempty(spec(k).alternatives)
        for n = 1:length(spec(k).alternatives)
            r = find_function(spec(k).alternatives{n},name);
            if ~isempty(r)
                return; end
        end
    elseif ~isempty(spec(k).children)
        r = find_function(spec(k).children,name);
        if ~isempty(r)
            return; end
    end
end