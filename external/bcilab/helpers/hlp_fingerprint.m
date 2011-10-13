function fp = hlp_fingerprint(data)
% Make a fingerprint (hash) of the given data structure.
% Fingerprint = hlp_fingerprint(Data)
%
% This includes all contents; however, large arrays (such as EEG.data) are only spot-checked.
%
% In:
%   Data        : some data structure
%
% Out:
%   Fingerprint : an integer that identifies the data
%
% Notes:
%   The fingerprint is not unique and identifies the dataset only with a certain (albeit high) probability.
%
%                                   Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                   2010-04-02

global bcilab;

try    
    if hlp_matlab_version < 707
        % save & override random state
        randstate = rand('state'); %#ok<*RAND>
        rand('state',5183);
    else
        % create a legacy-compatible RandStream for the fingerprinting
        if isfield(bcilab,'temp') && isfield(bcilab.temp,'randstream_fingerprint')
            reset(bcilab.temp.randstream_fingerprint,5183);
        else
            bcilab.temp.randstream_fingerprint = RandStream('swb2712','Seed',5183);
        end
    end
    
    % convert data into a string representation
    data = summarize(data);
    % make sure that it does not contain 0's
    data(data==0) = 'x';
    % obtain a hash code via Java (MATLAB does not support proper integer arithmetic...) 
    fp = typecast(int32(java.lang.String(data).hashCode()),'uint32');
    
catch
    e = lasterror; %#ok<LERR>
    env_handleerror(e);
end

if hlp_matlab_version < 707
    % restore random state
    rand('state',randstate); 
end



% get a recursive string summary of arbitrary data
function x = summarize(x)
if iscell(x)
    % check for standard homogeneous of cell arrays
    if cellfun('isclass',x,'double')
        x = ['cd' char(typecast([size(x) x{:}],'uint8'))];
    elseif cellfun('isclass',x,'char')
        x = ['cc' char(typecast(size(x),'uint8')) x{:}];
    else
        % recurse for heterogenous cell arrays
        tmp = cellfun(@summarize,x,'UniformOutput',false);
        x = ['cg' char(typecast(size(x),'uint8')) tmp{:}];
    end
elseif isnumeric(x)
    if ~isreal(x)
        x = [real(x) imag(x)]; end
    if issparse(x)
        x = [find(x) nonzeros(x)]; end
    if numel(x) <= 256
        % small matrices are hashed completely
        x = char(typecast([size(x) x(:)'],'uint8'));
    else        
        % large matrices are spot-checked
        global bcilab %#ok<TLEV>
        ne = numel(x);
        count = floor(256 + (ne-256)/1000);
        if hlp_matlab_version < 707
            indices = 1+floor((ne-1)*rand(1,count));
        else
            indices = 1+floor((ne-1)*rand(bcilab.temp.randstream_fingerprint,1,count));
        end 
        x = char(typecast([size(x) x(indices)],'uint8'));
    end
elseif isa(x,'function_handle')
    x = ['f(' char(x) ')'];
elseif isstruct(x)
	fn = fieldnames(x)';
    if numel(x) > length(fn)
        % iterate over struct fields
        x = cellfun(@(f) [f ':' summarize({x.(f)}) ';'],fn,'UniformOutput',false);
        x = ['s(' x{:} ')'];
    else
        % iterate over struct elements
        x = ['s(' summarize(fn) ':' summarize(struct2cell(x)) ')'];
    end
elseif ischar(x)
    x = x(:)';
elseif isobject(x)
    x = ['o(' char(class(x)) ':' summarize(struct(x)) ')'];
elseif islogical(x)
    x = summarize(uint8(x));
else
    warning('BCILAB:utl_fingerprint:unknown_type','utl_fingerprint: unknown type!');
    error; %#ok<LTARG>
end
