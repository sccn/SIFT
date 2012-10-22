function x = hlp_flattensearch(x,form)
% Flatten search() clauses in a nested data structure into a flat search() clause. 
% Result = hlp_flattensearch(Expression, Output-Form)
%
% Internal tool used by utl_gridsearch to enable the specification of search parameters using search() clauses.
%
% In:
%   Expression  : some data structure, usually an argument to utl_gridsearch, may or may not contain search clauses.
%
%   Output-Form : form of the output (default: 'search')
%                  * 'search': the output shall be a flattened search clause (or a plain value if no search)
%                  * 'cell': the output shall be a cell array of elements to search over
%
% Out:
%   Result      : a flattened search clause (or plain value), or a cell array of search possibilities.
% 
% See also:
%   search(), utl_gridsearch()
%
%                                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                       2010-06-29

if ~exist('form','var') || isempty(form)
    form = 'search'; end

% recursively factor search expressions out of a data structure, to give an overall search expression or data structure
if is_search(x)
    % search: nothing to do
elseif iscell(x)
    % cell array: create a cartesian product over cell-wise searches
    parts = {x};
    x = cellfun(@hlp_flattensearch,x,'UniformOutput',false);
    for c=vectorize(find(cellfun(@is_search,x))) 
        % replicate the parts once for each search item in the current clause
        partnum = length(parts);
        parts = repmat(parts,1,length(x{c}.parts));
        % and fill in the new item in the appropriate place        
        for j=1:length(parts)
            parts{j}{c} = x{c}.parts{ceil(j/partnum)}; end
    end
    x = search(parts);
elseif isfield(x,{'data','srate','chanlocs','event','epoch'})
    % dataset: do not descend further
elseif isstruct(x)
    % structure: create a cartesian product over field-wise searches
    if isscalar(x)
        % scalar structure
        parts = {x};
        fnames = fieldnames(x);
        x = structfun(@hlp_flattensearch,x,'UniformOutput',false);
        for f=fnames(structfun(@is_search,x))'
            % replicate the parts once for each search item in the current clause
            partnum = length(parts);
            parts = repmat(parts,1,length(x.(f{1}).parts));
            % and fill in the new item in the appropriate place
            for j=1:length(parts)
                parts{j}.(f{1}) = x.(f{1}).parts{ceil(j/partnum)}; end
        end
        x = search(parts);
    elseif ~isempty(x)
        % struct array (either with nested searches or a concatenation of search() expressions):
        % handle as a cell array
        x = hlp_flattensearch(arrayfun(@(s){s},x));
        % and re-concatenate the cell contents of each part of the search expression
        if is_search(x)
            for i=1:length(x.parts)
                x.parts{i} = reshape([x.parts{i}{:}],size(x.parts{i})); end
        else
            x = reshape([x{:}],size(x));
        end
    end
else
    % anything else: wrap into a search
    x = search({x});
end

% format the result
if is_search(x)
    switch form
        case 'cell'
            x = x.parts;
        case 'search'
            if isscalar(x.parts)
                x = x.parts{1}; end
        otherwise
            error('unsupported output form.');
    end
end

function x = vectorize(x)
x = x(:)';
