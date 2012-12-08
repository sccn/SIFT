function varargout = hlp_microcache(dom, f, varargin)
% Cache results of functions for repeated calls with the same arguments.
% Results = hlp_microcache(Domain, Function, Arguments...)
%
% This is a very lightweight mechanism to memoize results of functions that are often repeatedly called with
% the same arguments. and is designed for small-scale situations (i.e. the function is called only with a small variety of arguments,
% for example less than 100, and the arguments (as well as results) are not too big (e.g. no large matrices or datasets).
% If too many different arguments are supplied, the function "forgets" the oldest ones. The memory is lost after MATLAB is closed.
% Different places of a program can independently memoize results of their functions, by calling hlp_microcache with their own unique
% 'domain' identifier.
%
% In:
%   Domain    : arbitrary (MATLAB-conformant) string identifier of the cache 'domain'
%               may be used to keep separate matters separate
%
%   Function  : function handle (or lambda function) to compute a result from some arguments
%
%   Arguments : arguments to pass to the function
%
% Out:
%   Result    : output of the function for the given arguments
%
%
%
% Notes:
%   Only *referentially transparent* functions are allowed; if a function can give different outputs
%   for the same arguments, this can lead to subtle bugs; examples include functions that refer to global
%   state or lambda functions that refer to workspace variables. Can be fixed by turning all dependencies
%   of the function into arguments.
%
%   Care is required when lambda functions appear in the function's arguments (even deeply nested).
%   Every instance of such a lambda function is considered unique by current MATLABs, so when a lambda
%   term is re-expressed (instead of merely copied), hlp_microcache will actually re-compute the result,
%   except if the flag lambda_equality is set to true for the given Domain.
%
%
% Examples:
%   str = hlp_microcache('test',@num2str,1.5);   % ... twice as slow as str=numstr(1.5)
%   str = hlp_microcache('test',@num2str,1.5);   % ... 25% faster than str=numstr(1.5)
%   m = hlp_microcache('test',@magic,2000);      % ... same as m=magic(2000)
%   m = hlp_microcache('test',@magic,2000);      % ... 3 orders of magnitude faster than m=magic(2000)
%
%
% Advanced:
%  When called as hlp_microcache(Domain, 'option1', value1, 'option2', value2, ...), the options
%  control the caching policy for the respective domain. Possible options are:
%   'resort_freq' : re-sort cached argument sets by their usage frequency every resort_freq lookups (default: 10)
%                   note: the number is not per domain but per group of equally-sized argument sets.
%   'group_size'  : retain only the group_size most-recently used argument sets (default: 100)
%                   note: the number is not per domain but per group of equally-sized argument sets.
%   'max_item_size': maximum size of items that are cached, in bytes (larger ones are not cached) (default: 100000000)
%   'lambda_equality': whether lambda functions (among Arguments) with identical code are considered equal, ignoring the contents of their bound variables (default: false)
%
% Advanced Examples:
%  hlp_microcache('test', 'resort_freq',30, 'group_size',100);
%
%                                   Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                   2010-06-18

% mc ("microcache"): struct
%   domain-id --> domainpool
% domainpool: struct
%   config --> configuration
%   size-id --> sizepool
% sizepool: struct
%   outs: cell array of {f's outputs}
%   inps: cell array of [{char(f)},varargin,{nargout}]
%   frqs: double array of usage frequencies
%   luse: double array of last use cputime
%   lcnt: lookup count for this pool; used to schedule resorting (by usage frequency)
% configuration: struct
%   resort_freq: re-sort a sizepool every this many lookups (into that pool)
%   group_size: maximum number of entries per size-pool
%   max_item_size: maximum size of items that are cached
%   lambda_equality: whether lambda functions among the arguments are considered equal, ignoring the contents of their bound variables
persistent mc;

[varargout{1:nargout}] = f(varargin{:}); return; % KILL SWITCH
 
% is this a regular call?
if isa(f,'function_handle')    

    % get the current lookup time
    now = cputime;
    
    % compute the key structure
    cf = char(f);
    if cf(1) == '@'
        key_f = cf;
    else
        key_f = f;
    end
    
    try
        % check if we need the char version of the Function (lambda case) or not (regular case)
        
        if mc.(dom).config.lambda_equality
            % lambdas with identical function code are considered equal (even if bound variables differ)
            key = [{key_f},varargin,{nargout}];
            for k = find(cellfun('isclass',key(2:end),'function_handle'))
                ck = char(key{k+1});
                if ck(1) == '@'
                    key{k+1} = ck; end
            end
        else
            % lambdas are generally considered unique (except of copies of each other); the Function argument is a special case
            key = [{key_f},varargin,{nargout}];
        end
    catch
        % the lambda lookup failed...
        mc.(dom).config.lambda_equality = false;
        key = [{key_f},varargin,{nargout}];
    end
    
    % get the size id (sid) of the key (MATLAB keeps track of that for every object)
    keyinfo = whos('key');
    keysid = sprintf('s%.0f',keyinfo.bytes);
    
    try
        % retrieve the pool of size-equivalent objects
        sizepool = mc.(dom).(keysid);
        % search for the key in the pool (checking the most-frequently used keys first)
        for k=1:length(sizepool.inps)
            if isequalwithequalnans(key,sizepool.inps{k}) % (isequalwithequalnans() is neatly optimized)
                % found the key, deliver outputs
                varargout = sizepool.outs{k};
                % update the db record...
                sizepool.frqs(k) = sizepool.frqs(k)+1;
                sizepool.luse(k) = now;
                sizepool.lcnt = sizepool.lcnt+1;
                % resort by lookup frequency every resort_freq lookups
                if sizepool.lcnt > mc.(dom).config.resort_freq
                    [sizepool.frqs,inds] = sort(sizepool.frqs,'descend');
                    sizepool.inps = sizepool.inps(inds);
                    sizepool.outs = sizepool.outs(inds);
                    sizepool.luse = sizepool.luse(inds);
                    sizepool.lcnt = 0;
                end
                % write back
                mc.(dom).(keysid) = sizepool;
                return;
            end
        end
    catch
        % domain+keysid not yet in the cache: create appropriate structures (this is rare)
        sizepool = struct('inps',{{}}, 'outs',{{}}, 'frqs',{[]}, 'luse',{[]}, 'lcnt',{0});
    end
    
    if ~exist('varargout','var')
        % set up the default configuration for a domain, if it's not yet present
        if ~isfield(mc,dom) || ~isfield(mc.(dom),'config')
            mc.(dom).config = struct(); end
        if ~isfield(mc.(dom).config,'resort_freq')
            mc.(dom).config.resort_freq = 10; end
        if ~isfield(mc.(dom).config,'group_size')
            mc.(dom).config.group_size = 100; end
        if ~isfield(mc.(dom).config,'max_item_size')
            mc.(dom).config.max_item_size = 100000000; end
        % did not find the entry in the size pool: compute it
        [varargout{1:nargout}] = f(varargin{:});
        iteminfo = whos('varargout');
        if iteminfo.bytes <= mc.(dom).config.max_item_size
            % add to pool
            sizepool.luse(end+1) = now;
            sizepool.frqs(end+1) = 1;
            sizepool.inps{end+1} = key;
            sizepool.outs{end+1} = varargout;
            sizepool.lcnt = 0;
            % remove least-recently used entries if necessary
            while length(sizepool.inps) > mc.(dom).config.group_size
                [x,idx] = min(sizepool.luse); %#ok<ASGLU>
                sizepool.luse(idx) = [];
                sizepool.frqs(idx) = [];
                sizepool.inps(idx) = [];
                sizepool.outs(idx) = [];
            end
            % write back
            mc.(dom).(keysid) = sizepool;
        end
    end
    
else
    % f and what follows are name-value pairs that define the domain configuration
    varargin = [{f} varargin];
    for k=1:2:length(varargin)
        mc.(dom).config.(varargin{k}) = varargin{k+1}; end
end
