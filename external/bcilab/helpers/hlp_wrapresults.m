function a = hlp_wrapresults(f,varargin)
% Wraps all outputs produced by function f for the given arguments into a cell array.
% Results = hlp_wrapresults(Function, Args)
%
% note: it is not (currently) possible to efficiently determine the number of out-args for a varargout function;
%       in this case, at most 10 outputs are supported.
%
%						Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-03-28

% find out how many arguments f can maximally produce
len = nargout(f);
if len < 0
    % for varargout functions, we assume at most 10
    len = 10; end
% but possibly it produces fewer than that (for the given inputs), 
% so we may have to retry with fewer outargs
while len >= 1
    try
        % try to obtain len results from f()
        [a{1:len}] = f(varargin{:});
        return;
    catch
        % got an exception, check if it is outarg-related
        e = lasterror; %#ok<*LERR>
        if ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs','MATLAB:unassignedOutputs'}))
            % it isn't: rethrow
            rethrow(e); end
        % but if it was, we need to retry with fewer out-arguments
        len = len-1;
    end
end

% len = 0: f produces no outputs
f(varargin{:});
a = {};
