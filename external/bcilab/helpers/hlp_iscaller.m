function result = hlp_iscaller(func, level)
% Test whether some function is calling this function at some level(s) of indirection.
% Result = hlp_iscaller(Function)
%
% It can be specified what levels of nesting outside the function which runs hlp_iscaller are considered.
%
% In:
%   Function : function handle to a function to be tested
%   Level    : nesting level(s) that shall be tested (default: all)
%              level 1 is the function that invokes hlp_iscaller
%
% Out: 
%   Result   : whether the current code is called by the Caller
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-14

try
    throw; %#ok<LTARG> % fastest way to get an exception
catch
    e = lasterror; %#ok<LERR>
    if ~exist('level','var')
        level = 2:length(e.stack); end
    level = level+1;
    level(level < 1 | level > length(e.stack)) = [];
    result = any(strcmp(char(func),{e.stack(level).name}));
end