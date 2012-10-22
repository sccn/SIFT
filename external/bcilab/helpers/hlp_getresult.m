function res = hlp_getresult(idx,f,varargin)
% Returns the Result-Idx's output of the function, given the supplied arguments.
% Result = hlp_getresult(Result-Idx, Function, Arguments...)
%
%						Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-03-28

[tmp{1:idx}] = f(varargin{:});
res = tmp{idx};