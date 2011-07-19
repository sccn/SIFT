function arg_toworkspace(args,doyield)
% Copy the arguments in the given Struct into the workspace of the calling function. 
% arg_toworkspace(Struct,Yield)
%
% In:
%   Struct : an argument structure, as produced by arg_define
%
%   Yield : When true, this function reports ("yields") the field 'report_args' to the framework, given that it is non-empty.
%
% Implementation Note:
%   The Yield functionality is implemented by means of an exception that is recognized by arg_report().
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-09-24

global bcilab;

% check if we need to return the Struct to the framework
if exist('doyield','var') && doyield && ~isempty(args.report_args)
    bcilab.temp.report_args = args.report_args;
    error('BCILAB:arg:report_args','This (internal) exception is destined to be caught by arg_report(); please do not interfere with it.');
end

% place the variables in the caller's workspace
for fn=fieldnames(args)'
    assignin('caller',fn{1},args.(fn{1})); end
