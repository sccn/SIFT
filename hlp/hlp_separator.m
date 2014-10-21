function s = hlp_separator(len)
% return a horizontal 'separator' string of a desired length:
% e.g. '-------------------\n'
if nargin==0, len = 50; end
s = sprintf('%s\n',repmat('-',1,len));