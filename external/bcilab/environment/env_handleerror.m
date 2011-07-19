function env_handleerror(e,level)
% Displays a formatted error message for some error object, including a full stack trace.
% env_handleerror(Error, Indent)
%
% In:
%   Error   : error object, as received from lasterror or via a catch clause
%   Indent  : optional indentation level, in characters
%
% Example:
%   try
%     ...
%   catch e
%     env_handleerror(e);
%   end
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22

% compute the appropriate indentation level
if ~exist('level','var')
    level = '';
else
    level = repmat(' ',1,level);
end

try
    % display the message
    for message = explode(e.message,10)
        disp([level message{1}]); end
    disp([level 'occured in']);
    for i = 1:length(e.stack)
        try
            % opentoline is more or less undocumented functionality
            disp([level '  <a href="matlab:opentoline(''' e.stack(i).file ''',' num2str(e.stack(i).line) ')">' e.stack(i).name '</a>: ' num2str(e.stack(i).line)]);
        catch
            disp([level '  <a href="matlab:edit ' e.stack(i).file '">' e.stack(i).name '</a>: ' num2str(e.stack(i).line)]);
        end        
    end
catch,end