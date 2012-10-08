function help = parse_help(help)
% helper function for the arg* specifiers, to parse the help into a first and second part.

if ischar(help)
    % string: split at the end of the first sentence and put into a cell array
    [a,b] = regexp(help,'\.\s+[A-Z]','once');
    if ~isempty(a)
        help = {help(1:a-1), help(b-1:end)};
    else
        help = {help};
    end
elseif ~iscellstr(help)
    error('The help text must be a string.');
end

% remove trailing dot
if length(help{1}) > 1 && help{1}(end) == '.'
    help{1} = help{1}(1:end-1); end

% check for length limit
if length(help{1}) > 60
    % Note: The first sentence in the description is used in some GUIs which have a size limit;
    %       to prevent text from being cut off, please use a shorter wording in the first sentence.
    %
    %       Also note that if this sentence is not followed by a capital letter, the remaining part 
    %       is not considered separate.
    error(['The executive summary for the given argument is too long (' help{1} ')']); 
end
