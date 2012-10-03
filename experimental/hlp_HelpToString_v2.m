function helpstring = hlp_HelpToString(func,args,maxwrapcols)
% Automatically create Help Text from the arg specification for a given
% function. This will only generate help for arguments that are reportable
% (i.e., any  arg_norep()  arguments will be ignored)
%
% Inputs:
%
%   func:   the handle to the function for which to generate helptext
%   args:   cell array of mandatory argument inputs for the function
%   maxwrapcols: the max number of columns before wrapping. Use Inf for
%                no wrapping
% Output:
%
%   helpstring: formatted help text which can be pasted as a function
%               header.
%
% Example:
%
% >> helpstring = hlp_HelpToString(@vis_causalBrainMovie3D,{struct([]),struct('Determined_by_data',[],'freqs',[1:10],'erWinCenterTimes',[1:10])});
%
% See Also: arg_report()
%
% Author: Tim Mullen 2011, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if nargin<3
    maxwrapcols = 100;
end

% get the full argument spec for this function
rep = arg_report('rich',func,args);

helpstring = [{'Name:','Information'}; parseHelpText(rep,0,maxwrapcols)];

helpstring = cprintf([helpstring cell(size(helpstring,1), 1)],'-la',true,'-it',true,'-E','','-Lcs','-');



function output = parseHelpText(rep,indentlevel,maxwrapcols)
% This functon generates formatted help text for an argument (and its
% children). 
% rep:              the rich report for an argument
% indentlevel:      the depth level for this argument
% maxwrapcols:      the max number of columns before wrapping. Use Inf for
%                   no wrapping

% indent delimiter
delim = '    ';

if nargin < 3
    maxwrapcols = 100;
end
if nargin < 2
    indentlevel = 0;
end

output = {};


% iterate through all arguments and generate help text
for i=1:length(rep)
    
    if strcmpi(rep(i).names{1},'arg_selection')
        continue;
    end

    % extract the name of the entry
    name = rep(i).names{fastif(length(rep(i).names)==1,1,2)};

    % extract the help text for this entry
    % apply textwrapping if necessary
    helptext = {};
    rep(i).help = strtrim(rep(i).help);
    for j=1:length(rep(i).help)
        helptext = [helptext; textwrap([],rep(i).help(j),maxwrapcols)];
    end
    
    % add the data range and default value to helptext
    range = rep(i).range;
    def   = rep(i).value;
    if isempty(def), def='n/a'; end
    
    if ~isempty(range)
        if iscell(rep(i).range)
            range = ['{' ['''' strrep(cell2str(rep(i).range),',',''',''') ''''] '}'];
            helptext = [helptext; ...
                        sprintf('Possible values: %s', range); ...
                        sprintf('Default value  : ''%s''', def)];
        elseif isnumeric(rep(i).range)
            range = ['[' num2str(rep(i).range) ']'];
            helptext = [helptext; ...
                        sprintf('Input Range  : %s', range); ...
                        sprintf('Default value: %s', fastif(ischar(def),def,num2str(def)))];
        else
            range = char(rep(i).range);
            helptext = [helptext; ...
                        sprintf('Possible values: %s', range); ...
                        sprintf('Default value  : ''%s''', def)];
        end
        
    end
    
    % add the input type to help text
    type  = simplifyDataType(rep(i).type);
    helptext = [helptext; sprintf('Input Data Type: %s',type)];

    % construct help text for this entry
%     output = cprintf({[name ':'], arg.(name)},'-d','\t','-cla',true,'-O',@(x) cprintf(x,'-E','','-d','\t','-la',true), ...
%                          '-L',repmat('\t',[1 indentlevel]));

%     output = [output; {name,cprintf(char(helptext),'-cla',true,'-O',@(x) cprintf(x,'-E','','-d','\t','-la',true),'-L',repmat('\t',[1 indentlevel])) } ];
    
    % create blank leading "row"
    output = [output; {'',''}];
    
    % fill row with [Name:   Help Text]
    output = [output; {[name ':'],char(helptext)} ];
    
    % check if this has children and if so, create underline
    if ismember(char(rep(i).head),{'arg_sub','arg_subswitch','arg_subtoggle'})
        output = [output; {repmat('-',1,length(name)+1), ''}];
    end
    
    % output = sprintf('%s%s''%s'':\t%s\n\t\tType: %s\n\t\tRange: %s\n',output,prefix,argname,helptext,rep.type,range);

    % now parse children and insert as needed
    if ~isempty(rep(i).alternatives)
        for k=2:length(rep(i).alternatives)
            for child=1:length(rep(i).alternatives{k})
                ret    =  parseHelpText(rep(i).alternatives{k}(child),indentlevel,maxwrapcols);
                if ~isempty(ret)
                    for q=1:size(ret,1)
                        % insert indented child argument block
                        output = [output; {[repmat(delim,[1 indentlevel+1]) ret{q,1}] ret{q,2}}];
                    end
                end
            end
        end
    end

    if ~isempty(rep(i).children)
        % parse children help text
        for child=1:length(rep(i).children)
            ret    =  parseHelpText(rep(i).children(child),indentlevel,maxwrapcols);
            if ~isempty(ret)
                for q=1:size(ret,1)
                    % insert indented child argument block
                    output = [output; {[repmat(delim,[1 indentlevel+1]) ret{q,1}] ret{q,2}}];
                end
            end
        end
    end
    
end

    



function c = cell2str(C)

c = '';

for k=1:length(C)
    c = [c char(C{k}) ','];
end

c(end) = [];


% convert coded datatypes into more human-readable expressions
function plaintext = simplifyDataType(datatype)

    switch datatype
        case 'logical'
            plaintext = 'boolean';
        case 'char'
            plaintext = 'string';
        case 'int8'
            plaintext = 'integer (int8)';
        case 'uint8'
            plaintext = 'integer (uint8)';
        case 'int16' 
            plaintext = 'integer (int16)';
        case 'uint16'
            plaintext = 'integer (uint16)';
        case 'int32'
            plaintext = 'integer (int32)';
        case 'int64'
            plaintext = 'integer (int64)';
        case 'uint64'
            plaintext = 'integer (uint64)';
        case 'denserealsingle'
            plaintext = 'real number (single)';
        case 'denserealdouble'
            plaintext = 'real number (double)';
        case 'densecomplexsingle'
            plaintext = 'complex number (single)';
        case 'densecomplexdouble'
            plaintext = 'complex number (double)';
        case 'sparserealsingle'
            plaintext = 'real number (sparse single)';
        case 'sparserealdouble'
            plaintext = 'real number (sparse double)';
        case 'sparsecomplexsingle'
            plaintext = 'complex number (sparse single)';
        case 'sparsecomplexdouble'
            plaintext = 'complex number (double)';
        case 'cellstr'
            plaintext = 'cell array of strings (cellstr)';
        case 'object'
            plaintext = 'any matlab object';
        case 'expression'
            plaintext = 'any evaluable Matlab expression.';
        otherwise
            plaintext = datatype;
    end
            
            
            
            
            
            
%     
% 
% function wstring = wraptext(string,numchar)
%     wstring = '';
%     
%     for i=1:numchar:length(string)
%         chunk = string(i:i+numchar);
%         if ~strcmp(chunk(end),' ')
%             i=strfind(string(
%         end
%         
%         [tmp string] = strtok(string(1:numchar),' ');
%         wstring = [wstring tmp];
%     end
% 




% cprintf example:

% cprintf({'a','b';'c',{'a','r'}},'-d',':','-E','[]','-la',true,'-O',@(x) cprintf(x,'-d',':'))




%================

% 
% 
% 
% 
% % helpstring = hlp_HelpToString(func,args)
% %
% % Automatically create Help Text from the arg specification for a given
% % function.
% %
% % Inputs:
% %
% %   func:   the handle to the function for which to generate helptext
% %   args:   cell array of mandatory argument inputs for the function
% %   
% % Output:
% %
% %   helpstring: formatted help text which can be pasted as a function
% %               header.
% %
% % Example:
% %
% % >> helpstring = hlp_HelpToString(@vis_causalBrainMovie3D,{struct([]),struct('foo',[],'freqs',1:10,'erWinCenterTimes',1:10)});
% %
% % See Also: arg_report()
% %
% % Author: Tim Mullen and Christian Kothe, 2011, SCCN/INC, UCSD. 
% % Email:  tim@sccn.ucsd.edu
% 
% % This function is part of the Source Information Flow Toolbox (SIFT)
% %
% % This program is free software; you can redistribute it and/or modify
% % it under the terms of the GNU General Public License as published by
% % the Free Software Foundation; either version 3 of the License, or
% % (at your option) any later version.
% %
% % This program is distributed in the hope that it will be useful,
% % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% % GNU General Public License for more details.
% %
% % You should have received a copy of the GNU General Public License
% % along with this program; if not, write to the Free Software
% % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% function helpstring = hlp_HelpToString(func,args,maxwrapcols)
% 
% % TODO: use arg_report
% 
% if nargin<3
%     maxwrapcols = 100;
% end
% 
% rep = arg_report('rich',func,args);
% 
% helpstring = parseHelpText(rep,'',maxwrapcols);
% 
% 
% function output = parseHelpText(rep,prefix,maxwrapcols)
% 
% if nargin < 3
%     maxwrapcols = 100;
% end
% if nargin < 2
%     prefix = '';
% end
% 
% output = '';
% for i=1:length(rep)
%     if strcmpi(rep(i).names{1},'arg_selection')
%         continue;
%     end
%     
%     hlp = rep(i).help;
%     helptext = '';
%     for j=1:length(hlp)
%         txt = textwrap([],rep(i).help(j),maxwrapcols);
%         for line=1:length(txt)
%             helptext = [helptext txt{line} '\n'];
%         end
%     end
%     helptext(end) = []; % remove trailing CRLF
%     
%     range = rep(i).range;
%     if isempty(range)
%         range = '';
%     elseif iscell(range)
%         range = ['{' cell2str(range) '}'];
%     elseif isnumeric(range)
%         range = ['[' num2str(range) ']'];
%     else
%         range = char(range);
%     end
%     
%     output = sprintf('%s%s''%s'':\t%s\n\t\tType: %s\n\t\tRange: %s\n',output,prefix,rep(i).names{fastif(length(rep(i).names)==1,1,2)},helptext,rep(i).type,range);
%     for 
%     
%     if ~isempty(rep(i).alternatives)
%         for k=2:length(rep(i).alternatives)
%             output = [output parseHelpText(rep(i).alternatives{k},'---')];
%         end
%     end
%     
%     if ~isempty(rep(i).children)
%         % parse children help text
%         output = [output parseHelpText(rep(i).children,'       ')];
%     end
%     
% end
% 
% 
% 
% function c = cell2str(C)
% 
% c = '';
% 
% for k=1:length(C)
%     c = [c char(C{k}) ','];
% end
% 
% c(end) = [];
% 
% 
% % 
% % function wstring = wraptext(string,numchar)
% %     wstring = '';
% %     
% %     for i=1:numchar:length(string)
% %         chunk = string(i:i+numchar);
% %         if ~strcmp(chunk(end),' ')
% %             i=strfind(string(
% %         end
% %         
% %         [tmp string] = strtok(string(1:numchar),' ');
% %         wstring = [wstring tmp];
% %     end
% % 
